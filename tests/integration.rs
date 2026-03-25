use pavan::airfoil::NacaProfile;
use pavan::boundary;
use pavan::forces::{self, AeroForce};
use pavan::panel;
use pavan::vehicle::AeroBody;
use pavan::wind::WindField;
use pavan::*;

#[test]
fn full_flight_computation() {
    let body = AeroBody::light_aircraft();
    let forces = body.compute_forces(60.0, 2000.0, 5.0_f64.to_radians());
    assert!(forces.lift > 0.0);
    assert!(forces.drag > 0.0);
    assert!(
        forces.lift > forces.drag,
        "L should exceed D at moderate AoA"
    );
}

#[test]
fn atmosphere_consistency() {
    let t = atmosphere::standard_temperature(0.0);
    let p = atmosphere::standard_pressure(0.0);
    let rho = atmosphere::standard_density(0.0);
    // Ideal gas check: P = ρRT
    let p_check = rho * atmosphere::GAS_CONSTANT_AIR * t;
    assert!(
        (p - p_check).abs() < 10.0,
        "ideal gas law should hold at sea level"
    );
}

#[test]
fn airfoil_thickness_reasonable() {
    let profile = NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(100);
    let max_thickness: f64 = upper
        .iter()
        .zip(lower.iter())
        .map(|(u, l)| u.1 - l.1)
        .fold(0.0_f64, f64::max);
    assert!(
        (max_thickness - 0.12).abs() < 0.02,
        "NACA 0012 max thickness should be ~12%, got {max_thickness}"
    );
}

#[test]
fn boundary_layer_transition() {
    let re_lam = 100_000.0;
    let re_turb = 1_000_000.0;
    assert!(!boundary::is_turbulent(re_lam));
    assert!(boundary::is_turbulent(re_turb));
}

// --- Cross-module pipeline tests ---

#[test]
fn atmosphere_to_forces_pipeline() {
    // Full pipeline: altitude → density/temp → viscosity → Reynolds → forces
    let alt = 3000.0;
    let velocity = 70.0;
    let chord = 1.5;

    let rho = atmosphere::standard_density(alt);
    let t = atmosphere::standard_temperature(alt);
    let mu = forces::air_dynamic_viscosity(t);
    let re = forces::reynolds_number(rho, velocity, chord, mu);
    let q = atmosphere::dynamic_pressure(rho, velocity);

    assert!(rho > 0.0 && rho < 1.225);
    assert!(
        re > 1e6,
        "Re should be > 1M for aircraft conditions, got {re}"
    );
    assert!(q > 0.0);

    let cl = lift_coefficient_thin_airfoil(5.0_f64.to_radians());
    let cd = drag_coefficient(0.027, cl, 7.5, 0.8);
    let l = forces::lift(q, 16.2, cl);
    let d = forces::drag(q, 16.2, cd);
    assert!(l > d, "lift should exceed drag at 5° AoA");
}

#[test]
fn atmosphere_consistency_at_multiple_altitudes() {
    // Ideal gas law P = ρRT must hold at all altitudes
    for h in [0.0, 2000.0, 5000.0, 8000.0, 11_000.0, 15_000.0, 20_000.0] {
        let t = atmosphere::standard_temperature(h);
        let p = atmosphere::standard_pressure(h);
        let rho = atmosphere::standard_density(h);
        let p_check = rho * atmosphere::GAS_CONSTANT_AIR * t;
        let rel_err = (p - p_check).abs() / p;
        assert!(
            rel_err < 1e-6,
            "ideal gas law violated at {h}m: P={p}, ρRT={p_check}"
        );
    }
}

#[test]
fn boundary_layer_with_real_reynolds() {
    // Use forces module to compute real Reynolds, then check boundary layer
    let rho = atmosphere::standard_density(0.0);
    let t = atmosphere::standard_temperature(0.0);
    let mu = forces::air_dynamic_viscosity(t);
    let velocity = 50.0;
    let x = 1.0;

    let re = forces::reynolds_number(rho, velocity, x, mu);
    assert!(
        boundary::is_turbulent(re),
        "Re at sea level/50m/s should be turbulent"
    );

    let delta_lam = boundary::blasius_thickness(x, re);
    let delta_turb = boundary::turbulent_thickness(x, re);
    assert!(delta_turb > delta_lam);
    assert!(delta_turb < 0.05, "BL thickness should be reasonable at 1m");
}

#[test]
fn serde_round_trip_aero_force() {
    let f = forces::compute_aero_force(1.225, 100.0, 10.0, 0.5, 0.05, -0.1, 1.5);
    let json = serde_json::to_string(&f).expect("serialize AeroForce");
    let back: AeroForce = serde_json::from_str(&json).expect("deserialize AeroForce");
    assert!((back.lift - f.lift).abs() < f64::EPSILON);
    assert!((back.drag - f.drag).abs() < f64::EPSILON);
    assert!((back.moment - f.moment).abs() < f64::EPSILON);
}

#[test]
fn serde_round_trip_wind_field() {
    let w = WindField::from_speed_direction(15.0, 0.5);
    let json = serde_json::to_string(&w).expect("serialize WindField");
    let back: WindField = serde_json::from_str(&json).expect("deserialize WindField");
    for i in 0..3 {
        assert!((back.velocity[i] - w.velocity[i]).abs() < f64::EPSILON);
    }
}

#[test]
fn serde_round_trip_aero_body() {
    let body = AeroBody::light_aircraft();
    let json = serde_json::to_string(&body).expect("serialize AeroBody");
    let back: AeroBody = serde_json::from_str(&json).expect("deserialize AeroBody");
    assert!((back.reference_area - body.reference_area).abs() < f64::EPSILON);
    assert!((back.cd0 - body.cd0).abs() < f64::EPSILON);
    assert!((back.aspect_ratio - body.aspect_ratio).abs() < f64::EPSILON);
}

#[test]
fn wind_profile_altitude_interaction() {
    // Higher altitude → lower density, but wind profile is independent of density
    let v_ref = 10.0;
    let z_ref = 10.0;
    let z0 = 0.03;

    let v_100m = pavan::wind::log_wind_profile(v_ref, 100.0, z_ref, z0);
    let v_200m = pavan::wind::log_wind_profile(v_ref, 200.0, z_ref, z0);
    assert!(v_200m > v_100m, "wind should increase with height");

    // Dynamic pressure with wind at different altitudes
    let rho_low = atmosphere::standard_density(100.0);
    let rho_high = atmosphere::standard_density(5000.0);
    let q_low = atmosphere::dynamic_pressure(rho_low, v_100m);
    let q_high = atmosphere::dynamic_pressure(rho_high, v_200m);
    // Despite higher wind at altitude, lower density may reduce dynamic pressure
    assert!(q_low > 0.0 && q_high > 0.0);
}

#[test]
fn glider_vs_aircraft_full_envelope() {
    let glider = AeroBody::glider();
    let aircraft = AeroBody::light_aircraft();

    // Glider should have better L/D across a range of AoA
    for deg in [2.0_f64, 3.0, 5.0, 7.0] {
        let alpha = deg.to_radians();
        let fg = glider.compute_forces(30.0, 0.0, alpha);
        let fa = aircraft.compute_forces(30.0, 0.0, alpha);
        let ld_g = fg.lift / fg.drag;
        let ld_a = fa.lift / fa.drag;
        assert!(ld_g > ld_a, "glider should have better L/D at {deg}° AoA");
    }
}

// --- Panel method integration tests ---

#[test]
fn panel_method_airfoil_to_forces_pipeline() {
    let profile = NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(80);
    let panels = panel::panels_from_surface(&upper, &lower);

    let alpha = 5.0_f64.to_radians();
    let sol = panel::solve(&panels, alpha).expect("panel solve");

    let rho = atmosphere::standard_density(0.0);
    let f = forces::compute_aero_force(rho, 60.0, 16.2, sol.cl, sol.cd, sol.cm, 1.5);
    assert!(
        f.lift > 0.0,
        "panel method should produce positive lift at 5°"
    );
    assert!(f.drag >= 0.0);
}

#[test]
fn panel_method_vs_thin_airfoil_comparison() {
    let profile = NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(100);
    let panels = panel::panels_from_surface(&upper, &lower);

    let alpha = 3.0_f64.to_radians();
    let sol = panel::solve(&panels, alpha).expect("solve");
    let cl_thin = lift_coefficient_thin_airfoil(alpha);

    let rel_diff = (sol.cl - cl_thin).abs() / cl_thin;
    assert!(
        rel_diff < 0.3,
        "panel Cl={:.3} vs thin Cl={:.3}, diff={:.0}%",
        sol.cl,
        cl_thin,
        rel_diff * 100.0
    );
}

#[test]
fn panel_method_serde_round_trip() {
    let profile = NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(50);
    let panels = panel::panels_from_surface(&upper, &lower);
    let sol = panel::solve(&panels, 0.0).expect("solve");

    let json = serde_json::to_string(&sol).expect("serialize");
    let back: PanelSolution = serde_json::from_str(&json).expect("deserialize");
    assert!((back.cl - sol.cl).abs() < f64::EPSILON);
    assert_eq!(back.cp.len(), sol.cp.len());
}

#[test]
fn panel_method_alpha_sweep() {
    let profile = NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(80);
    let panels = panel::panels_from_surface(&upper, &lower);

    let alphas: Vec<f64> = (-5..=10).map(|d| (d as f64).to_radians()).collect();
    let results = panel::solve_multi(&panels, &alphas).expect("solve_multi");

    for i in 1..results.len() {
        assert!(
            results[i].cl > results[i - 1].cl,
            "Cl should increase monotonically with AoA"
        );
    }
}
