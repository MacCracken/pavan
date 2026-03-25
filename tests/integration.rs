use pavan::*;
use pavan::airfoil::NacaProfile;
use pavan::vehicle::AeroBody;
use pavan::boundary;

#[test]
fn full_flight_computation() {
    let body = AeroBody::light_aircraft();
    let forces = body.compute_forces(60.0, 2000.0, 5.0_f64.to_radians());
    assert!(forces.lift > 0.0);
    assert!(forces.drag > 0.0);
    assert!(forces.lift > forces.drag, "L should exceed D at moderate AoA");
}

#[test]
fn atmosphere_consistency() {
    let t = atmosphere::standard_temperature(0.0);
    let p = atmosphere::standard_pressure(0.0);
    let rho = atmosphere::standard_density(0.0);
    // Ideal gas check: P = ρRT
    let p_check = rho * atmosphere::GAS_CONSTANT_AIR * t;
    assert!((p - p_check).abs() < 10.0, "ideal gas law should hold at sea level");
}

#[test]
fn airfoil_thickness_reasonable() {
    let profile = NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(100);
    let max_thickness: f64 = upper.iter().zip(lower.iter())
        .map(|(u, l)| u.1 - l.1)
        .fold(0.0_f64, f64::max);
    assert!((max_thickness - 0.12).abs() < 0.02, "NACA 0012 max thickness should be ~12%, got {max_thickness}");
}

#[test]
fn boundary_layer_transition() {
    let re_lam = 100_000.0;
    let re_turb = 1_000_000.0;
    assert!(!boundary::is_turbulent(re_lam));
    assert!(boundary::is_turbulent(re_turb));
}
