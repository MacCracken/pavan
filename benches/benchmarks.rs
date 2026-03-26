use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

// --- Atmosphere ---

fn bench_standard_density(c: &mut Criterion) {
    c.bench_function("atmosphere/standard_density", |b| {
        b.iter(|| pavan::atmosphere::standard_density(black_box(5000.0)));
    });
}

fn bench_dynamic_pressure(c: &mut Criterion) {
    c.bench_function("atmosphere/dynamic_pressure", |b| {
        b.iter(|| pavan::atmosphere::dynamic_pressure(black_box(1.225), black_box(100.0)));
    });
}

fn bench_speed_of_sound(c: &mut Criterion) {
    c.bench_function("atmosphere/speed_of_sound", |b| {
        b.iter(|| pavan::atmosphere::speed_of_sound(black_box(288.15)));
    });
}

fn bench_standard_pressure(c: &mut Criterion) {
    c.bench_function("atmosphere/standard_pressure", |b| {
        b.iter(|| pavan::atmosphere::standard_pressure(black_box(5000.0)));
    });
}

fn bench_pressure_altitude(c: &mut Criterion) {
    c.bench_function("atmosphere/pressure_altitude", |b| {
        b.iter(|| pavan::atmosphere::pressure_altitude(black_box(54_048.0)));
    });
}

fn bench_mach_number(c: &mut Criterion) {
    c.bench_function("atmosphere/mach_number", |b| {
        b.iter(|| pavan::atmosphere::mach_number(black_box(250.0), black_box(288.15)));
    });
}

// --- Coefficients ---

fn bench_lift_coefficient(c: &mut Criterion) {
    c.bench_function("coefficients/lift_thin_airfoil", |b| {
        b.iter(|| pavan::coefficients::lift_coefficient_thin_airfoil(black_box(0.087)));
    });
}

fn bench_drag_coefficient(c: &mut Criterion) {
    c.bench_function("coefficients/drag_coefficient", |b| {
        b.iter(|| {
            pavan::coefficients::drag_coefficient(
                black_box(0.027),
                black_box(0.548),
                black_box(7.5),
                black_box(0.8),
            )
        });
    });
}

fn bench_max_lift_to_drag_ratio(c: &mut Criterion) {
    c.bench_function("coefficients/max_ld_ratio", |b| {
        b.iter(|| {
            pavan::coefficients::max_lift_to_drag_ratio(
                black_box(0.01),
                black_box(20.0),
                black_box(0.85),
            )
        });
    });
}

fn bench_cl_at_max_ld(c: &mut Criterion) {
    c.bench_function("coefficients/cl_at_max_ld", |b| {
        b.iter(|| {
            pavan::coefficients::cl_at_max_ld(black_box(0.02), black_box(8.0), black_box(0.8))
        });
    });
}

// --- Forces ---

fn bench_lift(c: &mut Criterion) {
    c.bench_function("forces/lift", |b| {
        b.iter(|| pavan::forces::lift(black_box(6125.0), black_box(16.2), black_box(0.548)));
    });
}

fn bench_drag(c: &mut Criterion) {
    c.bench_function("forces/drag", |b| {
        b.iter(|| pavan::forces::drag(black_box(6125.0), black_box(16.2), black_box(0.035)));
    });
}

fn bench_reynolds(c: &mut Criterion) {
    c.bench_function("forces/reynolds_number", |b| {
        b.iter(|| {
            pavan::forces::reynolds_number(
                black_box(1.225),
                black_box(50.0),
                black_box(1.0),
                black_box(1.8e-5),
            )
        });
    });
}

fn bench_air_dynamic_viscosity(c: &mut Criterion) {
    c.bench_function("forces/air_dynamic_viscosity", |b| {
        b.iter(|| pavan::forces::air_dynamic_viscosity(black_box(288.15)));
    });
}

fn bench_compute_aero_force(c: &mut Criterion) {
    c.bench_function("forces/compute_aero_force", |b| {
        b.iter(|| {
            pavan::forces::compute_aero_force(
                black_box(1.225),
                black_box(100.0),
                black_box(16.2),
                black_box(0.548),
                black_box(0.035),
                black_box(-0.1),
                black_box(1.5),
            )
        });
    });
}

// --- Airfoil ---

fn bench_naca_surface(c: &mut Criterion) {
    let profile = pavan::airfoil::NacaProfile::naca2412();
    c.bench_function("airfoil/naca2412_surface_100pts", |b| {
        b.iter(|| black_box(&profile).surface_coordinates(black_box(100)));
    });
}

// --- Boundary ---

fn bench_blasius_thickness(c: &mut Criterion) {
    c.bench_function("boundary/blasius_thickness", |b| {
        b.iter(|| pavan::boundary::blasius_thickness(black_box(1.0), black_box(500_000.0)));
    });
}

fn bench_turbulent_thickness(c: &mut Criterion) {
    c.bench_function("boundary/turbulent_thickness", |b| {
        b.iter(|| pavan::boundary::turbulent_thickness(black_box(1.0), black_box(1_000_000.0)));
    });
}

fn bench_skin_friction_laminar(c: &mut Criterion) {
    c.bench_function("boundary/skin_friction_laminar", |b| {
        b.iter(|| pavan::boundary::skin_friction_laminar(black_box(500_000.0)));
    });
}

fn bench_skin_friction_turbulent(c: &mut Criterion) {
    c.bench_function("boundary/skin_friction_turbulent", |b| {
        b.iter(|| pavan::boundary::skin_friction_turbulent(black_box(1_000_000.0)));
    });
}

// --- Wind ---

fn bench_log_wind_profile(c: &mut Criterion) {
    c.bench_function("wind/log_profile", |b| {
        b.iter(|| {
            pavan::wind::log_wind_profile(
                black_box(10.0),
                black_box(50.0),
                black_box(10.0),
                black_box(0.03),
            )
        });
    });
}

fn bench_power_law_wind_profile(c: &mut Criterion) {
    c.bench_function("wind/power_law_profile", |b| {
        b.iter(|| {
            pavan::wind::power_law_wind_profile(
                black_box(10.0),
                black_box(50.0),
                black_box(10.0),
                black_box(0.14),
            )
        });
    });
}

fn bench_wind_chill(c: &mut Criterion) {
    c.bench_function("wind/wind_chill", |b| {
        b.iter(|| pavan::wind::wind_chill(black_box(-10.0), black_box(30.0)));
    });
}

// --- Vehicle ---

fn bench_vehicle_forces(c: &mut Criterion) {
    let body = pavan::vehicle::AeroBody::light_aircraft();
    c.bench_function("vehicle/light_aircraft_forces", |b| {
        b.iter(|| body.compute_forces(black_box(60.0), black_box(2000.0), black_box(0.087)));
    });
}

// --- Panel ---

fn bench_panel_solve_50(c: &mut Criterion) {
    let profile = pavan::airfoil::NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(30);
    let panels = pavan::panel::panels_from_surface(&upper, &lower);
    c.bench_function("panel/naca0012_solve_50panels", |b| {
        b.iter(|| pavan::panel::solve(black_box(&panels), black_box(0.087)));
    });
}

fn bench_panel_solve_100(c: &mut Criterion) {
    let profile = pavan::airfoil::NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(60);
    let panels = pavan::panel::panels_from_surface(&upper, &lower);
    c.bench_function("panel/naca0012_solve_100panels", |b| {
        b.iter(|| pavan::panel::solve(black_box(&panels), black_box(0.087)));
    });
}

fn bench_panel_solve_200(c: &mut Criterion) {
    let profile = pavan::airfoil::NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(120);
    let panels = pavan::panel::panels_from_surface(&upper, &lower);
    c.bench_function("panel/naca0012_solve_200panels", |b| {
        b.iter(|| pavan::panel::solve(black_box(&panels), black_box(0.087)));
    });
}

fn bench_panel_solve_multi(c: &mut Criterion) {
    let profile = pavan::airfoil::NacaProfile::naca0012();
    let (upper, lower) = profile.surface_coordinates(60);
    let panels = pavan::panel::panels_from_surface(&upper, &lower);
    let alphas = [0.0, 0.035, 0.052, 0.087, 0.140];
    c.bench_function("panel/solve_multi_5_alphas_100p", |b| {
        b.iter(|| pavan::panel::solve_multi(black_box(&panels), black_box(&alphas)));
    });
}

fn bench_panel_cambered(c: &mut Criterion) {
    let profile = pavan::airfoil::NacaProfile::naca2412();
    let (upper, lower) = profile.surface_coordinates(60);
    let panels = pavan::panel::panels_from_surface(&upper, &lower);
    c.bench_function("panel/naca2412_solve_100panels", |b| {
        b.iter(|| pavan::panel::solve(black_box(&panels), black_box(0.087)));
    });
}

criterion_group!(
    benches,
    // atmosphere
    bench_standard_density,
    bench_dynamic_pressure,
    bench_speed_of_sound,
    bench_standard_pressure,
    bench_pressure_altitude,
    bench_mach_number,
    // coefficients
    bench_lift_coefficient,
    bench_drag_coefficient,
    bench_max_lift_to_drag_ratio,
    bench_cl_at_max_ld,
    // forces
    bench_lift,
    bench_drag,
    bench_reynolds,
    bench_air_dynamic_viscosity,
    bench_compute_aero_force,
    // airfoil
    bench_naca_surface,
    // boundary
    bench_blasius_thickness,
    bench_turbulent_thickness,
    bench_skin_friction_laminar,
    bench_skin_friction_turbulent,
    // wind
    bench_log_wind_profile,
    bench_power_law_wind_profile,
    bench_wind_chill,
    // vehicle
    bench_vehicle_forces,
    // panel
    bench_panel_solve_50,
    bench_panel_solve_100,
    bench_panel_solve_200,
    bench_panel_solve_multi,
    bench_panel_cambered,
    // vlm
    bench_vlm_rect_10x2,
    bench_vlm_rect_20x2,
    bench_vlm_rect_40x2,
    bench_vlm_tapered_20x2,
    bench_vlm_solve_multi,
    bench_vlm_generate_panels,
    // compressible
    bench_isentropic_ratios,
    bench_normal_shock,
    bench_oblique_shock_angle,
    bench_prandtl_meyer,
    bench_fanno,
    bench_mach_from_area_ratio,
    bench_mach_from_prandtl_meyer,
    // stability
    bench_neutral_point,
    bench_flap_effectiveness,
    bench_trim,
    bench_trim_at_speed,
    // propulsion
    bench_jet_thrust,
    bench_propeller_efficiency,
    bench_compute_thrust,
);

// --- Compressible benchmarks ---

fn bench_isentropic_ratios(c: &mut Criterion) {
    c.bench_function("compressible/isentropic_m2", |b| {
        b.iter(|| {
            let m = black_box(2.0);
            let g = black_box(1.4);
            (
                pavan::compressible::isentropic_temperature_ratio(m, g),
                pavan::compressible::isentropic_pressure_ratio(m, g),
                pavan::compressible::isentropic_density_ratio(m, g),
            )
        });
    });
}

fn bench_normal_shock(c: &mut Criterion) {
    c.bench_function("compressible/normal_shock_m2", |b| {
        b.iter(|| {
            let m = black_box(2.0);
            let g = black_box(1.4);
            (
                pavan::compressible::normal_shock_mach(m, g),
                pavan::compressible::normal_shock_pressure_ratio(m, g),
                pavan::compressible::normal_shock_total_pressure_ratio(m, g),
            )
        });
    });
}

fn bench_oblique_shock_angle(c: &mut Criterion) {
    c.bench_function("compressible/oblique_shock_angle", |b| {
        b.iter(|| {
            pavan::compressible::oblique_shock_angle(
                black_box(2.0),
                black_box(0.1745),
                black_box(1.4),
                black_box(false),
            )
        });
    });
}

fn bench_prandtl_meyer(c: &mut Criterion) {
    c.bench_function("compressible/prandtl_meyer_m2", |b| {
        b.iter(|| pavan::compressible::prandtl_meyer_angle(black_box(2.0), black_box(1.4)));
    });
}

fn bench_fanno(c: &mut Criterion) {
    c.bench_function("compressible/fanno_parameter_m05", |b| {
        b.iter(|| pavan::compressible::fanno_parameter(black_box(0.5), black_box(1.4)));
    });
}

fn bench_mach_from_area_ratio(c: &mut Criterion) {
    c.bench_function("compressible/mach_from_area_ratio", |b| {
        b.iter(|| {
            pavan::compressible::mach_from_area_ratio(
                black_box(1.6875),
                black_box(1.4),
                black_box(true),
            )
        });
    });
}

fn bench_mach_from_prandtl_meyer(c: &mut Criterion) {
    c.bench_function("compressible/mach_from_prandtl_meyer", |b| {
        b.iter(|| pavan::compressible::mach_from_prandtl_meyer(black_box(0.4604), black_box(1.4)));
    });
}

// --- Stability benchmarks ---

fn bench_neutral_point(c: &mut Criterion) {
    c.bench_function("stability/neutral_point", |b| {
        b.iter(|| {
            pavan::stability::neutral_point(black_box(5.0), black_box(-1.0), black_box(0.25))
        });
    });
}

fn bench_flap_effectiveness(c: &mut Criterion) {
    c.bench_function("stability/flap_effectiveness", |b| {
        b.iter(|| pavan::stability::flap_effectiveness(black_box(0.25)));
    });
}

fn bench_trim(c: &mut Criterion) {
    c.bench_function("stability/trim", |b| {
        b.iter(|| {
            pavan::stability::trim(
                black_box(0.5),
                black_box(0.05),
                black_box(-1.0),
                black_box(-1.5),
                black_box(5.0),
            )
        });
    });
}

fn bench_trim_at_speed(c: &mut Criterion) {
    c.bench_function("stability/trim_at_speed", |b| {
        b.iter(|| {
            pavan::stability::trim_at_speed(
                black_box(10000.0),
                black_box(60.0),
                black_box(0.0),
                black_box(16.2),
                black_box(0.05),
                black_box(-1.0),
                black_box(-1.5),
                black_box(5.0),
            )
        });
    });
}

// --- Propulsion benchmarks ---

fn bench_jet_thrust(c: &mut Criterion) {
    c.bench_function("propulsion/jet_thrust_at_altitude", |b| {
        b.iter(|| {
            pavan::propulsion::jet_thrust_at_altitude(
                black_box(50000.0),
                black_box(10000.0),
                black_box(0.7),
            )
        });
    });
}

fn bench_propeller_efficiency(c: &mut Criterion) {
    c.bench_function("propulsion/propeller_efficiency", |b| {
        b.iter(|| {
            pavan::propulsion::propeller_efficiency(black_box(0.6), black_box(0.85), black_box(5.0))
        });
    });
}

fn bench_compute_thrust(c: &mut Criterion) {
    let jet = pavan::propulsion::PropulsionType::Jet {
        thrust_sl: 50000.0,
        tsfc: 2e-5,
    };
    c.bench_function("propulsion/compute_thrust_jet", |b| {
        b.iter(|| {
            pavan::propulsion::compute_thrust(black_box(&jet), black_box(250.0), black_box(10000.0))
        });
    });
}

// --- VLM benchmarks ---

fn bench_vlm_rect_10x2(c: &mut Criterion) {
    let wing = pavan::vlm::WingGeometry::rectangular(6.0, 1.0, 5, 2);
    let panels = pavan::vlm::generate_panels(&wing);
    c.bench_function("vlm/rect_10x2_solve", |b| {
        b.iter(|| {
            pavan::vlm::solve(
                black_box(&panels),
                black_box(&wing),
                black_box(0.087),
                black_box(1.0),
            )
        });
    });
}

fn bench_vlm_rect_20x2(c: &mut Criterion) {
    let wing = pavan::vlm::WingGeometry::rectangular(6.0, 1.0, 10, 2);
    let panels = pavan::vlm::generate_panels(&wing);
    c.bench_function("vlm/rect_20x2_solve", |b| {
        b.iter(|| {
            pavan::vlm::solve(
                black_box(&panels),
                black_box(&wing),
                black_box(0.087),
                black_box(1.0),
            )
        });
    });
}

fn bench_vlm_rect_40x2(c: &mut Criterion) {
    let wing = pavan::vlm::WingGeometry::rectangular(6.0, 1.0, 20, 2);
    let panels = pavan::vlm::generate_panels(&wing);
    c.bench_function("vlm/rect_40x2_solve", |b| {
        b.iter(|| {
            pavan::vlm::solve(
                black_box(&panels),
                black_box(&wing),
                black_box(0.087),
                black_box(1.0),
            )
        });
    });
}

fn bench_vlm_tapered_20x2(c: &mut Criterion) {
    let wing = pavan::vlm::WingGeometry::tapered(6.0, 1.5, 0.5, 10, 2);
    let panels = pavan::vlm::generate_panels(&wing);
    c.bench_function("vlm/tapered_20x2_solve", |b| {
        b.iter(|| {
            pavan::vlm::solve(
                black_box(&panels),
                black_box(&wing),
                black_box(0.087),
                black_box(1.0),
            )
        });
    });
}

fn bench_vlm_solve_multi(c: &mut Criterion) {
    let wing = pavan::vlm::WingGeometry::rectangular(6.0, 1.0, 10, 2);
    let panels = pavan::vlm::generate_panels(&wing);
    let alphas = [0.0, 0.035, 0.052, 0.087, 0.140];
    c.bench_function("vlm/solve_multi_5_alphas_20x2", |b| {
        b.iter(|| {
            pavan::vlm::solve_multi(
                black_box(&panels),
                black_box(&wing),
                black_box(&alphas),
                black_box(1.0),
            )
        });
    });
}

fn bench_vlm_generate_panels(c: &mut Criterion) {
    let wing = pavan::vlm::WingGeometry::rectangular(6.0, 1.0, 20, 4);
    c.bench_function("vlm/generate_panels_40x4", |b| {
        b.iter(|| pavan::vlm::generate_panels(black_box(&wing)));
    });
}

// --- CFD benchmarks (feature-gated) ---

#[cfg(feature = "cfd")]
fn cfd_config_100x50() -> pavan::cfd::AirfoilCfdConfig {
    let mut cfg = pavan::cfd::AirfoilCfdConfig::default_for(50.0, 0.0);
    cfg.grid_nx = 100;
    cfg.grid_ny = 50;
    cfg.domain_size = (10.0, 5.0);
    cfg.dt = 0.001;
    cfg.pressure_iterations = 20;
    cfg.use_multigrid = false;
    cfg.surface_points = 60;
    cfg
}

#[cfg(feature = "cfd")]
fn bench_cfd_new(c: &mut Criterion) {
    let profile = pavan::airfoil::NacaProfile::naca0012();
    let cfg = cfd_config_100x50();
    c.bench_function("cfd/airfoil_new_100x50", |b| {
        b.iter(|| pavan::cfd::AirfoilCfd::new(black_box(&profile), black_box(&cfg)));
    });
}

#[cfg(feature = "cfd")]
fn bench_cfd_step(c: &mut Criterion) {
    let profile = pavan::airfoil::NacaProfile::naca0012();
    let cfg = cfd_config_100x50();
    c.bench_function("cfd/step_100x50", |b| {
        b.iter_batched(
            || pavan::cfd::AirfoilCfd::new(&profile, &cfg).expect("new"),
            |mut cfd| cfd.step().expect("step"),
            criterion::BatchSize::SmallInput,
        );
    });
}

#[cfg(feature = "cfd")]
fn bench_bl_mesh_params(c: &mut Criterion) {
    c.bench_function("cfd/bl_mesh_params", |b| {
        b.iter(|| {
            pavan::cfd::bl_mesh_params(
                black_box(1.0),
                black_box(50.0),
                black_box(0.0),
                black_box(1.0),
                black_box(1.2),
            )
        });
    });
}

#[cfg(feature = "cfd")]
criterion_group!(
    cfd_benches,
    bench_cfd_new,
    bench_cfd_step,
    bench_bl_mesh_params,
);
#[cfg(not(feature = "cfd"))]
criterion_main!(benches);

#[cfg(feature = "cfd")]
criterion_main!(benches, cfd_benches);
