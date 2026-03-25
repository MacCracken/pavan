use criterion::{Criterion, black_box, criterion_group, criterion_main};

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
);
criterion_main!(benches);
