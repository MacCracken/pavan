use criterion::{black_box, criterion_group, criterion_main, Criterion};

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

fn bench_lift_coefficient(c: &mut Criterion) {
    c.bench_function("coefficients/lift_thin_airfoil", |b| {
        b.iter(|| pavan::coefficients::lift_coefficient_thin_airfoil(black_box(0.087)));
    });
}

fn bench_reynolds(c: &mut Criterion) {
    c.bench_function("forces/reynolds_number", |b| {
        b.iter(|| pavan::forces::reynolds_number(black_box(1.225), black_box(50.0), black_box(1.0), black_box(1.8e-5)));
    });
}

fn bench_naca_surface(c: &mut Criterion) {
    let profile = pavan::airfoil::NacaProfile::naca2412();
    c.bench_function("airfoil/naca2412_surface_100pts", |b| {
        b.iter(|| black_box(&profile).surface_coordinates(black_box(100)));
    });
}

fn bench_vehicle_forces(c: &mut Criterion) {
    let body = pavan::vehicle::AeroBody::light_aircraft();
    c.bench_function("vehicle/light_aircraft_forces", |b| {
        b.iter(|| body.compute_forces(black_box(60.0), black_box(2000.0), black_box(0.087)));
    });
}

fn bench_log_wind_profile(c: &mut Criterion) {
    c.bench_function("wind/log_profile", |b| {
        b.iter(|| pavan::wind::log_wind_profile(black_box(10.0), black_box(50.0), black_box(10.0), black_box(0.03)));
    });
}

criterion_group!(
    benches,
    bench_standard_density,
    bench_dynamic_pressure,
    bench_lift_coefficient,
    bench_reynolds,
    bench_naca_surface,
    bench_vehicle_forces,
    bench_log_wind_profile,
);
criterion_main!(benches);
