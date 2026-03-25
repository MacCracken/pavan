use pavan::{atmosphere, coefficients, vehicle::AeroBody, wind};

fn main() {
    // ISA atmosphere at 3000m
    let alt = 3000.0;
    let t = atmosphere::standard_temperature(alt);
    let p = atmosphere::standard_pressure(alt);
    let rho = atmosphere::standard_density(alt);
    println!("Altitude: {alt}m → T={t:.1}K, P={p:.0}Pa, ρ={rho:.3}kg/m³");

    // Speed of sound and Mach
    let a = atmosphere::speed_of_sound(t);
    let v = 80.0; // 80 m/s
    let m = atmosphere::mach_number(v, t);
    println!("Speed of sound: {a:.1} m/s, V={v} m/s → Mach {m:.3}");

    // Thin airfoil at 5°
    let alpha = 5.0_f64.to_radians();
    let cl = coefficients::lift_coefficient_thin_airfoil(alpha);
    println!("Cl at 5° AoA: {cl:.3}");

    // Light aircraft forces
    let body = AeroBody::light_aircraft();
    let f = body.compute_forces(v, alt, alpha);
    println!(
        "Light aircraft at {v}m/s, {alt}m: L={:.0}N, D={:.0}N, L/D={:.1}",
        f.lift,
        f.drag,
        f.lift / f.drag
    );

    // Wind profile
    let v_ground = 10.0;
    let v_100m = wind::log_wind_profile(v_ground, 100.0, 10.0, 0.03);
    println!("Wind: {v_ground} m/s at 10m → {v_100m:.1} m/s at 100m");
}
