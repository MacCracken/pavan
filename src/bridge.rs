//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into pavan aerodynamics parameters and vice versa.
//!
//! Always available — takes primitive values (f64), no science crate deps.
//!
//! # Architecture
//!
//! ```text
//! badal   (weather)        ──┐
//! pravash (fluid dynamics)   ┼──> bridge ──> pavan aero parameters
//! goonj   (acoustics)       ┤
//! ushma   (thermodynamics)  ┘
//! ```

use crate::atmosphere;

// ── Badal bridges (weather/atmosphere) ─────────────────────────────────────

/// Convert altitude (m) to ISA density (kg/m³), pressure (Pa), and
/// temperature (K).
///
/// Returns `(density, pressure, temperature)`.
#[must_use]
pub fn altitude_to_flight_conditions(altitude_m: f64) -> (f64, f64, f64) {
    let t = atmosphere::standard_temperature(altitude_m);
    let p = atmosphere::standard_pressure(altitude_m);
    let rho = atmosphere::standard_density(altitude_m);
    (rho, p, t)
}

/// Convert atmospheric wind speed (m/s) and aircraft velocity (m/s)
/// to effective airspeed (m/s).
///
/// For headwind: effective = TAS + headwind.
/// `headwind_ms`: positive = headwind, negative = tailwind.
#[must_use]
#[inline]
pub fn wind_to_effective_airspeed(true_airspeed_ms: f64, headwind_ms: f64) -> f64 {
    true_airspeed_ms + headwind_ms
}

// ── Pravash bridges (fluid dynamics) ───────────────────────────────────────

/// Convert surface pressure coefficient (Cp) distribution to a net
/// pressure force per unit span (N/m).
///
/// Integrates Cp around the airfoil using the trapezoidal rule.
/// `cp_values`: Cp at each panel, `panel_lengths`: length of each panel,
/// `panel_normals_y`: y-component of each panel's outward normal.
/// `dynamic_pressure`: q = 0.5 × ρ × V².
#[must_use]
pub fn cp_to_lift_per_span(
    cp_values: &[f64],
    panel_lengths: &[f64],
    panel_normals_y: &[f64],
    dynamic_pressure: f64,
) -> f64 {
    let n = cp_values
        .len()
        .min(panel_lengths.len())
        .min(panel_normals_y.len());
    let mut lift = 0.0;
    for i in 0..n {
        lift -= cp_values[i] * panel_lengths[i] * panel_normals_y[i] * dynamic_pressure;
    }
    lift
}

/// Convert free-stream velocity and fluid properties to Reynolds number.
///
/// Re = ρ × V × L / μ
#[must_use]
#[inline]
pub fn flow_to_reynolds(
    velocity_ms: f64,
    density_kg_m3: f64,
    characteristic_length_m: f64,
    viscosity_pa_s: f64,
) -> f64 {
    if viscosity_pa_s <= 0.0 {
        return 0.0;
    }
    density_kg_m3 * velocity_ms.abs() * characteristic_length_m / viscosity_pa_s
}

// ── Goonj bridges (acoustics) ──────────────────────────────────────────────

/// Convert flow velocity (m/s) and turbulent intensity to estimated
/// aeroacoustic source power level (dB re 1pW).
///
/// Based on Lighthill's analogy: P ∝ ρ × V⁸ / c⁵.
/// `turbulent_intensity`: fraction (typically 0.01–0.1).
#[must_use]
pub fn flow_to_aeroacoustic_power_db(
    velocity_ms: f64,
    density_kg_m3: f64,
    speed_of_sound: f64,
    turbulent_intensity: f64,
    reference_length_m: f64,
) -> f64 {
    if speed_of_sound <= 0.0 || reference_length_m <= 0.0 {
        return 0.0;
    }
    let ti = turbulent_intensity.clamp(0.0, 1.0);
    let mach = velocity_ms.abs() / speed_of_sound;
    // Acoustic power ∝ ρ × V³ × L² × M⁵ × ti²
    let power = density_kg_m3
        * velocity_ms.abs().powi(3)
        * reference_length_m.powi(2)
        * mach.powi(5)
        * ti
        * ti;
    if power <= 0.0 {
        return 0.0;
    }
    // dB re 1 pW (10^-12 W)
    10.0 * (power / 1e-12).log10()
}

/// Convert Mach number to estimated shock noise factor (dB increase).
///
/// Shock-associated noise becomes significant above Mach 0.8.
/// Returns additional dB above broadband turbulence noise.
#[must_use]
#[inline]
pub fn mach_to_shock_noise_db(mach: f64) -> f64 {
    if mach <= 0.8 {
        return 0.0;
    }
    // Empirical: ~20 dB increase from Mach 0.8 to 1.5
    let excess = (mach - 0.8).min(0.7);
    (excess / 0.7) * 20.0
}

// ── Ushma bridges (thermodynamics) ─────────────────────────────────────────

/// Convert stagnation temperature rise (K) to aerodynamic heating rate (W/m²).
///
/// T_stag = T_static × (1 + (γ-1)/2 × M²)
/// Heat flux ≈ h × (T_stag - T_wall) where h is convective coefficient.
#[must_use]
#[inline]
pub fn stagnation_to_heat_flux(
    stag_temperature_k: f64,
    wall_temperature_k: f64,
    convection_coeff_w_m2_k: f64,
) -> f64 {
    convection_coeff_w_m2_k * (stag_temperature_k - wall_temperature_k)
}

/// Convert Mach number and static temperature to stagnation temperature (K).
///
/// T₀ = T × (1 + (γ-1)/2 × M²), with γ = 1.4 for air.
#[must_use]
#[inline]
pub fn mach_to_stagnation_temperature(mach: f64, static_temperature_k: f64) -> f64 {
    static_temperature_k * (1.0 + 0.2 * mach * mach)
}

/// Convert skin friction coefficient and dynamic pressure to surface
/// heat flux (W/m²) using Reynolds analogy.
///
/// q = Cf/2 × ρ × V × c_p × (T_aw - T_wall)
/// Simplified: St ≈ Cf/2 (Reynolds analogy factor = 1).
#[must_use]
#[inline]
pub fn skin_friction_to_heat_flux(
    cf: f64,
    dynamic_pressure_pa: f64,
    velocity_ms: f64,
    specific_heat_j_per_kg_k: f64,
    temp_diff_k: f64,
) -> f64 {
    if velocity_ms <= 0.0 {
        return 0.0;
    }
    let rho_v = 2.0 * dynamic_pressure_pa / velocity_ms;
    (cf / 2.0) * rho_v * specific_heat_j_per_kg_k * temp_diff_k
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Badal ──────────────────────────────────────────────────────────

    #[test]
    fn flight_conditions_sea_level() {
        let (rho, _, t) = altitude_to_flight_conditions(0.0);
        assert!((t - 288.15).abs() < 0.1);
        assert!((rho - 1.225).abs() < 0.01);
    }

    #[test]
    fn effective_airspeed_headwind() {
        let eas = wind_to_effective_airspeed(100.0, 20.0);
        assert!((eas - 120.0).abs() < 0.01);
    }

    #[test]
    fn effective_airspeed_tailwind() {
        let eas = wind_to_effective_airspeed(100.0, -15.0);
        assert!((eas - 85.0).abs() < 0.01);
    }

    // ── Pravash ────────────────────────────────────────────────────────

    #[test]
    fn reynolds_basic() {
        // Air at sea level, 10 m/s, 1m chord
        let re = flow_to_reynolds(10.0, 1.225, 1.0, 1.789e-5);
        assert!((re - 684_739.0).abs() < 100.0);
    }

    #[test]
    fn reynolds_zero_viscosity() {
        assert_eq!(flow_to_reynolds(10.0, 1.225, 1.0, 0.0), 0.0);
    }

    #[test]
    fn cp_to_lift_simple() {
        // Two panels: top has Cp=-1 (suction), bottom has Cp=+0.5 (pressure)
        let cp = [-(1.0), 0.5];
        let lengths = [1.0, 1.0];
        let normals_y = [1.0, -1.0]; // top points up, bottom points down
        let q = 100.0;
        let lift = cp_to_lift_per_span(&cp, &lengths, &normals_y, q);
        // lift = -(-1 × 1 × 1 × 100) - (0.5 × 1 × (-1) × 100) = 100 + 50 = 150
        assert!((lift - 150.0).abs() < 0.1);
    }

    // ── Goonj ──────────────────────────────────────────────────────────

    #[test]
    fn aeroacoustic_power_positive() {
        let db = flow_to_aeroacoustic_power_db(50.0, 1.225, 343.0, 0.05, 1.0);
        assert!(db > 0.0);
    }

    #[test]
    fn aeroacoustic_power_zero_velocity() {
        let db = flow_to_aeroacoustic_power_db(0.0, 1.225, 343.0, 0.05, 1.0);
        assert_eq!(db, 0.0);
    }

    #[test]
    fn shock_noise_subsonic() {
        assert_eq!(mach_to_shock_noise_db(0.5), 0.0);
    }

    #[test]
    fn shock_noise_transonic() {
        let db = mach_to_shock_noise_db(1.2);
        assert!(db > 0.0 && db < 20.0);
    }

    // ── Ushma ──────────────────────────────────────────────────────────

    #[test]
    fn stagnation_temp_mach2() {
        // M=2, T=220K → T₀ = 220 × (1 + 0.2×4) = 220 × 1.8 = 396K
        let t0 = mach_to_stagnation_temperature(2.0, 220.0);
        assert!((t0 - 396.0).abs() < 0.1);
    }

    #[test]
    fn stagnation_heat_flux_basic() {
        let q = stagnation_to_heat_flux(500.0, 300.0, 100.0);
        assert!((q - 20_000.0).abs() < 0.1);
    }

    #[test]
    fn skin_friction_heat() {
        let q = skin_friction_to_heat_flux(0.003, 1000.0, 100.0, 1005.0, 50.0);
        assert!(q > 0.0);
    }

    #[test]
    fn skin_friction_zero_velocity() {
        assert_eq!(
            skin_friction_to_heat_flux(0.003, 1000.0, 0.0, 1005.0, 50.0),
            0.0
        );
    }
}
