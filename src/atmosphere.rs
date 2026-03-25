use std::f64::consts::PI;

/// Sea level temperature (K) — ISA standard.
pub const SEA_LEVEL_TEMPERATURE: f64 = 288.15;
/// Sea level pressure (Pa) — ISA standard.
pub const SEA_LEVEL_PRESSURE: f64 = 101_325.0;
/// Sea level density (kg/m³) — ISA standard.
pub const SEA_LEVEL_DENSITY: f64 = 1.225;
/// Temperature lapse rate (K/m) — troposphere.
pub const LAPSE_RATE: f64 = 0.0065;
/// Specific gas constant for dry air (J/(kg·K)).
pub const GAS_CONSTANT_AIR: f64 = 287.058;
/// Ratio of specific heats for air (γ = cp/cv).
pub const GAMMA: f64 = 1.4;
/// Gravitational acceleration (m/s²).
pub const G: f64 = 9.80665;
/// Tropopause altitude (m).
pub const TROPOPAUSE: f64 = 11_000.0;

/// ISA standard temperature at altitude (troposphere, h < 11000m).
///
/// T = 288.15 - 0.0065 × h
#[must_use]
#[inline]
pub fn standard_temperature(altitude_m: f64) -> f64 {
    if altitude_m <= TROPOPAUSE {
        SEA_LEVEL_TEMPERATURE - LAPSE_RATE * altitude_m
    } else {
        // Tropopause: constant 216.65 K
        216.65
    }
}

/// ISA standard pressure at altitude (troposphere).
///
/// P = P₀ × (T / T₀)^(g / (L × R))
#[must_use]
pub fn standard_pressure(altitude_m: f64) -> f64 {
    if altitude_m <= TROPOPAUSE {
        let temp_ratio = standard_temperature(altitude_m) / SEA_LEVEL_TEMPERATURE;
        SEA_LEVEL_PRESSURE * temp_ratio.powf(G / (LAPSE_RATE * GAS_CONSTANT_AIR))
    } else {
        // Above tropopause: exponential decay
        let p_tropo = standard_pressure(TROPOPAUSE);
        let t_tropo = 216.65;
        p_tropo * (-(G * (altitude_m - TROPOPAUSE)) / (GAS_CONSTANT_AIR * t_tropo)).exp()
    }
}

/// ISA standard density at altitude from ideal gas law.
///
/// ρ = P / (R × T)
#[must_use]
#[inline]
pub fn standard_density(altitude_m: f64) -> f64 {
    let t = standard_temperature(altitude_m);
    let p = standard_pressure(altitude_m);
    if t <= 0.0 { return 0.0; }
    p / (GAS_CONSTANT_AIR * t)
}

/// Dynamic pressure: q = ½ρV²
#[must_use]
#[inline]
pub fn dynamic_pressure(density: f64, velocity: f64) -> f64 {
    0.5 * density * velocity * velocity
}

/// Speed of sound in air: a = √(γRT)
#[must_use]
#[inline]
pub fn speed_of_sound(temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 { return 0.0; }
    (GAMMA * GAS_CONSTANT_AIR * temperature_k).sqrt()
}

/// Mach number: M = V / a
#[must_use]
#[inline]
pub fn mach_number(velocity: f64, temperature_k: f64) -> f64 {
    let a = speed_of_sound(temperature_k);
    if a <= 0.0 { return 0.0; }
    velocity / a
}

/// Pressure altitude from given pressure (inverse of standard_pressure, troposphere).
#[must_use]
pub fn pressure_altitude(pressure_pa: f64) -> f64 {
    if pressure_pa <= 0.0 { return 0.0; }
    let exponent = LAPSE_RATE * GAS_CONSTANT_AIR / G;
    SEA_LEVEL_TEMPERATURE / LAPSE_RATE * (1.0 - (pressure_pa / SEA_LEVEL_PRESSURE).powf(exponent))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sea_level_temperature() {
        let t = standard_temperature(0.0);
        assert!((t - 288.15).abs() < 0.01, "sea level T should be 288.15 K, got {t}");
    }

    #[test]
    fn temperature_at_11km() {
        let t = standard_temperature(11_000.0);
        assert!((t - 216.65).abs() < 0.01, "tropopause T should be 216.65 K, got {t}");
    }

    #[test]
    fn temperature_above_tropopause() {
        let t = standard_temperature(15_000.0);
        assert!((t - 216.65).abs() < 0.01, "above tropopause T should be constant 216.65 K");
    }

    #[test]
    fn sea_level_pressure_value() {
        let p = standard_pressure(0.0);
        assert!((p - 101_325.0).abs() < 1.0, "sea level P should be 101325 Pa, got {p}");
    }

    #[test]
    fn pressure_decreases_with_altitude() {
        let p0 = standard_pressure(0.0);
        let p5k = standard_pressure(5000.0);
        let p10k = standard_pressure(10_000.0);
        assert!(p5k < p0);
        assert!(p10k < p5k);
    }

    #[test]
    fn sea_level_density_value() {
        let rho = standard_density(0.0);
        assert!((rho - 1.225).abs() < 0.01, "sea level density should be ~1.225, got {rho}");
    }

    #[test]
    fn dynamic_pressure_basic() {
        // 1.225 kg/m³, 100 m/s → q = 6125 Pa
        let q = dynamic_pressure(1.225, 100.0);
        assert!((q - 6125.0).abs() < 1.0, "q should be ~6125 Pa, got {q}");
    }

    #[test]
    fn speed_of_sound_sea_level() {
        let a = speed_of_sound(288.15);
        assert!((a - 340.3).abs() < 0.5, "speed of sound at sea level should be ~340 m/s, got {a}");
    }

    #[test]
    fn mach_one_at_sea_level() {
        let t = standard_temperature(0.0);
        let a = speed_of_sound(t);
        let m = mach_number(a, t);
        assert!((m - 1.0).abs() < 0.001, "Mach should be 1.0 at speed of sound, got {m}");
    }

    #[test]
    fn subsonic_mach() {
        let t = standard_temperature(0.0);
        let m = mach_number(100.0, t);
        assert!(m < 1.0, "100 m/s at sea level should be subsonic, got M={m}");
    }

    #[test]
    fn pressure_altitude_roundtrip() {
        let h = 5000.0;
        let p = standard_pressure(h);
        let h_back = pressure_altitude(p);
        assert!((h_back - h).abs() < 10.0, "pressure altitude roundtrip should be close, got {h_back}");
    }
}
