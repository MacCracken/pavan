use serde::{Deserialize, Serialize};

/// NWS wind chill upper temperature limit (°C).
const WIND_CHILL_TEMP_LIMIT: f64 = 10.0;
/// NWS wind chill minimum wind speed (km/h).
const WIND_CHILL_WIND_LIMIT: f64 = 4.8;
/// NWS wind chill wind speed exponent.
const WIND_CHILL_EXPONENT: f64 = 0.16;

/// A uniform wind field.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct WindField {
    /// Wind velocity vector (m/s) in [x, y, z].
    pub velocity: [f64; 3],
}

impl WindField {
    /// Create a wind field from speed and direction.
    ///
    /// Direction in radians: 0 = +x, π/2 = +y.
    #[must_use]
    #[inline]
    pub fn from_speed_direction(speed: f64, direction_rad: f64) -> Self {
        Self {
            velocity: [
                speed * direction_rad.cos(),
                speed * direction_rad.sin(),
                0.0,
            ],
        }
    }

    /// Wind speed magnitude.
    #[must_use]
    #[inline]
    pub fn speed(&self) -> f64 {
        (self.velocity[0] * self.velocity[0]
            + self.velocity[1] * self.velocity[1]
            + self.velocity[2] * self.velocity[2])
            .sqrt()
    }
}

/// Logarithmic wind profile: velocity at height z.
///
/// V(z) = V_ref × ln(z / z0) / ln(z_ref / z0)
///
/// z0 = surface roughness length (m). Typical: grass=0.03, urban=1.0, sea=0.001.
#[must_use]
#[inline]
pub fn log_wind_profile(v_ref: f64, z: f64, z_ref: f64, roughness_length: f64) -> f64 {
    if z <= roughness_length || z_ref <= roughness_length || roughness_length <= 0.0 {
        return 0.0;
    }
    v_ref * (z / roughness_length).ln() / (z_ref / roughness_length).ln()
}

/// Power law wind profile (simpler alternative to log profile).
///
/// V(z) = V_ref × (z / z_ref)^α
///
/// α = Hellman exponent (typical: 0.14 for open terrain, 0.30 for urban).
#[must_use]
#[inline]
pub fn power_law_wind_profile(v_ref: f64, z: f64, z_ref: f64, alpha: f64) -> f64 {
    if z <= 0.0 || z_ref <= 0.0 {
        return 0.0;
    }
    v_ref * (z / z_ref).powf(alpha)
}

/// Wind chill temperature (NWS formula, valid for T ≤ 10°C and V ≥ 4.8 km/h).
#[must_use]
#[inline]
pub fn wind_chill(temp_celsius: f64, wind_speed_kmh: f64) -> f64 {
    if temp_celsius > WIND_CHILL_TEMP_LIMIT || wind_speed_kmh < WIND_CHILL_WIND_LIMIT {
        return temp_celsius;
    }
    let v_exp = wind_speed_kmh.powf(WIND_CHILL_EXPONENT);
    13.12 + 0.6215 * temp_celsius - 11.37 * v_exp + 0.3965 * temp_celsius * v_exp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wind_field_speed() {
        let w = WindField {
            velocity: [3.0, 4.0, 0.0],
        };
        assert!((w.speed() - 5.0).abs() < 0.001);
    }

    #[test]
    fn wind_from_speed_direction() {
        let w = WindField::from_speed_direction(10.0, 0.0);
        assert!((w.velocity[0] - 10.0).abs() < 0.001);
        assert!(w.velocity[1].abs() < 0.001);
    }

    #[test]
    fn log_profile_increases_with_height() {
        let v10 = log_wind_profile(10.0, 10.0, 10.0, 0.03);
        let v50 = log_wind_profile(10.0, 50.0, 10.0, 0.03);
        assert!(v50 > v10, "wind should increase with height");
    }

    #[test]
    fn log_profile_at_reference_height() {
        let v = log_wind_profile(10.0, 10.0, 10.0, 0.03);
        assert!(
            (v - 10.0).abs() < 0.01,
            "at reference height, speed should equal reference speed"
        );
    }

    #[test]
    fn power_law_at_reference() {
        let v = power_law_wind_profile(10.0, 10.0, 10.0, 0.14);
        assert!((v - 10.0).abs() < 0.001);
    }

    #[test]
    fn wind_chill_makes_it_colder() {
        let wc = wind_chill(-10.0, 30.0);
        assert!(
            wc < -10.0,
            "wind chill should be colder than actual temp, got {wc}"
        );
    }

    #[test]
    fn wind_chill_no_effect_warm() {
        let wc = wind_chill(15.0, 30.0);
        assert!((wc - 15.0).abs() < 0.01, "no wind chill above 10°C");
    }

    // --- Edge cases ---

    #[test]
    fn log_profile_zero_roughness_guard() {
        let v = log_wind_profile(10.0, 50.0, 10.0, 0.0);
        assert_eq!(v, 0.0, "zero roughness should return 0");
    }

    #[test]
    fn log_profile_height_below_roughness_guard() {
        let v = log_wind_profile(10.0, 0.01, 10.0, 0.03);
        assert_eq!(v, 0.0, "height below roughness should return 0");
    }

    #[test]
    fn power_law_zero_height_guard() {
        let v = power_law_wind_profile(10.0, 0.0, 10.0, 0.14);
        assert_eq!(v, 0.0, "zero height should return 0");
    }

    #[test]
    fn power_law_negative_height_guard() {
        let v = power_law_wind_profile(10.0, -5.0, 10.0, 0.14);
        assert_eq!(v, 0.0, "negative height should return 0");
    }

    #[test]
    fn power_law_increases_with_height() {
        let v10 = power_law_wind_profile(10.0, 10.0, 10.0, 0.14);
        let v50 = power_law_wind_profile(10.0, 50.0, 10.0, 0.14);
        assert!(v50 > v10, "wind should increase with height");
    }

    #[test]
    fn wind_chill_low_wind_no_effect() {
        let wc = wind_chill(-10.0, 3.0);
        assert!(
            (wc - (-10.0)).abs() < 0.01,
            "wind below 4.8 km/h should have no chill effect"
        );
    }

    #[test]
    fn wind_chill_known_value() {
        // NWS: -10°C at 30 km/h → approximately -20°C
        let wc = wind_chill(-10.0, 30.0);
        assert!(
            wc > -25.0 && wc < -15.0,
            "wind chill at -10°C/30kmh should be ~-20°C, got {wc}"
        );
    }

    #[test]
    fn wind_field_3d_speed() {
        let w = WindField {
            velocity: [1.0, 2.0, 2.0],
        };
        assert!(
            (w.speed() - 3.0).abs() < 0.001,
            "3D magnitude of (1,2,2) should be 3"
        );
    }

    #[test]
    fn wind_from_speed_direction_90_degrees() {
        let w = WindField::from_speed_direction(10.0, std::f64::consts::FRAC_PI_2);
        assert!(w.velocity[0].abs() < 0.001, "x should be ~0 at 90°");
        assert!(
            (w.velocity[1] - 10.0).abs() < 0.001,
            "y should be ~10 at 90°"
        );
    }
}
