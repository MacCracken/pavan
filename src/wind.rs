use serde::{Deserialize, Serialize};

/// A uniform wind field.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindField {
    /// Wind velocity vector (m/s) in [x, y, z].
    pub velocity: [f64; 3],
}

impl WindField {
    /// Create a wind field from speed and direction.
    ///
    /// Direction in radians: 0 = +x, π/2 = +y.
    #[must_use]
    pub fn from_speed_direction(speed: f64, direction_rad: f64) -> Self {
        Self {
            velocity: [speed * direction_rad.cos(), speed * direction_rad.sin(), 0.0],
        }
    }

    /// Wind speed magnitude.
    #[must_use]
    #[inline]
    pub fn speed(&self) -> f64 {
        (self.velocity[0] * self.velocity[0] + self.velocity[1] * self.velocity[1] + self.velocity[2] * self.velocity[2]).sqrt()
    }
}

/// Logarithmic wind profile: velocity at height z.
///
/// V(z) = V_ref × ln(z / z0) / ln(z_ref / z0)
///
/// z0 = surface roughness length (m). Typical: grass=0.03, urban=1.0, sea=0.001.
#[must_use]
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
    if z <= 0.0 || z_ref <= 0.0 { return 0.0; }
    v_ref * (z / z_ref).powf(alpha)
}

/// Wind chill temperature (NWS formula, valid for T ≤ 10°C and V ≥ 4.8 km/h).
#[must_use]
pub fn wind_chill(temp_celsius: f64, wind_speed_kmh: f64) -> f64 {
    if temp_celsius > 10.0 || wind_speed_kmh < 4.8 {
        return temp_celsius;
    }
    let v016 = wind_speed_kmh.powf(0.16);
    13.12 + 0.6215 * temp_celsius - 11.37 * v016 + 0.3965 * temp_celsius * v016
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wind_field_speed() {
        let w = WindField { velocity: [3.0, 4.0, 0.0] };
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
        assert!((v - 10.0).abs() < 0.01, "at reference height, speed should equal reference speed");
    }

    #[test]
    fn power_law_at_reference() {
        let v = power_law_wind_profile(10.0, 10.0, 10.0, 0.14);
        assert!((v - 10.0).abs() < 0.001);
    }

    #[test]
    fn wind_chill_makes_it_colder() {
        let wc = wind_chill(-10.0, 30.0);
        assert!(wc < -10.0, "wind chill should be colder than actual temp, got {wc}");
    }

    #[test]
    fn wind_chill_no_effect_warm() {
        let wc = wind_chill(15.0, 30.0);
        assert!((wc - 15.0).abs() < 0.01, "no wind chill above 10°C");
    }
}
