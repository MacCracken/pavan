use serde::{Deserialize, Serialize};

/// Surface coordinate list: Vec of (x, y) points normalized to chord = 1.
pub type SurfacePoints = Vec<(f64, f64)>;

/// NACA 4-digit thickness distribution coefficients.
const NACA_A0: f64 = 0.2969;
const NACA_A1: f64 = 0.1260;
const NACA_A2: f64 = 0.3516;
const NACA_A3: f64 = 0.2843;
const NACA_A4: f64 = 0.1015;

/// NACA 4-digit airfoil profile.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct NacaProfile {
    /// Maximum camber as fraction of chord (first digit / 100).
    pub max_camber: f64,
    /// Location of max camber as fraction of chord (second digit / 10).
    pub camber_position: f64,
    /// Maximum thickness as fraction of chord (last two digits / 100).
    pub max_thickness: f64,
}

impl NacaProfile {
    /// Create from NACA 4 digits (e.g., 2412 → camber=0.02, pos=0.4, thickness=0.12).
    #[must_use]
    #[inline]
    pub fn from_digits(d1: u8, d2: u8, d3: u8, d4: u8) -> Self {
        Self {
            max_camber: d1 as f64 / 100.0,
            camber_position: d2 as f64 / 10.0,
            max_thickness: (d3 * 10 + d4) as f64 / 100.0,
        }
    }

    /// NACA 0012 — symmetric, 12% thickness.
    #[must_use]
    pub fn naca0012() -> Self {
        Self::from_digits(0, 0, 1, 2)
    }

    /// NACA 2412 — 2% camber at 40% chord, 12% thickness.
    #[must_use]
    pub fn naca2412() -> Self {
        Self::from_digits(2, 4, 1, 2)
    }

    /// NACA 4415 — 4% camber at 40% chord, 15% thickness.
    #[must_use]
    pub fn naca4415() -> Self {
        Self::from_digits(4, 4, 1, 5)
    }

    /// Is this a symmetric airfoil (zero camber)?
    #[must_use]
    #[inline]
    pub fn is_symmetric(&self) -> bool {
        self.max_camber.abs() < f64::EPSILON
    }

    /// Generate upper and lower surface coordinates.
    ///
    /// Returns (upper_points, lower_points), each as Vec<(x, y)> normalized to chord = 1.
    #[must_use]
    pub fn surface_coordinates(&self, num_points: usize) -> (SurfacePoints, SurfacePoints) {
        let mut upper = Vec::with_capacity(num_points);
        let mut lower = Vec::with_capacity(num_points);
        let t = self.max_thickness;

        for i in 0..num_points {
            let x = i as f64 / (num_points - 1).max(1) as f64;

            // Thickness distribution (NACA formula)
            let yt = 5.0
                * t
                * (NACA_A0 * x.sqrt() - NACA_A1 * x - NACA_A2 * x * x + NACA_A3 * x * x * x
                    - NACA_A4 * x * x * x * x);

            if self.is_symmetric() {
                upper.push((x, yt));
                lower.push((x, -yt));
            } else {
                let (yc, dyc_dx) = self.camber_at(x);
                let theta = dyc_dx.atan();
                upper.push((x - yt * theta.sin(), yc + yt * theta.cos()));
                lower.push((x + yt * theta.sin(), yc - yt * theta.cos()));
            }
        }

        (upper, lower)
    }

    /// Camber line height and slope at position x (0–1).
    fn camber_at(&self, x: f64) -> (f64, f64) {
        let m = self.max_camber;
        let p = self.camber_position;

        if m.abs() < f64::EPSILON || p.abs() < f64::EPSILON {
            return (0.0, 0.0);
        }

        if x < p {
            let yc = m / (p * p) * (2.0 * p * x - x * x);
            let dyc = 2.0 * m / (p * p) * (p - x);
            (yc, dyc)
        } else {
            let one_minus_p = 1.0 - p;
            let yc = m / (one_minus_p * one_minus_p) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x);
            let dyc = 2.0 * m / (one_minus_p * one_minus_p) * (p - x);
            (yc, dyc)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json;

    #[test]
    fn naca0012_is_symmetric() {
        assert!(NacaProfile::naca0012().is_symmetric());
    }

    #[test]
    fn naca2412_is_not_symmetric() {
        assert!(!NacaProfile::naca2412().is_symmetric());
    }

    #[test]
    fn naca0012_symmetric_surface() {
        let profile = NacaProfile::naca0012();
        let (upper, lower) = profile.surface_coordinates(50);
        // Symmetric: upper y = -lower y at each x
        for (u, l) in upper.iter().zip(lower.iter()) {
            assert!(
                (u.1 + l.1).abs() < 1e-10,
                "symmetric airfoil should have y_upper = -y_lower"
            );
        }
    }

    #[test]
    fn naca2412_has_camber() {
        let profile = NacaProfile::naca2412();
        let (upper, _lower) = profile.surface_coordinates(50);
        // At ~40% chord, upper surface should be above symmetric position
        let mid = &upper[20]; // ~40% of 50 points
        assert!(
            mid.1 > 0.0,
            "cambered airfoil should have positive y on upper surface"
        );
    }

    #[test]
    fn surface_starts_and_ends_near_zero() {
        let profile = NacaProfile::naca0012();
        let (upper, _) = profile.surface_coordinates(100);
        assert!(upper[0].1.abs() < 0.01, "leading edge should be near y=0");
        assert!(
            upper.last().unwrap().1.abs() < 0.01,
            "trailing edge should be near y=0"
        );
    }

    #[test]
    fn from_digits_values() {
        let p = NacaProfile::from_digits(2, 4, 1, 2);
        assert!((p.max_camber - 0.02).abs() < f64::EPSILON);
        assert!((p.camber_position - 0.4).abs() < f64::EPSILON);
        assert!((p.max_thickness - 0.12).abs() < f64::EPSILON);
    }

    #[test]
    fn thickness_at_30_percent() {
        let profile = NacaProfile::naca0012();
        let (upper, lower) = profile.surface_coordinates(100);
        let idx = 30; // 30% chord
        let thickness = upper[idx].1 - lower[idx].1;
        assert!(
            thickness > 0.0 && thickness < 0.15,
            "thickness at 30% should be reasonable, got {thickness}"
        );
    }

    // --- Edge cases ---

    #[test]
    fn naca4415_properties() {
        let p = NacaProfile::naca4415();
        assert!((p.max_camber - 0.04).abs() < f64::EPSILON);
        assert!((p.camber_position - 0.4).abs() < f64::EPSILON);
        assert!((p.max_thickness - 0.15).abs() < f64::EPSILON);
        assert!(!p.is_symmetric());
    }

    #[test]
    fn surface_coordinates_two_points() {
        // Minimum useful: 2 points (leading + trailing edge)
        let profile = NacaProfile::naca0012();
        let (upper, lower) = profile.surface_coordinates(2);
        assert_eq!(upper.len(), 2);
        assert_eq!(lower.len(), 2);
    }

    #[test]
    fn cambered_upper_above_lower() {
        let profile = NacaProfile::naca2412();
        let (upper, lower) = profile.surface_coordinates(50);
        // Mid-chord points: upper should generally be above lower
        for i in 1..49 {
            assert!(
                upper[i].1 > lower[i].1,
                "upper surface should be above lower at point {i}"
            );
        }
    }

    #[test]
    fn naca4415_thicker_than_0012() {
        let thin = NacaProfile::naca0012();
        let thick = NacaProfile::naca4415();
        let (u_thin, l_thin) = thin.surface_coordinates(100);
        let (u_thick, l_thick) = thick.surface_coordinates(100);
        let max_t_thin: f64 = u_thin
            .iter()
            .zip(l_thin.iter())
            .map(|(u, l)| u.1 - l.1)
            .fold(0.0_f64, f64::max);
        let max_t_thick: f64 = u_thick
            .iter()
            .zip(l_thick.iter())
            .map(|(u, l)| u.1 - l.1)
            .fold(0.0_f64, f64::max);
        assert!(
            max_t_thick > max_t_thin,
            "NACA 4415 should be thicker than NACA 0012"
        );
    }

    #[test]
    fn camber_at_zero_camber_position() {
        // NacaProfile with zero camber position → camber_at returns (0,0)
        let p = NacaProfile::from_digits(2, 0, 1, 2);
        let (upper, lower) = p.surface_coordinates(10);
        // Should still produce valid coordinates
        assert!(upper.len() == 10);
        assert!(lower.len() == 10);
    }

    #[test]
    fn serde_round_trip() {
        let profile = NacaProfile::naca2412();
        let json = serde_json::to_string(&profile).expect("serialize");
        let back: NacaProfile = serde_json::from_str(&json).expect("deserialize");
        assert!((back.max_camber - profile.max_camber).abs() < f64::EPSILON);
        assert!((back.camber_position - profile.camber_position).abs() < f64::EPSILON);
        assert!((back.max_thickness - profile.max_thickness).abs() < f64::EPSILON);
    }
}
