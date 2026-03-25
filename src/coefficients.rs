use std::f64::consts::PI;

/// Lift coefficient from thin airfoil theory (small angles of attack).
///
/// Cl = 2π × sin(α) ≈ 2πα for small α
#[must_use]
#[inline]
pub fn lift_coefficient_thin_airfoil(angle_of_attack_rad: f64) -> f64 {
    2.0 * PI * angle_of_attack_rad.sin()
}

/// Total drag coefficient (parasitic + induced).
///
/// Cd = Cd0 + Cl² / (π × e × AR)
///
/// Where Cd0 = zero-lift drag, e = Oswald efficiency (typically 0.7-0.85), AR = aspect ratio.
#[must_use]
#[inline]
pub fn drag_coefficient(cd0: f64, cl: f64, aspect_ratio: f64, oswald_efficiency: f64) -> f64 {
    if aspect_ratio <= 0.0 || oswald_efficiency <= 0.0 {
        return cd0;
    }
    cd0 + cl * cl / (PI * oswald_efficiency * aspect_ratio)
}

/// Lift-to-drag ratio: L/D = Cl / Cd
#[must_use]
#[inline]
pub fn lift_to_drag_ratio(cl: f64, cd: f64) -> f64 {
    if cd.abs() < f64::EPSILON { return 0.0; }
    cl / cd
}

/// Maximum lift-to-drag ratio for a given configuration.
///
/// (L/D)_max = ½ × √(π × e × AR / Cd0)
#[must_use]
#[inline]
pub fn max_lift_to_drag_ratio(cd0: f64, aspect_ratio: f64, oswald_efficiency: f64) -> f64 {
    if cd0 <= 0.0 { return 0.0; }
    0.5 * (PI * oswald_efficiency * aspect_ratio / cd0).sqrt()
}

/// Cl at maximum L/D.
///
/// Cl_opt = √(π × e × AR × Cd0)
#[must_use]
#[inline]
pub fn cl_at_max_ld(cd0: f64, aspect_ratio: f64, oswald_efficiency: f64) -> f64 {
    (PI * oswald_efficiency * aspect_ratio * cd0).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thin_airfoil_at_5_degrees() {
        let alpha = 5.0_f64.to_radians();
        let cl = lift_coefficient_thin_airfoil(alpha);
        assert!((cl - 0.548).abs() < 0.01, "Cl at 5° should be ~0.548, got {cl}");
    }

    #[test]
    fn thin_airfoil_at_zero() {
        let cl = lift_coefficient_thin_airfoil(0.0);
        assert!(cl.abs() < 0.001, "Cl at 0° should be ~0");
    }

    #[test]
    fn thin_airfoil_at_10_degrees() {
        let alpha = 10.0_f64.to_radians();
        let cl = lift_coefficient_thin_airfoil(alpha);
        assert!((cl - 1.09).abs() < 0.02, "Cl at 10° should be ~1.09, got {cl}");
    }

    #[test]
    fn drag_increases_with_lift() {
        let cd_low = drag_coefficient(0.02, 0.2, 8.0, 0.8);
        let cd_high = drag_coefficient(0.02, 1.0, 8.0, 0.8);
        assert!(cd_high > cd_low, "drag should increase with lift");
    }

    #[test]
    fn drag_zero_lift_equals_parasitic() {
        let cd = drag_coefficient(0.02, 0.0, 8.0, 0.8);
        assert!((cd - 0.02).abs() < f64::EPSILON, "zero-lift drag should equal Cd0");
    }

    #[test]
    fn higher_aspect_ratio_less_induced_drag() {
        let cd_low_ar = drag_coefficient(0.02, 0.5, 4.0, 0.8);
        let cd_high_ar = drag_coefficient(0.02, 0.5, 12.0, 0.8);
        assert!(cd_high_ar < cd_low_ar, "higher AR should have less induced drag");
    }

    #[test]
    fn lift_to_drag_ratio_basic() {
        let ld = lift_to_drag_ratio(1.0, 0.05);
        assert!((ld - 20.0).abs() < 0.01);
    }

    #[test]
    fn max_ld_ratio_typical_glider() {
        // Cd0=0.01, AR=20, e=0.85 → (L/D)max ≈ 36.5
        let ld_max = max_lift_to_drag_ratio(0.01, 20.0, 0.85);
        assert!(ld_max > 30.0 && ld_max < 40.0, "glider L/D max should be ~36, got {ld_max}");
    }
}
