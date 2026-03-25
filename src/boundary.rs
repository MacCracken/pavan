/// Blasius laminar boundary layer thickness.
///
/// δ = 5.0 × x / √Re_x
#[must_use]
#[inline]
pub fn blasius_thickness(x: f64, reynolds_x: f64) -> f64 {
    if reynolds_x <= 0.0 { return 0.0; }
    5.0 * x / reynolds_x.sqrt()
}

/// Turbulent boundary layer thickness (1/7th power law).
///
/// δ = 0.37 × x / Re_x^0.2
#[must_use]
#[inline]
pub fn turbulent_thickness(x: f64, reynolds_x: f64) -> f64 {
    if reynolds_x <= 0.0 { return 0.0; }
    0.37 * x / reynolds_x.powf(0.2)
}

/// Critical Reynolds number for boundary layer transition.
pub const TRANSITION_REYNOLDS: f64 = 500_000.0;

/// Check if flow is turbulent at given Reynolds number.
#[must_use]
#[inline]
pub fn is_turbulent(reynolds: f64) -> bool {
    reynolds > TRANSITION_REYNOLDS
}

/// Skin friction coefficient for laminar flow (Blasius).
///
/// Cf = 0.664 / √Re_x
#[must_use]
#[inline]
pub fn skin_friction_laminar(reynolds_x: f64) -> f64 {
    if reynolds_x <= 0.0 { return 0.0; }
    0.664 / reynolds_x.sqrt()
}

/// Skin friction coefficient for turbulent flow (Schlichting).
///
/// Cf = 0.027 / Re_x^(1/7)
#[must_use]
#[inline]
pub fn skin_friction_turbulent(reynolds_x: f64) -> f64 {
    if reynolds_x <= 0.0 { return 0.0; }
    0.027 / reynolds_x.powf(1.0 / 7.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blasius_basic() {
        // x=1m, Re=100000 → δ = 5/316 ≈ 0.0158m
        let delta = blasius_thickness(1.0, 100_000.0);
        assert!((delta - 0.0158).abs() < 0.001, "Blasius thickness should be ~0.016m, got {delta}");
    }

    #[test]
    fn turbulent_thicker_than_laminar() {
        let re = 1_000_000.0;
        let lam = blasius_thickness(1.0, re);
        let turb = turbulent_thickness(1.0, re);
        assert!(turb > lam, "turbulent BL should be thicker than laminar at same Re");
    }

    #[test]
    fn transition_check() {
        assert!(!is_turbulent(100_000.0));
        assert!(is_turbulent(1_000_000.0));
    }

    #[test]
    fn skin_friction_laminar_decreases_with_re() {
        let cf_low = skin_friction_laminar(10_000.0);
        let cf_high = skin_friction_laminar(1_000_000.0);
        assert!(cf_high < cf_low);
    }

    #[test]
    fn turbulent_friction_higher_than_laminar() {
        let re = 1_000_000.0;
        let cf_lam = skin_friction_laminar(re);
        let cf_turb = skin_friction_turbulent(re);
        assert!(cf_turb > cf_lam, "turbulent friction should be higher");
    }

    #[test]
    fn zero_reynolds_safe() {
        assert_eq!(blasius_thickness(1.0, 0.0), 0.0);
        assert_eq!(turbulent_thickness(1.0, 0.0), 0.0);
        assert_eq!(skin_friction_laminar(0.0), 0.0);
    }
}
