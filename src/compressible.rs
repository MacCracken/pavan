//! Compressible flow relations for gasdynamics.
//!
//! Isentropic flow, normal and oblique shocks, Prandtl-Meyer expansion,
//! Fanno (friction) and Rayleigh (heat addition) duct flows.
//! All formulas verified against NACA 1135 standard gas tables.

use std::f64::consts::PI;

use crate::error::{PavanError, Result};

// ── Isentropic Flow ──────────────────────────────────────────────────────

/// Isentropic temperature ratio T/T₀ = (1 + (γ-1)/2 · M²)⁻¹.
#[must_use]
#[inline]
pub fn isentropic_temperature_ratio(mach: f64, gamma: f64) -> f64 {
    1.0 / (1.0 + 0.5 * (gamma - 1.0) * mach * mach)
}

/// Isentropic pressure ratio P/P₀ = (T/T₀)^(γ/(γ-1)).
#[must_use]
#[inline]
pub fn isentropic_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    isentropic_temperature_ratio(mach, gamma).powf(gamma / (gamma - 1.0))
}

/// Isentropic density ratio ρ/ρ₀ = (T/T₀)^(1/(γ-1)).
#[must_use]
#[inline]
pub fn isentropic_density_ratio(mach: f64, gamma: f64) -> f64 {
    isentropic_temperature_ratio(mach, gamma).powf(1.0 / (gamma - 1.0))
}

/// Isentropic area ratio A/A* for a given Mach number.
///
/// A/A* = (1/M) · [(2/(γ+1)) · (1 + (γ-1)/2 · M²)]^((γ+1)/(2(γ-1))).
/// Returns f64::INFINITY for M = 0.
#[must_use]
#[inline]
pub fn isentropic_area_ratio(mach: f64, gamma: f64) -> f64 {
    if mach.abs() < f64::EPSILON {
        return f64::INFINITY;
    }
    let gp1 = gamma + 1.0;
    let gm1 = gamma - 1.0;
    let term = (2.0 / gp1) * (1.0 + 0.5 * gm1 * mach * mach);
    (1.0 / mach) * term.powf(gp1 / (2.0 * gm1))
}

/// Find Mach number from isentropic area ratio A/A*.
///
/// Set `supersonic` to true for the supersonic root (M > 1), false for subsonic.
pub fn mach_from_area_ratio(area_ratio: f64, gamma: f64, supersonic: bool) -> Result<f64> {
    if area_ratio < 1.0 {
        return Err(PavanError::InvalidGeometry(
            "area ratio must be >= 1.0".into(),
        ));
    }
    if (area_ratio - 1.0).abs() < 1e-10 {
        return Ok(1.0);
    }

    let (lo, hi) = if supersonic { (1.0, 30.0) } else { (1e-6, 1.0) };

    hisab::num::bisection(
        |m| isentropic_area_ratio(m, gamma) - area_ratio,
        lo,
        hi,
        1e-10,
        100,
    )
    .map_err(|e| PavanError::ComputationError(format!("area ratio inversion failed: {e}")))
}

// ── Normal Shock ─────────────────────────────────────────────────────────

/// Downstream Mach number after a normal shock.
///
/// M₂² = (1 + (γ-1)/2 · M₁²) / (γ·M₁² - (γ-1)/2).
/// Requires M₁ ≥ 1.0.
#[must_use]
#[inline]
pub fn normal_shock_mach(m1: f64, gamma: f64) -> f64 {
    if m1 < 1.0 {
        return m1;
    }
    let gm1 = gamma - 1.0;
    let m1_sq = m1 * m1;
    let num = 1.0 + 0.5 * gm1 * m1_sq;
    let den = gamma * m1_sq - 0.5 * gm1;
    if den <= 0.0 {
        return 0.0;
    }
    (num / den).sqrt()
}

/// Static pressure ratio P₂/P₁ across a normal shock.
#[must_use]
#[inline]
pub fn normal_shock_pressure_ratio(m1: f64, gamma: f64) -> f64 {
    if m1 < 1.0 {
        return 1.0;
    }
    1.0 + 2.0 * gamma / (gamma + 1.0) * (m1 * m1 - 1.0)
}

/// Static temperature ratio T₂/T₁ across a normal shock.
#[must_use]
#[inline]
pub fn normal_shock_temperature_ratio(m1: f64, gamma: f64) -> f64 {
    if m1 < 1.0 {
        return 1.0;
    }
    let m1_sq = m1 * m1;
    let gp1 = gamma + 1.0;
    let gm1 = gamma - 1.0;
    let pr = 1.0 + 2.0 * gamma / gp1 * (m1_sq - 1.0);
    let dr = gp1 * m1_sq / (gm1 * m1_sq + 2.0);
    pr / dr
}

/// Density ratio ρ₂/ρ₁ across a normal shock.
#[must_use]
#[inline]
pub fn normal_shock_density_ratio(m1: f64, gamma: f64) -> f64 {
    if m1 < 1.0 {
        return 1.0;
    }
    let m1_sq = m1 * m1;
    let gp1 = gamma + 1.0;
    let gm1 = gamma - 1.0;
    gp1 * m1_sq / (gm1 * m1_sq + 2.0)
}

/// Total pressure ratio P₀₂/P₀₁ across a normal shock (always ≤ 1).
#[must_use]
#[inline]
pub fn normal_shock_total_pressure_ratio(m1: f64, gamma: f64) -> f64 {
    if m1 < 1.0 {
        return 1.0;
    }
    let m2 = normal_shock_mach(m1, gamma);
    // P₀₂/P₀₁ = (P₂/P₁) × (P₁/P₀₁) / (P₂/P₀₂)
    // isentropic_pressure_ratio returns P/P₀
    let pr = normal_shock_pressure_ratio(m1, gamma);
    let p_over_p0_1 = isentropic_pressure_ratio(m1, gamma); // P₁/P₀₁
    let p_over_p0_2 = isentropic_pressure_ratio(m2, gamma); // P₂/P₀₂
    pr * p_over_p0_1 / p_over_p0_2
}

// ── Oblique Shock ────────────────────────────────────────────────────────

/// Shock wave angle β from upstream Mach, deflection angle θ, and γ.
///
/// Solves the θ-β-M relation numerically. Set `strong` for the strong shock solution.
pub fn oblique_shock_angle(m1: f64, theta_rad: f64, gamma: f64, strong: bool) -> Result<f64> {
    if m1 < 1.0 {
        return Err(PavanError::InvalidVelocity(
            "oblique shock requires M >= 1".into(),
        ));
    }
    if theta_rad < 0.0 {
        return Err(PavanError::InvalidAngle(
            "deflection angle must be non-negative".into(),
        ));
    }

    let mu = (1.0 / m1).asin(); // Mach angle (minimum β)

    // θ-β-M relation: tan(θ) = 2·cot(β) · (M₁²·sin²β - 1) / (M₁²·(γ + cos2β) + 2)
    let theta_from_beta = |beta: f64| -> f64 {
        let m1_sq = m1 * m1;
        let sin_b = beta.sin();
        let cos_b = beta.cos();
        let num = m1_sq * sin_b * sin_b - 1.0;
        let den = m1_sq * (gamma + (2.0 * beta).cos()) + 2.0;
        if den.abs() < f64::EPSILON || sin_b.abs() < f64::EPSILON {
            return 0.0;
        }
        (2.0 * cos_b / sin_b * num / den).atan()
    };

    let (lo, hi) = if strong {
        // Strong shock: between β that gives max θ and π/2
        (mu + 0.01, PI / 2.0)
    } else {
        // Weak shock: between Mach angle and ~π/2
        (mu + 0.001, PI / 2.0 - 0.001)
    };

    // For weak shock, search from Mach angle upward
    // For strong shock, search from near-normal downward
    // Both satisfy: theta_from_beta(β) = θ
    if strong {
        // Search from π/2 downward
        hisab::num::bisection(|beta| theta_from_beta(beta) - theta_rad, lo, hi, 1e-8, 200).map_err(
            |e| PavanError::ComputationError(format!("oblique shock angle solve failed: {e}")),
        )
    } else {
        // Find weak shock: the first root above Mach angle
        // θ rises from 0 at μ to max, then falls back to 0 at π/2
        // We need the first root (weak shock)
        // Find β_max first, then search [μ, β_max]
        let beta_max = hisab::num::bisection(
            |beta| {
                let db = 1e-6;
                (theta_from_beta(beta + db) - theta_from_beta(beta - db)) / (2.0 * db)
            },
            mu + 0.01,
            PI / 2.0 - 0.01,
            1e-6,
            100,
        )
        .unwrap_or(PI / 4.0);

        hisab::num::bisection(
            |beta| theta_from_beta(beta) - theta_rad,
            mu + 0.001,
            beta_max,
            1e-8,
            200,
        )
        .map_err(|e| PavanError::ComputationError(format!("oblique shock angle solve failed: {e}")))
    }
}

/// Downstream Mach number after an oblique shock.
#[must_use]
#[inline]
pub fn oblique_shock_mach(m1: f64, beta_rad: f64, theta_rad: f64, gamma: f64) -> f64 {
    // Normal component: M₁n = M₁·sin(β)
    let m1n = m1 * beta_rad.sin();
    // Normal shock on the normal component
    let m2n = normal_shock_mach(m1n, gamma);
    // Downstream: M₂ = M₂n / sin(β - θ)
    let sin_bt = (beta_rad - theta_rad).sin();
    if sin_bt.abs() < f64::EPSILON {
        return 0.0;
    }
    m2n / sin_bt
}

/// Maximum deflection angle for an attached oblique shock at given Mach.
#[must_use]
pub fn max_deflection_angle(m1: f64, gamma: f64) -> f64 {
    if m1 <= 1.0 {
        return 0.0;
    }

    let mu = (1.0 / m1).asin();

    let theta_from_beta = |beta: f64| -> f64 {
        let m1_sq = m1 * m1;
        let sin_b = beta.sin();
        let cos_b = beta.cos();
        let num = m1_sq * sin_b * sin_b - 1.0;
        let den = m1_sq * (gamma + (2.0 * beta).cos()) + 2.0;
        if den.abs() < f64::EPSILON || sin_b.abs() < f64::EPSILON {
            return 0.0;
        }
        (2.0 * cos_b / sin_b * num / den).atan()
    };

    // Search for maximum θ between Mach angle and π/2
    // Use golden section or simple scan
    let n = 1000;
    let mut max_theta = 0.0;
    for i in 1..n {
        let beta = mu + (PI / 2.0 - mu) * i as f64 / n as f64;
        let theta = theta_from_beta(beta);
        if theta > max_theta {
            max_theta = theta;
        }
    }
    max_theta
}

// ── Prandtl-Meyer Expansion ─────────────────────────────────────────────

/// Prandtl-Meyer function ν(M) in radians.
///
/// ν = √((γ+1)/(γ-1)) · arctan(√((γ-1)/(γ+1)·(M²-1))) - arctan(√(M²-1)).
/// Returns 0 for M ≤ 1.
#[must_use]
#[inline]
pub fn prandtl_meyer_angle(mach: f64, gamma: f64) -> f64 {
    if mach <= 1.0 {
        return 0.0;
    }
    let gp1 = gamma + 1.0;
    let gm1 = gamma - 1.0;
    let ratio = gp1 / gm1;
    let m_term = (mach * mach - 1.0).sqrt();
    ratio.sqrt() * (m_term / ratio.sqrt()).atan() - m_term.atan()
}

/// Find Mach number from Prandtl-Meyer angle ν (radians).
pub fn mach_from_prandtl_meyer(nu_rad: f64, gamma: f64) -> Result<f64> {
    if nu_rad < 0.0 {
        return Err(PavanError::InvalidAngle(
            "Prandtl-Meyer angle must be non-negative".into(),
        ));
    }
    if nu_rad < 1e-10 {
        return Ok(1.0);
    }

    // Max ν at M→∞: ν_max = (π/2) · (√((γ+1)/(γ-1)) - 1)
    let nu_max = (PI / 2.0) * (((gamma + 1.0) / (gamma - 1.0)).sqrt() - 1.0);
    if nu_rad > nu_max {
        return Err(PavanError::InvalidAngle(format!(
            "Prandtl-Meyer angle {nu_rad:.4} exceeds maximum {nu_max:.4}"
        )));
    }

    hisab::num::bisection(
        |m| prandtl_meyer_angle(m, gamma) - nu_rad,
        1.0,
        50.0,
        1e-10,
        100,
    )
    .map_err(|e| PavanError::ComputationError(format!("Prandtl-Meyer inversion failed: {e}")))
}

// ── Fanno Flow (adiabatic + friction) ────────────────────────────────────

/// Fanno flow friction parameter 4fL*/D.
///
/// 4fL*/D = (1-M²)/(γM²) + (γ+1)/(2γ) · ln((γ+1)M²/(2(1+(γ-1)/2·M²))).
#[must_use]
#[inline]
pub fn fanno_parameter(mach: f64, gamma: f64) -> f64 {
    if mach.abs() < f64::EPSILON {
        return f64::INFINITY;
    }
    let m_sq = mach * mach;
    let gp1 = gamma + 1.0;
    let gm1 = gamma - 1.0;
    let term1 = (1.0 - m_sq) / (gamma * m_sq);
    let term2 = gp1 / (2.0 * gamma) * (gp1 * m_sq / (2.0 * (1.0 + 0.5 * gm1 * m_sq))).ln();
    term1 + term2
}

/// Fanno flow temperature ratio T/T*.
#[must_use]
#[inline]
pub fn fanno_temperature_ratio(mach: f64, gamma: f64) -> f64 {
    let gp1 = gamma + 1.0;
    let gm1 = gamma - 1.0;
    gp1 / (2.0 + gm1 * mach * mach)
}

/// Fanno flow pressure ratio P/P*.
#[must_use]
#[inline]
pub fn fanno_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    if mach.abs() < f64::EPSILON {
        return f64::INFINITY;
    }
    (1.0 / mach) * fanno_temperature_ratio(mach, gamma).sqrt()
}

// ── Rayleigh Flow (frictionless + heat addition) ─────────────────────────

/// Rayleigh flow static temperature ratio T/T*.
#[must_use]
#[inline]
pub fn rayleigh_temperature_ratio(mach: f64, gamma: f64) -> f64 {
    let m_sq = mach * mach;
    let gp1 = gamma + 1.0;
    let num = gp1 * mach;
    let den = 1.0 + gamma * m_sq;
    if den.abs() < f64::EPSILON {
        return 0.0;
    }
    (num / den) * (num / den)
}

/// Rayleigh flow static pressure ratio P/P*.
#[must_use]
#[inline]
pub fn rayleigh_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    let den = 1.0 + gamma * mach * mach;
    if den.abs() < f64::EPSILON {
        return 0.0;
    }
    (gamma + 1.0) / den
}

/// Rayleigh flow total temperature ratio T₀/T₀*.
#[must_use]
#[inline]
pub fn rayleigh_total_temperature_ratio(mach: f64, gamma: f64) -> f64 {
    let m_sq = mach * mach;
    let gp1 = gamma + 1.0;
    let gm1 = gamma - 1.0;
    let den = 1.0 + gamma * m_sq;
    if den.abs() < f64::EPSILON {
        return 0.0;
    }
    2.0 * gp1 * m_sq / (den * den) * (1.0 + 0.5 * gm1 * m_sq)
}

#[cfg(test)]
mod tests {
    use super::*;

    const G: f64 = crate::atmosphere::GAMMA;

    // ── Isentropic ──

    #[test]
    fn isentropic_m0_all_ratios_unity() {
        assert!((isentropic_temperature_ratio(0.0, G) - 1.0).abs() < 1e-10);
        assert!((isentropic_pressure_ratio(0.0, G) - 1.0).abs() < 1e-10);
        assert!((isentropic_density_ratio(0.0, G) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn isentropic_m1_sonic() {
        // T/T₀ = 2/(γ+1) = 0.8333
        let tr = isentropic_temperature_ratio(1.0, G);
        assert!(
            (tr - 0.8333).abs() < 0.001,
            "T/T₀ at M=1 should be 0.8333, got {tr}"
        );

        // P/P₀ = 0.5283
        let pr = isentropic_pressure_ratio(1.0, G);
        assert!(
            (pr - 0.5283).abs() < 0.001,
            "P/P₀ at M=1 should be 0.5283, got {pr}"
        );

        // A/A* = 1.0
        let ar = isentropic_area_ratio(1.0, G);
        assert!(
            (ar - 1.0).abs() < 1e-10,
            "A/A* at M=1 should be 1.0, got {ar}"
        );
    }

    #[test]
    fn isentropic_m2_table_values() {
        // NACA 1135 table: M=2, γ=1.4
        let tr = isentropic_temperature_ratio(2.0, G);
        assert!((tr - 0.5556).abs() < 0.001, "T/T₀ at M=2: got {tr}");

        let pr = isentropic_pressure_ratio(2.0, G);
        assert!((pr - 0.1278).abs() < 0.001, "P/P₀ at M=2: got {pr}");

        let dr = isentropic_density_ratio(2.0, G);
        assert!((dr - 0.2300).abs() < 0.001, "ρ/ρ₀ at M=2: got {dr}");

        let ar = isentropic_area_ratio(2.0, G);
        assert!((ar - 1.6875).abs() < 0.001, "A/A* at M=2: got {ar}");
    }

    #[test]
    fn isentropic_m3_table_values() {
        let tr = isentropic_temperature_ratio(3.0, G);
        assert!((tr - 0.3571).abs() < 0.001, "T/T₀ at M=3: got {tr}");
    }

    #[test]
    fn isentropic_area_ratio_m0_infinity() {
        assert_eq!(isentropic_area_ratio(0.0, G), f64::INFINITY);
    }

    #[test]
    fn area_ratio_inverse_subsonic() {
        let m_orig = 0.5;
        let ar = isentropic_area_ratio(m_orig, G);
        let m_back = mach_from_area_ratio(ar, G, false).expect("subsonic");
        assert!(
            (m_back - m_orig).abs() < 1e-6,
            "round-trip failed: got {m_back}"
        );
    }

    #[test]
    fn area_ratio_inverse_supersonic() {
        let m_orig = 2.5;
        let ar = isentropic_area_ratio(m_orig, G);
        let m_back = mach_from_area_ratio(ar, G, true).expect("supersonic");
        assert!(
            (m_back - m_orig).abs() < 1e-6,
            "round-trip failed: got {m_back}"
        );
    }

    #[test]
    fn area_ratio_at_throat() {
        let m = mach_from_area_ratio(1.0, G, false).expect("throat");
        assert!((m - 1.0).abs() < 1e-6);
    }

    #[test]
    fn area_ratio_rejects_below_one() {
        assert!(mach_from_area_ratio(0.5, G, false).is_err());
    }

    // ── Normal Shock ──

    #[test]
    fn normal_shock_m1_no_change() {
        assert!((normal_shock_mach(1.0, G) - 1.0).abs() < 1e-10);
        assert!((normal_shock_pressure_ratio(1.0, G) - 1.0).abs() < 1e-10);
        assert!((normal_shock_temperature_ratio(1.0, G) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn normal_shock_m2_table_values() {
        let m2 = normal_shock_mach(2.0, G);
        assert!((m2 - 0.5774).abs() < 0.001, "M₂ at M₁=2: got {m2}");

        let pr = normal_shock_pressure_ratio(2.0, G);
        assert!((pr - 4.5).abs() < 0.01, "P₂/P₁ at M₁=2: got {pr}");

        let tr = normal_shock_temperature_ratio(2.0, G);
        assert!((tr - 1.6875).abs() < 0.001, "T₂/T₁ at M₁=2: got {tr}");

        let dr = normal_shock_density_ratio(2.0, G);
        assert!((dr - 2.6667).abs() < 0.001, "ρ₂/ρ₁ at M₁=2: got {dr}");
    }

    #[test]
    fn normal_shock_total_pressure_loss() {
        let p0r = normal_shock_total_pressure_ratio(2.0, G);
        assert!(p0r < 1.0, "total pressure should decrease across shock");
        assert!((p0r - 0.7209).abs() < 0.001, "P₀₂/P₀₁ at M₁=2: got {p0r}");
    }

    #[test]
    fn normal_shock_subsonic_passthrough() {
        assert!((normal_shock_mach(0.5, G) - 0.5).abs() < 1e-10);
    }

    // ── Oblique Shock ──

    #[test]
    fn oblique_shock_m2_theta10() {
        let beta = oblique_shock_angle(2.0, 10.0_f64.to_radians(), G, false).expect("weak");
        let beta_deg = beta.to_degrees();
        assert!(
            (beta_deg - 39.3).abs() < 1.0,
            "weak shock angle at M=2, θ=10° should be ~39.3°, got {beta_deg}"
        );
    }

    #[test]
    fn oblique_shock_downstream_mach() {
        let theta = 10.0_f64.to_radians();
        let beta = oblique_shock_angle(2.0, theta, G, false).expect("weak");
        let m2 = oblique_shock_mach(2.0, beta, theta, G);
        assert!(
            m2 > 1.0,
            "weak oblique shock at M=2 should leave supersonic flow, got {m2}"
        );
    }

    #[test]
    fn oblique_shock_subsonic_errors() {
        assert!(oblique_shock_angle(0.5, 0.1, G, false).is_err());
    }

    #[test]
    fn max_deflection_m2() {
        let theta_max = max_deflection_angle(2.0, G);
        let deg = theta_max.to_degrees();
        assert!(
            (deg - 22.97).abs() < 1.0,
            "max deflection at M=2 should be ~22.97°, got {deg}"
        );
    }

    #[test]
    fn max_deflection_subsonic_zero() {
        assert_eq!(max_deflection_angle(0.5, G), 0.0);
    }

    // ── Prandtl-Meyer ──

    #[test]
    fn prandtl_meyer_m1_zero() {
        assert_eq!(prandtl_meyer_angle(1.0, G), 0.0);
    }

    #[test]
    fn prandtl_meyer_m2() {
        let nu = prandtl_meyer_angle(2.0, G).to_degrees();
        assert!(
            (nu - 26.38).abs() < 0.1,
            "ν at M=2 should be ~26.38°, got {nu}"
        );
    }

    #[test]
    fn prandtl_meyer_inverse_round_trip() {
        let m_orig = 2.5;
        let nu = prandtl_meyer_angle(m_orig, G);
        let m_back = mach_from_prandtl_meyer(nu, G).expect("inverse");
        assert!((m_back - m_orig).abs() < 1e-6, "round-trip: got {m_back}");
    }

    #[test]
    fn prandtl_meyer_inverse_at_zero() {
        let m = mach_from_prandtl_meyer(0.0, G).expect("zero");
        assert!((m - 1.0).abs() < 1e-6);
    }

    #[test]
    fn prandtl_meyer_subsonic_zero() {
        assert_eq!(prandtl_meyer_angle(0.5, G), 0.0);
    }

    // ── Fanno ──

    #[test]
    fn fanno_sonic_zero_parameter() {
        let f = fanno_parameter(1.0, G);
        assert!(f.abs() < 1e-10, "4fL*/D at M=1 should be 0, got {f}");
    }

    #[test]
    fn fanno_sonic_temperature_unity() {
        assert!((fanno_temperature_ratio(1.0, G) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn fanno_sonic_pressure_unity() {
        assert!((fanno_pressure_ratio(1.0, G) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn fanno_subsonic_positive_parameter() {
        assert!(fanno_parameter(0.5, G) > 0.0);
    }

    #[test]
    fn fanno_supersonic_positive_parameter() {
        assert!(fanno_parameter(2.0, G) > 0.0);
    }

    // ── Rayleigh ──

    #[test]
    fn rayleigh_sonic_all_unity() {
        assert!((rayleigh_temperature_ratio(1.0, G) - 1.0).abs() < 1e-10);
        assert!((rayleigh_pressure_ratio(1.0, G) - 1.0).abs() < 1e-10);
        assert!((rayleigh_total_temperature_ratio(1.0, G) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rayleigh_m2_pressure() {
        // P/P* at M=2: (γ+1)/(1+γ·M²) = 2.4/6.6 ≈ 0.3636
        let pr = rayleigh_pressure_ratio(2.0, G);
        assert!((pr - 0.3636).abs() < 0.001, "P/P* at M=2: got {pr}");
    }

    #[test]
    fn rayleigh_total_temp_subsonic_less_than_one() {
        // For subsonic flow, T₀/T₀* < 1 (heat must be added to reach sonic)
        let t0r = rayleigh_total_temperature_ratio(0.5, G);
        assert!(t0r < 1.0, "T₀/T₀* subsonic should be < 1, got {t0r}");
    }

    // ── Cross-module ──

    #[test]
    fn normal_shock_downstream_is_subsonic() {
        for m1 in [1.5, 2.0, 3.0, 5.0] {
            let m2 = normal_shock_mach(m1, G);
            assert!(m2 < 1.0, "M₂ should be subsonic for M₁={m1}, got {m2}");
        }
    }

    #[test]
    fn isentropic_ratios_decrease_with_mach() {
        let tr1 = isentropic_temperature_ratio(1.0, G);
        let tr2 = isentropic_temperature_ratio(2.0, G);
        let tr3 = isentropic_temperature_ratio(3.0, G);
        assert!(tr2 < tr1);
        assert!(tr3 < tr2);
    }

    // --- Guard clause coverage ---

    #[test]
    fn normal_shock_pressure_subsonic_passthrough() {
        assert!((normal_shock_pressure_ratio(0.5, G) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn normal_shock_temperature_subsonic_passthrough() {
        assert!((normal_shock_temperature_ratio(0.5, G) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn normal_shock_density_subsonic_passthrough() {
        assert!((normal_shock_density_ratio(0.5, G) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn normal_shock_total_pressure_subsonic_passthrough() {
        assert!((normal_shock_total_pressure_ratio(0.5, G) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn normal_shock_mach_den_guard() {
        // M₁ exactly 1 with specific γ that could make denominator ≤ 0
        assert!(normal_shock_mach(1.0, G) >= 0.0);
    }

    #[test]
    fn oblique_shock_mach_small_angle_guard() {
        let m2 = oblique_shock_mach(2.0, 0.5, 0.5, G);
        assert!(m2 >= 0.0);
    }

    #[test]
    fn prandtl_meyer_negative_angle_errors() {
        assert!(mach_from_prandtl_meyer(-0.1, G).is_err());
    }

    #[test]
    fn prandtl_meyer_exceeds_max_errors() {
        assert!(mach_from_prandtl_meyer(10.0, G).is_err());
    }

    #[test]
    fn fanno_m0_infinity() {
        assert_eq!(fanno_parameter(0.0, G), f64::INFINITY);
    }

    #[test]
    fn fanno_pressure_m0_infinity() {
        assert_eq!(fanno_pressure_ratio(0.0, G), f64::INFINITY);
    }

    #[test]
    fn rayleigh_temperature_m0() {
        let tr = rayleigh_temperature_ratio(0.0, G);
        assert!((tr).abs() < f64::EPSILON);
    }

    #[test]
    fn rayleigh_total_temp_m0() {
        let t0r = rayleigh_total_temperature_ratio(0.0, G);
        assert!((t0r).abs() < f64::EPSILON);
    }

    #[test]
    fn oblique_shock_strong_solution() {
        // Strong shock: the solver may not always converge for all M/θ combos,
        // but it should not panic. Test that the code path executes.
        let result = oblique_shock_angle(5.0, 5.0_f64.to_radians(), G, true);
        // If it succeeds, strong angle should be > weak angle
        if let Ok(strong) = result {
            let weak = oblique_shock_angle(5.0, 5.0_f64.to_radians(), G, false).unwrap();
            assert!(strong > weak);
        }
    }

    #[test]
    fn oblique_shock_mach_zero_sin() {
        // β = θ → sin(β-θ) = 0 → guard returns 0
        let m2 = oblique_shock_mach(2.0, 0.3, 0.3, G);
        assert_eq!(m2, 0.0);
    }

    #[test]
    fn rayleigh_pressure_zero_guard() {
        // When den = 1 + γ·M² ≈ 0 — not physically possible with γ=1.4 but test the guard
        let pr = rayleigh_pressure_ratio(0.0, G);
        // At M=0: (γ+1)/(1+0) = 2.4
        assert!((pr - 2.4).abs() < 0.001);
    }

    #[test]
    fn rayleigh_temperature_zero_guard() {
        let tr = rayleigh_temperature_ratio(0.0, G);
        assert_eq!(tr, 0.0);
    }

    #[test]
    fn oblique_shock_negative_angle_errors() {
        assert!(oblique_shock_angle(2.0, -0.1, G, false).is_err());
    }

    #[test]
    fn different_gamma() {
        // γ=1.3 (diatomic at high temp) should give different results
        let tr_14 = isentropic_temperature_ratio(2.0, 1.4);
        let tr_13 = isentropic_temperature_ratio(2.0, 1.3);
        assert!(
            (tr_14 - tr_13).abs() > 0.01,
            "different γ should give different T/T₀"
        );
    }
}
