//! Longitudinal static stability, control surfaces, and trim.
//!
//! Neutral point, static margin, elevator effectiveness, flap effectiveness,
//! and trim calculation for conventional aircraft configurations.

use serde::{Deserialize, Serialize};

use crate::atmosphere;
use crate::error::{PavanError, Result};

/// Longitudinal stability parameters for an aircraft configuration.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[non_exhaustive]
pub struct LongitudinalStability {
    /// Neutral point location as fraction of MAC (from leading edge).
    pub neutral_point: f64,
    /// Static margin: (x_np - x_cg) / MAC. Positive = stable.
    pub static_margin: f64,
    /// Pitch moment derivative dCm/dα (1/rad). Negative = stable.
    pub cm_alpha: f64,
    /// Lift curve slope dCl/dα (1/rad).
    pub cl_alpha: f64,
}

/// Control surface definition.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[non_exhaustive]
pub struct ControlSurface {
    /// Flap chord as fraction of local chord (e.g., 0.25 = 25%).
    pub chord_fraction: f64,
    /// Inboard edge as fraction of semi-span.
    pub span_start: f64,
    /// Outboard edge as fraction of semi-span.
    pub span_end: f64,
    /// Maximum deflection angle (rad).
    pub max_deflection_rad: f64,
}

/// Trim solution at a given flight condition.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[non_exhaustive]
pub struct TrimSolution {
    /// Trimmed angle of attack (rad).
    pub alpha_rad: f64,
    /// Elevator deflection for trim (rad).
    pub elevator_rad: f64,
    /// Lift coefficient at trim.
    pub cl: f64,
    /// Drag coefficient at trim (if available, 0 otherwise).
    pub cd: f64,
    /// Residual moment coefficient (should be ~0 for successful trim).
    pub cm_residual: f64,
}

/// Neutral point location as fraction of MAC.
///
/// x_np = x_ac - Cm_α / CL_α.
/// Where x_ac is the aerodynamic center (typically 0.25 for thin airfoils).
#[must_use]
#[inline]
pub fn neutral_point(cl_alpha: f64, cm_alpha: f64, x_ac: f64) -> f64 {
    if cl_alpha.abs() < f64::EPSILON {
        return x_ac;
    }
    x_ac - cm_alpha / cl_alpha
}

/// Static margin: (x_np - x_cg) / MAC. Positive = stable.
#[must_use]
#[inline]
pub fn static_margin(x_np: f64, x_cg: f64) -> f64 {
    x_np - x_cg
}

/// Aircraft Cm_α including wing and horizontal tail contributions.
///
/// Cm_α = CL_α_w · (x_cg - x_ac_w) - η_t · V_H · CL_α_t · (1 - dε/dα)
///
/// Where:
/// - CL_α_w: wing lift curve slope (1/rad)
/// - x_cg, x_ac_w: CG and wing AC as fractions of MAC
/// - η_t: tail efficiency (dynamic pressure ratio, typically 0.9)
/// - V_H: horizontal tail volume coefficient (S_t · l_t / (S_w · MAC))
/// - CL_α_t: tail lift curve slope
/// - dε/dα: downwash gradient at tail (typically 0.3-0.5)
#[must_use]
#[inline]
pub fn cm_alpha(
    cl_alpha_wing: f64,
    x_ac_wing: f64,
    x_cg: f64,
    tail_volume: f64,
    cl_alpha_tail: f64,
    tail_efficiency: f64,
    downwash_gradient: f64,
) -> f64 {
    cl_alpha_wing * (x_cg - x_ac_wing)
        - tail_efficiency * tail_volume * cl_alpha_tail * (1.0 - downwash_gradient)
}

/// Flap effectiveness factor τ from thin airfoil theory.
///
/// τ = 1 - (θ_f - sin(θ_f)) / π, where θ_f = arccos(1 - 2·cf/c).
/// cf/c is the flap chord fraction (0 to 1).
#[must_use]
#[inline]
pub fn flap_effectiveness(chord_fraction: f64) -> f64 {
    let cf = chord_fraction.clamp(0.0, 1.0);
    if cf < f64::EPSILON {
        return 0.0;
    }
    if (cf - 1.0).abs() < f64::EPSILON {
        return 1.0;
    }
    let theta_f = (2.0 * cf - 1.0).acos();
    1.0 - (theta_f - theta_f.sin()) / std::f64::consts::PI
}

/// Elevator pitch moment effectiveness dCm/dδe.
///
/// dCm/dδe = -η_t · V_H · CL_α_t · τ(cf/c)
///
/// Negative value means trailing-edge-down elevator produces nose-down moment (correct).
#[must_use]
#[inline]
pub fn elevator_effectiveness(
    cl_alpha_tail: f64,
    elevator_chord_fraction: f64,
    tail_volume: f64,
    tail_efficiency: f64,
) -> f64 {
    let tau = flap_effectiveness(elevator_chord_fraction);
    -tail_efficiency * tail_volume * cl_alpha_tail * tau
}

/// Solve for trim (α and δe) given required CL.
///
/// Linear system:
/// - CL = CL_α · α → α = CL_req / CL_α
/// - Cm = Cm_0 + Cm_α · α + Cm_δe · δe = 0 → δe = -(Cm_0 + Cm_α · α) / Cm_δe
pub fn trim(
    cl_required: f64,
    cm_0: f64,
    cm_alpha_val: f64,
    cm_delta_e: f64,
    cl_alpha_val: f64,
) -> Result<TrimSolution> {
    if cl_alpha_val.abs() < f64::EPSILON {
        return Err(PavanError::ComputationError(
            "zero lift curve slope: cannot trim".into(),
        ));
    }
    if cm_delta_e.abs() < f64::EPSILON {
        return Err(PavanError::ComputationError(
            "zero elevator effectiveness: cannot trim".into(),
        ));
    }

    let alpha = cl_required / cl_alpha_val;
    let delta_e = -(cm_0 + cm_alpha_val * alpha) / cm_delta_e;
    let cm_residual = cm_0 + cm_alpha_val * alpha + cm_delta_e * delta_e;

    Ok(TrimSolution {
        alpha_rad: alpha,
        elevator_rad: delta_e,
        cl: cl_required,
        cd: 0.0,
        cm_residual,
    })
}

/// Solve for trim at a given flight speed, altitude, and weight.
///
/// CL_required = W / (q · S), then solves the trim system.
#[allow(clippy::too_many_arguments)]
pub fn trim_at_speed(
    weight_n: f64,
    velocity: f64,
    altitude: f64,
    wing_area: f64,
    cm_0: f64,
    cm_alpha_val: f64,
    cm_delta_e: f64,
    cl_alpha_val: f64,
) -> Result<TrimSolution> {
    if velocity <= 0.0 {
        return Err(PavanError::InvalidVelocity(
            "velocity must be positive".into(),
        ));
    }
    if wing_area <= 0.0 {
        return Err(PavanError::InvalidGeometry(
            "wing area must be positive".into(),
        ));
    }

    let rho = atmosphere::standard_density(altitude);
    let q = atmosphere::dynamic_pressure(rho, velocity);
    if q <= 0.0 {
        return Err(PavanError::ComputationError("zero dynamic pressure".into()));
    }

    let cl_required = weight_n / (q * wing_area);
    trim(cl_required, cm_0, cm_alpha_val, cm_delta_e, cl_alpha_val)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Neutral point & static margin ──

    #[test]
    fn neutral_point_basic() {
        // Wing AC at 0.25, Cm_α = -0.5, CL_α = 5.0
        // x_np = 0.25 - (-0.5)/5.0 = 0.25 + 0.1 = 0.35
        let xnp = neutral_point(5.0, -0.5, 0.25);
        assert!((xnp - 0.35).abs() < 0.001, "x_np should be 0.35, got {xnp}");
    }

    #[test]
    fn neutral_point_zero_cl_alpha() {
        let xnp = neutral_point(0.0, -0.5, 0.25);
        assert!((xnp - 0.25).abs() < f64::EPSILON);
    }

    #[test]
    fn static_margin_stable() {
        let sm = static_margin(0.35, 0.25);
        assert!(sm > 0.0, "CG ahead of NP should be stable");
        assert!((sm - 0.10).abs() < 0.001);
    }

    #[test]
    fn static_margin_unstable() {
        let sm = static_margin(0.25, 0.35);
        assert!(sm < 0.0, "CG behind NP should be unstable");
    }

    #[test]
    fn static_margin_neutral() {
        let sm = static_margin(0.25, 0.25);
        assert!(sm.abs() < f64::EPSILON);
    }

    // ── Cm_alpha ──

    #[test]
    fn cm_alpha_with_tail() {
        // Wing: CL_α=5.0, x_ac=0.25, x_cg=0.30
        // Tail: V_H=0.7, CL_α_t=4.0, η=0.9, dε/dα=0.4
        let cma = cm_alpha(5.0, 0.25, 0.30, 0.7, 4.0, 0.9, 0.4);
        // Wing: 5.0*(0.30-0.25) = 0.25 (destabilizing, CG behind AC)
        // Tail: -0.9*0.7*4.0*(1-0.4) = -1.512 (stabilizing)
        // Total: 0.25 - 1.512 = -1.262
        assert!(
            (cma - (-1.262)).abs() < 0.01,
            "Cm_α should be ~-1.262, got {cma}"
        );
        assert!(cma < 0.0, "with tail, Cm_α should be negative (stable)");
    }

    #[test]
    fn cm_alpha_wing_only() {
        // No tail (V_H = 0): Cm_α = CL_α * (x_cg - x_ac)
        let cma = cm_alpha(5.0, 0.25, 0.30, 0.0, 0.0, 0.0, 0.0);
        assert!(cma > 0.0, "wing only with CG behind AC should be unstable");
    }

    // ── Flap effectiveness ──

    #[test]
    fn flap_effectiveness_25_percent() {
        // Thin airfoil theory: θ_f = arccos(-0.5) = 2π/3, τ ≈ 0.609
        let tau = flap_effectiveness(0.25);
        assert!(
            (tau - 0.609).abs() < 0.02,
            "τ at cf/c=0.25 should be ~0.609, got {tau}"
        );
    }

    #[test]
    fn flap_effectiveness_50_percent() {
        // θ_f = arccos(0) = π/2, τ = 1 - (π/2 - 1)/π ≈ 0.818
        let tau = flap_effectiveness(0.50);
        assert!(
            (tau - 0.818).abs() < 0.02,
            "τ at cf/c=0.50 should be ~0.818, got {tau}"
        );
    }

    #[test]
    fn flap_effectiveness_zero() {
        assert_eq!(flap_effectiveness(0.0), 0.0);
    }

    #[test]
    fn flap_effectiveness_full() {
        assert!((flap_effectiveness(1.0) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn flap_effectiveness_monotonic() {
        let t1 = flap_effectiveness(0.1);
        let t2 = flap_effectiveness(0.3);
        let t3 = flap_effectiveness(0.5);
        assert!(t2 > t1);
        assert!(t3 > t2);
    }

    // ── Elevator effectiveness ──

    #[test]
    fn elevator_effectiveness_basic() {
        let cm_de = elevator_effectiveness(4.0, 0.25, 0.7, 0.9);
        assert!(cm_de < 0.0, "elevator should produce negative dCm/dδe");
    }

    // ── Trim ──

    #[test]
    fn trim_basic() {
        let sol = trim(0.5, 0.05, -1.0, -1.5, 5.0).expect("trim");
        assert!(
            (sol.cl - 0.5).abs() < f64::EPSILON,
            "trimmed CL should match required"
        );
        assert!(
            sol.cm_residual.abs() < 1e-10,
            "moment should be zero at trim"
        );
        // α = 0.5/5.0 = 0.1 rad
        assert!((sol.alpha_rad - 0.1).abs() < 1e-10);
    }

    #[test]
    fn trim_zero_cl_alpha_errors() {
        assert!(trim(0.5, 0.0, -1.0, -1.5, 0.0).is_err());
    }

    #[test]
    fn trim_zero_elevator_errors() {
        assert!(trim(0.5, 0.0, -1.0, 0.0, 5.0).is_err());
    }

    #[test]
    fn trim_at_speed_cessna_class() {
        // W=10000N, V=60m/s, sea level, S=16.2m²
        let sol = trim_at_speed(10000.0, 60.0, 0.0, 16.2, 0.05, -1.0, -1.5, 5.0).expect("trim");
        assert!(
            sol.alpha_rad > 0.0 && sol.alpha_rad < 0.3,
            "α should be reasonable"
        );
        assert!(sol.cm_residual.abs() < 1e-10);
    }

    #[test]
    fn trim_at_speed_zero_velocity_errors() {
        assert!(trim_at_speed(10000.0, 0.0, 0.0, 16.2, 0.0, -1.0, -1.5, 5.0).is_err());
    }

    #[test]
    fn trim_at_speed_zero_area_errors() {
        assert!(trim_at_speed(10000.0, 60.0, 0.0, 0.0, 0.0, -1.0, -1.5, 5.0).is_err());
    }

    // ── Serde ──

    #[test]
    fn serde_trim_solution() {
        let sol = trim(0.5, 0.05, -1.0, -1.5, 5.0).expect("trim");
        let json = serde_json::to_string(&sol).expect("serialize");
        let back: TrimSolution = serde_json::from_str(&json).expect("deserialize");
        assert!((back.alpha_rad - sol.alpha_rad).abs() < f64::EPSILON);
    }

    #[test]
    fn serde_longitudinal_stability() {
        let ls = LongitudinalStability {
            neutral_point: 0.35,
            static_margin: 0.10,
            cm_alpha: -1.0,
            cl_alpha: 5.0,
        };
        let json = serde_json::to_string(&ls).expect("serialize");
        let back: LongitudinalStability = serde_json::from_str(&json).expect("deserialize");
        assert!((back.neutral_point - ls.neutral_point).abs() < f64::EPSILON);
    }

    #[test]
    fn serde_control_surface() {
        let cs = ControlSurface {
            chord_fraction: 0.25,
            span_start: 0.3,
            span_end: 0.9,
            max_deflection_rad: 25.0_f64.to_radians(),
        };
        let json = serde_json::to_string(&cs).expect("serialize");
        let back: ControlSurface = serde_json::from_str(&json).expect("deserialize");
        assert!((back.chord_fraction - cs.chord_fraction).abs() < f64::EPSILON);
    }
}
