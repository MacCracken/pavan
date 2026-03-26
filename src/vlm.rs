//! Vortex Lattice Method (VLM) for 3D finite wings.
//!
//! Models a wing as a lattice of horseshoe vortices, solving for circulation
//! distribution via Biot-Savart influence coefficients. Produces span loading,
//! lift coefficient, induced drag, and Oswald efficiency factor.

use std::f64::consts::PI;

use hisab::DVec3;
use serde::{Deserialize, Serialize};

use crate::error::{PavanError, Result};

/// Vortex core cutoff to prevent singularity in Biot-Savart.
const VORTEX_CORE_EPS: f64 = 1e-8;

/// Far-field truncation length for semi-infinite trailing legs (chord-lengths).
const TRAILING_LEG_FACTOR: f64 = 50.0;

/// Minimum number of panels for a meaningful solve.
const MIN_VLM_PANELS: usize = 2;

/// Wing planform definition for VLM analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct WingGeometry {
    /// Full wingspan (m), tip to tip.
    pub span: f64,
    /// Root chord length (m).
    pub root_chord: f64,
    /// Tip chord length (m). Equal to root_chord for rectangular wings.
    pub tip_chord: f64,
    /// Leading edge sweep angle (rad). 0 for unswept.
    pub sweep_le_rad: f64,
    /// Dihedral angle (rad). 0 for flat wing.
    pub dihedral_rad: f64,
    /// Geometric twist at tip relative to root (rad). Negative = washout.
    pub twist_tip_rad: f64,
    /// Spanwise panels per semi-span.
    pub n_span: usize,
    /// Chordwise panels.
    pub n_chord: usize,
}

impl WingGeometry {
    /// Rectangular wing (no taper, no sweep, no twist).
    #[must_use]
    #[inline]
    pub fn rectangular(span: f64, chord: f64, n_span: usize, n_chord: usize) -> Self {
        Self {
            span,
            root_chord: chord,
            tip_chord: chord,
            sweep_le_rad: 0.0,
            dihedral_rad: 0.0,
            twist_tip_rad: 0.0,
            n_span,
            n_chord,
        }
    }

    /// Tapered wing (no sweep, no twist).
    #[must_use]
    #[inline]
    pub fn tapered(
        span: f64,
        root_chord: f64,
        tip_chord: f64,
        n_span: usize,
        n_chord: usize,
    ) -> Self {
        Self {
            span,
            root_chord,
            tip_chord,
            sweep_le_rad: 0.0,
            dihedral_rad: 0.0,
            twist_tip_rad: 0.0,
            n_span,
            n_chord,
        }
    }

    /// Wing reference area: S = span × (root + tip) / 2.
    #[must_use]
    #[inline]
    pub fn reference_area(&self) -> f64 {
        self.span * (self.root_chord + self.tip_chord) * 0.5
    }

    /// Aspect ratio: AR = span² / S.
    #[must_use]
    #[inline]
    pub fn aspect_ratio(&self) -> f64 {
        let s = self.reference_area();
        if s <= 0.0 {
            return 0.0;
        }
        self.span * self.span / s
    }

    /// Total number of panels (both semi-spans × chordwise).
    #[must_use]
    #[inline]
    pub fn total_panels(&self) -> usize {
        2 * self.n_span * self.n_chord
    }
}

/// A single VLM panel with horseshoe vortex.
#[derive(Debug, Clone)]
pub struct VlmPanel {
    /// Panel corners: [LE-left, LE-right, TE-right, TE-left].
    pub corners: [DVec3; 4],
    /// Bound vortex endpoints (at 1/4 chord): [left, right].
    pub bound: [DVec3; 2],
    /// Control point at 3/4 chord, mid-span.
    pub control_point: DVec3,
    /// Outward unit normal.
    pub normal: DVec3,
    /// Panel area (m²).
    pub area: f64,
    /// Local chord length (m).
    pub chord: f64,
    /// Spanwise width (m).
    pub dy: f64,
}

/// Result of a VLM solve.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct VlmSolution {
    /// Total lift coefficient.
    pub cl: f64,
    /// Induced drag coefficient.
    pub cdi: f64,
    /// Pitching moment coefficient about root leading edge.
    pub cm: f64,
    /// Sectional lift coefficient per spanwise strip.
    pub span_cl: Vec<f64>,
    /// Circulation strength per panel (m²/s).
    pub circulation: Vec<f64>,
    /// Oswald span efficiency factor (1.0 = elliptic loading).
    pub oswald_efficiency: f64,
}

/// Generate VLM panels from wing planform geometry.
///
/// Creates a lattice of panels across both semi-spans, with bound vortices
/// at 1/4 chord and control points at 3/4 chord (Pistolesi's theorem).
#[must_use]
pub fn generate_panels(wing: &WingGeometry) -> Vec<VlmPanel> {
    let n_span_total = 2 * wing.n_span;
    let n_chord = wing.n_chord;
    let n_total = n_span_total * n_chord;
    let mut panels = Vec::with_capacity(n_total);

    let half_span = wing.span * 0.5;

    for i_span in 0..n_span_total {
        // Spanwise position: -half_span to +half_span
        let eta_left = (i_span as f64 / n_span_total as f64) * 2.0 - 1.0;
        let eta_right = ((i_span + 1) as f64 / n_span_total as f64) * 2.0 - 1.0;
        let y_left = eta_left * half_span;
        let y_right = eta_right * half_span;

        // Interpolate chord and twist at each spanwise station
        let frac_left = eta_left.abs();
        let frac_right = eta_right.abs();
        let chord_left = wing.root_chord + (wing.tip_chord - wing.root_chord) * frac_left;
        let chord_right = wing.root_chord + (wing.tip_chord - wing.root_chord) * frac_right;

        // Leading edge x-offset from sweep
        let x_le_left = frac_left * half_span * wing.sweep_le_rad.tan();
        let x_le_right = frac_right * half_span * wing.sweep_le_rad.tan();

        // Dihedral z-offset
        let z_left = y_left.abs() * wing.dihedral_rad.sin();
        let z_right = y_right.abs() * wing.dihedral_rad.sin();

        // Twist (linear from root to tip)
        let twist_left = wing.twist_tip_rad * frac_left;
        let twist_right = wing.twist_tip_rad * frac_right;

        for i_chord in 0..n_chord {
            let xi_le = i_chord as f64 / n_chord as f64;
            let xi_te = (i_chord + 1) as f64 / n_chord as f64;

            // Panel corners (x downstream, y spanwise, z up)
            let le_left = DVec3::new(
                x_le_left + xi_le * chord_left,
                y_left,
                z_left - xi_le * chord_left * twist_left.sin(),
            );
            let le_right = DVec3::new(
                x_le_right + xi_le * chord_right,
                y_right,
                z_right - xi_le * chord_right * twist_right.sin(),
            );
            let te_right = DVec3::new(
                x_le_right + xi_te * chord_right,
                y_right,
                z_right - xi_te * chord_right * twist_right.sin(),
            );
            let te_left = DVec3::new(
                x_le_left + xi_te * chord_left,
                y_left,
                z_left - xi_te * chord_left * twist_left.sin(),
            );

            // Bound vortex at 1/4 chord of this panel
            let quarter = 0.25;
            let bound_frac = xi_le + quarter * (xi_te - xi_le);
            let bound_left = DVec3::new(
                x_le_left + bound_frac * chord_left,
                y_left,
                z_left - bound_frac * chord_left * twist_left.sin(),
            );
            let bound_right = DVec3::new(
                x_le_right + bound_frac * chord_right,
                y_right,
                z_right - bound_frac * chord_right * twist_right.sin(),
            );

            // Control point at 3/4 chord, mid-span
            let three_quarter = 0.75;
            let cp_frac = xi_le + three_quarter * (xi_te - xi_le);
            let y_mid = 0.5 * (y_left + y_right);
            let frac_mid = y_mid.abs() / half_span;
            let chord_mid = wing.root_chord + (wing.tip_chord - wing.root_chord) * frac_mid;
            let x_le_mid = frac_mid * half_span * wing.sweep_le_rad.tan();
            let z_mid = y_mid.abs() * wing.dihedral_rad.sin();
            let twist_mid = wing.twist_tip_rad * frac_mid;
            let cp = DVec3::new(
                x_le_mid + cp_frac * chord_mid,
                y_mid,
                z_mid - cp_frac * chord_mid * twist_mid.sin(),
            );

            // Normal from cross product of diagonals (order gives +Z for flat wing)
            let diag1 = te_right - le_left;
            let diag2 = te_left - le_right;
            let normal = diag2.cross(diag1);
            let normal_len = normal.length();
            let normal = if normal_len > f64::EPSILON {
                normal / normal_len
            } else {
                tracing::warn!(
                    "degenerate VLM panel at span index {i_span}, chord index {i_chord}: zero normal"
                );
                DVec3::Z
            };

            // Area = 0.5 * |diag1 × diag2|
            let area = 0.5 * normal_len;
            let local_chord = 0.5 * (chord_left + chord_right) * (xi_te - xi_le);
            let dy = (y_right - y_left).abs();

            panels.push(VlmPanel {
                corners: [le_left, le_right, te_right, te_left],
                bound: [bound_left, bound_right],
                control_point: cp,
                normal,
                area,
                chord: local_chord,
                dy,
            });
        }
    }

    panels
}

/// Velocity induced by a finite vortex segment a→b at point p (unit circulation).
///
/// Uses the Biot-Savart law with vortex core cutoff for numerical stability.
#[must_use]
#[inline]
pub fn biot_savart_segment(p: DVec3, a: DVec3, b: DVec3) -> DVec3 {
    let r1 = p - a;
    let r2 = p - b;
    let r0 = b - a;

    let cross = r1.cross(r2);
    let cross_sq = cross.length_squared();

    // Core cutoff: avoid singularity when p is on the vortex line
    if cross_sq < VORTEX_CORE_EPS {
        return DVec3::ZERO;
    }

    let r1_len = r1.length();
    let r2_len = r2.length();

    if r1_len < VORTEX_CORE_EPS || r2_len < VORTEX_CORE_EPS {
        return DVec3::ZERO;
    }

    // Biot-Savart: V = Γ/(4π) · (r1×r2)/|r1×r2|² · r0·(r1/|r1| - r2/|r2|)
    let dot_factor = r0.dot(r1 / r1_len - r2 / r2_len);

    (1.0 / (4.0 * PI)) * cross * (dot_factor / cross_sq)
}

/// Induced velocity from a horseshoe vortex (unit circulation) at a control point.
///
/// Horseshoe = bound segment + two trailing semi-infinite legs extending downstream (+x).
#[inline]
fn horseshoe_influence(panel: &VlmPanel, control: DVec3, far_field: f64) -> DVec3 {
    let a = panel.bound[0]; // left end of bound vortex
    let b = panel.bound[1]; // right end of bound vortex

    // Trailing leg direction: downstream (+x)
    let trail = DVec3::new(far_field, 0.0, 0.0);

    // Left trailing leg: from far upstream-left to bound-left (semi-infinite, approximated)
    let left_far = a + trail;
    let v_left = biot_savart_segment(control, left_far, a);

    // Bound segment: from left to right
    let v_bound = biot_savart_segment(control, a, b);

    // Right trailing leg: from bound-right to far downstream-right
    let right_far = b + trail;
    let v_right = biot_savart_segment(control, b, right_far);

    v_left + v_bound + v_right
}

/// Post-process VLM solution: compute CL, CDi, Cm, span loading, Oswald efficiency.
#[allow(clippy::needless_range_loop)]
fn vlm_post_process(
    panels: &[VlmPanel],
    gamma: Vec<f64>,
    wing: &WingGeometry,
    v_inf: f64,
) -> VlmSolution {
    let s_ref = wing.reference_area();
    let q_inf = 0.5 * v_inf * v_inf;
    let ar = wing.aspect_ratio();
    let half_span = wing.span * 0.5;
    let n_span_total = 2 * wing.n_span;
    let n_chord = wing.n_chord;

    // Aggregate circulation per spanwise strip
    let mut strip_gamma = vec![0.0; n_span_total];
    for i_span in 0..n_span_total {
        for i_chord in 0..n_chord {
            strip_gamma[i_span] += gamma[i_span * n_chord + i_chord];
        }
    }

    // Lift and moment via Kutta-Joukowski
    let mut cl_total = 0.0;
    let mut cm_total = 0.0;
    let mut span_cl = Vec::with_capacity(n_span_total);

    for i_span in 0..n_span_total {
        let idx0 = i_span * n_chord;
        let dy = panels[idx0].dy;
        let local_chord: f64 = (0..n_chord).map(|ic| panels[idx0 + ic].chord).sum();

        let cl_section = if local_chord > 0.0 && q_inf > 0.0 {
            strip_gamma[i_span] * v_inf / (q_inf * local_chord)
        } else {
            0.0
        };
        span_cl.push(cl_section);

        let dcl = strip_gamma[i_span] * dy * v_inf / (q_inf * s_ref);
        cl_total += dcl;

        let x_cp = panels[idx0].control_point.x;
        let mac = 0.5 * (wing.root_chord + wing.tip_chord);
        cm_total -= x_cp * dcl / mac;
    }

    // Induced drag via Trefftz plane trailing vortex sheet
    let mut trail_gamma = Vec::with_capacity(n_span_total + 1);
    trail_gamma.push(strip_gamma[0]);
    for j in 0..n_span_total - 1 {
        trail_gamma.push(strip_gamma[j + 1] - strip_gamma[j]);
    }
    trail_gamma.push(-strip_gamma[n_span_total - 1]);

    let mut trail_y = Vec::with_capacity(n_span_total + 1);
    for j in 0..=n_span_total {
        let eta = (j as f64 / n_span_total as f64) * 2.0 - 1.0;
        trail_y.push(eta * half_span);
    }

    let mut cdi_total = 0.0;
    for i_span in 0..n_span_total {
        let idx0 = i_span * n_chord;
        let y_i = panels[idx0].control_point.y;
        let dy = panels[idx0].dy;

        let mut w_i = 0.0;
        for j in 0..trail_gamma.len() {
            let dist = y_i - trail_y[j];
            if dist.abs() > 1e-10 {
                w_i -= trail_gamma[j] / (4.0 * PI * dist);
            }
        }
        cdi_total -= strip_gamma[i_span] * w_i * dy / (q_inf * s_ref);
    }

    let oswald = if cdi_total.abs() > f64::EPSILON && ar > 0.0 {
        cl_total * cl_total / (PI * ar * cdi_total)
    } else {
        1.0
    };

    VlmSolution {
        cl: cl_total,
        cdi: cdi_total,
        cm: cm_total,
        span_cl,
        circulation: gamma,
        oswald_efficiency: oswald,
    }
}

/// Solve the VLM for a single angle of attack.
///
/// Panels should be generated from [`generate_panels`]. Freestream velocity
/// `v_inf` in m/s; angle of attack in radians.
#[allow(clippy::needless_range_loop)]
pub fn solve(
    panels: &[VlmPanel],
    wing: &WingGeometry,
    alpha_rad: f64,
    v_inf: f64,
) -> Result<VlmSolution> {
    let n = panels.len();
    if n < MIN_VLM_PANELS {
        return Err(PavanError::InvalidGeometry(format!(
            "need at least {MIN_VLM_PANELS} VLM panels, got {n}"
        )));
    }
    if v_inf <= 0.0 {
        return Err(PavanError::InvalidVelocity(
            "freestream velocity must be positive".into(),
        ));
    }

    let far_field = wing.root_chord * TRAILING_LEG_FACTOR;

    // Freestream velocity vector
    let v_free = DVec3::new(v_inf * alpha_rad.cos(), 0.0, v_inf * alpha_rad.sin());

    // Build influence matrix: A[i][j] = V_induced_by_horseshoe_j at control_i · normal_i
    let mut matrix: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    let mut rhs: Vec<f64> = vec![0.0; n];

    for i in 0..n {
        for j in 0..n {
            let v_ind = horseshoe_influence(&panels[j], panels[i].control_point, far_field);
            matrix[i][j] = v_ind.dot(panels[i].normal);
        }
        rhs[i] = -v_free.dot(panels[i].normal);
    }

    // Solve for circulation strengths
    let (lu, pivot) = hisab::num::lu_decompose(&matrix)
        .map_err(|e| PavanError::ComputationError(format!("VLM LU decomposition failed: {e}")))?;
    let gamma = hisab::num::lu_solve(&lu, &pivot, &rhs)
        .map_err(|e| PavanError::ComputationError(format!("VLM LU solve failed: {e}")))?;

    Ok(vlm_post_process(panels, gamma, wing, v_inf))
}

/// Solve VLM for multiple angles of attack efficiently.
///
/// The influence matrix is LU-decomposed once and reused for each angle.
#[allow(clippy::needless_range_loop)]
pub fn solve_multi(
    panels: &[VlmPanel],
    wing: &WingGeometry,
    alphas: &[f64],
    v_inf: f64,
) -> Result<Vec<VlmSolution>> {
    let n = panels.len();
    if n < MIN_VLM_PANELS {
        return Err(PavanError::InvalidGeometry(format!(
            "need at least {MIN_VLM_PANELS} VLM panels, got {n}"
        )));
    }
    if v_inf <= 0.0 {
        return Err(PavanError::InvalidVelocity(
            "freestream velocity must be positive".into(),
        ));
    }

    let far_field = wing.root_chord * TRAILING_LEG_FACTOR;

    // Build influence matrix (angle-independent)
    let mut matrix: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let v_ind = horseshoe_influence(&panels[j], panels[i].control_point, far_field);
            matrix[i][j] = v_ind.dot(panels[i].normal);
        }
    }

    let (lu, pivot) = hisab::num::lu_decompose(&matrix)
        .map_err(|e| PavanError::ComputationError(format!("VLM LU decomposition failed: {e}")))?;

    let mut results = Vec::with_capacity(alphas.len());

    for &alpha in alphas {
        let v_free = DVec3::new(v_inf * alpha.cos(), 0.0, v_inf * alpha.sin());

        let rhs: Vec<f64> = panels.iter().map(|p| -v_free.dot(p.normal)).collect();

        let gamma = hisab::num::lu_solve(&lu, &pivot, &rhs)
            .map_err(|e| PavanError::ComputationError(format!("VLM LU solve failed: {e}")))?;

        results.push(vlm_post_process(panels, gamma, wing, v_inf));
    }

    Ok(results)
}

/// Compute Oswald span efficiency from VLM results.
///
/// e = CL² / (π · AR · CDi). Returns 1.0 for elliptic loading.
#[must_use]
#[inline]
pub fn oswald_efficiency(cl: f64, cdi: f64, aspect_ratio: f64) -> f64 {
    if cdi.abs() < f64::EPSILON || aspect_ratio <= 0.0 {
        return 1.0;
    }
    cl * cl / (PI * aspect_ratio * cdi)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rect_wing(n_span: usize, n_chord: usize) -> (WingGeometry, Vec<VlmPanel>) {
        // AR=6 rectangular wing: span=6m, chord=1m
        let wing = WingGeometry::rectangular(6.0, 1.0, n_span, n_chord);
        let panels = generate_panels(&wing);
        (wing, panels)
    }

    // --- Panel generation ---

    #[test]
    fn panel_count_correct() {
        let wing = WingGeometry::rectangular(6.0, 1.0, 10, 4);
        let panels = generate_panels(&wing);
        assert_eq!(panels.len(), 2 * 10 * 4);
    }

    #[test]
    fn panels_have_positive_area() {
        let (_, panels) = rect_wing(10, 2);
        for (i, p) in panels.iter().enumerate() {
            assert!(p.area > 0.0, "panel {i} has non-positive area {}", p.area);
        }
    }

    #[test]
    fn panel_normals_are_unit() {
        let (_, panels) = rect_wing(10, 2);
        for (i, p) in panels.iter().enumerate() {
            let mag = p.normal.length();
            assert!(
                (mag - 1.0).abs() < 1e-10,
                "panel {i} normal not unit: {mag}"
            );
        }
    }

    #[test]
    fn panel_normals_point_upward() {
        // For a flat wing in the XY plane, normals should point +Z
        let (_, panels) = rect_wing(10, 2);
        for (i, p) in panels.iter().enumerate() {
            assert!(
                p.normal.z > 0.5,
                "panel {i} normal should point up, got z={}",
                p.normal.z
            );
        }
    }

    #[test]
    fn control_points_within_wing() {
        let (wing, panels) = rect_wing(10, 2);
        let half_span = wing.span * 0.5;
        for p in &panels {
            assert!(p.control_point.y.abs() <= half_span + 0.01);
            assert!(p.control_point.x >= -0.01);
            assert!(p.control_point.x <= wing.root_chord + 0.01);
        }
    }

    // --- Biot-Savart ---

    #[test]
    fn biot_savart_known_geometry() {
        // Point at (0.5,1,0) from segment along x-axis (0,0,0)→(1,0,0)
        // Right-hand rule: dl(+x) × r(+y) = +z
        let v = biot_savart_segment(
            DVec3::new(0.5, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
        );
        assert!(
            v.z > 0.0,
            "Biot-Savart should induce +z velocity, got {v:?}"
        );
        assert!(v.x.abs() < 1e-10, "x component should be ~0");
        assert!(v.y.abs() < 1e-10, "y component should be ~0");
    }

    #[test]
    fn biot_savart_point_on_line_returns_zero() {
        // Point on the vortex line itself → core cutoff returns zero
        let v = biot_savart_segment(
            DVec3::new(0.5, 0.0, 0.0),
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
        );
        assert_eq!(v, DVec3::ZERO);
    }

    // --- Solve correctness ---

    #[test]
    fn rectangular_wing_ar6_at_5deg() {
        let (wing, panels) = rect_wing(10, 2);
        let sol = solve(&panels, &wing, 5.0_f64.to_radians(), 1.0).expect("solve");
        // Finite wing CL < 2D Cl (0.55). For AR=6, CL ≈ 0.35-0.50
        assert!(
            sol.cl > 0.2 && sol.cl < 0.6,
            "AR=6 wing at 5° should have CL in [0.2, 0.6], got {}",
            sol.cl
        );
    }

    #[test]
    fn positive_induced_drag() {
        let (wing, panels) = rect_wing(10, 2);
        let sol = solve(&panels, &wing, 5.0_f64.to_radians(), 1.0).expect("solve");
        assert!(
            sol.cdi > 0.0,
            "induced drag should be positive, got {}",
            sol.cdi
        );
    }

    #[test]
    fn symmetric_wing_zero_aoa_zero_lift() {
        let (wing, panels) = rect_wing(8, 2);
        let sol = solve(&panels, &wing, 0.0, 1.0).expect("solve");
        assert!(
            sol.cl.abs() < 0.01,
            "symmetric wing at 0° should have CL ≈ 0, got {}",
            sol.cl
        );
    }

    #[test]
    fn lift_increases_with_aoa() {
        let (wing, panels) = rect_wing(8, 2);
        let sol2 = solve(&panels, &wing, 2.0_f64.to_radians(), 1.0).expect("solve");
        let sol5 = solve(&panels, &wing, 5.0_f64.to_radians(), 1.0).expect("solve");
        let sol8 = solve(&panels, &wing, 8.0_f64.to_radians(), 1.0).expect("solve");
        assert!(sol5.cl > sol2.cl);
        assert!(sol8.cl > sol5.cl);
    }

    #[test]
    fn oswald_efficiency_reasonable() {
        let (wing, panels) = rect_wing(10, 2);
        let sol = solve(&panels, &wing, 5.0_f64.to_radians(), 1.0).expect("solve");
        // Rectangular wing: e ≈ 0.7-1.0
        assert!(
            sol.oswald_efficiency > 0.5 && sol.oswald_efficiency < 1.5,
            "Oswald e should be reasonable, got {}",
            sol.oswald_efficiency
        );
    }

    #[test]
    fn span_loading_has_correct_length() {
        let (wing, panels) = rect_wing(10, 2);
        let sol = solve(&panels, &wing, 5.0_f64.to_radians(), 1.0).expect("solve");
        assert_eq!(sol.span_cl.len(), 2 * 10);
    }

    #[test]
    fn tapered_wing_positive_lift() {
        let wing = WingGeometry::tapered(6.0, 1.5, 0.5, 10, 2);
        let panels = generate_panels(&wing);
        let sol = solve(&panels, &wing, 5.0_f64.to_radians(), 1.0).expect("solve");
        assert!(sol.cl > 0.1, "tapered wing should produce lift");
    }

    #[test]
    fn convergence_with_panel_count() {
        let wing_lo = WingGeometry::rectangular(6.0, 1.0, 5, 1);
        let wing_hi = WingGeometry::rectangular(6.0, 1.0, 20, 2);
        let panels_lo = generate_panels(&wing_lo);
        let panels_hi = generate_panels(&wing_hi);

        let sol_lo = solve(&panels_lo, &wing_lo, 5.0_f64.to_radians(), 1.0).expect("solve");
        let sol_hi = solve(&panels_hi, &wing_hi, 5.0_f64.to_radians(), 1.0).expect("solve");

        // Both should give roughly similar CL
        let diff = (sol_lo.cl - sol_hi.cl).abs();
        assert!(
            diff < 0.15,
            "CL should converge: lo={}, hi={}, diff={}",
            sol_lo.cl,
            sol_hi.cl,
            diff
        );
    }

    // --- solve_multi ---

    #[test]
    fn solve_multi_matches_individual() {
        let (wing, panels) = rect_wing(6, 2);
        let alphas = [0.0, 3.0_f64.to_radians(), 5.0_f64.to_radians()];
        let multi = solve_multi(&panels, &wing, &alphas, 1.0).expect("solve_multi");

        for (i, &alpha) in alphas.iter().enumerate() {
            let single = solve(&panels, &wing, alpha, 1.0).expect("solve");
            assert!(
                (multi[i].cl - single.cl).abs() < 1e-10,
                "CL mismatch at alpha[{i}]"
            );
            assert!(
                (multi[i].cdi - single.cdi).abs() < 1e-10,
                "CDi mismatch at alpha[{i}]"
            );
            assert!(
                (multi[i].cm - single.cm).abs() < 1e-10,
                "Cm mismatch at alpha[{i}]"
            );
        }
    }

    #[test]
    fn solve_multi_empty_alphas() {
        let (wing, panels) = rect_wing(6, 2);
        let results = solve_multi(&panels, &wing, &[], 1.0).expect("solve_multi");
        assert!(results.is_empty());
    }

    // --- Error cases ---

    #[test]
    fn too_few_panels_errors() {
        let wing = WingGeometry::rectangular(6.0, 1.0, 1, 1);
        let panels = vec![]; // empty
        assert!(solve(&panels, &wing, 0.0, 1.0).is_err());
    }

    #[test]
    fn zero_velocity_errors() {
        let (wing, panels) = rect_wing(6, 2);
        assert!(solve(&panels, &wing, 0.0, 0.0).is_err());
    }

    // --- Geometry presets ---

    #[test]
    fn rectangular_reference_area() {
        let wing = WingGeometry::rectangular(6.0, 1.0, 10, 2);
        assert!((wing.reference_area() - 6.0).abs() < f64::EPSILON);
    }

    #[test]
    fn rectangular_aspect_ratio() {
        let wing = WingGeometry::rectangular(6.0, 1.0, 10, 2);
        assert!((wing.aspect_ratio() - 6.0).abs() < f64::EPSILON);
    }

    #[test]
    fn tapered_reference_area() {
        let wing = WingGeometry::tapered(10.0, 2.0, 1.0, 10, 2);
        // S = 10 × (2 + 1) / 2 = 15
        assert!((wing.reference_area() - 15.0).abs() < f64::EPSILON);
    }

    #[test]
    fn oswald_efficiency_fn() {
        // CL=0.5, CDi=0.005, AR=8 → e = 0.25/(π·8·0.005) = 0.25/0.1257 ≈ 1.99
        // That's high — just testing the formula works
        let e = oswald_efficiency(0.5, 0.0133, 6.0);
        assert!(e > 0.5 && e < 1.5, "e should be reasonable, got {e}");
    }

    // --- Guard clause coverage ---

    #[test]
    fn total_panels_count() {
        let wing = WingGeometry::rectangular(6.0, 1.0, 10, 3);
        assert_eq!(wing.total_panels(), 60);
    }

    #[test]
    fn zero_span_area_and_ar() {
        let wing = WingGeometry::rectangular(0.0, 1.0, 5, 2);
        assert_eq!(wing.reference_area(), 0.0);
        assert_eq!(wing.aspect_ratio(), 0.0);
    }

    #[test]
    fn solve_zero_velocity_errors() {
        let (wing, panels) = rect_wing(6, 2);
        assert!(solve(&panels, &wing, 0.0, 0.0).is_err());
    }

    #[test]
    fn solve_multi_zero_velocity_errors() {
        let (wing, panels) = rect_wing(6, 2);
        assert!(solve_multi(&panels, &wing, &[0.0], 0.0).is_err());
    }

    #[test]
    fn solve_multi_too_few_panels_errors() {
        let wing = WingGeometry::rectangular(6.0, 1.0, 1, 1);
        assert!(solve_multi(&[], &wing, &[0.0], 1.0).is_err());
    }

    #[test]
    fn biot_savart_coincident_endpoints() {
        // a == b → zero-length segment
        let v = biot_savart_segment(
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(0.0, 0.0, 0.0),
        );
        assert_eq!(v, DVec3::ZERO);
    }

    #[test]
    fn biot_savart_point_near_endpoint() {
        let v = biot_savart_segment(
            DVec3::new(1e-15, 0.0, 0.0),
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
        );
        // Should return zero due to core cutoff
        assert_eq!(v, DVec3::ZERO);
    }

    #[test]
    fn wing_geometry_serde() {
        let wing = WingGeometry::rectangular(6.0, 1.0, 10, 2);
        let json = serde_json::to_string(&wing).expect("serialize");
        let back: WingGeometry = serde_json::from_str(&json).expect("deserialize");
        assert!((back.span - wing.span).abs() < f64::EPSILON);
    }

    #[test]
    fn vlm_solution_serde() {
        let (wing, panels) = rect_wing(6, 2);
        let sol = solve(&panels, &wing, 5.0_f64.to_radians(), 1.0).expect("solve");
        let json = serde_json::to_string(&sol).expect("serialize");
        let back: VlmSolution = serde_json::from_str(&json).expect("deserialize");
        assert!((back.cl - sol.cl).abs() < f64::EPSILON);
    }
}
