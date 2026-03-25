//! Hess-Smith 2D panel method for potential flow around airfoils.
//!
//! Solves for pressure distribution, lift, drag, and moment coefficients
//! using constant-strength source panels with a uniform vortex distribution
//! and the Kutta condition at the trailing edge.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::airfoil::SurfacePoints;
use crate::error::{PavanError, Result};

/// A single panel element between two boundary points.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct Panel {
    /// Start point (x, y).
    pub start: (f64, f64),
    /// End point (x, y).
    pub end: (f64, f64),
    /// Midpoint (control point).
    pub center: (f64, f64),
    /// Outward unit normal.
    pub normal: (f64, f64),
    /// Unit tangent (start → end direction).
    pub tangent: (f64, f64),
    /// Panel length.
    pub length: f64,
}

/// Result of a panel method solve.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct PanelSolution {
    /// Lift coefficient.
    pub cl: f64,
    /// Pressure drag coefficient.
    pub cd: f64,
    /// Pitch moment coefficient about quarter-chord.
    pub cm: f64,
    /// Pressure coefficient at each panel center.
    pub cp: Vec<f64>,
    /// Total vortex strength (circulation).
    pub gamma: f64,
}

/// Minimum panels required for a meaningful solve.
const MIN_PANELS: usize = 4;

/// Build panels from upper and lower surface coordinates.
///
/// Merges into a single closed counterclockwise loop:
/// trailing edge (upper) → upper surface (reversed) → leading edge → lower surface → trailing edge (lower).
///
/// CCW ordering ensures outward normals (ty, -tx) point away from the airfoil body.
#[must_use]
pub fn panels_from_surface(upper: &SurfacePoints, lower: &SurfacePoints) -> Vec<Panel> {
    if upper.len() < 2 || lower.len() < 2 {
        return Vec::new();
    }

    // Build closed CCW loop: upper reversed (TE→LE) then lower (LE→TE)
    let mut points: Vec<(f64, f64)> = Vec::with_capacity(upper.len() + lower.len());

    // Upper surface reversed: TE → LE
    for &pt in upper.iter().rev() {
        points.push(pt);
    }
    // Lower surface forward: LE → TE (skip first point, it duplicates the LE)
    for &pt in lower.iter().skip(1) {
        points.push(pt);
    }

    let n = points.len();
    let mut panels = Vec::with_capacity(n - 1);

    // Open trailing edge: no closing panel between last and first points
    for i in 0..n - 1 {
        let (x0, y0) = points[i];
        let (x1, y1) = points[i + 1];
        let dx = x1 - x0;
        let dy = y1 - y0;
        let len = (dx * dx + dy * dy).sqrt();

        if len < f64::EPSILON {
            continue;
        }

        let tx = dx / len;
        let ty = dy / len;
        // Outward normal: for clockwise ordering, normal points outward (right of tangent)
        let nx = ty;
        let ny = -tx;

        panels.push(Panel {
            start: (x0, y0),
            end: (x1, y1),
            center: (0.5 * (x0 + x1), 0.5 * (y0 + y1)),
            normal: (nx, ny),
            tangent: (tx, ty),
            length: len,
        });
    }

    panels
}

/// Compute the source and vortex influence coefficients for panel j at control point i.
///
/// Returns (A_n, A_t, B_n, B_t):
/// - A_n: normal velocity from unit source on panel j
/// - A_t: tangential velocity from unit source on panel j
/// - B_n: normal velocity from unit vortex on panel j
/// - B_t: tangential velocity from unit vortex on panel j
#[inline]
fn panel_influence(panels: &[Panel], i: usize, j: usize) -> (f64, f64, f64, f64) {
    if i == j {
        // Self-influence: source contributes ±0.5 to normal, vortex ±0.5 to tangential
        return (0.5, 0.0, 0.0, 0.5);
    }

    let pi = &panels[i];
    let pj = &panels[j];

    // Transform control point i into panel j's local coordinates
    let dx = pi.center.0 - pj.start.0;
    let dy = pi.center.1 - pj.start.1;

    // Local coordinates aligned with panel j
    let xi = dx * pj.tangent.0 + dy * pj.tangent.1;
    let eta = dx * pj.normal.0 + dy * pj.normal.1;
    let sj = pj.length;

    // Distances from panel j endpoints to control point i
    let r1_sq = xi * xi + eta * eta;
    let r2_sq = (xi - sj) * (xi - sj) + eta * eta;

    // Guard against degenerate geometry
    if r1_sq < f64::EPSILON || r2_sq < f64::EPSILON {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let ln_ratio = 0.5 * (r1_sq / r2_sq).ln();
    let theta = eta.atan2(xi - sj) - eta.atan2(xi);

    // Source influence in panel j's local frame
    let u_s = (1.0 / (2.0 * PI)) * ln_ratio; // tangential
    let v_s = (1.0 / (2.0 * PI)) * theta; // normal

    // Vortex influence in panel j's local frame (rotated 90° from source)
    let u_v = (1.0 / (2.0 * PI)) * theta;
    let v_v = -(1.0 / (2.0 * PI)) * ln_ratio;

    // Rotate back to global frame and project onto panel i's normal and tangent
    // Source velocity in global: (u_s * tj.tx + v_s * tj.nx, u_s * tj.ty + v_s * tj.ny)
    let vx_s = u_s * pj.tangent.0 + v_s * pj.normal.0;
    let vy_s = u_s * pj.tangent.1 + v_s * pj.normal.1;
    let vx_v = u_v * pj.tangent.0 + v_v * pj.normal.0;
    let vy_v = u_v * pj.tangent.1 + v_v * pj.normal.1;

    // Project onto panel i's normal and tangent
    let a_n = vx_s * pi.normal.0 + vy_s * pi.normal.1;
    let a_t = vx_s * pi.tangent.0 + vy_s * pi.tangent.1;
    let b_n = vx_v * pi.normal.0 + vy_v * pi.normal.1;
    let b_t = vx_v * pi.tangent.0 + vy_v * pi.tangent.1;

    (a_n, a_t, b_n, b_t)
}

/// Solve the Hess-Smith panel method for a single angle of attack.
///
/// Panels should be generated from [`panels_from_surface`]. The angle of attack
/// is in radians. Freestream velocity is normalized to 1.
///
/// Returns a [`PanelSolution`] with Cl, Cd, Cm, Cp distribution, and circulation.
pub fn solve(panels: &[Panel], alpha_rad: f64) -> Result<PanelSolution> {
    let n = panels.len();
    if n < MIN_PANELS {
        return Err(PavanError::InvalidGeometry(format!(
            "need at least {MIN_PANELS} panels, got {n}"
        )));
    }

    // Freestream components
    let cos_a = alpha_rad.cos();
    let sin_a = alpha_rad.sin();

    // Build (N+1) x (N+1) augmented system:
    //   Rows 0..N: zero normal velocity at each panel (no-penetration)
    //   Row N: Kutta condition (equal pressure at TE panels)
    // Columns 0..N: source strengths σ_j
    // Column N: vortex strength γ

    let size = n + 1;
    let mut matrix: Vec<Vec<f64>> = vec![vec![0.0; size]; size];
    let mut rhs: Vec<f64> = vec![0.0; size];

    // Store tangential influence for velocity computation later
    let mut a_t_store: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    let mut b_t_store: Vec<Vec<f64>> = vec![vec![0.0; n]; n];

    for i in 0..n {
        let mut b_n_sum = 0.0;

        for j in 0..n {
            let (a_n, a_t, b_n, b_t) = panel_influence(panels, i, j);

            // Source influence on normal velocity
            matrix[i][j] = a_n;
            // Vortex influence on normal velocity (all panels share same γ)
            b_n_sum += b_n;

            a_t_store[i][j] = a_t;
            b_t_store[i][j] = b_t;
        }

        // Vortex column
        matrix[i][n] = b_n_sum;

        // RHS: negative freestream normal component
        rhs[i] = -(cos_a * panels[i].normal.0 + sin_a * panels[i].normal.1);
    }

    // Kutta condition: sum of tangential velocities at first and last panels = 0
    // V_t(0) + V_t(N-1) = 0
    for j in 0..n {
        matrix[n][j] = a_t_store[0][j] + a_t_store[n - 1][j];
    }
    let b_t_sum: f64 = b_t_store[0]
        .iter()
        .zip(b_t_store[n - 1].iter())
        .map(|(a, b)| a + b)
        .sum();
    matrix[n][n] = b_t_sum;
    rhs[n] = -(cos_a * (panels[0].tangent.0 + panels[n - 1].tangent.0)
        + sin_a * (panels[0].tangent.1 + panels[n - 1].tangent.1));

    // Solve linear system
    let (lu, pivot) = hisab::num::lu_decompose(&matrix)
        .map_err(|e| PavanError::ComputationError(format!("LU decomposition failed: {e}")))?;
    let solution = hisab::num::lu_solve(&lu, &pivot, &rhs)
        .map_err(|e| PavanError::ComputationError(format!("LU solve failed: {e}")))?;

    let gamma = solution[n];

    // Compute tangential velocity and Cp at each panel
    let mut cp = Vec::with_capacity(n);
    let mut cl = 0.0;
    let mut cd = 0.0;
    let mut cm = 0.0;

    for i in 0..n {
        // Tangential velocity = freestream tangential + source + vortex contributions
        let mut vt = cos_a * panels[i].tangent.0 + sin_a * panels[i].tangent.1;

        for j in 0..n {
            vt += solution[j] * a_t_store[i][j];
            vt += gamma * b_t_store[i][j];
        }

        let cp_i = 1.0 - vt * vt;
        cp.push(cp_i);

        // Integrate forces: Cl = -∮ Cp·sin(θ)·ds, Cd = -∮ Cp·cos(θ)·ds (normalized)
        // Using panel geometry directly:
        let dx = panels[i].end.0 - panels[i].start.0;
        let dy = panels[i].end.1 - panels[i].start.1;

        // Force coefficients from pressure integration
        // Lift direction: perpendicular to freestream
        // In body frame: Cn = -∮ Cp·dx (normal force), Ca = ∮ Cp·dy (axial force)
        // Then: Cl = Cn·cos(α) - Ca·sin(α), Cd = Cn·sin(α) + Ca·cos(α)
        let cn_i = cp_i * dx; // normal force contribution (per unit chord)
        let ca_i = -cp_i * dy; // axial force contribution

        cl += cn_i;
        cd += ca_i;

        // Moment about quarter-chord: Cm = -∮ Cp·(x - 0.25)·dx (nose-up positive)
        let xc = panels[i].center.0;
        cm -= cp_i * (xc - 0.25) * dx;
    }

    // Rotate from body to wind frame
    let cl_wind = cl * cos_a - cd * sin_a;
    let cd_wind = cl * sin_a + cd * cos_a;

    Ok(PanelSolution {
        cl: cl_wind,
        cd: cd_wind,
        cm,
        cp,
        gamma,
    })
}

/// Solve the panel method for multiple angles of attack efficiently.
///
/// The influence matrix is decomposed once (LU) and reused for each angle,
/// since only the RHS changes with angle of attack.
pub fn solve_multi(panels: &[Panel], alphas: &[f64]) -> Result<Vec<PanelSolution>> {
    let n = panels.len();
    if n < MIN_PANELS {
        return Err(PavanError::InvalidGeometry(format!(
            "need at least {MIN_PANELS} panels, got {n}"
        )));
    }

    let size = n + 1;
    let mut matrix: Vec<Vec<f64>> = vec![vec![0.0; size]; size];

    // Store all influence coefficients
    let mut a_t_store: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    let mut b_t_store: Vec<Vec<f64>> = vec![vec![0.0; n]; n];

    for i in 0..n {
        let mut b_n_sum = 0.0;

        for j in 0..n {
            let (a_n, a_t, b_n, b_t) = panel_influence(panels, i, j);
            matrix[i][j] = a_n;
            b_n_sum += b_n;
            a_t_store[i][j] = a_t;
            b_t_store[i][j] = b_t;
        }

        matrix[i][n] = b_n_sum;
    }

    // Kutta condition row
    for j in 0..n {
        matrix[n][j] = a_t_store[0][j] + a_t_store[n - 1][j];
    }
    let b_t_sum: f64 = b_t_store[0]
        .iter()
        .zip(b_t_store[n - 1].iter())
        .map(|(a, b)| a + b)
        .sum();
    matrix[n][n] = b_t_sum;

    // LU decompose once
    let (lu, pivot) = hisab::num::lu_decompose(&matrix)
        .map_err(|e| PavanError::ComputationError(format!("LU decomposition failed: {e}")))?;

    let mut results = Vec::with_capacity(alphas.len());

    for &alpha in alphas {
        let cos_a = alpha.cos();
        let sin_a = alpha.sin();

        let mut rhs: Vec<f64> = panels
            .iter()
            .map(|p| -(cos_a * p.normal.0 + sin_a * p.normal.1))
            .collect();
        rhs.push(
            -(cos_a * (panels[0].tangent.0 + panels[n - 1].tangent.0)
                + sin_a * (panels[0].tangent.1 + panels[n - 1].tangent.1)),
        );

        let solution = hisab::num::lu_solve(&lu, &pivot, &rhs)
            .map_err(|e| PavanError::ComputationError(format!("LU solve failed: {e}")))?;

        let gamma = solution[n];

        let mut cp = Vec::with_capacity(n);
        let mut cl = 0.0;
        let mut cd = 0.0;
        let mut cm = 0.0;

        for i in 0..n {
            let mut vt = cos_a * panels[i].tangent.0 + sin_a * panels[i].tangent.1;
            for j in 0..n {
                vt += solution[j] * a_t_store[i][j];
                vt += gamma * b_t_store[i][j];
            }

            let cp_i = 1.0 - vt * vt;
            cp.push(cp_i);

            let dx = panels[i].end.0 - panels[i].start.0;
            let dy = panels[i].end.1 - panels[i].start.1;
            let cn_i = cp_i * dx;
            let ca_i = -cp_i * dy;
            cl += cn_i;
            cd += ca_i;

            // Moment about quarter-chord: Cm = -∮ Cp·(x - 0.25)·dx (nose-up positive)
            let xc = panels[i].center.0;
            cm -= cp_i * (xc - 0.25) * dx;
        }

        let cl_wind = cl * cos_a - cd * sin_a;
        let cd_wind = cl * sin_a + cd * cos_a;

        results.push(PanelSolution {
            cl: cl_wind,
            cd: cd_wind,
            cm,
            cp,
            gamma,
        });
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::airfoil::NacaProfile;

    fn naca0012_panels(n: usize) -> Vec<Panel> {
        let profile = NacaProfile::naca0012();
        let (upper, lower) = profile.surface_coordinates(n);
        panels_from_surface(&upper, &lower)
    }

    fn naca2412_panels(n: usize) -> Vec<Panel> {
        let profile = NacaProfile::naca2412();
        let (upper, lower) = profile.surface_coordinates(n);
        panels_from_surface(&upper, &lower)
    }

    // --- Panel construction ---

    #[test]
    fn panels_from_surface_count() {
        let panels = naca0012_panels(50);
        // ~2*(50-1) panels (upper + lower, minus shared LE point, plus closing panel)
        assert!(
            panels.len() > 80,
            "should have ~98 panels, got {}",
            panels.len()
        );
    }

    #[test]
    fn panels_have_positive_length() {
        let panels = naca0012_panels(50);
        for (i, p) in panels.iter().enumerate() {
            assert!(
                p.length > 0.0,
                "panel {i} has non-positive length {}",
                p.length
            );
        }
    }

    #[test]
    fn panel_normals_are_unit() {
        let panels = naca0012_panels(50);
        for (i, p) in panels.iter().enumerate() {
            let mag = (p.normal.0 * p.normal.0 + p.normal.1 * p.normal.1).sqrt();
            assert!(
                (mag - 1.0).abs() < 1e-10,
                "panel {i} normal not unit: mag={mag}"
            );
        }
    }

    #[test]
    fn panel_tangents_are_unit() {
        let panels = naca0012_panels(50);
        for (i, p) in panels.iter().enumerate() {
            let mag = (p.tangent.0 * p.tangent.0 + p.tangent.1 * p.tangent.1).sqrt();
            assert!(
                (mag - 1.0).abs() < 1e-10,
                "panel {i} tangent not unit: mag={mag}"
            );
        }
    }

    #[test]
    fn empty_surface_returns_empty() {
        let panels = panels_from_surface(&vec![], &vec![]);
        assert!(panels.is_empty());
    }

    #[test]
    fn single_point_surface_returns_empty() {
        let panels = panels_from_surface(&vec![(0.0, 0.0)], &vec![(0.0, 0.0)]);
        assert!(panels.is_empty());
    }

    // --- Solve correctness ---

    #[test]
    fn symmetric_airfoil_zero_aoa_zero_lift() {
        let panels = naca0012_panels(80);
        let sol = solve(&panels, 0.0).expect("solve should succeed");
        assert!(
            sol.cl.abs() < 0.02,
            "NACA 0012 at 0° should have Cl ≈ 0, got {}",
            sol.cl
        );
    }

    #[test]
    fn symmetric_airfoil_5deg() {
        let panels = naca0012_panels(100);
        let alpha = 5.0_f64.to_radians();
        let sol = solve(&panels, alpha).expect("solve");
        // Thin airfoil theory: Cl ≈ 0.55. Panel method should be close.
        assert!(
            (sol.cl - 0.55).abs() < 0.10,
            "NACA 0012 at 5° Cl should be ~0.55, got {}",
            sol.cl
        );
    }

    #[test]
    fn symmetric_airfoil_10deg() {
        let panels = naca0012_panels(100);
        let alpha = 10.0_f64.to_radians();
        let sol = solve(&panels, alpha).expect("solve");
        assert!(
            (sol.cl - 1.09).abs() < 0.15,
            "NACA 0012 at 10° Cl should be ~1.09, got {}",
            sol.cl
        );
    }

    #[test]
    fn cambered_airfoil_zero_aoa_positive_lift() {
        let panels = naca2412_panels(100);
        let sol = solve(&panels, 0.0).expect("solve");
        assert!(
            sol.cl > 0.05,
            "NACA 2412 at 0° should have positive Cl (camber effect), got {}",
            sol.cl
        );
    }

    #[test]
    fn negative_aoa_negative_lift() {
        let panels = naca0012_panels(80);
        let alpha = (-5.0_f64).to_radians();
        let sol = solve(&panels, alpha).expect("solve");
        assert!(
            sol.cl < -0.1,
            "negative AoA should produce negative Cl, got {}",
            sol.cl
        );
    }

    #[test]
    fn lift_increases_with_aoa() {
        let panels = naca0012_panels(80);
        let sol_2 = solve(&panels, 2.0_f64.to_radians()).expect("solve");
        let sol_5 = solve(&panels, 5.0_f64.to_radians()).expect("solve");
        let sol_8 = solve(&panels, 8.0_f64.to_radians()).expect("solve");
        assert!(sol_5.cl > sol_2.cl, "Cl should increase with AoA");
        assert!(sol_8.cl > sol_5.cl, "Cl should increase with AoA");
    }

    #[test]
    fn cp_distribution_has_correct_length() {
        let panels = naca0012_panels(50);
        let sol = solve(&panels, 0.0).expect("solve");
        assert_eq!(sol.cp.len(), panels.len());
    }

    #[test]
    fn symmetric_cp_at_zero_aoa() {
        let panels = naca0012_panels(80);
        let sol = solve(&panels, 0.0).expect("solve");
        // For symmetric airfoil at 0°, Cp should be roughly symmetric
        let n = sol.cp.len();
        let half = n / 2;
        let mut diff_sum = 0.0;
        for i in 0..half {
            diff_sum += (sol.cp[i] - sol.cp[n - 1 - i]).abs();
        }
        let avg_diff = diff_sum / half as f64;
        assert!(
            avg_diff < 0.1,
            "Cp should be approximately symmetric, avg diff = {avg_diff}"
        );
    }

    #[test]
    fn convergence_with_panel_count() {
        let sol_50 = solve(&naca0012_panels(30), 5.0_f64.to_radians()).expect("solve");
        let sol_100 = solve(&naca0012_panels(60), 5.0_f64.to_radians()).expect("solve");
        let sol_200 = solve(&naca0012_panels(120), 5.0_f64.to_radians()).expect("solve");

        // Difference should decrease with more panels
        let diff_1 = (sol_50.cl - sol_100.cl).abs();
        let diff_2 = (sol_100.cl - sol_200.cl).abs();
        assert!(
            diff_2 < diff_1,
            "Cl should converge: diff(50,100)={diff_1}, diff(100,200)={diff_2}"
        );
    }

    // --- solve_multi ---

    #[test]
    fn solve_multi_matches_individual() {
        let panels = naca0012_panels(60);
        let alphas = [0.0, 3.0_f64.to_radians(), 5.0_f64.to_radians()];
        let multi = solve_multi(&panels, &alphas).expect("solve_multi");

        for (i, &alpha) in alphas.iter().enumerate() {
            let single = solve(&panels, alpha).expect("solve");
            assert!(
                (multi[i].cl - single.cl).abs() < 1e-10,
                "solve_multi[{i}] Cl={} != solve Cl={}",
                multi[i].cl,
                single.cl
            );
            assert!(
                (multi[i].cd - single.cd).abs() < 1e-10,
                "solve_multi[{i}] Cd={} != solve Cd={}",
                multi[i].cd,
                single.cd
            );
            assert!(
                (multi[i].cm - single.cm).abs() < 1e-10,
                "solve_multi[{i}] Cm={} != solve Cm={}",
                multi[i].cm,
                single.cm
            );
        }
    }

    #[test]
    fn solve_multi_empty_alphas() {
        let panels = naca0012_panels(50);
        let results = solve_multi(&panels, &[]).expect("solve_multi");
        assert!(results.is_empty());
    }

    // --- Error cases ---

    #[test]
    fn too_few_panels_errors() {
        let panels = vec![Panel {
            start: (0.0, 0.0),
            end: (1.0, 0.0),
            center: (0.5, 0.0),
            normal: (0.0, 1.0),
            tangent: (1.0, 0.0),
            length: 1.0,
        }];
        assert!(solve(&panels, 0.0).is_err());
    }

    // --- Moment ---

    #[test]
    fn symmetric_zero_aoa_zero_moment() {
        let panels = naca0012_panels(80);
        let sol = solve(&panels, 0.0).expect("solve");
        assert!(
            sol.cm.abs() < 0.02,
            "NACA 0012 at 0° should have Cm ≈ 0, got {}",
            sol.cm
        );
    }

    #[test]
    fn cambered_airfoil_negative_moment() {
        let panels = naca2412_panels(100);
        let sol = solve(&panels, 0.0).expect("solve");
        // Cambered airfoils typically have negative (nose-down) Cm about c/4
        assert!(
            sol.cm < 0.0,
            "NACA 2412 should have negative Cm about c/4, got {}",
            sol.cm
        );
    }
}
