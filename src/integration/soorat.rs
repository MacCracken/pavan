//! Soorat integration — visualization data structures for aerodynamic analysis.
//!
//! Provides structured types that soorat can render: panel meshes with pressure,
//! airfoil profiles, flow fields, boundary layer ribbons, and wake/vortex lines.

use serde::{Deserialize, Serialize};

// ── Panel mesh visualization ───────────────────────────────────────────────

/// Panel mesh with pressure coefficients for colored surface rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PanelMeshVisualization {
    /// Panel endpoints `[x, y]` — the boundary of the airfoil.
    pub points: Vec<[f64; 2]>,
    /// Pressure coefficient at each panel center.
    pub cp: Vec<f64>,
    /// Lift coefficient from the solution.
    pub cl: f64,
    /// Drag coefficient from the solution.
    pub cd: f64,
    /// Minimum Cp (most negative = strongest suction).
    pub cp_min: f64,
    /// Maximum Cp.
    pub cp_max: f64,
}

impl PanelMeshVisualization {
    /// Create from pavan panels and a solution.
    #[must_use]
    pub fn from_panels_and_solution(
        panels: &[crate::panel::Panel],
        solution: &crate::panel::PanelSolution,
    ) -> Self {
        let mut points: Vec<[f64; 2]> = panels.iter().map(|p| [p.start.0, p.start.1]).collect();
        if let Some(last) = panels.last() {
            points.push([last.end.0, last.end.1]);
        }

        let cp_min = solution.cp.iter().cloned().fold(f64::MAX, f64::min);
        let cp_max = solution.cp.iter().cloned().fold(f64::MIN, f64::max);

        Self {
            points,
            cp: solution.cp.clone(),
            cl: solution.cl,
            cd: solution.cd,
            cp_min,
            cp_max,
        }
    }
}

// ── Airfoil profile ────────────────────────────────────────────────────────

/// Airfoil profile coordinates for line/surface rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AirfoilProfile {
    /// Upper surface points `[x, y]` from LE to TE.
    pub upper: Vec<[f64; 2]>,
    /// Lower surface points `[x, y]` from LE to TE.
    pub lower: Vec<[f64; 2]>,
    /// Camber line points `[x, y]`.
    pub camber: Vec<[f64; 2]>,
    /// NACA designation (e.g. "2412").
    pub designation: String,
    /// Chord length (1.0 for normalized).
    pub chord: f64,
}

impl AirfoilProfile {
    /// Create from a NACA profile with a given number of surface points.
    #[must_use]
    pub fn from_naca(profile: &crate::airfoil::NacaProfile, num_points: usize) -> Self {
        let (upper_pts, lower_pts) = profile.surface_coordinates(num_points);
        let upper: Vec<[f64; 2]> = upper_pts.iter().map(|&(x, y)| [x, y]).collect();
        let lower: Vec<[f64; 2]> = lower_pts.iter().map(|&(x, y)| [x, y]).collect();

        // Camber line is the average of upper and lower at each x
        let camber: Vec<[f64; 2]> = upper
            .iter()
            .zip(lower.iter())
            .map(|(u, l)| [u[0], (u[1] + l[1]) * 0.5])
            .collect();

        let d1 = (profile.max_camber * 100.0).round() as u32;
        let d2 = (profile.camber_position * 10.0).round() as u32;
        let d34 = (profile.max_thickness * 100.0).round() as u32;

        Self {
            upper,
            lower,
            camber,
            designation: format!("{d1}{d2}{d34:02}"),
            chord: 1.0,
        }
    }
}

// ── Flow field ─────────────────────────────────────────────────────────────

/// A 2D flow field for heatmap/streamline/arrow rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FlowField2D {
    /// Velocity vectors `[u, v]` at each grid point (m/s).
    /// Flattened row-major: `velocities[y * nx + x]`.
    pub velocities: Vec<[f64; 2]>,
    /// Scalar field (pressure, Cp, or speed) at each grid point.
    pub scalar: Vec<f64>,
    /// Grid dimensions (nx, ny).
    pub dimensions: [usize; 2],
    /// World-space origin `[x, y]` in metres.
    pub origin: [f64; 2],
    /// Grid spacing in metres.
    pub spacing: f64,
    /// Name of the scalar field (e.g. "pressure", "speed", "Cp").
    pub scalar_name: String,
    /// Maximum velocity magnitude (for normalization).
    pub max_speed: f64,
}

impl FlowField2D {
    /// Create a uniform freestream flow field.
    #[must_use]
    pub fn uniform(
        nx: usize,
        ny: usize,
        origin: [f64; 2],
        spacing: f64,
        freestream_u: f64,
        freestream_v: f64,
    ) -> Self {
        let count = nx * ny;
        let speed = (freestream_u * freestream_u + freestream_v * freestream_v).sqrt();
        Self {
            velocities: vec![[freestream_u, freestream_v]; count],
            scalar: vec![speed; count],
            dimensions: [nx, ny],
            origin,
            spacing,
            scalar_name: "speed".to_string(),
            max_speed: speed,
        }
    }
}

// ── Boundary layer profile ─────────────────────────────────────────────────

/// Boundary layer profile data for ribbon rendering along a surface.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoundaryLayerProfile {
    /// Surface x-positions where the profile is sampled.
    pub x_stations: Vec<f64>,
    /// Boundary layer thickness at each station (m).
    pub thickness: Vec<f64>,
    /// Whether each station is laminar (true) or turbulent (false).
    pub is_laminar: Vec<bool>,
    /// Skin friction coefficient at each station.
    pub cf: Vec<f64>,
}

impl BoundaryLayerProfile {
    /// Compute boundary layer profile along a flat plate.
    ///
    /// `chord`: plate length (m). `re_chord`: Reynolds number based on chord.
    /// `num_stations`: number of sample points along the chord.
    #[must_use]
    pub fn flat_plate(chord: f64, re_chord: f64, num_stations: usize) -> Self {
        let n = num_stations.max(2);
        let mut x_stations = Vec::with_capacity(n);
        let mut thickness = Vec::with_capacity(n);
        let mut is_laminar = Vec::with_capacity(n);
        let mut cf = Vec::with_capacity(n);

        for i in 0..n {
            let t = (i as f64 + 0.5) / n as f64; // avoid x=0
            let x = t * chord;
            let re_x = re_chord * t;

            x_stations.push(x);
            let lam = !crate::boundary::is_turbulent(re_x);
            is_laminar.push(lam);

            if lam {
                thickness.push(crate::boundary::blasius_thickness(x, re_x));
                cf.push(crate::boundary::skin_friction_laminar(re_x));
            } else {
                thickness.push(crate::boundary::turbulent_thickness(x, re_x));
                cf.push(crate::boundary::skin_friction_turbulent(re_x));
            }
        }

        Self {
            x_stations,
            thickness,
            is_laminar,
            cf,
        }
    }
}

// ── Wake/vortex visualization ──────────────────────────────────────────────

/// Vortex filament data for line rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VortexVisualization {
    /// Vortex filament segments.
    pub filaments: Vec<VortexFilament>,
}

/// A single vortex filament segment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VortexFilament {
    /// Start point `[x, y, z]`.
    pub start: [f64; 3],
    /// End point `[x, y, z]`.
    pub end: [f64; 3],
    /// Circulation strength (m²/s).
    pub gamma: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn airfoil_profile_naca0012() {
        let naca = crate::airfoil::NacaProfile::naca0012();
        let prof = AirfoilProfile::from_naca(&naca, 20);
        assert_eq!(prof.upper.len(), 20);
        assert_eq!(prof.lower.len(), 20);
        assert_eq!(prof.designation, "0012");
        // Symmetric → camber line at y≈0
        for c in &prof.camber {
            assert!(c[1].abs() < 0.01, "symmetric airfoil camber should be ~0");
        }
    }

    #[test]
    fn airfoil_profile_cambered() {
        let naca = crate::airfoil::NacaProfile::from_digits(2, 4, 1, 2);
        let prof = AirfoilProfile::from_naca(&naca, 30);
        assert_eq!(prof.designation, "2412");
        // Cambered → some camber line points above zero
        let has_camber = prof.camber.iter().any(|c| c[1] > 0.001);
        assert!(has_camber);
    }

    #[test]
    fn flow_field_uniform() {
        let field = FlowField2D::uniform(5, 5, [0.0, 0.0], 0.1, 10.0, 0.0);
        assert_eq!(field.velocities.len(), 25);
        assert!((field.max_speed - 10.0).abs() < 0.01);
    }

    #[test]
    fn boundary_layer_flat_plate() {
        // Re = 1e6 over 1m chord
        let bl = BoundaryLayerProfile::flat_plate(1.0, 1e6, 10);
        assert_eq!(bl.x_stations.len(), 10);
        assert_eq!(bl.thickness.len(), 10);
        // Thickness should grow along the plate
        assert!(bl.thickness.last().unwrap() > bl.thickness.first().unwrap());
        // Should have laminar start, turbulent end
        assert!(bl.is_laminar[0]);
        assert!(!bl.is_laminar.last().unwrap());
    }

    #[test]
    fn vortex_viz_serializes() {
        let viz = VortexVisualization {
            filaments: vec![VortexFilament {
                start: [0.0, 0.0, 0.0],
                end: [1.0, 0.0, 0.0],
                gamma: 5.0,
            }],
        };
        let json = serde_json::to_string(&viz);
        assert!(json.is_ok());
    }

    #[test]
    fn panel_mesh_from_panels() {
        let naca = crate::airfoil::NacaProfile::naca0012();
        let (upper, lower) = naca.surface_coordinates(20);
        let panels = crate::panel::panels_from_surface(&upper, &lower);
        let solution = crate::panel::solve(&panels, 0.0).unwrap();
        let viz = PanelMeshVisualization::from_panels_and_solution(&panels, &solution);
        assert!(!viz.points.is_empty());
        assert!(!viz.cp.is_empty());
        assert!(viz.cp_min <= viz.cp_max);
    }

    #[test]
    fn flow_field_serializes() {
        let field = FlowField2D::uniform(2, 2, [0.0; 2], 1.0, 5.0, 3.0);
        let json = serde_json::to_string(&field);
        assert!(json.is_ok());
    }
}
