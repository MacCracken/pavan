//! CFD integration with pravash Navier-Stokes solver.
//!
//! Bridges pavan's aerodynamic types with pravash's grid-based fluid solver
//! for viscous flow simulation around airfoils. Feature-gated behind `cfd`.

use serde::{Deserialize, Serialize};

use pravash::grid::{BoundaryCondition, FluidGrid, GridConfig};

use crate::airfoil::NacaProfile;
use crate::atmosphere;
use crate::error::{PavanError, Result};
use crate::forces;
use crate::panel::{Panel, PanelSolution};

/// Configuration for a 2D CFD simulation around an airfoil.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct AirfoilCfdConfig {
    /// Grid cells in x (streamwise).
    pub grid_nx: usize,
    /// Grid cells in y (crossflow).
    pub grid_ny: usize,
    /// Physical domain size (width, height) in meters.
    pub domain_size: (f64, f64),
    /// Timestep in seconds. 0.0 for auto-CFL each step.
    pub dt: f64,
    /// Freestream velocity (m/s).
    pub freestream_velocity: f64,
    /// Altitude (m) — used to derive density and viscosity from ISA.
    pub altitude: f64,
    /// Angle of attack (radians).
    pub angle_of_attack_rad: f64,
    /// Number of pressure Poisson iterations per step.
    pub pressure_iterations: usize,
    /// Use multigrid pressure solver (faster for large grids).
    pub use_multigrid: bool,
    /// Number of surface coordinate points per side for rasterization.
    pub surface_points: usize,
}

impl AirfoilCfdConfig {
    /// Sensible defaults for a given freestream velocity and altitude.
    ///
    /// Grid is 200×100, domain is 10×5 chord-lengths, auto-CFL timestep.
    #[must_use]
    pub fn default_for(velocity: f64, altitude: f64) -> Self {
        Self {
            grid_nx: 200,
            grid_ny: 100,
            domain_size: (10.0, 5.0),
            dt: 0.0,
            freestream_velocity: velocity,
            altitude,
            angle_of_attack_rad: 0.0,
            pressure_iterations: 40,
            use_multigrid: true,
            surface_points: 80,
        }
    }
}

/// Snapshot of a CFD solution at a given timestep.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct CfdSnapshot {
    /// Number of steps completed.
    pub step: usize,
    /// Simulation time (seconds).
    pub time: f64,
    /// Lift coefficient.
    pub cl: f64,
    /// Drag coefficient (pressure).
    pub cd: f64,
    /// Pitch moment coefficient about quarter-chord.
    pub cm: f64,
    /// Maximum flow speed in the domain.
    pub max_speed: f64,
    /// Total kinetic energy in the domain.
    pub kinetic_energy: f64,
}

/// Manages a 2D airfoil CFD simulation wrapping [`pravash::grid::FluidGrid`].
pub struct AirfoilCfd {
    grid: FluidGrid,
    grid_config: GridConfig,
    /// Grid cells marked as solid (airfoil interior).
    solid: Vec<bool>,
    /// Surface boundary cells with outward normals for force integration.
    surface_cells: Vec<(usize, usize)>,
    surface_normals: Vec<(f64, f64)>,
    /// Freestream velocity components (vx, vy).
    freestream: (f64, f64),
    /// Dynamic pressure for coefficient normalization.
    q_inf: f64,
    /// Chord length in grid cells.
    chord_cells: usize,
    /// Grid cell size (m).
    dx: f64,
    /// ISA viscosity for auto-CFL.
    viscosity: f64,
    /// Auto-CFL mode.
    auto_dt: bool,
    /// Current step count.
    steps: usize,
    /// Current simulation time.
    time: f64,
}

impl AirfoilCfd {
    /// Create a new CFD simulation for the given airfoil and configuration.
    pub fn new(profile: &NacaProfile, config: &AirfoilCfdConfig) -> Result<Self> {
        if config.freestream_velocity <= 0.0 {
            return Err(PavanError::InvalidVelocity(
                "freestream velocity must be positive".into(),
            ));
        }
        if config.grid_nx < 20 || config.grid_ny < 10 {
            return Err(PavanError::InvalidGeometry(
                "grid must be at least 20×10".into(),
            ));
        }

        let dx = config.domain_size.0 / config.grid_nx as f64;

        // ISA properties at altitude
        let rho = atmosphere::standard_density(config.altitude);
        let temp = atmosphere::standard_temperature(config.altitude);
        let mu = forces::air_dynamic_viscosity(temp);
        let viscosity = mu / rho; // kinematic viscosity

        // Freestream components
        let cos_a = config.angle_of_attack_rad.cos();
        let sin_a = config.angle_of_attack_rad.sin();
        let vx_inf = config.freestream_velocity * cos_a;
        let vy_inf = config.freestream_velocity * sin_a;
        let q_inf = atmosphere::dynamic_pressure(rho, config.freestream_velocity);

        // Create grid
        let mut grid = FluidGrid::new(config.grid_nx, config.grid_ny, dx).map_err(|e| {
            PavanError::ComputationError(format!("failed to create fluid grid: {e}"))
        })?;

        // Initialize freestream everywhere
        for val in grid.vx.iter_mut() {
            *val = vx_inf;
        }
        for val in grid.vy.iter_mut() {
            *val = vy_inf;
        }
        for val in grid.density.iter_mut() {
            *val = rho;
        }

        // Rasterize airfoil onto grid
        let chord_cells = (config.grid_nx as f64 * 0.1).max(10.0) as usize; // 10% of grid width
        let x_offset = config.grid_nx / 4; // airfoil starts at 25% of domain
        let y_center = config.grid_ny / 2;

        let (upper, lower) = profile.surface_coordinates(config.surface_points);

        let mut solid = vec![false; config.grid_nx * config.grid_ny];
        let mut surface_cells = Vec::new();
        let mut surface_normals = Vec::new();

        // Mark solid cells inside the airfoil
        Self::rasterize_airfoil(
            &upper,
            &lower,
            chord_cells,
            x_offset,
            y_center,
            config.grid_nx,
            config.grid_ny,
            &mut solid,
            &mut surface_cells,
            &mut surface_normals,
        );

        // Zero velocity in solid cells
        for &(cx, cy) in &surface_cells {
            let idx = grid.idx(cx, cy);
            grid.vx[idx] = 0.0;
            grid.vy[idx] = 0.0;
        }

        // Grid solver config
        let auto_dt = config.dt <= 0.0;
        let dt = if auto_dt {
            grid.cfl_dt(viscosity, 0.3)
        } else {
            config.dt
        };

        let mut grid_config = GridConfig::default();
        grid_config.dt = dt;
        grid_config.viscosity = viscosity;
        grid_config.diffusion_iterations = 20;
        grid_config.pressure_iterations = config.pressure_iterations;
        grid_config.boundary = BoundaryCondition::FreeSlip;
        grid_config.vorticity_confinement = 0.0;
        grid_config.buoyancy_alpha = 0.0;
        grid_config.ambient_density = rho;
        grid_config.use_maccormack = true;
        grid_config.use_bfecc = false;
        grid_config.smagorinsky_cs = 0.0;
        grid_config.use_multigrid = config.use_multigrid;

        Ok(Self {
            grid,
            grid_config,
            solid,
            surface_cells,
            surface_normals,
            freestream: (vx_inf, vy_inf),
            q_inf,
            chord_cells,
            dx,
            viscosity,
            auto_dt,
            steps: 0,
            time: 0.0,
        })
    }

    /// Advance the simulation by one timestep.
    pub fn step(&mut self) -> Result<CfdSnapshot> {
        // Auto-CFL timestep
        if self.auto_dt {
            let dt = self.grid.cfl_dt(self.viscosity, 0.3).max(1e-8);
            self.grid_config.dt = dt;
        }

        // Step the grid solver
        self.grid
            .step(&self.grid_config)
            .map_err(|e| PavanError::ComputationError(format!("grid step failed: {e}")))?;

        // Enforce boundary conditions
        self.enforce_boundaries();

        self.steps += 1;
        self.time += self.grid_config.dt;

        Ok(self.snapshot())
    }

    /// Run multiple timesteps, returning a snapshot after each.
    pub fn run(&mut self, num_steps: usize) -> Result<Vec<CfdSnapshot>> {
        let mut history = Vec::with_capacity(num_steps);
        for _ in 0..num_steps {
            history.push(self.step()?);
        }
        Ok(history)
    }

    /// Current solution snapshot without advancing.
    #[must_use]
    pub fn snapshot(&self) -> CfdSnapshot {
        let (cl, cd, cm) = self.integrate_surface_forces();
        CfdSnapshot {
            step: self.steps,
            time: self.time,
            cl,
            cd,
            cm,
            max_speed: self.grid.max_speed(),
            kinetic_energy: self.grid.total_kinetic_energy(),
        }
    }

    /// Access the pressure field (row-major, nx × ny).
    #[must_use]
    pub fn pressure_field(&self) -> &[f64] {
        &self.grid.pressure
    }

    /// Access the velocity fields (vx, vy), each row-major nx × ny.
    #[must_use]
    pub fn velocity_field(&self) -> (&[f64], &[f64]) {
        (&self.grid.vx, &self.grid.vy)
    }

    /// Number of cells on the airfoil surface boundary.
    #[must_use]
    pub fn surface_cell_count(&self) -> usize {
        self.surface_cells.len()
    }

    /// Grid dimensions (nx, ny).
    #[must_use]
    pub fn grid_size(&self) -> (usize, usize) {
        (self.grid.nx, self.grid.ny)
    }

    // ── Internal ──────────────────────────────────────────────────────────

    /// Enforce no-slip on airfoil surface and freestream at domain boundaries.
    fn enforce_boundaries(&mut self) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;

        // No-slip on solid cells
        for y in 0..ny {
            for x in 0..nx {
                let idx = y * nx + x;
                if self.solid[idx] {
                    self.grid.vx[idx] = 0.0;
                    self.grid.vy[idx] = 0.0;
                }
            }
        }

        // Freestream at inflow (left column)
        for y in 0..ny {
            let idx = y * nx;
            self.grid.vx[idx] = self.freestream.0;
            self.grid.vy[idx] = self.freestream.1;
        }

        // Outflow (right column): zero-gradient (copy from interior)
        for y in 0..ny {
            let idx_out = y * nx + (nx - 1);
            let idx_in = y * nx + (nx - 2);
            self.grid.vx[idx_out] = self.grid.vx[idx_in];
            self.grid.vy[idx_out] = self.grid.vy[idx_in];
        }
    }

    /// Integrate pressure around the airfoil surface to get Cl, Cd, Cm.
    fn integrate_surface_forces(&self) -> (f64, f64, f64) {
        if self.surface_cells.is_empty() || self.q_inf <= 0.0 {
            return (0.0, 0.0, 0.0);
        }

        let nx = self.grid.nx;
        let x_offset = nx / 4;
        let cos_a = self.freestream.0 / (self.freestream.0.hypot(self.freestream.1)).max(1e-20);
        let sin_a = self.freestream.1 / (self.freestream.0.hypot(self.freestream.1)).max(1e-20);

        let mut cn = 0.0; // normal to chord (body frame)
        let mut ca = 0.0; // along chord (body frame)
        let mut cm = 0.0;

        let chord_m = self.chord_cells as f64 * self.dx;

        for (i, &(cx, cy)) in self.surface_cells.iter().enumerate() {
            let idx = cy * nx + cx;
            let p = self.grid.pressure[idx];
            let (nx_i, ny_i) = self.surface_normals[i];

            // Force contribution from this cell
            let fx = p * nx_i * self.dx;
            let fy = p * ny_i * self.dx;

            cn += fy;
            ca += fx;

            // Moment about quarter-chord
            let x_body = (cx as f64 - x_offset as f64) * self.dx;
            cm -= fy * (x_body - 0.25 * chord_m);
        }

        // Normalize by dynamic pressure × chord
        let norm = self.q_inf * chord_m;
        if norm <= 0.0 {
            return (0.0, 0.0, 0.0);
        }

        let cn_norm = cn / norm;
        let ca_norm = ca / norm;
        let cm_norm = cm / (norm * chord_m);

        // Rotate to wind frame
        let cl = cn_norm * cos_a - ca_norm * sin_a;
        let cd = cn_norm * sin_a + ca_norm * cos_a;

        (cl, cd, cm_norm)
    }

    /// Rasterize an airfoil onto the grid, marking solid cells and extracting surface cells.
    #[allow(clippy::too_many_arguments)]
    fn rasterize_airfoil(
        upper: &[(f64, f64)],
        lower: &[(f64, f64)],
        chord_cells: usize,
        x_offset: usize,
        y_center: usize,
        nx: usize,
        ny: usize,
        solid: &mut [bool],
        surface_cells: &mut Vec<(usize, usize)>,
        surface_normals: &mut Vec<(f64, f64)>,
    ) {
        // For each x column in the chord, find upper and lower y bounds
        for gx in 0..chord_cells {
            let x_frac = gx as f64 / chord_cells.max(1) as f64;

            // Interpolate upper and lower y at this x
            let y_upper = Self::interpolate_y(upper, x_frac);
            let y_lower = Self::interpolate_y(lower, x_frac);

            let gy_upper = y_center as f64 + y_upper * chord_cells as f64;
            let gy_lower = y_center as f64 + y_lower * chord_cells as f64;

            let gy_lo = (gy_lower.floor() as usize).max(1).min(ny - 2);
            let gy_hi = (gy_upper.ceil() as usize).max(1).min(ny - 2);

            let cx = x_offset + gx;
            if cx >= nx - 1 {
                continue;
            }

            // Mark interior as solid
            for gy in gy_lo..=gy_hi {
                let idx = gy * nx + cx;
                solid[idx] = true;
            }

            // Upper surface cell
            if gy_hi < ny - 1 {
                surface_cells.push((cx, gy_hi));
                // Approximate normal: pointing upward
                let dy_dx = Self::slope(upper, x_frac);
                let len = (1.0 + dy_dx * dy_dx).sqrt();
                surface_normals.push((-dy_dx / len, 1.0 / len));
            }

            // Lower surface cell
            if gy_lo > 0 {
                surface_cells.push((cx, gy_lo));
                let dy_dx = Self::slope(lower, x_frac);
                let len = (1.0 + dy_dx * dy_dx).sqrt();
                surface_normals.push((dy_dx / len, -1.0 / len));
            }
        }
    }

    /// Interpolate y value from surface points at a given x fraction.
    fn interpolate_y(points: &[(f64, f64)], x: f64) -> f64 {
        if points.is_empty() {
            return 0.0;
        }
        if points.len() == 1 {
            return points[0].1;
        }

        // Linear search (points are sorted by x: 0→1)
        for i in 0..points.len() - 1 {
            let (x0, y0) = points[i];
            let (x1, y1) = points[i + 1];
            if x >= x0 && x <= x1 {
                if (x1 - x0).abs() < f64::EPSILON {
                    return y0;
                }
                let t = (x - x0) / (x1 - x0);
                return y0 + t * (y1 - y0);
            }
        }

        points.last().map_or(0.0, |p| p.1)
    }

    /// Approximate surface slope dy/dx at a given x fraction.
    fn slope(points: &[(f64, f64)], x: f64) -> f64 {
        let h = 0.01;
        let y1 = Self::interpolate_y(points, (x + h).min(1.0));
        let y0 = Self::interpolate_y(points, (x - h).max(0.0));
        (y1 - y0) / (2.0 * h)
    }
}

/// Seed CFD velocity/pressure fields from a panel method inviscid solution.
///
/// This "warm start" accelerates CFD convergence by initializing with the
/// potential flow solution rather than uniform freestream.
pub fn init_from_panel(
    cfd: &mut AirfoilCfd,
    panel_sol: &PanelSolution,
    panels: &[Panel],
) -> Result<()> {
    if panels.is_empty() {
        return Err(PavanError::InvalidGeometry(
            "no panels for initialization".into(),
        ));
    }

    let nx = cfd.grid.nx;
    let ny = cfd.grid.ny;
    let x_offset = nx / 4;

    // Seed pressure field from Cp distribution
    // Map panel Cp values onto nearby grid cells
    for (i, panel) in panels.iter().enumerate() {
        if i >= panel_sol.cp.len() {
            break;
        }
        let cp = panel_sol.cp[i];

        let gx = x_offset as f64 + panel.center.0 * cfd.chord_cells as f64;
        let gy = (ny / 2) as f64 + panel.center.1 * cfd.chord_cells as f64;

        let ix = gx.round() as usize;
        let iy = gy.round() as usize;

        if ix < nx && iy < ny {
            let idx = iy * nx + ix;
            if !cfd.solid[idx] {
                cfd.grid.pressure[idx] = cp * cfd.q_inf;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_config() -> AirfoilCfdConfig {
        AirfoilCfdConfig {
            grid_nx: 60,
            grid_ny: 30,
            domain_size: (6.0, 3.0),
            dt: 0.001,
            freestream_velocity: 50.0,
            altitude: 0.0,
            angle_of_attack_rad: 0.0,
            pressure_iterations: 20,
            use_multigrid: false,
            surface_points: 40,
        }
    }

    #[test]
    fn default_for_produces_valid_config() {
        let cfg = AirfoilCfdConfig::default_for(60.0, 0.0);
        assert_eq!(cfg.grid_nx, 200);
        assert_eq!(cfg.grid_ny, 100);
        assert!((cfg.freestream_velocity - 60.0).abs() < f64::EPSILON);
    }

    #[test]
    fn new_creates_grid_successfully() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let cfd = AirfoilCfd::new(&profile, &cfg);
        assert!(
            cfd.is_ok(),
            "AirfoilCfd::new should succeed: {:?}",
            cfd.err()
        );
    }

    #[test]
    fn new_rejects_zero_velocity() {
        let profile = NacaProfile::naca0012();
        let mut cfg = default_config();
        cfg.freestream_velocity = 0.0;
        assert!(AirfoilCfd::new(&profile, &cfg).is_err());
    }

    #[test]
    fn new_rejects_tiny_grid() {
        let profile = NacaProfile::naca0012();
        let mut cfg = default_config();
        cfg.grid_nx = 5;
        cfg.grid_ny = 5;
        assert!(AirfoilCfd::new(&profile, &cfg).is_err());
    }

    #[test]
    fn surface_cells_are_populated() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let cfd = AirfoilCfd::new(&profile, &cfg).expect("new");
        assert!(
            cfd.surface_cell_count() > 5,
            "should have surface cells, got {}",
            cfd.surface_cell_count()
        );
    }

    #[test]
    fn grid_size_matches_config() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let cfd = AirfoilCfd::new(&profile, &cfg).expect("new");
        assert_eq!(cfd.grid_size(), (60, 30));
    }

    #[test]
    fn step_produces_valid_snapshot() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let mut cfd = AirfoilCfd::new(&profile, &cfg).expect("new");
        let snap = cfd.step().expect("step");
        assert_eq!(snap.step, 1);
        assert!(snap.time > 0.0);
        assert!(snap.cl.is_finite());
        assert!(snap.cd.is_finite());
        assert!(snap.max_speed >= 0.0);
    }

    #[test]
    fn run_returns_correct_count() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let mut cfd = AirfoilCfd::new(&profile, &cfg).expect("new");
        let history = cfd.run(5).expect("run");
        assert_eq!(history.len(), 5);
        assert_eq!(history[4].step, 5);
    }

    #[test]
    fn snapshot_matches_step() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let mut cfd = AirfoilCfd::new(&profile, &cfg).expect("new");
        let step_snap = cfd.step().expect("step");
        let snap = cfd.snapshot();
        assert!((step_snap.cl - snap.cl).abs() < f64::EPSILON);
        assert_eq!(step_snap.step, snap.step);
    }

    #[test]
    fn pressure_field_has_correct_size() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let cfd = AirfoilCfd::new(&profile, &cfg).expect("new");
        assert_eq!(cfd.pressure_field().len(), 60 * 30);
    }

    #[test]
    fn velocity_field_has_correct_size() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let cfd = AirfoilCfd::new(&profile, &cfg).expect("new");
        let (vx, vy) = cfd.velocity_field();
        assert_eq!(vx.len(), 60 * 30);
        assert_eq!(vy.len(), 60 * 30);
    }

    #[test]
    fn freestream_at_angle() {
        let profile = NacaProfile::naca0012();
        let mut cfg = default_config();
        cfg.angle_of_attack_rad = 5.0_f64.to_radians();
        let cfd = AirfoilCfd::new(&profile, &cfg).expect("new");
        let expected_vx = 50.0 * 5.0_f64.to_radians().cos();
        let expected_vy = 50.0 * 5.0_f64.to_radians().sin();
        assert!((cfd.freestream.0 - expected_vx).abs() < 0.01);
        assert!((cfd.freestream.1 - expected_vy).abs() < 0.01);
    }

    #[test]
    fn auto_cfl_timestep() {
        let profile = NacaProfile::naca0012();
        let mut cfg = default_config();
        cfg.dt = 0.0; // auto
        let mut cfd = AirfoilCfd::new(&profile, &cfg).expect("new");
        let snap = cfd.step().expect("step");
        assert!(snap.time > 0.0, "auto-CFL should produce positive timestep");
    }

    #[test]
    fn init_from_panel_succeeds() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let mut cfd = AirfoilCfd::new(&profile, &cfg).expect("new");

        let (upper, lower) = profile.surface_coordinates(60);
        let panels = crate::panel::panels_from_surface(&upper, &lower);
        let sol = crate::panel::solve(&panels, 0.0).expect("panel solve");

        let result = init_from_panel(&mut cfd, &sol, &panels);
        assert!(result.is_ok());
    }

    #[test]
    fn init_from_panel_rejects_empty() {
        let profile = NacaProfile::naca0012();
        let cfg = default_config();
        let mut cfd = AirfoilCfd::new(&profile, &cfg).expect("new");

        let sol = PanelSolution {
            cl: 0.0,
            cd: 0.0,
            cm: 0.0,
            cp: vec![],
            gamma: 0.0,
        };
        assert!(init_from_panel(&mut cfd, &sol, &[]).is_err());
    }

    #[test]
    fn serde_round_trip_snapshot() {
        let snap = CfdSnapshot {
            step: 10,
            time: 0.01,
            cl: 0.5,
            cd: 0.03,
            cm: -0.05,
            max_speed: 55.0,
            kinetic_energy: 1000.0,
        };
        let json = serde_json::to_string(&snap).expect("serialize");
        let back: CfdSnapshot = serde_json::from_str(&json).expect("deserialize");
        assert!((back.cl - snap.cl).abs() < f64::EPSILON);
        assert_eq!(back.step, snap.step);
    }

    #[test]
    fn serde_round_trip_config() {
        let cfg = AirfoilCfdConfig::default_for(60.0, 2000.0);
        let json = serde_json::to_string(&cfg).expect("serialize");
        let back: AirfoilCfdConfig = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(back.grid_nx, cfg.grid_nx);
        assert!((back.altitude - cfg.altitude).abs() < f64::EPSILON);
    }
}
