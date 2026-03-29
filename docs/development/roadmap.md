# Pavan Roadmap

## v1.0.0 — Released 2026-03-25

All criteria met:
- API frozen ✓
- Zero unwrap/panic in library code ✓
- 94.27% test coverage (276 tests) ✓
- ISA values verified against standard tables ✓
- 45 benchmarks baselined ✓
- 3 downstream consumers defined (kiran/joshua, impetus, badal)

---

## Cross-Crate Bridges

- [ ] **`bridge.rs` module** — primitive-value conversions for cross-crate aerodynamics
- [ ] **badal bridge**: altitude (m) → ISA atmosphere (density, pressure, temperature, viscosity)
- [ ] **pravash bridge**: surface pressure coefficient → fluid boundary condition; lift/drag coefficients → external forces
- [ ] **goonj bridge**: flow velocity (m/s), turbulent intensity → aeroacoustic source power; Mach number → shock noise level
- [ ] **ushma bridge**: stagnation temperature (K) → aerodynamic heating rate; skin friction → surface heat flux

---

## Soorat Integration (rendering visualization)

- [ ] **`integration/soorat.rs` module** — feature-gated `soorat-compat`
- [ ] **Panel mesh visualization**: panel vertices, normals, and pressure coefficients for colored surface rendering
- [ ] **Flow field**: velocity/pressure scalar and vector fields on regular grids for heatmap/streamline rendering
- [ ] **Airfoil profile**: NACA/custom airfoil coordinates for line/surface rendering
- [ ] **Boundary layer profile**: velocity and temperature profiles along surfaces for ribbon rendering
- [ ] **Wake/vortex visualization**: vortex filament positions and strengths for line rendering
