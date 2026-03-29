# Changelog

## [1.1.0]

### Added
- **bridge** — cross-crate primitive-value bridges for badal (ISA flight conditions, effective airspeed), pravash (Cp to lift, Reynolds number), goonj (aeroacoustic power level, shock noise), ushma (stagnation temperature, heat flux, skin friction to heat)
- **integration/soorat** — feature-gated `soorat-compat` module with visualization data structures: `PanelMeshVisualization` (pressure-colored panels from solution), `AirfoilProfile` (upper/lower/camber from NACA), `FlowField2D` (velocity + scalar grid), `BoundaryLayerProfile` (thickness/transition along surface), `VortexVisualization` (filament segments with circulation)

### Updated
- hisab 1.1.0 -> 1.3.0, pravash 1.0.0 -> 1.2.0, zerocopy 0.8.47 -> 0.8.48

## [1.0.0] - 2026-03-25

API frozen. 15 modules, 276 tests (94.27% coverage), 45 benchmarks.

### Added
- **panel** — 2D Hess-Smith panel method (Cp distribution, Cl, Cd, Cm, solve/solve_multi)
- **vlm** — 3D Vortex Lattice Method (horseshoe vortices, Biot-Savart, span loading, induced drag, Oswald efficiency)
- **compressible** — isentropic flow relations, normal/oblique shock, Prandtl-Meyer expansion, Fanno/Rayleigh duct flow
- **stability** — longitudinal static stability (neutral point, static margin), flap/elevator effectiveness, trim solver
- **propulsion** — jet thrust lapse with altitude, TSFC, propeller efficiency (advance ratio model), Froude efficiency
- **cfd** — pravash Navier-Stokes integration (feature-gated `cfd`), airfoil rasterization, pressure/velocity field access, panel warm-start, boundary layer mesh parameters
- `body_to_wind()` utility in forces.rs for body→wind frame rotation

### Changed
- Dependencies updated: criterion 0.5 → 0.8, pravash 0.24 → 1.0
- `#[non_exhaustive]` added to all public structs (AeroForce, NacaProfile, AeroBody, WindField, Panel, PanelSolution, VlmSolution, WingGeometry, etc.)
- `#[inline]` added to all hot-path functions
- Magic numbers extracted to named constants (NACA coefficients, Sutherland's law, wind chill, vortex core cutoff)
- `SurfacePoints` type alias for airfoil coordinate output
- `TROPOPAUSE_TEMPERATURE` constant added to atmosphere
- ISA atmosphere docs updated for stratosphere support
- Panel/VLM post-processing extracted into shared helpers (reduced duplication)
- Re-exports expanded in lib.rs
- Degenerate geometry now logs tracing warnings instead of silent fallback

### Fixed
- Unused `PI` import in atmosphere.rs
- Unused `forces` import in examples/basic.rs
- Panel method Cm sign in solve_multi (was `+=`, corrected to `-=`)

## [0.1.0] - 2026-03-25

Initial scaffold with real aerodynamics implementations.

### Modules
- **atmosphere** — ISA standard model (troposphere + stratosphere), dynamic pressure, Mach, speed of sound
- **airfoil** — NACA 4-digit profile generation (camber + thickness), surface coordinates, presets
- **coefficients** — thin airfoil Cl, total Cd (parasitic + induced), L/D ratio, max L/D
- **forces** — lift, drag, Reynolds number, Sutherland viscosity, full force computation
- **boundary** — Blasius laminar, turbulent 1/7th power, transition detection, skin friction
- **wind** — WindField, logarithmic and power law profiles, wind chill
- **vehicle** — AeroBody (light aircraft, glider presets), altitude-dependent forces
