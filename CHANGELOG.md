# Changelog

## [Unreleased]

### Added
- **panel** — 2D Hess-Smith panel method (Cp distribution, Cl, Cd, Cm, solve/solve_multi)
- **vlm** — 3D Vortex Lattice Method (horseshoe vortices, Biot-Savart, span loading, induced drag, Oswald efficiency)
- **cfd** — pravash Navier-Stokes integration (feature-gated `cfd`), airfoil rasterization, pressure/velocity field access, panel warm-start

### Changed
- Dependencies updated: criterion 0.5 → 0.8, pravash 0.24 → 1.0
- `#[non_exhaustive]` added to all public structs (AeroForce, NacaProfile, AeroBody, WindField)
- `#[inline]` added to all hot-path functions
- Magic numbers extracted to named constants (NACA coefficients, Sutherland's law, wind chill)
- `SurfacePoints` type alias for airfoil coordinate output
- `TROPOPAUSE_TEMPERATURE` constant added
- ISA atmosphere docs updated for stratosphere support
- Re-exports expanded: Panel, PanelSolution, NacaProfile, AeroBody, WindField, WingGeometry, VlmSolution

### Fixed
- Unused `PI` import in atmosphere.rs
- Unused `forces` import in examples/basic.rs

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
