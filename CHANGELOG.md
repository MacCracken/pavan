# Changelog

## [0.1.0] - 2026-03-25

Initial scaffold with real aerodynamics implementations.

### Modules
- **atmosphere** — ISA standard model (troposphere + tropopause), dynamic pressure, Mach, speed of sound
- **airfoil** — NACA 4-digit profile generation (camber + thickness), surface coordinates, presets
- **coefficients** — thin airfoil Cl, total Cd (parasitic + induced), L/D ratio, max L/D
- **forces** — lift, drag, Reynolds number, Sutherland viscosity, full force computation
- **boundary** — Blasius laminar, turbulent 1/7th power, transition detection, skin friction
- **wind** — WindField, logarithmic and power law profiles, wind chill
- **vehicle** — AeroBody (light aircraft, glider presets), altitude-dependent forces
