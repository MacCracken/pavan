# Pavan

**Pavan** (पवन — Sanskrit for "wind, breeze") — aerodynamics engine for the [AGNOS](https://github.com/MacCracken/agnosticos) ecosystem.

## Features

- **ISA Atmosphere** — temperature, pressure, density vs altitude, speed of sound, Mach number
- **NACA Airfoils** — 4-digit profile generation with camber and thickness distribution
- **Coefficients** — thin airfoil lift (Cl=2πα), induced drag, L/D ratio optimization
- **Forces** — lift, drag, Reynolds number, Sutherland viscosity, full force computation
- **Boundary Layer** — Blasius laminar, turbulent (1/7th power), transition detection, skin friction
- **Wind** — uniform fields, logarithmic profile, power law profile, wind chill
- **Vehicle** — AeroBody presets (light aircraft, glider), altitude-dependent force computation

## Quick Start

```rust
use pavan::{atmosphere, coefficients, forces};

let rho = atmosphere::standard_density(0.0);  // sea level: 1.225 kg/m³
let q = atmosphere::dynamic_pressure(rho, 100.0);  // 6125 Pa at 100 m/s

let alpha = 5.0_f64.to_radians();
let cl = coefficients::lift_coefficient_thin_airfoil(alpha);  // ~0.548
let l = forces::lift(q, 16.0, cl);  // lift on 16 m² wing
```

## License

GPL-3.0
