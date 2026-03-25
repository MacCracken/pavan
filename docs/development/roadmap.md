# Pavan Roadmap

## Status

**v0.1.0** — Initial scaffold with real aerodynamics implementations.

## Future Features (demand-gated)

### Compressible Flow
- Prandtl-Meyer expansion waves
- Normal/oblique shock relations
- Isentropic flow tables
- Fanno/Rayleigh flow

### Panel Methods
- 2D panel method (Hess-Smith)
- Vortex lattice method for 3D wings
- Induced drag from circulation distribution

### Stability & Control
- Longitudinal static stability (neutral point)
- Control surface effectiveness (elevator, aileron)
- Trim calculation

### Propulsion
- Thrust-specific fuel consumption
- Propeller efficiency model
- Jet engine thrust lapse with altitude

### CFD Integration (pravash)
- Couple with pravash Navier-Stokes solver
- Pressure/velocity field visualization data for soorat
- Boundary layer mesh generation

## v1.0.0 Criteria
- API frozen
- Zero unwrap/panic in library code
- 90%+ test coverage
- ISA values verified against standard tables
- Benchmark history with golden numbers
- 3+ downstream consumers
