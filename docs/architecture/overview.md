# Pavan Architecture

## Module Map

```
pavan
├── error.rs          — PavanError (5 variants)
├── atmosphere.rs     — ISA standard model, dynamic pressure, Mach, speed of sound
├── airfoil.rs        — NACA 4-digit profiles, surface coordinates, camber/thickness
├── coefficients.rs   — Cl (thin airfoil), Cd (parasitic + induced), L/D optimization
├── forces.rs         — lift, drag, Reynolds, Sutherland viscosity, AeroForce, body_to_wind
├── boundary.rs       — Blasius/turbulent BL thickness, transition, skin friction
├── wind.rs           — WindField, log/power law profiles, wind chill
├── vehicle.rs        — AeroBody presets, altitude-dependent force computation
├── panel.rs          — 2D Hess-Smith panel method (Cp, Cl, Cd, Cm)
├── vlm.rs            — 3D Vortex Lattice Method (span loading, CDi, Oswald e)
├── compressible.rs   — isentropic flow, shocks, Prandtl-Meyer, Fanno, Rayleigh
├── stability.rs      — neutral point, static margin, flap/elevator, trim
├── propulsion.rs     — jet thrust lapse, TSFC, propeller efficiency, Froude
├── cfd.rs            — pravash N-S integration, BL mesh params (feature: cfd)
└── logging.rs        — tracing init (feature: logging)
```

## Consumers
- **kiran/joshua** — flight sim, projectile physics, wind effects on NPCs
- **impetus** — aerodynamic forces on rigid bodies
- **badal** — atmospheric state feeds aerodynamics
