# Pavan Architecture

## Module Map

```
pavan
├── error.rs        — PavanError (5 variants)
├── atmosphere.rs   — ISA standard model, dynamic pressure, Mach, speed of sound
├── airfoil.rs      — NACA 4-digit profiles, surface coordinates, camber/thickness
├── coefficients.rs — Cl (thin airfoil), Cd (parasitic + induced), L/D optimization
├── forces.rs       — lift, drag, Reynolds, Sutherland viscosity, AeroForce
├── boundary.rs     — Blasius/turbulent BL thickness, transition, skin friction
├── wind.rs         — WindField, log/power law profiles, wind chill
└── vehicle.rs      — AeroBody presets, altitude-dependent force computation
```

## Consumers
- **kiran/joshua** — flight sim, projectile physics, wind effects on NPCs
- **impetus** — aerodynamic forces on rigid bodies
- **badal** — atmospheric state feeds aerodynamics
