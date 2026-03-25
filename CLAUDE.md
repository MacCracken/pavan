# Pavan — Claude Code Instructions

## Project Identity

**Pavan** (Sanskrit: पवन — wind, breeze) — Aerodynamics engine for AGNOS

- **Type**: Flat library crate
- **License**: GPL-3.0
- **MSRV**: 1.89
- **Version**: SemVer 0.1.0

## Consumers

kiran/joshua (flight sim, projectiles, wind effects), impetus (aero forces on bodies), badal (atmospheric feed)

## Development Process

### P(-1): Scaffold Hardening (before any new features)
1. Test + benchmark sweep
2. Cleanliness check: fmt, clippy, audit, deny
3. Baseline benchmarks
4. Audit + refactor
5. Additional tests/benchmarks
6. Prove wins with benchmarks
7. Repeat if heavy

### Key Principles
- Never skip benchmarks
- Own the stack (hisab types)
- `#[non_exhaustive]` on enums, `#[must_use]` on pure fns, `#[inline]` on hot paths
- ISA atmosphere values must match standard tables
- All physics formulas must have correctness tests with known engineering values

## DO NOT
- Commit/push — user handles git
- NEVER use `gh` CLI
- No unwrap/panic in library code
- No skip benchmarks
