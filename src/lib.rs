//! # Pavan
//!
//! **Pavan** (पवन — Sanskrit for "wind, breeze") — aerodynamics engine for the AGNOS ecosystem.
//!
//! Provides ISA atmosphere model, NACA airfoil generation, lift/drag coefficients,
//! aerodynamic forces, boundary layer analysis, and wind field modeling.
//! Built on [`hisab`] for math.

pub mod error;
pub mod atmosphere;
pub mod airfoil;
pub mod coefficients;
pub mod forces;
pub mod boundary;
pub mod wind;
pub mod vehicle;

#[cfg(feature = "logging")]
pub mod logging;

pub use error::{PavanError, Result};
pub use atmosphere::{standard_temperature, standard_pressure, standard_density, dynamic_pressure, speed_of_sound, mach_number};
pub use coefficients::{lift_coefficient_thin_airfoil, drag_coefficient};
pub use forces::{lift, drag, reynolds_number, AeroForce};
