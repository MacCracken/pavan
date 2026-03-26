//! # Pavan
//!
//! **Pavan** (पवन — Sanskrit for "wind, breeze") — aerodynamics engine for the AGNOS ecosystem.
//!
//! Provides ISA atmosphere model, NACA airfoil generation, lift/drag coefficients,
//! aerodynamic forces, boundary layer analysis, and wind field modeling.
//! Built on [`hisab`] for math.

pub mod airfoil;
pub mod atmosphere;
pub mod boundary;
pub mod coefficients;
pub mod compressible;
pub mod error;
pub mod forces;
pub mod panel;
pub mod propulsion;
pub mod stability;
pub mod vehicle;
pub mod vlm;
pub mod wind;

#[cfg(feature = "cfd")]
pub mod cfd;

#[cfg(feature = "logging")]
pub mod logging;

pub use atmosphere::{
    dynamic_pressure, mach_number, speed_of_sound, standard_density, standard_pressure,
    standard_temperature,
};
pub use coefficients::{drag_coefficient, lift_coefficient_thin_airfoil};
pub use error::{PavanError, Result};
pub use forces::{AeroForce, body_to_wind, drag, lift, reynolds_number};

pub use airfoil::{NacaProfile, SurfacePoints};
pub use panel::{Panel, PanelSolution};
pub use vehicle::AeroBody;
pub use vlm::{VlmSolution, WingGeometry};
pub use wind::WindField;

#[cfg(feature = "cfd")]
pub use cfd::{AirfoilCfd, AirfoilCfdConfig, CfdSnapshot};
