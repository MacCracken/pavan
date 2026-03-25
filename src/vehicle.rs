use crate::atmosphere;
use crate::coefficients;
use crate::forces::AeroForce;
use serde::{Deserialize, Serialize};

/// An aerodynamic body with reference properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct AeroBody {
    /// Reference area (m²) — typically wing area for aircraft.
    pub reference_area: f64,
    /// Zero-lift drag coefficient.
    pub cd0: f64,
    /// Wing aspect ratio (span² / area).
    pub aspect_ratio: f64,
    /// Oswald efficiency factor (0.7–0.85 typical).
    pub oswald_efficiency: f64,
    /// Chord length (m) — for moment computation.
    pub chord: f64,
}

impl AeroBody {
    /// Typical light aircraft (Cessna 172 class).
    #[must_use]
    #[inline]
    pub fn light_aircraft() -> Self {
        Self {
            reference_area: 16.2,
            cd0: 0.027,
            aspect_ratio: 7.5,
            oswald_efficiency: 0.8,
            chord: 1.5,
        }
    }

    /// High-performance glider.
    #[must_use]
    #[inline]
    pub fn glider() -> Self {
        Self {
            reference_area: 9.5,
            cd0: 0.01,
            aspect_ratio: 23.0,
            oswald_efficiency: 0.9,
            chord: 0.8,
        }
    }

    /// Compute aerodynamic forces at given conditions.
    #[must_use]
    #[inline]
    pub fn compute_forces(
        &self,
        velocity: f64,
        altitude: f64,
        angle_of_attack_rad: f64,
    ) -> AeroForce {
        let rho = atmosphere::standard_density(altitude);
        let q = atmosphere::dynamic_pressure(rho, velocity);
        let cl = coefficients::lift_coefficient_thin_airfoil(angle_of_attack_rad);
        let cd =
            coefficients::drag_coefficient(self.cd0, cl, self.aspect_ratio, self.oswald_efficiency);

        AeroForce {
            lift: q * self.reference_area * cl,
            drag: q * self.reference_area * cd,
            moment: 0.0, // simplified: no pitch moment
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn light_aircraft_generates_lift() {
        let body = AeroBody::light_aircraft();
        let forces = body.compute_forces(60.0, 0.0, 5.0_f64.to_radians());
        assert!(forces.lift > 0.0, "should generate lift at 5° AoA");
        assert!(forces.drag > 0.0, "should have drag");
    }

    #[test]
    fn zero_aoa_no_lift() {
        let body = AeroBody::light_aircraft();
        let forces = body.compute_forces(60.0, 0.0, 0.0);
        assert!(
            forces.lift.abs() < 1.0,
            "zero AoA should produce near-zero lift"
        );
    }

    #[test]
    fn glider_better_ld_ratio() {
        let aircraft = AeroBody::light_aircraft();
        let glider = AeroBody::glider();
        let alpha = 3.0_f64.to_radians();

        let f_aircraft = aircraft.compute_forces(40.0, 0.0, alpha);
        let f_glider = glider.compute_forces(40.0, 0.0, alpha);

        let ld_aircraft = f_aircraft.lift / f_aircraft.drag;
        let ld_glider = f_glider.lift / f_glider.drag;
        assert!(
            ld_glider > ld_aircraft,
            "glider should have better L/D ratio"
        );
    }

    #[test]
    fn forces_decrease_with_altitude() {
        let body = AeroBody::light_aircraft();
        let alpha = 5.0_f64.to_radians();
        let f_low = body.compute_forces(60.0, 0.0, alpha);
        let f_high = body.compute_forces(60.0, 5000.0, alpha);
        assert!(
            f_high.lift < f_low.lift,
            "lift should decrease with altitude (lower density)"
        );
    }

    // --- Edge cases ---

    #[test]
    fn moment_is_zero_simplified() {
        let body = AeroBody::light_aircraft();
        let f = body.compute_forces(60.0, 0.0, 5.0_f64.to_radians());
        assert_eq!(
            f.moment, 0.0,
            "vehicle compute_forces should produce zero moment (simplified)"
        );
    }

    #[test]
    fn zero_velocity_zero_forces() {
        let body = AeroBody::light_aircraft();
        let f = body.compute_forces(0.0, 0.0, 5.0_f64.to_radians());
        assert_eq!(f.lift, 0.0);
        assert_eq!(f.drag, 0.0);
    }

    #[test]
    fn glider_preset_values() {
        let g = AeroBody::glider();
        assert!((g.reference_area - 9.5).abs() < f64::EPSILON);
        assert!((g.cd0 - 0.01).abs() < f64::EPSILON);
        assert!((g.aspect_ratio - 23.0).abs() < f64::EPSILON);
        assert!((g.oswald_efficiency - 0.9).abs() < f64::EPSILON);
        assert!((g.chord - 0.8).abs() < f64::EPSILON);
    }

    #[test]
    fn light_aircraft_preset_values() {
        let a = AeroBody::light_aircraft();
        assert!((a.reference_area - 16.2).abs() < f64::EPSILON);
        assert!((a.cd0 - 0.027).abs() < f64::EPSILON);
        assert!((a.aspect_ratio - 7.5).abs() < f64::EPSILON);
    }

    #[test]
    fn drag_always_positive_at_nonzero_velocity() {
        let body = AeroBody::light_aircraft();
        // Even at zero AoA, parasitic drag exists
        let f = body.compute_forces(60.0, 0.0, 0.0);
        assert!(f.drag > 0.0, "parasitic drag should exist even at zero AoA");
    }
}
