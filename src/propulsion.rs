//! Propulsion models for thrust, fuel consumption, and efficiency.
//!
//! Thrust-specific fuel consumption (TSFC), propeller efficiency,
//! and jet engine thrust lapse with altitude.

use serde::{Deserialize, Serialize};

use crate::atmosphere;
use crate::error::{PavanError, Result};

/// Propulsion system type with performance parameters.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[non_exhaustive]
pub enum PropulsionType {
    /// Jet engine with TSFC and sea-level static thrust.
    Jet {
        /// Sea-level static thrust (N).
        thrust_sl: f64,
        /// Thrust-specific fuel consumption (kg/(N·s)).
        tsfc: f64,
    },
    /// Propeller with efficiency model.
    Propeller {
        /// Shaft power (W).
        power: f64,
        /// Maximum propeller efficiency (typically 0.7-0.85).
        eta_max: f64,
        /// Propeller diameter (m).
        diameter: f64,
    },
}

/// Thrust available at a given flight condition.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[non_exhaustive]
pub struct ThrustResult {
    /// Available thrust (N).
    pub thrust: f64,
    /// Fuel flow rate (kg/s).
    pub fuel_flow: f64,
    /// Propulsive efficiency (0-1).
    pub efficiency: f64,
}

// ── Jet Engine ───────────────────────────────────────────────────────────

/// Jet engine thrust lapse with altitude (simple density ratio model).
///
/// T = T_sl × (ρ/ρ₀)^n, where n ≈ 0.7 for turbofan, 1.0 for turbojet.
#[must_use]
#[inline]
pub fn jet_thrust_at_altitude(thrust_sl: f64, altitude: f64, lapse_exponent: f64) -> f64 {
    let rho = atmosphere::standard_density(altitude);
    let sigma = rho / atmosphere::SEA_LEVEL_DENSITY;
    if sigma <= 0.0 {
        return 0.0;
    }
    thrust_sl * sigma.powf(lapse_exponent)
}

/// Fuel flow rate for a jet engine: ṁ_f = TSFC × T.
#[must_use]
#[inline]
pub fn jet_fuel_flow(tsfc: f64, thrust: f64) -> f64 {
    tsfc * thrust
}

/// Specific range for jet aircraft: SR = V / (TSFC × T) = V × (L/D) / (TSFC × W).
///
/// Returns distance per unit fuel (m/kg).
#[must_use]
#[inline]
pub fn jet_specific_range(velocity: f64, tsfc: f64, thrust: f64) -> f64 {
    let fuel_flow = jet_fuel_flow(tsfc, thrust);
    if fuel_flow <= 0.0 {
        return 0.0;
    }
    velocity / fuel_flow
}

// ── Propeller ────────────────────────────────────────────────────────────

/// Propeller advance ratio J = V / (n × D).
///
/// V in m/s, rps in rev/s, diameter in m.
#[must_use]
#[inline]
pub fn advance_ratio(velocity: f64, rps: f64, diameter: f64) -> f64 {
    let nd = rps * diameter;
    if nd.abs() < f64::EPSILON {
        return 0.0;
    }
    velocity / nd
}

/// Simple propeller efficiency model: η = η_max × (1 - e^(-k·J)).
///
/// Rises from 0 at static (J=0) to η_max at cruise. k controls the sharpness
/// of the transition (typically k ≈ 4-6).
#[must_use]
#[inline]
pub fn propeller_efficiency(advance_ratio_j: f64, eta_max: f64, k: f64) -> f64 {
    if advance_ratio_j <= 0.0 {
        return 0.0;
    }
    eta_max * (1.0 - (-k * advance_ratio_j).exp())
}

/// Propeller thrust from shaft power: T = η × P / V.
///
/// Returns 0 at static (V=0) since this model doesn't cover static thrust.
#[must_use]
#[inline]
pub fn propeller_thrust(power: f64, efficiency: f64, velocity: f64) -> f64 {
    if velocity <= 0.0 {
        return 0.0;
    }
    efficiency * power / velocity
}

/// Ideal (Froude) propeller efficiency: η_ideal = 2 / (1 + √(1 + T/(q·A_disc))).
///
/// Upper bound on propeller efficiency from momentum theory.
#[must_use]
#[inline]
pub fn froude_efficiency(thrust: f64, velocity: f64, disc_area: f64) -> f64 {
    if velocity <= 0.0 || disc_area <= 0.0 {
        return 0.0;
    }
    let q = 0.5 * velocity * velocity; // per unit density
    let tc = thrust / (q * disc_area);
    if tc < 0.0 {
        return 0.0;
    }
    2.0 / (1.0 + (1.0 + tc).sqrt())
}

// ── Combined ─────────────────────────────────────────────────────────────

/// Compute thrust and fuel flow for a propulsion system at given conditions.
pub fn compute_thrust(
    propulsion: &PropulsionType,
    velocity: f64,
    altitude: f64,
) -> Result<ThrustResult> {
    match *propulsion {
        PropulsionType::Jet { thrust_sl, tsfc } => {
            if thrust_sl <= 0.0 {
                return Err(PavanError::InvalidGeometry(
                    "sea-level thrust must be positive".into(),
                ));
            }
            // Turbofan lapse: n=0.7
            let thrust = jet_thrust_at_altitude(thrust_sl, altitude, 0.7);
            let fuel_flow = jet_fuel_flow(tsfc, thrust);
            Ok(ThrustResult {
                thrust,
                fuel_flow,
                efficiency: 1.0, // jet efficiency tracked separately
            })
        }
        PropulsionType::Propeller {
            power,
            eta_max,
            diameter,
        } => {
            if power <= 0.0 {
                return Err(PavanError::InvalidGeometry(
                    "shaft power must be positive".into(),
                ));
            }
            // Assume constant RPM, compute J and η
            // Typical cruise RPM for small prop: ~2500 RPM = 41.7 rps
            let rps = 41.7;
            let j = advance_ratio(velocity, rps, diameter);
            let eta = propeller_efficiency(j, eta_max, 5.0);
            let thrust = propeller_thrust(power, eta, velocity);
            // Simple fuel model: SFC for piston ≈ 8.5e-8 kg/(W·s)
            let sfc_piston = 8.5e-8;
            let fuel_flow = sfc_piston * power;
            Ok(ThrustResult {
                thrust,
                fuel_flow,
                efficiency: eta,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Jet ──

    // --- Guard clause coverage ---

    #[test]
    fn jet_thrust_extreme_altitude_zero() {
        let t = jet_thrust_at_altitude(50000.0, 100_000.0, 1.0);
        assert!(t >= 0.0);
    }

    #[test]
    fn advance_ratio_zero_diameter() {
        assert_eq!(advance_ratio(50.0, 40.0, 0.0), 0.0);
    }

    #[test]
    fn propeller_thrust_zero_eta() {
        assert_eq!(propeller_thrust(150_000.0, 0.0, 50.0), 0.0);
    }

    #[test]
    fn propeller_efficiency_negative_j() {
        assert_eq!(propeller_efficiency(-1.0, 0.85, 5.0), 0.0);
    }

    #[test]
    fn propeller_thrust_zero_velocity() {
        assert_eq!(propeller_thrust(150_000.0, 0.8, 0.0), 0.0);
    }

    #[test]
    fn froude_efficiency_zero_area() {
        assert_eq!(froude_efficiency(2400.0, 50.0, 0.0), 0.0);
    }

    #[test]
    fn jet_thrust_sea_level() {
        let t = jet_thrust_at_altitude(50000.0, 0.0, 0.7);
        assert!(
            (t - 50000.0).abs() < 1.0,
            "sea level thrust should equal T_sl"
        );
    }

    #[test]
    fn jet_thrust_decreases_with_altitude() {
        let t0 = jet_thrust_at_altitude(50000.0, 0.0, 0.7);
        let t5 = jet_thrust_at_altitude(50000.0, 5000.0, 0.7);
        let t10 = jet_thrust_at_altitude(50000.0, 10000.0, 0.7);
        assert!(t5 < t0);
        assert!(t10 < t5);
    }

    #[test]
    fn jet_fuel_flow_proportional() {
        let ff1 = jet_fuel_flow(2e-5, 40000.0);
        let ff2 = jet_fuel_flow(2e-5, 20000.0);
        assert!((ff1 / ff2 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn jet_specific_range_basic() {
        // V=250m/s, TSFC=2e-5, T=40000N → SR = 250/(2e-5*40000) = 250/0.8 = 312.5 m/kg
        let sr = jet_specific_range(250.0, 2e-5, 40000.0);
        assert!((sr - 312.5).abs() < 0.1);
    }

    #[test]
    fn jet_specific_range_zero_thrust() {
        assert_eq!(jet_specific_range(250.0, 2e-5, 0.0), 0.0);
    }

    #[test]
    fn turbojet_lapse_steeper() {
        let t_fan = jet_thrust_at_altitude(50000.0, 10000.0, 0.7);
        let t_jet = jet_thrust_at_altitude(50000.0, 10000.0, 1.0);
        assert!(
            t_jet < t_fan,
            "turbojet (n=1) should lose more thrust than turbofan (n=0.7)"
        );
    }

    // ── Propeller ──

    #[test]
    fn advance_ratio_basic() {
        // V=50m/s, n=40rps, D=2m → J = 50/(40*2) = 0.625
        let j = advance_ratio(50.0, 40.0, 2.0);
        assert!((j - 0.625).abs() < 1e-10);
    }

    #[test]
    fn advance_ratio_static() {
        assert_eq!(advance_ratio(0.0, 40.0, 2.0), 0.0);
    }

    #[test]
    fn propeller_efficiency_increases_with_j() {
        let e1 = propeller_efficiency(0.2, 0.85, 5.0);
        let e2 = propeller_efficiency(0.5, 0.85, 5.0);
        let e3 = propeller_efficiency(1.0, 0.85, 5.0);
        assert!(e2 > e1);
        assert!(e3 > e2);
    }

    #[test]
    fn propeller_efficiency_static_zero() {
        assert_eq!(propeller_efficiency(0.0, 0.85, 5.0), 0.0);
    }

    #[test]
    fn propeller_efficiency_max_limit() {
        let e = propeller_efficiency(10.0, 0.85, 5.0);
        assert!(
            (e - 0.85).abs() < 0.01,
            "at high J, η should approach η_max"
        );
    }

    #[test]
    fn propeller_thrust_basic() {
        // P=150kW, η=0.8, V=50m/s → T = 0.8*150000/50 = 2400N
        let t = propeller_thrust(150_000.0, 0.8, 50.0);
        assert!((t - 2400.0).abs() < 0.1);
    }

    #[test]
    fn propeller_thrust_static_zero() {
        assert_eq!(propeller_thrust(150_000.0, 0.8, 0.0), 0.0);
    }

    #[test]
    fn froude_efficiency_basic() {
        // Upper bound on efficiency
        let eta = froude_efficiency(2400.0, 50.0, std::f64::consts::PI);
        assert!(
            eta > 0.5 && eta < 1.0,
            "Froude efficiency should be 0.5-1.0, got {eta}"
        );
    }

    #[test]
    fn froude_efficiency_static_zero() {
        assert_eq!(froude_efficiency(2400.0, 0.0, std::f64::consts::PI), 0.0);
    }

    // ── Combined ──

    #[test]
    fn compute_thrust_jet() {
        let jet = PropulsionType::Jet {
            thrust_sl: 50000.0,
            tsfc: 2e-5,
        };
        let result = compute_thrust(&jet, 250.0, 10000.0).expect("compute");
        assert!(result.thrust > 0.0 && result.thrust < 50000.0);
        assert!(result.fuel_flow > 0.0);
    }

    #[test]
    fn compute_thrust_propeller() {
        let prop = PropulsionType::Propeller {
            power: 150_000.0,
            eta_max: 0.85,
            diameter: 2.0,
        };
        let result = compute_thrust(&prop, 50.0, 0.0).expect("compute");
        assert!(result.thrust > 0.0);
        assert!(result.efficiency > 0.0 && result.efficiency <= 0.85);
    }

    #[test]
    fn compute_thrust_jet_zero_thrust_errors() {
        let jet = PropulsionType::Jet {
            thrust_sl: 0.0,
            tsfc: 2e-5,
        };
        assert!(compute_thrust(&jet, 250.0, 0.0).is_err());
    }

    #[test]
    fn compute_thrust_prop_zero_power_errors() {
        let prop = PropulsionType::Propeller {
            power: 0.0,
            eta_max: 0.85,
            diameter: 2.0,
        };
        assert!(compute_thrust(&prop, 50.0, 0.0).is_err());
    }

    // ── Serde ──

    #[test]
    fn serde_propulsion_type() {
        let jet = PropulsionType::Jet {
            thrust_sl: 50000.0,
            tsfc: 2e-5,
        };
        let json = serde_json::to_string(&jet).expect("serialize");
        let back: PropulsionType = serde_json::from_str(&json).expect("deserialize");
        match back {
            PropulsionType::Jet { thrust_sl, .. } => {
                assert!((thrust_sl - 50000.0).abs() < f64::EPSILON)
            }
            _ => panic!("wrong variant"),
        }
    }

    #[test]
    fn serde_thrust_result() {
        let jet = PropulsionType::Jet {
            thrust_sl: 50000.0,
            tsfc: 2e-5,
        };
        let result = compute_thrust(&jet, 250.0, 0.0).expect("compute");
        let json = serde_json::to_string(&result).expect("serialize");
        let back: ThrustResult = serde_json::from_str(&json).expect("deserialize");
        assert!((back.thrust - result.thrust).abs() < f64::EPSILON);
    }
}
