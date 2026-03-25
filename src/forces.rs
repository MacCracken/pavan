use serde::{Deserialize, Serialize};

/// Aerodynamic forces acting on a body.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AeroForce {
    pub lift: f64,
    pub drag: f64,
    pub moment: f64,
}

/// Lift force: L = q × S × Cl
#[must_use]
#[inline]
pub fn lift(dynamic_pressure: f64, reference_area: f64, cl: f64) -> f64 {
    dynamic_pressure * reference_area * cl
}

/// Drag force: D = q × S × Cd
#[must_use]
#[inline]
pub fn drag(dynamic_pressure: f64, reference_area: f64, cd: f64) -> f64 {
    dynamic_pressure * reference_area * cd
}

/// Reynolds number: Re = ρVL / μ
///
/// density in kg/m³, velocity in m/s, length in m, viscosity in Pa·s.
#[must_use]
#[inline]
pub fn reynolds_number(density: f64, velocity: f64, length: f64, dynamic_viscosity: f64) -> f64 {
    if dynamic_viscosity <= 0.0 { return 0.0; }
    density * velocity * length / dynamic_viscosity
}

/// Dynamic viscosity of air at temperature T (Sutherland's law).
///
/// μ = μ_ref × (T/T_ref)^(3/2) × (T_ref + S) / (T + S)
///
/// Where μ_ref = 1.716e-5 Pa·s, T_ref = 273.15 K, S = 110.4 K.
#[must_use]
#[inline]
pub fn air_dynamic_viscosity(temperature_k: f64) -> f64 {
    let mu_ref = 1.716e-5;
    let t_ref = 273.15;
    let s = 110.4;
    if temperature_k <= 0.0 { return mu_ref; }
    mu_ref * (temperature_k / t_ref).powf(1.5) * (t_ref + s) / (temperature_k + s)
}

/// Compute full aerodynamic force from flight conditions.
#[must_use]
pub fn compute_aero_force(
    density: f64,
    velocity: f64,
    reference_area: f64,
    cl: f64,
    cd: f64,
    cm: f64,
    chord: f64,
) -> AeroForce {
    let q = crate::atmosphere::dynamic_pressure(density, velocity);
    AeroForce {
        lift: q * reference_area * cl,
        drag: q * reference_area * cd,
        moment: q * reference_area * chord * cm,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lift_basic() {
        // q=6125, S=10, Cl=0.5 → L=30625 N
        let l = lift(6125.0, 10.0, 0.5);
        assert!((l - 30625.0).abs() < 1.0, "lift should be ~30625 N, got {l}");
    }

    #[test]
    fn drag_basic() {
        let d = drag(6125.0, 10.0, 0.05);
        assert!((d - 3062.5).abs() < 1.0);
    }

    #[test]
    fn zero_velocity_zero_forces() {
        let q = crate::atmosphere::dynamic_pressure(1.225, 0.0);
        assert_eq!(lift(q, 10.0, 0.5), 0.0);
        assert_eq!(drag(q, 10.0, 0.05), 0.0);
    }

    #[test]
    fn reynolds_at_sea_level() {
        // 1.225 kg/m³, 50 m/s, 1m chord, μ≈1.8e-5 → Re ≈ 3.4M
        let re = reynolds_number(1.225, 50.0, 1.0, 1.8e-5);
        assert!(re > 3e6 && re < 4e6, "Re should be ~3.4M, got {re}");
    }

    #[test]
    fn reynolds_increases_with_velocity() {
        let re_slow = reynolds_number(1.225, 10.0, 1.0, 1.8e-5);
        let re_fast = reynolds_number(1.225, 100.0, 1.0, 1.8e-5);
        assert!(re_fast > re_slow);
    }

    #[test]
    fn air_viscosity_at_sea_level() {
        let mu = air_dynamic_viscosity(288.15);
        assert!((mu - 1.79e-5).abs() < 0.1e-5, "viscosity at sea level should be ~1.79e-5, got {mu}");
    }

    #[test]
    fn viscosity_increases_with_temperature() {
        let mu_cold = air_dynamic_viscosity(250.0);
        let mu_hot = air_dynamic_viscosity(350.0);
        assert!(mu_hot > mu_cold, "viscosity should increase with temperature");
    }

    #[test]
    fn aero_force_computation() {
        let f = compute_aero_force(1.225, 100.0, 10.0, 0.5, 0.05, -0.1, 1.5);
        assert!(f.lift > 0.0);
        assert!(f.drag > 0.0);
        assert!(f.moment < 0.0); // negative Cm → nose-down moment
    }
}
