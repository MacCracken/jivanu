//! Pharmacokinetics — drug absorption, distribution, and elimination.
//!
//! Models for predicting drug concentration over time in biological
//! compartments. Used alongside the [`resistance`](crate::resistance) module
//! to connect antibiotic dosing to bactericidal effect.

use serde::{Deserialize, Serialize};

use crate::error::{JivanuError, Result, validate_non_negative, validate_positive};

/// Elimination half-life from an elimination rate constant.
///
/// `t½ = ln(2) / k_e`
///
/// # Errors
///
/// Returns error if `k_e` is non-positive.
#[inline]
#[must_use = "returns the half-life without side effects"]
pub fn half_life(k_e: f64) -> Result<f64> {
    validate_positive(k_e, "k_e")?;
    Ok(core::f64::consts::LN_2 / k_e)
}

/// Elimination rate constant from a half-life.
///
/// `k_e = ln(2) / t½`
///
/// # Errors
///
/// Returns error if `t_half` is non-positive.
#[inline]
#[must_use = "returns the elimination rate constant without side effects"]
pub fn elimination_rate(t_half: f64) -> Result<f64> {
    validate_positive(t_half, "t_half")?;
    Ok(core::f64::consts::LN_2 / t_half)
}

/// One-compartment IV bolus: plasma concentration at time `t`.
///
/// `C(t) = (dose / V_d) × e^(-k_e × t)`
///
/// Models instantaneous drug distribution after intravenous injection.
///
/// # Arguments
///
/// - `dose` — administered dose (mg)
/// - `v_d` — volume of distribution (L)
/// - `k_e` — elimination rate constant (1/hr)
/// - `t` — time after administration (hr)
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the plasma concentration without side effects"]
pub fn iv_bolus_concentration(dose: f64, v_d: f64, k_e: f64, t: f64) -> Result<f64> {
    validate_positive(dose, "dose")?;
    validate_positive(v_d, "v_d")?;
    validate_positive(k_e, "k_e")?;
    validate_non_negative(t, "t")?;
    let c0 = dose / v_d;
    Ok(c0 * (-k_e * t).exp())
}

/// One-compartment oral administration: plasma concentration at time `t`.
///
/// `C(t) = (F × dose × k_a) / (V_d × (k_a - k_e)) × (e^(-k_e×t) - e^(-k_a×t))`
///
/// Bateman equation for first-order absorption and elimination.
///
/// # Arguments
///
/// - `dose` — administered dose (mg)
/// - `bioavailability` — fraction absorbed (F, 0–1)
/// - `v_d` — volume of distribution (L)
/// - `k_a` — absorption rate constant (1/hr)
/// - `k_e` — elimination rate constant (1/hr)
/// - `t` — time after administration (hr)
///
/// # Errors
///
/// Returns error if parameters are invalid or `k_a == k_e`.
#[inline]
#[must_use = "returns the plasma concentration without side effects"]
pub fn oral_concentration(
    dose: f64,
    bioavailability: f64,
    v_d: f64,
    k_a: f64,
    k_e: f64,
    t: f64,
) -> Result<f64> {
    validate_positive(dose, "dose")?;
    validate_positive(bioavailability, "bioavailability")?;
    if bioavailability > 1.0 {
        return Err(JivanuError::ComputationError(
            "bioavailability must be in (0, 1]".into(),
        ));
    }
    validate_positive(v_d, "v_d")?;
    validate_positive(k_a, "k_a")?;
    validate_positive(k_e, "k_e")?;
    validate_non_negative(t, "t")?;

    if (k_a - k_e).abs() < 1e-15 {
        return Err(JivanuError::ComputationError(
            "k_a must differ from k_e (flip-flop kinetics not supported)".into(),
        ));
    }

    let numerator = bioavailability * dose * k_a;
    let denominator = v_d * (k_a - k_e);
    Ok(numerator / denominator * ((-k_e * t).exp() - (-k_a * t).exp()))
}

/// Time to maximum concentration (Tmax) for oral one-compartment model.
///
/// `T_max = ln(k_a / k_e) / (k_a - k_e)`
///
/// # Errors
///
/// Returns error if rate constants are invalid or equal.
#[inline]
#[must_use = "returns the Tmax without side effects"]
pub fn oral_tmax(k_a: f64, k_e: f64) -> Result<f64> {
    validate_positive(k_a, "k_a")?;
    validate_positive(k_e, "k_e")?;
    if (k_a - k_e).abs() < 1e-15 {
        return Err(JivanuError::ComputationError(
            "k_a must differ from k_e".into(),
        ));
    }
    Ok((k_a / k_e).ln() / (k_a - k_e))
}

/// Maximum concentration (Cmax) for oral one-compartment model.
///
/// Evaluates `oral_concentration` at `T_max`.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the Cmax without side effects"]
pub fn oral_cmax(dose: f64, bioavailability: f64, v_d: f64, k_a: f64, k_e: f64) -> Result<f64> {
    let tmax = oral_tmax(k_a, k_e)?;
    oral_concentration(dose, bioavailability, v_d, k_a, k_e, tmax)
}

/// Area under the concentration-time curve (AUC) via the trapezoidal rule.
///
/// Computes the definite integral of concentration over time from paired
/// `(time, concentration)` data points. Data must be sorted by time.
///
/// AUC is the primary measure of total drug exposure.
///
/// # Errors
///
/// Returns error if fewer than 2 data points or times are not increasing.
#[must_use = "returns the AUC without side effects"]
pub fn auc_trapezoidal(times: &[f64], concentrations: &[f64]) -> Result<f64> {
    if times.len() != concentrations.len() {
        return Err(JivanuError::ComputationError(
            "times and concentrations must have the same length".into(),
        ));
    }
    if times.len() < 2 {
        return Err(JivanuError::ComputationError(
            "at least 2 data points required for AUC".into(),
        ));
    }
    let mut auc = 0.0;
    for i in 1..times.len() {
        let dt = times[i] - times[i - 1];
        if dt < 0.0 {
            return Err(JivanuError::ComputationError(
                "times must be monotonically increasing".into(),
            ));
        }
        auc += 0.5 * (concentrations[i - 1] + concentrations[i]) * dt;
    }
    Ok(auc)
}

/// Analytical AUC (0→∞) for one-compartment IV bolus.
///
/// `AUC = dose / (V_d × k_e) = C_0 / k_e`
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the AUC without side effects"]
pub fn auc_iv_bolus(dose: f64, v_d: f64, k_e: f64) -> Result<f64> {
    validate_positive(dose, "dose")?;
    validate_positive(v_d, "v_d")?;
    validate_positive(k_e, "k_e")?;
    Ok(dose / (v_d * k_e))
}

/// Two-compartment IV bolus model state.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TwoCompartmentState {
    /// Drug amount in the central compartment (mg).
    pub central: f64,
    /// Drug amount in the peripheral compartment (mg).
    pub peripheral: f64,
}

/// Parameters for the two-compartment IV bolus model.
///
/// ```text
/// dA_c/dt = -(k_e + k_12) × A_c + k_21 × A_p
/// dA_p/dt =  k_12 × A_c - k_21 × A_p
/// ```
///
/// Plasma concentration is `A_c / V_c`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TwoCompartmentParams {
    /// Volume of the central compartment (L).
    pub v_c: f64,
    /// Elimination rate constant from central compartment (1/hr).
    pub k_e: f64,
    /// Transfer rate central → peripheral (1/hr).
    pub k_12: f64,
    /// Transfer rate peripheral → central (1/hr).
    pub k_21: f64,
}

/// One step of the two-compartment IV bolus model.
///
/// When the `hisab` feature is enabled, uses RK4 integration.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the new state without side effects"]
pub fn two_compartment_step(
    state: &TwoCompartmentState,
    params: &TwoCompartmentParams,
    dt: f64,
) -> Result<TwoCompartmentState> {
    validate_non_negative(state.central, "central")?;
    validate_non_negative(state.peripheral, "peripheral")?;
    validate_positive(params.v_c, "v_c")?;
    validate_positive(params.k_e, "k_e")?;
    validate_non_negative(params.k_12, "k_12")?;
    validate_non_negative(params.k_21, "k_21")?;
    validate_positive(dt, "dt")?;

    two_compartment_step_inner(state, params, dt)
}

#[cfg(feature = "hisab")]
fn two_compartment_step_inner(
    state: &TwoCompartmentState,
    params: &TwoCompartmentParams,
    dt: f64,
) -> Result<TwoCompartmentState> {
    let TwoCompartmentParams {
        k_e, k_12, k_21, ..
    } = *params;
    let y0 = [state.central, state.peripheral];
    let result = hisab::num::rk4(
        |_t, y, dy| {
            dy[0] = -(k_e + k_12) * y[0] + k_21 * y[1];
            dy[1] = k_12 * y[0] - k_21 * y[1];
        },
        0.0,
        &y0,
        dt,
        1,
    )
    .map_err(|e| JivanuError::SimulationFailed(e.to_string()))?;
    Ok(TwoCompartmentState {
        central: result[0].max(0.0),
        peripheral: result[1].max(0.0),
    })
}

#[cfg(not(feature = "hisab"))]
fn two_compartment_step_inner(
    state: &TwoCompartmentState,
    params: &TwoCompartmentParams,
    dt: f64,
) -> Result<TwoCompartmentState> {
    let ac = state.central;
    let ap = state.peripheral;
    let dac = -(params.k_e + params.k_12) * ac + params.k_21 * ap;
    let dap = params.k_12 * ac - params.k_21 * ap;
    Ok(TwoCompartmentState {
        central: (ac + dac * dt).max(0.0),
        peripheral: (ap + dap * dt).max(0.0),
    })
}

/// Plasma concentration from a two-compartment state.
#[inline]
#[must_use]
pub fn plasma_concentration(state: &TwoCompartmentState, v_c: f64) -> f64 {
    if v_c <= 0.0 {
        return 0.0;
    }
    state.central / v_c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_half_life() {
        let t = half_life(core::f64::consts::LN_2).unwrap();
        assert!((t - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_elimination_rate() {
        let k = elimination_rate(1.0).unwrap();
        assert!((k - core::f64::consts::LN_2).abs() < 1e-10);
    }

    #[test]
    fn test_half_life_elimination_rate_roundtrip() {
        let k = 0.35;
        let t = half_life(k).unwrap();
        let k_back = elimination_rate(t).unwrap();
        assert!((k - k_back).abs() < 1e-10);
    }

    #[test]
    fn test_iv_bolus_at_t0() {
        let c = iv_bolus_concentration(500.0, 50.0, 0.1, 0.0).unwrap();
        assert!((c - 10.0).abs() < 1e-10); // C0 = dose/Vd = 10
    }

    #[test]
    fn test_iv_bolus_decay() {
        let c0 = iv_bolus_concentration(500.0, 50.0, 0.1, 0.0).unwrap();
        let c1 = iv_bolus_concentration(500.0, 50.0, 0.1, 1.0).unwrap();
        assert!(c1 < c0);
        assert!(c1 > 0.0);
    }

    #[test]
    fn test_iv_bolus_at_half_life() {
        // At t = t½, concentration should be C0/2
        let k_e = 0.1;
        let t_half = half_life(k_e).unwrap();
        let c0 = iv_bolus_concentration(500.0, 50.0, k_e, 0.0).unwrap();
        let c_half = iv_bolus_concentration(500.0, 50.0, k_e, t_half).unwrap();
        assert!((c_half - c0 / 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_oral_concentration_at_t0() {
        let c = oral_concentration(500.0, 1.0, 50.0, 1.0, 0.1, 0.0).unwrap();
        assert!((c - 0.0).abs() < 1e-10); // No absorption yet
    }

    #[test]
    fn test_oral_concentration_rises_then_falls() {
        let c1 = oral_concentration(500.0, 1.0, 50.0, 1.0, 0.1, 1.0).unwrap();
        let c_peak = oral_concentration(500.0, 1.0, 50.0, 1.0, 0.1, 3.0).unwrap();
        let c_late = oral_concentration(500.0, 1.0, 50.0, 1.0, 0.1, 20.0).unwrap();
        assert!(c1 > 0.0);
        assert!(c_late < c_peak);
    }

    #[test]
    fn test_oral_tmax() {
        let tmax = oral_tmax(1.0, 0.1).unwrap();
        assert!(tmax > 0.0);
        // Tmax = ln(1.0/0.1) / (1.0 - 0.1) = ln(10) / 0.9 ≈ 2.558
        assert!((tmax - 10.0_f64.ln() / 0.9).abs() < 1e-10);
    }

    #[test]
    fn test_oral_tmax_equal_rates_error() {
        assert!(oral_tmax(0.5, 0.5).is_err());
    }

    #[test]
    fn test_oral_cmax() {
        let cmax = oral_cmax(500.0, 1.0, 50.0, 1.0, 0.1).unwrap();
        assert!(cmax > 0.0);
        // Cmax should be less than C0 of equivalent IV bolus
        let c0_iv = 500.0 / 50.0;
        assert!(cmax < c0_iv);
    }

    #[test]
    fn test_auc_trapezoidal_rectangle() {
        // Constant concentration of 10 from t=0 to t=5
        let auc = auc_trapezoidal(&[0.0, 5.0], &[10.0, 10.0]).unwrap();
        assert!((auc - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_auc_trapezoidal_triangle() {
        // Linear decay from 10 to 0 over t=0..10
        let auc = auc_trapezoidal(&[0.0, 10.0], &[10.0, 0.0]).unwrap();
        assert!((auc - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_auc_trapezoidal_multipoint() {
        let times = [0.0, 1.0, 2.0, 3.0];
        let concs = [10.0, 8.0, 4.0, 1.0];
        let auc = auc_trapezoidal(&times, &concs).unwrap();
        // Trapezoids: (10+8)/2*1 + (8+4)/2*1 + (4+1)/2*1 = 9 + 6 + 2.5 = 17.5
        assert!((auc - 17.5).abs() < 1e-10);
    }

    #[test]
    fn test_auc_trapezoidal_too_few_points() {
        assert!(auc_trapezoidal(&[0.0], &[10.0]).is_err());
    }

    #[test]
    fn test_auc_trapezoidal_non_monotonic() {
        assert!(auc_trapezoidal(&[0.0, 2.0, 1.0], &[10.0, 5.0, 3.0]).is_err());
    }

    #[test]
    fn test_auc_iv_bolus_analytical() {
        // AUC = dose / (Vd * ke) = 500 / (50 * 0.1) = 100
        let auc = auc_iv_bolus(500.0, 50.0, 0.1).unwrap();
        assert!((auc - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_two_compartment_step_elimination() {
        let state = TwoCompartmentState {
            central: 500.0,
            peripheral: 0.0,
        };
        let params = TwoCompartmentParams {
            v_c: 50.0,
            k_e: 0.1,
            k_12: 0.2,
            k_21: 0.1,
        };
        let next = two_compartment_step(&state, &params, 0.1).unwrap();
        // Central should decrease (elimination + distribution)
        assert!(next.central < 500.0);
        // Peripheral should increase (distribution in)
        assert!(next.peripheral > 0.0);
    }

    #[test]
    fn test_two_compartment_mass_decreases() {
        // Total drug mass should decrease due to elimination
        let state = TwoCompartmentState {
            central: 500.0,
            peripheral: 100.0,
        };
        let params = TwoCompartmentParams {
            v_c: 50.0,
            k_e: 0.1,
            k_12: 0.2,
            k_21: 0.1,
        };
        let next = two_compartment_step(&state, &params, 0.1).unwrap();
        let total_before = state.central + state.peripheral;
        let total_after = next.central + next.peripheral;
        assert!(total_after < total_before);
    }

    #[test]
    fn test_two_compartment_no_transfer() {
        // k_12 = k_21 = 0 → behaves like one-compartment
        let state = TwoCompartmentState {
            central: 500.0,
            peripheral: 0.0,
        };
        let params = TwoCompartmentParams {
            v_c: 50.0,
            k_e: 0.1,
            k_12: 0.0,
            k_21: 0.0,
        };
        let next = two_compartment_step(&state, &params, 0.1).unwrap();
        assert!(next.central < 500.0);
        assert!((next.peripheral - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_plasma_concentration() {
        let state = TwoCompartmentState {
            central: 500.0,
            peripheral: 100.0,
        };
        let c = plasma_concentration(&state, 50.0);
        assert!((c - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_two_compartment_state_serde_roundtrip() {
        let state = TwoCompartmentState {
            central: 500.0,
            peripheral: 100.0,
        };
        let json = serde_json::to_string(&state).unwrap();
        let back: TwoCompartmentState = serde_json::from_str(&json).unwrap();
        assert!((state.central - back.central).abs() < 1e-10);
    }

    #[test]
    fn test_two_compartment_params_serde_roundtrip() {
        let params = TwoCompartmentParams {
            v_c: 50.0,
            k_e: 0.1,
            k_12: 0.2,
            k_21: 0.1,
        };
        let json = serde_json::to_string(&params).unwrap();
        let back: TwoCompartmentParams = serde_json::from_str(&json).unwrap();
        assert!((params.k_e - back.k_e).abs() < 1e-10);
    }
}
