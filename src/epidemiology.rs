//! Epidemiology — SIR/SEIR models, R0, herd immunity.
//!
//! When the `hisab` feature is enabled, compartmental model stepping uses
//! fourth-order Runge-Kutta (RK4) integration via [`hisab::num::rk4`] for
//! higher accuracy and Neumaier-compensated accumulation. Without the
//! feature, a forward-Euler integrator is used as a fallback.

use serde::{Deserialize, Serialize};

use crate::error::{Result, validate_non_negative, validate_positive};

/// SIR compartmental model state.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SirState {
    /// Susceptible fraction.
    pub s: f64,
    /// Infected fraction.
    pub i: f64,
    /// Recovered fraction.
    pub r: f64,
}

/// One step of the SIR model.
///
/// When the `hisab` feature is enabled, uses fourth-order Runge-Kutta (RK4)
/// integration for higher accuracy. Otherwise falls back to forward Euler.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the new SIR state without side effects"]
pub fn sir_step(s: f64, i: f64, r: f64, beta: f64, gamma: f64, dt: f64) -> Result<SirState> {
    validate_non_negative(s, "s")?;
    validate_non_negative(i, "i")?;
    validate_non_negative(r, "r")?;
    validate_positive(beta, "beta")?;
    validate_positive(gamma, "gamma")?;
    validate_positive(dt, "dt")?;

    sir_step_inner(s, i, r, beta, gamma, dt)
}

#[cfg(feature = "hisab")]
fn sir_step_inner(s: f64, i: f64, r: f64, beta: f64, gamma: f64, dt: f64) -> Result<SirState> {
    let y0 = [s, i, r];
    let result = hisab::num::rk4(
        |_t, y, dy| {
            dy[0] = -beta * y[0] * y[1];
            dy[1] = beta * y[0] * y[1] - gamma * y[1];
            dy[2] = gamma * y[1];
        },
        0.0,
        &y0,
        dt,
        1,
    )
    .map_err(|e| crate::error::JivanuError::SimulationFailed(e.to_string()))?;
    Ok(SirState {
        s: result[0].max(0.0),
        i: result[1].max(0.0),
        r: result[2].max(0.0),
    })
}

#[cfg(not(feature = "hisab"))]
fn sir_step_inner(s: f64, i: f64, r: f64, beta: f64, gamma: f64, dt: f64) -> Result<SirState> {
    let ds = -beta * s * i;
    let di = beta * s * i - gamma * i;
    let dr = gamma * i;
    Ok(SirState {
        s: (s + ds * dt).max(0.0),
        i: (i + di * dt).max(0.0),
        r: (r + dr * dt).max(0.0),
    })
}

/// SEIR model state (with exposed compartment).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SeirState {
    /// Susceptible fraction.
    pub s: f64,
    /// Exposed (infected but not yet infectious) fraction.
    pub e: f64,
    /// Infected (and infectious) fraction.
    pub i: f64,
    /// Recovered fraction.
    pub r: f64,
}

/// Parameters for the SEIR model.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SeirParams {
    /// Transmission rate.
    pub beta: f64,
    /// Rate at which exposed become infectious (1/incubation period).
    pub sigma: f64,
    /// Recovery rate.
    pub gamma: f64,
    /// Time step.
    pub dt: f64,
}

/// One step of the SEIR model.
///
/// `sigma` is the rate at which exposed become infectious (1/incubation period).
///
/// When the `hisab` feature is enabled, uses fourth-order Runge-Kutta (RK4)
/// integration for higher accuracy. Otherwise falls back to forward Euler.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the new SEIR state without side effects"]
pub fn seir_step(state: &SeirState, params: &SeirParams) -> Result<SeirState> {
    let SeirState { s, e, i, r } = *state;
    let SeirParams {
        beta,
        sigma,
        gamma,
        dt,
    } = *params;
    validate_non_negative(s, "s")?;
    validate_non_negative(e, "e")?;
    validate_non_negative(i, "i")?;
    validate_non_negative(r, "r")?;
    validate_positive(beta, "beta")?;
    validate_positive(sigma, "sigma")?;
    validate_positive(gamma, "gamma")?;
    validate_positive(dt, "dt")?;

    seir_step_inner(s, e, i, r, beta, sigma, gamma, dt)
}

#[cfg(feature = "hisab")]
#[allow(clippy::too_many_arguments)]
fn seir_step_inner(
    s: f64,
    e: f64,
    i: f64,
    r: f64,
    beta: f64,
    sigma: f64,
    gamma: f64,
    dt: f64,
) -> Result<SeirState> {
    let y0 = [s, e, i, r];
    let result = hisab::num::rk4(
        |_t, y, dy| {
            dy[0] = -beta * y[0] * y[2];
            dy[1] = beta * y[0] * y[2] - sigma * y[1];
            dy[2] = sigma * y[1] - gamma * y[2];
            dy[3] = gamma * y[2];
        },
        0.0,
        &y0,
        dt,
        1,
    )
    .map_err(|e| crate::error::JivanuError::SimulationFailed(e.to_string()))?;
    Ok(SeirState {
        s: result[0].max(0.0),
        e: result[1].max(0.0),
        i: result[2].max(0.0),
        r: result[3].max(0.0),
    })
}

#[cfg(not(feature = "hisab"))]
#[allow(clippy::too_many_arguments)]
fn seir_step_inner(
    s: f64,
    e: f64,
    i: f64,
    r: f64,
    beta: f64,
    sigma: f64,
    gamma: f64,
    dt: f64,
) -> Result<SeirState> {
    let ds = -beta * s * i;
    let de = beta * s * i - sigma * e;
    let di = sigma * e - gamma * i;
    let dr = gamma * i;
    Ok(SeirState {
        s: (s + ds * dt).max(0.0),
        e: (e + de * dt).max(0.0),
        i: (i + di * dt).max(0.0),
        r: (r + dr * dt).max(0.0),
    })
}

/// Basic reproduction number.
///
/// `R0 = beta / gamma`
///
/// # Errors
///
/// Returns error if gamma is non-positive.
#[inline]
#[must_use = "returns R0 without side effects"]
pub fn r_naught(beta: f64, gamma: f64) -> Result<f64> {
    validate_positive(beta, "beta")?;
    validate_positive(gamma, "gamma")?;
    Ok(beta / gamma)
}

/// Herd immunity threshold.
///
/// `H = 1 - 1/R0`
///
/// # Errors
///
/// Returns error if R0 is non-positive.
#[inline]
#[must_use = "returns the herd immunity threshold without side effects"]
pub fn herd_immunity_threshold(r0: f64) -> Result<f64> {
    validate_positive(r0, "r0")?;
    Ok(1.0 - 1.0 / r0)
}

/// Case fatality rate.
///
/// `CFR = deaths / cases`
///
/// # Errors
///
/// Returns error if cases is zero.
#[inline]
#[must_use = "returns the CFR without side effects"]
pub fn case_fatality_rate(deaths: u64, cases: u64) -> Result<f64> {
    if cases == 0 {
        return Err(crate::error::JivanuError::ComputationError(
            "cases must be > 0".into(),
        ));
    }
    Ok(deaths as f64 / cases as f64)
}

/// Run a full SIR trajectory.
///
/// When the `hisab` feature is enabled, uses [`hisab::num::rk4_trajectory`]
/// for fourth-order Runge-Kutta integration with Neumaier-compensated
/// accumulation across all steps.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the trajectory without side effects"]
pub fn sir_trajectory(
    s0: f64,
    i0: f64,
    r0: f64,
    beta: f64,
    gamma: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<SirState>> {
    validate_non_negative(s0, "s0")?;
    validate_non_negative(i0, "i0")?;
    validate_non_negative(r0, "r0")?;
    validate_positive(beta, "beta")?;
    validate_positive(gamma, "gamma")?;
    validate_positive(dt, "dt")?;

    sir_trajectory_inner(s0, i0, r0, beta, gamma, dt, steps)
}

#[cfg(feature = "hisab")]
fn sir_trajectory_inner(
    s0: f64,
    i0: f64,
    r0: f64,
    beta: f64,
    gamma: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<SirState>> {
    let y0 = [s0, i0, r0];
    let t_end = dt * steps as f64;
    let raw = hisab::num::rk4_trajectory(
        |_t, y, dy| {
            dy[0] = -beta * y[0] * y[1];
            dy[1] = beta * y[0] * y[1] - gamma * y[1];
            dy[2] = gamma * y[1];
        },
        0.0,
        &y0,
        t_end,
        steps,
    )
    .map_err(|e| crate::error::JivanuError::SimulationFailed(e.to_string()))?;
    Ok(raw
        .into_iter()
        .map(|(_t, y)| SirState {
            s: y[0].max(0.0),
            i: y[1].max(0.0),
            r: y[2].max(0.0),
        })
        .collect())
}

#[cfg(not(feature = "hisab"))]
fn sir_trajectory_inner(
    s0: f64,
    i0: f64,
    r0: f64,
    beta: f64,
    gamma: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<SirState>> {
    let mut traj = Vec::with_capacity(steps + 1);
    let mut state = SirState {
        s: s0,
        i: i0,
        r: r0,
    };
    traj.push(state);
    for _ in 0..steps {
        state = sir_step(state.s, state.i, state.r, beta, gamma, dt)?;
        traj.push(state);
    }
    Ok(traj)
}

/// SIRS model state (SIR with waning immunity).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SirsState {
    /// Susceptible fraction.
    pub s: f64,
    /// Infected fraction.
    pub i: f64,
    /// Recovered fraction.
    pub r: f64,
}

/// Parameters for the SIRS model with optional vaccination.
///
/// ```text
/// dS/dt = -β·S·I + δ·R - ν·S
/// dI/dt =  β·S·I - γ·I
/// dR/dt =  γ·I - δ·R + ν·S
/// ```
///
/// - `delta` — rate of immunity waning (1/duration of immunity)
/// - `nu` — vaccination rate (fraction of susceptibles vaccinated per time unit)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SirsParams {
    /// Transmission rate.
    pub beta: f64,
    /// Recovery rate.
    pub gamma: f64,
    /// Waning immunity rate (0 = permanent immunity = standard SIR).
    pub delta: f64,
    /// Vaccination rate of susceptibles (0 = no vaccination).
    pub nu: f64,
    /// Time step.
    pub dt: f64,
}

/// One step of the SIRS model with optional vaccination.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the new SIRS state without side effects"]
pub fn sirs_step(state: &SirsState, params: &SirsParams) -> Result<SirsState> {
    validate_non_negative(state.s, "s")?;
    validate_non_negative(state.i, "i")?;
    validate_non_negative(state.r, "r")?;
    validate_positive(params.beta, "beta")?;
    validate_positive(params.gamma, "gamma")?;
    validate_non_negative(params.delta, "delta")?;
    validate_non_negative(params.nu, "nu")?;
    validate_positive(params.dt, "dt")?;

    sirs_step_inner(state, params)
}

#[cfg(feature = "hisab")]
fn sirs_step_inner(state: &SirsState, params: &SirsParams) -> Result<SirsState> {
    let SirsParams {
        beta,
        gamma,
        delta,
        nu,
        dt,
    } = *params;
    let y0 = [state.s, state.i, state.r];
    let result = hisab::num::rk4(
        |_t, y, dy| {
            dy[0] = -beta * y[0] * y[1] + delta * y[2] - nu * y[0];
            dy[1] = beta * y[0] * y[1] - gamma * y[1];
            dy[2] = gamma * y[1] - delta * y[2] + nu * y[0];
        },
        0.0,
        &y0,
        dt,
        1,
    )
    .map_err(|e| crate::error::JivanuError::SimulationFailed(e.to_string()))?;
    Ok(SirsState {
        s: result[0].max(0.0),
        i: result[1].max(0.0),
        r: result[2].max(0.0),
    })
}

#[cfg(not(feature = "hisab"))]
fn sirs_step_inner(state: &SirsState, params: &SirsParams) -> Result<SirsState> {
    let SirsState { s, i, r } = *state;
    let SirsParams {
        beta,
        gamma,
        delta,
        nu,
        dt,
    } = *params;
    let ds = -beta * s * i + delta * r - nu * s;
    let di = beta * s * i - gamma * i;
    let dr = gamma * i - delta * r + nu * s;
    Ok(SirsState {
        s: (s + ds * dt).max(0.0),
        i: (i + di * dt).max(0.0),
        r: (r + dr * dt).max(0.0),
    })
}

/// Effective reproduction number accounting for vaccination.
///
/// `R_eff = R0 × (1 - vaccination_coverage)`
///
/// When `R_eff < 1`, the disease cannot sustain transmission.
///
/// # Errors
///
/// Returns error if R0 is non-positive or coverage is outside [0, 1].
#[inline]
#[must_use = "returns R_eff without side effects"]
pub fn effective_r(r0: f64, vaccination_coverage: f64) -> Result<f64> {
    validate_positive(r0, "r0")?;
    if !(0.0..=1.0).contains(&vaccination_coverage) {
        return Err(crate::error::JivanuError::ComputationError(
            "vaccination_coverage must be in [0, 1]".into(),
        ));
    }
    Ok(r0 * (1.0 - vaccination_coverage))
}

/// Critical vaccination coverage to achieve herd immunity.
///
/// `V_c = 1 - 1/R0`
///
/// Same formula as herd immunity threshold, but framed as the
/// vaccination target needed to prevent epidemics.
///
/// # Errors
///
/// Returns error if R0 is non-positive.
#[inline]
#[must_use = "returns the critical vaccination coverage without side effects"]
pub fn critical_vaccination_coverage(r0: f64) -> Result<f64> {
    validate_positive(r0, "r0")?;
    Ok((1.0 - 1.0 / r0).max(0.0))
}

/// Run a full SEIR trajectory.
///
/// When the `hisab` feature is enabled, uses RK4 integration.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the trajectory without side effects"]
pub fn seir_trajectory(
    state0: &SeirState,
    params: &SeirParams,
    steps: usize,
) -> Result<Vec<SeirState>> {
    validate_non_negative(state0.s, "s")?;
    validate_non_negative(state0.e, "e")?;
    validate_non_negative(state0.i, "i")?;
    validate_non_negative(state0.r, "r")?;
    validate_positive(params.beta, "beta")?;
    validate_positive(params.sigma, "sigma")?;
    validate_positive(params.gamma, "gamma")?;
    validate_positive(params.dt, "dt")?;

    seir_trajectory_inner(state0, params, steps)
}

#[cfg(feature = "hisab")]
fn seir_trajectory_inner(
    state0: &SeirState,
    params: &SeirParams,
    steps: usize,
) -> Result<Vec<SeirState>> {
    let SeirParams {
        beta,
        sigma,
        gamma,
        dt,
    } = *params;
    let y0 = [state0.s, state0.e, state0.i, state0.r];
    let t_end = dt * steps as f64;
    let raw = hisab::num::rk4_trajectory(
        |_t, y, dy| {
            dy[0] = -beta * y[0] * y[2];
            dy[1] = beta * y[0] * y[2] - sigma * y[1];
            dy[2] = sigma * y[1] - gamma * y[2];
            dy[3] = gamma * y[2];
        },
        0.0,
        &y0,
        t_end,
        steps,
    )
    .map_err(|e| crate::error::JivanuError::SimulationFailed(e.to_string()))?;
    Ok(raw
        .into_iter()
        .map(|(_t, y)| SeirState {
            s: y[0].max(0.0),
            e: y[1].max(0.0),
            i: y[2].max(0.0),
            r: y[3].max(0.0),
        })
        .collect())
}

#[cfg(not(feature = "hisab"))]
fn seir_trajectory_inner(
    state0: &SeirState,
    params: &SeirParams,
    steps: usize,
) -> Result<Vec<SeirState>> {
    let mut traj = Vec::with_capacity(steps + 1);
    let mut state = *state0;
    traj.push(state);
    for _ in 0..steps {
        state = seir_step(&state, params)?;
        traj.push(state);
    }
    Ok(traj)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_r_naught() {
        let r0 = r_naught(0.5, 0.2).unwrap();
        assert!((r0 - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_herd_immunity_r0_2_5() {
        let h = herd_immunity_threshold(2.5).unwrap();
        assert!((h - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_herd_immunity_r0_3() {
        let h = herd_immunity_threshold(3.0).unwrap();
        assert!((h - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_sir_conservation() {
        let state = sir_step(0.99, 0.01, 0.0, 0.5, 0.1, 0.01).unwrap();
        let total = state.s + state.i + state.r;
        assert!((total - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_sir_declining_r0_lt_1() {
        let state = sir_step(0.99, 0.01, 0.0, 0.1, 0.5, 0.1).unwrap();
        assert!(state.i < 0.01);
    }

    #[test]
    fn test_seir_step() {
        let state = SeirState {
            s: 0.98,
            e: 0.01,
            i: 0.01,
            r: 0.0,
        };
        let params = SeirParams {
            beta: 0.5,
            sigma: 0.2,
            gamma: 0.1,
            dt: 0.01,
        };
        let new_state = seir_step(&state, &params).unwrap();
        let total = new_state.s + new_state.e + new_state.i + new_state.r;
        assert!((total - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_case_fatality_rate() {
        let cfr = case_fatality_rate(10, 1000).unwrap();
        assert!((cfr - 0.01).abs() < 1e-10);
    }

    #[test]
    fn test_case_fatality_rate_zero_cases() {
        assert!(case_fatality_rate(0, 0).is_err());
    }

    #[test]
    fn test_sir_trajectory() {
        let traj = sir_trajectory(0.99, 0.01, 0.0, 0.5, 0.1, 0.1, 10).unwrap();
        assert_eq!(traj.len(), 11);
    }

    #[test]
    fn test_sir_state_serde_roundtrip() {
        let state = SirState {
            s: 0.9,
            i: 0.05,
            r: 0.05,
        };
        let json = serde_json::to_string(&state).unwrap();
        let back: SirState = serde_json::from_str(&json).unwrap();
        assert!((state.s - back.s).abs() < 1e-10);
    }

    #[test]
    fn test_seir_state_serde_roundtrip() {
        let state = SeirState {
            s: 0.9,
            e: 0.03,
            i: 0.05,
            r: 0.02,
        };
        let json = serde_json::to_string(&state).unwrap();
        let back: SeirState = serde_json::from_str(&json).unwrap();
        assert!((state.e - back.e).abs() < 1e-10);
    }

    #[test]
    fn test_seir_params_serde_roundtrip() {
        let params = SeirParams {
            beta: 0.5,
            sigma: 0.2,
            gamma: 0.1,
            dt: 0.01,
        };
        let json = serde_json::to_string(&params).unwrap();
        let back: SeirParams = serde_json::from_str(&json).unwrap();
        assert!((params.sigma - back.sigma).abs() < 1e-10);
    }

    /// With RK4, conservation should hold to near machine precision even
    /// with a coarse time step. With Euler, it drifts.
    #[test]
    fn test_sir_trajectory_conservation_tight() {
        let traj = sir_trajectory(0.99, 0.01, 0.0, 0.5, 0.1, 0.1, 100).unwrap();
        for (idx, state) in traj.iter().enumerate() {
            let total = state.s + state.i + state.r;
            // RK4 (hisab): < 1e-10 drift. Euler: < 1e-3 drift at dt=0.1.
            #[cfg(feature = "hisab")]
            assert!(
                (total - 1.0).abs() < 1e-10,
                "step {idx}: S+I+R = {total}, drift = {}",
                (total - 1.0).abs()
            );
            #[cfg(not(feature = "hisab"))]
            assert!(
                (total - 1.0).abs() < 1e-3,
                "step {idx}: S+I+R = {total}, drift = {}",
                (total - 1.0).abs()
            );
        }
    }

    /// SEIR conservation over many steps.
    #[test]
    fn test_seir_trajectory_conservation() {
        let mut state = SeirState {
            s: 0.97,
            e: 0.02,
            i: 0.01,
            r: 0.0,
        };
        let params = SeirParams {
            beta: 0.5,
            sigma: 0.2,
            gamma: 0.1,
            dt: 0.1,
        };
        for step in 0..100 {
            state = seir_step(&state, &params).unwrap();
            let total = state.s + state.e + state.i + state.r;
            #[cfg(feature = "hisab")]
            assert!(
                (total - 1.0).abs() < 1e-10,
                "step {step}: S+E+I+R = {total}"
            );
            #[cfg(not(feature = "hisab"))]
            assert!((total - 1.0).abs() < 1e-3, "step {step}: S+E+I+R = {total}");
        }
    }

    // --- SIRS + vaccination tests ---

    #[test]
    fn test_sirs_waning_immunity() {
        // With waning immunity (delta > 0), recovered return to susceptible
        let state = SirsState {
            s: 0.1,
            i: 0.0,
            r: 0.9,
        };
        let params = SirsParams {
            beta: 0.5,
            gamma: 0.1,
            delta: 0.05, // immunity wanes
            nu: 0.0,
            dt: 1.0,
        };
        let next = sirs_step(&state, &params).unwrap();
        assert!(next.s > state.s, "susceptibles should increase from waning");
        assert!(next.r < state.r, "recovered should decrease from waning");
    }

    #[test]
    fn test_sirs_vaccination() {
        // Vaccination moves susceptibles to recovered
        let state = SirsState {
            s: 0.9,
            i: 0.0,
            r: 0.1,
        };
        let params = SirsParams {
            beta: 0.5,
            gamma: 0.1,
            delta: 0.0,
            nu: 0.1, // 10% vaccination rate
            dt: 0.1,
        };
        let next = sirs_step(&state, &params).unwrap();
        assert!(
            next.s < state.s,
            "susceptibles should decrease from vaccination"
        );
        assert!(
            next.r > state.r,
            "recovered should increase from vaccination"
        );
    }

    #[test]
    fn test_sirs_conservation() {
        let state = SirsState {
            s: 0.8,
            i: 0.1,
            r: 0.1,
        };
        let params = SirsParams {
            beta: 0.5,
            gamma: 0.1,
            delta: 0.05,
            nu: 0.02,
            dt: 0.01,
        };
        let next = sirs_step(&state, &params).unwrap();
        let total = next.s + next.i + next.r;
        assert!((total - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_sirs_delta_zero_is_sir() {
        // SIRS with delta=0, nu=0 should behave like SIR
        let sirs_state = SirsState {
            s: 0.99,
            i: 0.01,
            r: 0.0,
        };
        let sirs_params = SirsParams {
            beta: 0.5,
            gamma: 0.1,
            delta: 0.0,
            nu: 0.0,
            dt: 0.01,
        };
        let sir_result = sir_step(0.99, 0.01, 0.0, 0.5, 0.1, 0.01).unwrap();
        let sirs_result = sirs_step(&sirs_state, &sirs_params).unwrap();
        assert!((sir_result.s - sirs_result.s).abs() < 1e-6);
        assert!((sir_result.i - sirs_result.i).abs() < 1e-6);
    }

    #[test]
    fn test_effective_r() {
        // R0=3, 60% vaccinated → R_eff = 1.2
        let r_eff = effective_r(3.0, 0.6).unwrap();
        assert!((r_eff - 1.2).abs() < 1e-10);
    }

    #[test]
    fn test_effective_r_full_coverage() {
        let r_eff = effective_r(3.0, 1.0).unwrap();
        assert!((r_eff - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_effective_r_invalid_coverage() {
        assert!(effective_r(3.0, 1.5).is_err());
        assert!(effective_r(3.0, -0.1).is_err());
    }

    #[test]
    fn test_critical_vaccination_coverage() {
        // R0=3 → Vc = 1 - 1/3 = 0.667
        let vc = critical_vaccination_coverage(3.0).unwrap();
        assert!((vc - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_critical_vaccination_r0_lt_1() {
        // R0 < 1 → no vaccination needed, Vc = 0
        let vc = critical_vaccination_coverage(0.5).unwrap();
        assert!((vc - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_seir_trajectory_length() {
        let state = SeirState {
            s: 0.97,
            e: 0.02,
            i: 0.01,
            r: 0.0,
        };
        let params = SeirParams {
            beta: 0.5,
            sigma: 0.2,
            gamma: 0.1,
            dt: 0.1,
        };
        let traj = seir_trajectory(&state, &params, 10).unwrap();
        assert_eq!(traj.len(), 11);
    }

    #[test]
    fn test_seir_trajectory_starts_at_initial() {
        let state = SeirState {
            s: 0.97,
            e: 0.02,
            i: 0.01,
            r: 0.0,
        };
        let params = SeirParams {
            beta: 0.5,
            sigma: 0.2,
            gamma: 0.1,
            dt: 0.1,
        };
        let traj = seir_trajectory(&state, &params, 5).unwrap();
        assert!((traj[0].s - 0.97).abs() < 1e-10);
        assert!((traj[0].e - 0.02).abs() < 1e-10);
    }

    #[test]
    fn test_sirs_state_serde_roundtrip() {
        let state = SirsState {
            s: 0.8,
            i: 0.1,
            r: 0.1,
        };
        let json = serde_json::to_string(&state).unwrap();
        let back: SirsState = serde_json::from_str(&json).unwrap();
        assert!((state.s - back.s).abs() < 1e-10);
    }

    #[test]
    fn test_sirs_params_serde_roundtrip() {
        let params = SirsParams {
            beta: 0.5,
            gamma: 0.1,
            delta: 0.05,
            nu: 0.02,
            dt: 0.01,
        };
        let json = serde_json::to_string(&params).unwrap();
        let back: SirsParams = serde_json::from_str(&json).unwrap();
        assert!((params.delta - back.delta).abs() < 1e-10);
    }
}
