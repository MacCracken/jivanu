//! Microbial growth — exponential, logistic, Monod kinetics.

use serde::{Deserialize, Serialize};

use crate::error::{validate_finite, validate_non_negative, validate_positive, JivanuError, Result};

/// Exponential growth: `N = N0 * e^(r * t)`.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the population without side effects"]
pub fn exponential_growth(n0: f64, rate: f64, time: f64) -> Result<f64> {
    validate_positive(n0, "n0")?;
    validate_finite(rate, "rate")?;
    validate_non_negative(time, "time")?;
    Ok(n0 * (rate * time).exp())
}

/// Logistic growth: `N = K / (1 + ((K - N0) / N0) * e^(-r * t))`.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the population without side effects"]
pub fn logistic_growth(n0: f64, rate: f64, capacity: f64, time: f64) -> Result<f64> {
    validate_positive(n0, "n0")?;
    validate_finite(rate, "rate")?;
    validate_positive(capacity, "capacity")?;
    validate_non_negative(time, "time")?;
    let ratio = (capacity - n0) / n0;
    Ok(capacity / (1.0 + ratio * (-rate * time).exp()))
}

/// Doubling time from growth rate.
///
/// `t_d = ln(2) / r`
///
/// # Errors
///
/// Returns error if rate is non-positive.
#[inline]
#[must_use = "returns the doubling time without side effects"]
pub fn doubling_time(rate: f64) -> Result<f64> {
    validate_positive(rate, "rate")?;
    Ok(core::f64::consts::LN_2 / rate)
}

/// Monod kinetics: specific growth rate as a function of substrate concentration.
///
/// `mu = mu_max * S / (K_s + S)`
///
/// At `S = K_s`, `mu = mu_max / 2` (half-saturation).
/// At `S >> K_s`, `mu ≈ mu_max`.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the specific growth rate without side effects"]
pub fn monod_kinetics(substrate: f64, mu_max: f64, k_s: f64) -> Result<f64> {
    validate_non_negative(substrate, "substrate")?;
    validate_positive(mu_max, "mu_max")?;
    validate_positive(k_s, "k_s")?;
    Ok(mu_max * substrate / (k_s + substrate))
}

/// Growth phase of a microbial culture.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[non_exhaustive]
pub enum GrowthPhase {
    /// Lag phase: adaptation, no significant growth.
    Lag,
    /// Exponential (log) phase: rapid, constant-rate growth.
    Exponential,
    /// Stationary phase: growth rate equals death rate.
    Stationary,
    /// Death (decline) phase: cell death exceeds growth.
    Death,
}

/// Chemostat steady-state biomass concentration.
///
/// `X = Y * (S_in - S_ss)` where `S_ss = K_s * D / (mu_max - D)`
///
/// # Errors
///
/// Returns error if dilution rate exceeds mu_max (washout).
#[must_use = "returns the steady-state biomass without side effects"]
pub fn chemostat_steady_state(
    dilution_rate: f64,
    mu_max: f64,
    k_s: f64,
    substrate_feed: f64,
    yield_coeff: f64,
) -> Result<(f64, f64)> {
    validate_positive(dilution_rate, "dilution_rate")?;
    validate_positive(mu_max, "mu_max")?;
    validate_positive(k_s, "k_s")?;
    validate_non_negative(substrate_feed, "substrate_feed")?;
    validate_positive(yield_coeff, "yield_coeff")?;

    if dilution_rate >= mu_max {
        return Err(JivanuError::SimulationFailed(
            "dilution rate exceeds mu_max: washout condition".into(),
        ));
    }

    let s_ss = k_s * dilution_rate / (mu_max - dilution_rate);
    let x_ss = yield_coeff * (substrate_feed - s_ss);
    Ok((x_ss.max(0.0), s_ss))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exponential_growth_t0() {
        let n = exponential_growth(100.0, 0.5, 0.0).unwrap();
        assert!((n - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_exponential_growth_doubling() {
        // At t = ln(2)/r, population should double.
        let r = core::f64::consts::LN_2; // r = ln(2), so doubling time = 1
        let n = exponential_growth(100.0, r, 1.0).unwrap();
        assert!((n - 200.0).abs() < 1e-6);
    }

    #[test]
    fn test_logistic_growth_approaches_capacity() {
        let n = logistic_growth(10.0, 0.5, 1000.0, 100.0).unwrap();
        assert!((n - 1000.0).abs() < 1.0); // should be very close to K
    }

    #[test]
    fn test_logistic_growth_at_t0() {
        let n = logistic_growth(100.0, 0.5, 1000.0, 0.0).unwrap();
        assert!((n - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_doubling_time() {
        let td = doubling_time(core::f64::consts::LN_2).unwrap();
        assert!((td - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_doubling_time_at_rate_0_693() {
        // mu = 0.693/hr → doubling time ≈ 1 hour
        let td = doubling_time(0.693).unwrap();
        assert!((td - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_monod_at_high_substrate() {
        // S >> K_s: mu ≈ mu_max
        let mu = monod_kinetics(1000.0, 0.5, 1.0).unwrap();
        assert!((mu - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_monod_at_half_saturation() {
        // S = K_s: mu = mu_max / 2
        let mu = monod_kinetics(1.0, 0.5, 1.0).unwrap();
        assert!((mu - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_monod_at_zero_substrate() {
        let mu = monod_kinetics(0.0, 0.5, 1.0).unwrap();
        assert!((mu - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_growth_phase_ordering() {
        assert!(GrowthPhase::Lag < GrowthPhase::Exponential);
        assert!(GrowthPhase::Exponential < GrowthPhase::Stationary);
        assert!(GrowthPhase::Stationary < GrowthPhase::Death);
    }

    #[test]
    fn test_chemostat_steady_state() {
        let (x, s) = chemostat_steady_state(0.2, 0.5, 0.1, 10.0, 0.5).unwrap();
        assert!(x > 0.0);
        assert!(s > 0.0);
        assert!(s < 10.0); // substrate consumed
    }

    #[test]
    fn test_chemostat_washout() {
        // D >= mu_max → washout
        assert!(chemostat_steady_state(0.6, 0.5, 0.1, 10.0, 0.5).is_err());
    }

    #[test]
    fn test_growth_phase_serde_roundtrip() {
        let phase = GrowthPhase::Exponential;
        let json = serde_json::to_string(&phase).unwrap();
        let back: GrowthPhase = serde_json::from_str(&json).unwrap();
        assert_eq!(phase, back);
    }
}
