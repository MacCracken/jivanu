//! Microbial growth — exponential, logistic, Monod kinetics.

use serde::{Deserialize, Serialize};

use crate::error::{
    JivanuError, Result, validate_finite, validate_non_negative, validate_positive,
};

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

/// State of a two-strain competition system.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TwoStrainState {
    /// Population density of strain 1.
    pub n1: f64,
    /// Population density of strain 2.
    pub n2: f64,
}

/// Parameters for two-strain Lotka-Volterra competition.
///
/// Each strain has its own growth rate and carrying capacity. The competition
/// coefficients `alpha12` and `alpha21` describe the inhibitory effect of
/// one strain on the other:
///
/// ```text
/// dN1/dt = r1 * N1 * (1 - (N1 + α12*N2) / K1)
/// dN2/dt = r2 * N2 * (1 - (N2 + α21*N1) / K2)
/// ```
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CompetitionParams {
    /// Intrinsic growth rate of strain 1.
    pub r1: f64,
    /// Intrinsic growth rate of strain 2.
    pub r2: f64,
    /// Carrying capacity of strain 1.
    pub k1: f64,
    /// Carrying capacity of strain 2.
    pub k2: f64,
    /// Effect of strain 2 on strain 1 (competition coefficient).
    pub alpha12: f64,
    /// Effect of strain 1 on strain 2 (competition coefficient).
    pub alpha21: f64,
}

/// Outcome of two-strain competition at equilibrium.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum CompetitionOutcome {
    /// Strain 1 excludes strain 2.
    Strain1Wins,
    /// Strain 2 excludes strain 1.
    Strain2Wins,
    /// Stable coexistence (each strain inhibits itself more than the other).
    Coexistence,
    /// Unstable — outcome depends on initial conditions (priority effects).
    Unstable,
}

/// Predict the equilibrium outcome of two-strain Lotka-Volterra competition.
///
/// Uses the competitive exclusion principle:
/// - Coexistence when `α12 < K1/K2` AND `α21 < K2/K1`
/// - Strain 1 wins when `α12 < K1/K2` AND `α21 > K2/K1`
/// - Strain 2 wins when `α12 > K1/K2` AND `α21 < K2/K1`
/// - Unstable when `α12 > K1/K2` AND `α21 > K2/K1`
///
/// # Errors
///
/// Returns error if carrying capacities are non-positive or competition
/// coefficients are negative.
#[inline]
#[must_use = "returns the competition outcome without side effects"]
pub fn competition_outcome(params: &CompetitionParams) -> Result<CompetitionOutcome> {
    validate_positive(params.k1, "k1")?;
    validate_positive(params.k2, "k2")?;
    validate_non_negative(params.alpha12, "alpha12")?;
    validate_non_negative(params.alpha21, "alpha21")?;

    let ratio_12 = params.k1 / params.k2;
    let ratio_21 = params.k2 / params.k1;

    let s1_tolerates_s2 = params.alpha12 < ratio_12;
    let s2_tolerates_s1 = params.alpha21 < ratio_21;

    Ok(match (s1_tolerates_s2, s2_tolerates_s1) {
        (true, true) => CompetitionOutcome::Coexistence,
        (true, false) => CompetitionOutcome::Strain1Wins,
        (false, true) => CompetitionOutcome::Strain2Wins,
        (false, false) => CompetitionOutcome::Unstable,
    })
}

/// One step of two-strain Lotka-Volterra competition.
///
/// When the `hisab` feature is enabled, uses fourth-order Runge-Kutta.
/// Otherwise uses forward Euler.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the new state without side effects"]
pub fn competition_step(
    state: &TwoStrainState,
    params: &CompetitionParams,
    dt: f64,
) -> Result<TwoStrainState> {
    validate_non_negative(state.n1, "n1")?;
    validate_non_negative(state.n2, "n2")?;
    validate_positive(params.r1, "r1")?;
    validate_positive(params.r2, "r2")?;
    validate_positive(params.k1, "k1")?;
    validate_positive(params.k2, "k2")?;
    validate_non_negative(params.alpha12, "alpha12")?;
    validate_non_negative(params.alpha21, "alpha21")?;
    validate_positive(dt, "dt")?;

    competition_step_inner(state, params, dt)
}

#[cfg(feature = "hisab")]
fn competition_step_inner(
    state: &TwoStrainState,
    params: &CompetitionParams,
    dt: f64,
) -> Result<TwoStrainState> {
    let CompetitionParams {
        r1,
        r2,
        k1,
        k2,
        alpha12,
        alpha21,
    } = *params;
    let y0 = [state.n1, state.n2];
    let result = hisab::num::rk4(
        |_t, y, dy| {
            dy[0] = r1 * y[0] * (1.0 - (y[0] + alpha12 * y[1]) / k1);
            dy[1] = r2 * y[1] * (1.0 - (y[1] + alpha21 * y[0]) / k2);
        },
        0.0,
        &y0,
        dt,
        1,
    )
    .map_err(|e| JivanuError::SimulationFailed(e.to_string()))?;
    Ok(TwoStrainState {
        n1: result[0].max(0.0),
        n2: result[1].max(0.0),
    })
}

#[cfg(not(feature = "hisab"))]
fn competition_step_inner(
    state: &TwoStrainState,
    params: &CompetitionParams,
    dt: f64,
) -> Result<TwoStrainState> {
    let n1 = state.n1;
    let n2 = state.n2;
    let dn1 = params.r1 * n1 * (1.0 - (n1 + params.alpha12 * n2) / params.k1);
    let dn2 = params.r2 * n2 * (1.0 - (n2 + params.alpha21 * n1) / params.k2);
    Ok(TwoStrainState {
        n1: (n1 + dn1 * dt).max(0.0),
        n2: (n2 + dn2 * dt).max(0.0),
    })
}

/// Baranyi-Roberts growth model with explicit lag phase.
///
/// Uses the commonly applied three-parameter form:
/// `N(t) = N_max / (1 + ((N_max/N_0) - 1) × e^(-μ_max × (t - lag)))`
/// for `t > lag`, and `N(t) = N_0` for `t ≤ lag`.
///
/// Reference: Baranyi & Roberts (1994) Int J Food Microbiol 23:277-294.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the population without side effects"]
pub fn baranyi_roberts(n0: f64, n_max: f64, mu_max: f64, lag: f64, t: f64) -> Result<f64> {
    validate_positive(n0, "n0")?;
    validate_positive(n_max, "n_max")?;
    validate_positive(mu_max, "mu_max")?;
    validate_non_negative(lag, "lag")?;
    validate_non_negative(t, "t")?;

    if t <= lag {
        return Ok(n0);
    }
    let t_eff = t - lag;
    let ratio = (n_max - n0) / n0;
    Ok(n_max / (1.0 + ratio * (-mu_max * t_eff).exp()))
}

/// Gompertz growth model.
///
/// `ln(N/N_0) = A × exp(-exp((μ_max × e / A) × (lag - t) + 1))`
///
/// where `A = ln(N_max / N_0)`.
///
/// Returns the population at time t.
///
/// Reference: Zwietering et al. (1990) Appl Environ Microbiol 56:1875-1881.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the population without side effects"]
pub fn gompertz(n0: f64, n_max: f64, mu_max: f64, lag: f64, t: f64) -> Result<f64> {
    validate_positive(n0, "n0")?;
    validate_positive(n_max, "n_max")?;
    validate_positive(mu_max, "mu_max")?;
    validate_non_negative(lag, "lag")?;
    validate_non_negative(t, "t")?;

    let a = (n_max / n0).ln();
    let e = core::f64::consts::E;
    let exponent = (mu_max * e / a) * (lag - t) + 1.0;
    let ln_ratio = a * (-exponent.exp()).exp();
    Ok(n0 * ln_ratio.exp())
}

/// Temperature growth rate modifier (Ratkowsky square-root model).
///
/// `√μ = b × (T - T_min)` → `μ = b² × (T - T_min)²`
///
/// Valid for `T > T_min`. Returns 0 at or below `T_min`.
///
/// Reference: Ratkowsky et al. (1982) J Bacteriol 149:1-5.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the growth rate without side effects"]
pub fn ratkowsky_temperature(b: f64, t_celsius: f64, t_min: f64) -> Result<f64> {
    validate_positive(b, "b")?;
    validate_finite(t_celsius, "t_celsius")?;
    validate_finite(t_min, "t_min")?;

    if t_celsius <= t_min {
        return Ok(0.0);
    }
    let diff = t_celsius - t_min;
    Ok(b * b * diff * diff)
}

/// Cardinal temperature model (Rosso).
///
/// Returns a growth modifier in [0, 1] based on temperature.
/// Growth is zero outside [T_min, T_max], optimal at T_opt.
///
/// Reference: Rosso et al. (1993) J Theor Biol 162:447-463.
///
/// # Errors
///
/// Returns error if cardinal temperatures are not ordered T_min < T_opt < T_max.
#[must_use = "returns the growth rate modifier without side effects"]
pub fn cardinal_temperature(t_celsius: f64, t_min: f64, t_opt: f64, t_max: f64) -> Result<f64> {
    validate_finite(t_celsius, "t_celsius")?;
    validate_finite(t_min, "t_min")?;
    validate_finite(t_opt, "t_opt")?;
    validate_finite(t_max, "t_max")?;

    if t_min >= t_opt || t_opt >= t_max {
        return Err(JivanuError::ComputationError(
            "cardinal temperatures must satisfy T_min < T_opt < T_max".into(),
        ));
    }

    if t_celsius <= t_min || t_celsius >= t_max {
        return Ok(0.0);
    }

    let num = (t_celsius - t_max) * (t_celsius - t_min).powi(2);
    let denom = (t_opt - t_min)
        * ((t_opt - t_min) * (t_celsius - t_opt)
            - (t_opt - t_max) * (t_opt + t_min - 2.0 * t_celsius));
    Ok((num / denom).clamp(0.0, 1.0))
}

/// Cardinal pH model.
///
/// Returns a growth modifier in [0, 1] based on pH.
/// Growth is zero outside [pH_min, pH_max], optimal at pH_opt.
///
/// Reference: Rosso et al. (1995) Appl Environ Microbiol.
///
/// # Errors
///
/// Returns error if cardinal pH values are not ordered.
#[must_use = "returns the pH growth modifier without side effects"]
pub fn cardinal_ph(ph: f64, ph_min: f64, ph_opt: f64, ph_max: f64) -> Result<f64> {
    validate_finite(ph, "ph")?;
    validate_finite(ph_min, "ph_min")?;
    validate_finite(ph_opt, "ph_opt")?;
    validate_finite(ph_max, "ph_max")?;

    if ph_min >= ph_opt || ph_opt >= ph_max {
        return Err(JivanuError::ComputationError(
            "cardinal pH must satisfy pH_min < pH_opt < pH_max".into(),
        ));
    }

    if ph <= ph_min || ph >= ph_max {
        return Ok(0.0);
    }

    let num = (ph - ph_min) * (ph - ph_max);
    let denom = num - (ph - ph_opt).powi(2);
    if denom.abs() < 1e-15 {
        return Ok(0.0);
    }
    Ok((num / denom).clamp(0.0, 1.0))
}

/// Herbert-Pirt maintenance energy: specific substrate consumption rate.
///
/// `q_s = μ / Y_true + m_s`
///
/// Reference: Pirt (1965) Proc R Soc Lond B 163:224-231.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the specific substrate consumption rate without side effects"]
pub fn maintenance_substrate_rate(mu: f64, y_true: f64, m_s: f64) -> Result<f64> {
    validate_non_negative(mu, "mu")?;
    validate_positive(y_true, "y_true")?;
    validate_non_negative(m_s, "m_s")?;
    Ok(mu / y_true + m_s)
}

/// Apparent yield accounting for maintenance energy.
///
/// `Y_app = μ / (μ/Y_true + m_s)`
///
/// # Errors
///
/// Returns error if μ is zero or parameters are invalid.
#[inline]
#[must_use = "returns the apparent yield without side effects"]
pub fn apparent_yield(mu: f64, y_true: f64, m_s: f64) -> Result<f64> {
    validate_positive(mu, "mu")?;
    validate_positive(y_true, "y_true")?;
    validate_non_negative(m_s, "m_s")?;
    let q_s = mu / y_true + m_s;
    Ok(mu / q_s)
}

/// Tilman R* (minimum resource requirement) for a species.
///
/// `R* = K_s × D / (μ_max - D)`
///
/// The species with the lowest R* wins competition for a single limiting
/// resource in a chemostat. This is the resource-ratio theory of competition.
///
/// Reference: Tilman (1982) "Resource Competition and Community Structure".
///
/// # Errors
///
/// Returns error if dilution rate exceeds μ_max.
#[inline]
#[must_use = "returns the R* value without side effects"]
pub fn tilman_r_star(mu_max: f64, k_s: f64, dilution_rate: f64) -> Result<f64> {
    validate_positive(mu_max, "mu_max")?;
    validate_positive(k_s, "k_s")?;
    validate_positive(dilution_rate, "dilution_rate")?;
    if dilution_rate >= mu_max {
        return Err(JivanuError::SimulationFailed(
            "dilution rate must be less than mu_max".into(),
        ));
    }
    Ok(k_s * dilution_rate / (mu_max - dilution_rate))
}

/// Predict the winner of two-species resource competition (Tilman R* theory).
///
/// The species with the lower R* excludes the other when competing for
/// a single limiting resource. Returns `true` if species 1 wins.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns whether species 1 wins without side effects"]
pub fn r_star_competition_winner(
    mu_max_1: f64,
    k_s_1: f64,
    mu_max_2: f64,
    k_s_2: f64,
    dilution_rate: f64,
) -> Result<bool> {
    let r1 = tilman_r_star(mu_max_1, k_s_1, dilution_rate)?;
    let r2 = tilman_r_star(mu_max_2, k_s_2, dilution_rate)?;
    Ok(r1 < r2)
}

/// Cross-feeding growth rate: species B grows on the byproduct of species A.
///
/// Models syntrophy where species A produces a metabolite at rate proportional
/// to its growth, and species B consumes it via Monod kinetics.
///
/// `byproduct_conc = yield_ab × N_a × μ_a`
/// `μ_b = μ_max_b × byproduct_conc / (K_s_b + byproduct_conc)`
///
/// Returns the specific growth rate of species B.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the cross-feeding growth rate without side effects"]
pub fn cross_feeding_growth(
    n_a: f64,
    mu_a: f64,
    yield_ab: f64,
    mu_max_b: f64,
    k_s_b: f64,
) -> Result<f64> {
    validate_non_negative(n_a, "n_a")?;
    validate_non_negative(mu_a, "mu_a")?;
    validate_non_negative(yield_ab, "yield_ab")?;
    validate_positive(mu_max_b, "mu_max_b")?;
    validate_positive(k_s_b, "k_s_b")?;

    let byproduct = yield_ab * n_a * mu_a;
    monod_kinetics(byproduct, mu_max_b, k_s_b)
}

/// Multi-species chemostat: steady-state substrate concentration.
///
/// With N species competing for a single substrate in a chemostat,
/// at steady state only the species with the lowest R* persists
/// (competitive exclusion). This returns the steady-state substrate
/// level, which equals the winner's R*.
///
/// Returns `(winning_species_index, steady_state_substrate)`.
///
/// # Errors
///
/// Returns error if any species would wash out at the given dilution rate.
#[must_use = "returns the chemostat outcome without side effects"]
pub fn multi_species_chemostat(
    mu_max: &[f64],
    k_s: &[f64],
    dilution_rate: f64,
) -> Result<(usize, f64)> {
    if mu_max.len() != k_s.len() || mu_max.is_empty() {
        return Err(JivanuError::ComputationError(
            "mu_max and k_s must have the same non-zero length".into(),
        ));
    }
    validate_positive(dilution_rate, "dilution_rate")?;

    let mut best_idx = 0;
    let mut best_r_star = f64::MAX;

    for i in 0..mu_max.len() {
        let r = tilman_r_star(mu_max[i], k_s[i], dilution_rate)?;
        if r < best_r_star {
            best_r_star = r;
            best_idx = i;
        }
    }

    Ok((best_idx, best_r_star))
}

/// Stochastic birth-death event from the Gillespie algorithm.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct GillespieEvent {
    /// Time at which the event occurred.
    pub time: f64,
    /// Population after the event.
    pub population: u64,
}

/// Run a Gillespie stochastic simulation for a simple birth-death process.
///
/// Implements the direct method (Gillespie, 1977) for a population with:
/// - Birth rate: `birth_rate × N`
/// - Death rate: `death_rate × N`
///
/// Produces exact stochastic trajectories for discrete populations.
/// Uses a simple linear congruential RNG seeded by `seed`.
///
/// # Arguments
///
/// - `n0` — initial population count
/// - `birth_rate` — per-capita birth rate (1/time)
/// - `death_rate` — per-capita death rate (1/time)
/// - `t_end` — simulation end time
/// - `seed` — RNG seed for reproducibility
///
/// # Errors
///
/// Returns error if rates are negative or `t_end` is non-positive.
#[must_use = "returns the stochastic trajectory without side effects"]
pub fn gillespie_birth_death(
    n0: u64,
    birth_rate: f64,
    death_rate: f64,
    t_end: f64,
    seed: u64,
) -> Result<Vec<GillespieEvent>> {
    validate_non_negative(birth_rate, "birth_rate")?;
    validate_non_negative(death_rate, "death_rate")?;
    validate_positive(t_end, "t_end")?;

    let mut events = Vec::new();
    let mut t = 0.0;
    let mut n = n0;
    let mut rng_state = seed.wrapping_add(1); // avoid zero seed

    events.push(GillespieEvent {
        time: t,
        population: n,
    });

    while t < t_end && n > 0 {
        let total_rate = (birth_rate + death_rate) * n as f64;
        if total_rate <= 0.0 {
            break;
        }

        // Generate two uniform random numbers via LCG
        rng_state = rng_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        let u1 = (rng_state >> 33) as f64 / (1u64 << 31) as f64;
        rng_state = rng_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        let u2 = (rng_state >> 33) as f64 / (1u64 << 31) as f64;

        // Time to next event (exponential distribution)
        let u1_clamped = u1.max(1e-15); // avoid ln(0)
        let dt = -(u1_clamped.ln()) / total_rate;
        t += dt;

        if t > t_end {
            break;
        }

        // Determine event type
        let birth_prob = birth_rate * n as f64 / total_rate;
        if u2 < birth_prob {
            n += 1; // birth
        } else {
            n -= 1; // death
        }

        events.push(GillespieEvent {
            time: t,
            population: n,
        });
    }

    Ok(events)
}

/// N-strain Lotka-Volterra competition step.
///
/// Generalizes two-strain competition to an arbitrary number of strains
/// using an interaction matrix:
///
/// ```text
/// dN_i/dt = r_i * N_i * (1 - Σ_j(α_ij * N_j) / K_i)
/// ```
///
/// where `alpha[i][j]` is the competitive effect of strain j on strain i
/// (with `alpha[i][i] = 1.0` by convention, meaning self-limitation).
///
/// When the `hisab` feature is enabled, uses RK4 integration.
///
/// # Arguments
///
/// - `populations` — current density of each strain
/// - `growth_rates` — intrinsic growth rate of each strain
/// - `carrying_capacities` — carrying capacity of each strain
/// - `alpha` — interaction matrix (row-major, `n × n`)
/// - `dt` — time step
///
/// # Errors
///
/// Returns error if dimensions are inconsistent or parameters are invalid.
#[must_use = "returns the new populations without side effects"]
pub fn n_strain_competition_step(
    populations: &[f64],
    growth_rates: &[f64],
    carrying_capacities: &[f64],
    alpha: &[Vec<f64>],
    dt: f64,
) -> Result<Vec<f64>> {
    let n = populations.len();
    if growth_rates.len() != n || carrying_capacities.len() != n || alpha.len() != n {
        return Err(JivanuError::ComputationError(
            "all input arrays must have the same length".into(),
        ));
    }
    for (i, row) in alpha.iter().enumerate() {
        if row.len() != n {
            return Err(JivanuError::ComputationError(format!(
                "alpha row {i} has length {}, expected {n}",
                row.len()
            )));
        }
    }
    for i in 0..n {
        validate_non_negative(populations[i], "population")?;
        validate_positive(growth_rates[i], "growth_rate")?;
        validate_positive(carrying_capacities[i], "carrying_capacity")?;
        for val in &alpha[i] {
            validate_non_negative(*val, "alpha element")?;
        }
    }
    validate_positive(dt, "dt")?;

    n_strain_step_inner(populations, growth_rates, carrying_capacities, alpha, dt, n)
}

#[cfg(feature = "hisab")]
fn n_strain_step_inner(
    populations: &[f64],
    growth_rates: &[f64],
    carrying_capacities: &[f64],
    alpha: &[Vec<f64>],
    dt: f64,
    n: usize,
) -> Result<Vec<f64>> {
    let r = growth_rates.to_vec();
    let k = carrying_capacities.to_vec();
    let a: Vec<Vec<f64>> = alpha.to_vec();
    let result = hisab::num::rk4(
        move |_t, y, dy| {
            for i in 0..n {
                let mut sum = 0.0;
                for j in 0..n {
                    sum += a[i][j] * y[j];
                }
                dy[i] = r[i] * y[i] * (1.0 - sum / k[i]);
            }
        },
        0.0,
        populations,
        dt,
        1,
    )
    .map_err(|e| JivanuError::SimulationFailed(e.to_string()))?;
    Ok(result.into_iter().map(|v| v.max(0.0)).collect())
}

#[cfg(not(feature = "hisab"))]
fn n_strain_step_inner(
    populations: &[f64],
    growth_rates: &[f64],
    carrying_capacities: &[f64],
    alpha: &[Vec<f64>],
    dt: f64,
    n: usize,
) -> Result<Vec<f64>> {
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            sum += alpha[i][j] * populations[j];
        }
        let dn = growth_rates[i] * populations[i] * (1.0 - sum / carrying_capacities[i]);
        result.push((populations[i] + dn * dt).max(0.0));
    }
    Ok(result)
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

    // --- Two-strain competition tests ---

    #[test]
    fn test_competition_outcome_coexistence() {
        // Both alphas small: each inhibits itself more than the other
        let params = CompetitionParams {
            r1: 0.5,
            r2: 0.4,
            k1: 1000.0,
            k2: 800.0,
            alpha12: 0.5, // < K1/K2 = 1.25
            alpha21: 0.5, // < K2/K1 = 0.8
        };
        assert_eq!(
            competition_outcome(&params).unwrap(),
            CompetitionOutcome::Coexistence
        );
    }

    #[test]
    fn test_competition_outcome_strain1_wins() {
        let params = CompetitionParams {
            r1: 0.5,
            r2: 0.4,
            k1: 1000.0,
            k2: 800.0,
            alpha12: 0.5, // < K1/K2 = 1.25 → strain 1 tolerates
            alpha21: 1.5, // > K2/K1 = 0.8  → strain 2 cannot tolerate
        };
        assert_eq!(
            competition_outcome(&params).unwrap(),
            CompetitionOutcome::Strain1Wins
        );
    }

    #[test]
    fn test_competition_outcome_strain2_wins() {
        let params = CompetitionParams {
            r1: 0.5,
            r2: 0.4,
            k1: 1000.0,
            k2: 800.0,
            alpha12: 2.0, // > K1/K2 = 1.25 → strain 1 cannot tolerate
            alpha21: 0.5, // < K2/K1 = 0.8  → strain 2 tolerates
        };
        assert_eq!(
            competition_outcome(&params).unwrap(),
            CompetitionOutcome::Strain2Wins
        );
    }

    #[test]
    fn test_competition_outcome_unstable() {
        let params = CompetitionParams {
            r1: 0.5,
            r2: 0.4,
            k1: 1000.0,
            k2: 800.0,
            alpha12: 2.0, // > K1/K2
            alpha21: 1.5, // > K2/K1
        };
        assert_eq!(
            competition_outcome(&params).unwrap(),
            CompetitionOutcome::Unstable
        );
    }

    #[test]
    fn test_competition_step_populations_grow() {
        let state = TwoStrainState { n1: 10.0, n2: 10.0 };
        let params = CompetitionParams {
            r1: 0.5,
            r2: 0.4,
            k1: 1000.0,
            k2: 800.0,
            alpha12: 0.5,
            alpha21: 0.5,
        };
        let next = competition_step(&state, &params, 0.1).unwrap();
        assert!(next.n1 > state.n1);
        assert!(next.n2 > state.n2);
    }

    #[test]
    fn test_competition_step_at_capacity() {
        // Strain 1 at K, strain 2 absent → no growth
        let state = TwoStrainState {
            n1: 1000.0,
            n2: 0.0,
        };
        let params = CompetitionParams {
            r1: 0.5,
            r2: 0.4,
            k1: 1000.0,
            k2: 800.0,
            alpha12: 0.5,
            alpha21: 0.5,
        };
        let next = competition_step(&state, &params, 0.1).unwrap();
        assert!((next.n1 - 1000.0).abs() < 1.0);
        assert!((next.n2 - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_competition_exclusion_long_run() {
        // Strain 1 wins scenario: run many steps, strain 2 should decline
        let params = CompetitionParams {
            r1: 0.5,
            r2: 0.4,
            k1: 1000.0,
            k2: 800.0,
            alpha12: 0.5,
            alpha21: 1.5,
        };
        let mut state = TwoStrainState {
            n1: 100.0,
            n2: 100.0,
        };
        for _ in 0..10000 {
            state = competition_step(&state, &params, 0.01).unwrap();
        }
        // Strain 1 should dominate near K1, strain 2 near 0
        assert!(state.n1 > 900.0, "n1 = {}", state.n1);
        assert!(state.n2 < 10.0, "n2 = {}", state.n2);
    }

    #[test]
    fn test_two_strain_state_serde_roundtrip() {
        let state = TwoStrainState {
            n1: 500.0,
            n2: 300.0,
        };
        let json = serde_json::to_string(&state).unwrap();
        let back: TwoStrainState = serde_json::from_str(&json).unwrap();
        assert!((state.n1 - back.n1).abs() < 1e-10);
        assert!((state.n2 - back.n2).abs() < 1e-10);
    }

    #[test]
    fn test_competition_outcome_serde_roundtrip() {
        let co = CompetitionOutcome::Coexistence;
        let json = serde_json::to_string(&co).unwrap();
        let back: CompetitionOutcome = serde_json::from_str(&json).unwrap();
        assert_eq!(co, back);
    }

    // --- N-strain competition tests ---

    #[test]
    fn test_n_strain_single_species_logistic() {
        // 1 strain with alpha=1.0 should behave like logistic growth
        let pops = [100.0];
        let rates = [0.5];
        let caps = [1000.0];
        let alpha = [vec![1.0]];
        let next = n_strain_competition_step(&pops, &rates, &caps, &alpha, 0.1).unwrap();
        assert!(next[0] > 100.0); // should grow
        assert!(next[0] < 1000.0); // below capacity
    }

    #[test]
    fn test_n_strain_matches_two_strain() {
        // N-strain with n=2 should match two-strain API
        let state = TwoStrainState {
            n1: 100.0,
            n2: 50.0,
        };
        let params = CompetitionParams {
            r1: 0.5,
            r2: 0.4,
            k1: 1000.0,
            k2: 800.0,
            alpha12: 0.6,
            alpha21: 0.3,
        };
        let two = competition_step(&state, &params, 0.1).unwrap();
        let n_result = n_strain_competition_step(
            &[100.0, 50.0],
            &[0.5, 0.4],
            &[1000.0, 800.0],
            &[vec![1.0, 0.6], vec![0.3, 1.0]],
            0.1,
        )
        .unwrap();
        assert!((two.n1 - n_result[0]).abs() < 1e-6);
        assert!((two.n2 - n_result[1]).abs() < 1e-6);
    }

    #[test]
    fn test_n_strain_three_species() {
        let pops = [100.0, 100.0, 100.0];
        let rates = [0.5, 0.4, 0.3];
        let caps = [1000.0, 800.0, 600.0];
        let alpha = [
            vec![1.0, 0.5, 0.3],
            vec![0.4, 1.0, 0.2],
            vec![0.3, 0.3, 1.0],
        ];
        let next = n_strain_competition_step(&pops, &rates, &caps, &alpha, 0.1).unwrap();
        assert_eq!(next.len(), 3);
        for p in &next {
            assert!(*p > 0.0);
        }
    }

    #[test]
    fn test_n_strain_dimension_mismatch() {
        assert!(
            n_strain_competition_step(
                &[100.0, 50.0],
                &[0.5],
                &[1000.0, 800.0],
                &[vec![1.0, 0.5], vec![0.5, 1.0]],
                0.1,
            )
            .is_err()
        );
    }

    #[test]
    fn test_n_strain_alpha_row_mismatch() {
        assert!(
            n_strain_competition_step(
                &[100.0, 50.0],
                &[0.5, 0.4],
                &[1000.0, 800.0],
                &[vec![1.0, 0.5], vec![0.5]], // row 1 too short
                0.1,
            )
            .is_err()
        );
    }

    #[test]
    fn test_competition_params_serde_roundtrip() {
        let params = CompetitionParams {
            r1: 0.5,
            r2: 0.4,
            k1: 1000.0,
            k2: 800.0,
            alpha12: 0.5,
            alpha21: 0.5,
        };
        let json = serde_json::to_string(&params).unwrap();
        let back: CompetitionParams = serde_json::from_str(&json).unwrap();
        assert!((params.r1 - back.r1).abs() < 1e-10);
    }

    // --- Gillespie stochastic tests ---

    #[test]
    fn test_gillespie_pure_birth() {
        // Only births: population should grow
        let events = gillespie_birth_death(10, 1.0, 0.0, 1.0, 42).unwrap();
        assert!(events.len() >= 2);
        assert_eq!(events[0].population, 10);
        assert!(events.last().unwrap().population > 10);
    }

    #[test]
    fn test_gillespie_pure_death() {
        // Only deaths: population should decline toward 0
        let events = gillespie_birth_death(100, 0.0, 1.0, 10.0, 42).unwrap();
        assert!(events.last().unwrap().population < 100);
    }

    #[test]
    fn test_gillespie_extinction() {
        // High death rate, small population: should go extinct
        let events = gillespie_birth_death(5, 0.0, 10.0, 10.0, 42).unwrap();
        assert_eq!(events.last().unwrap().population, 0);
    }

    #[test]
    fn test_gillespie_reproducible() {
        // Same seed → same trajectory
        let a = gillespie_birth_death(50, 0.5, 0.1, 2.0, 123).unwrap();
        let b = gillespie_birth_death(50, 0.5, 0.1, 2.0, 123).unwrap();
        assert_eq!(a.len(), b.len());
        for (ea, eb) in a.iter().zip(b.iter()) {
            assert_eq!(ea.population, eb.population);
            assert!((ea.time - eb.time).abs() < 1e-10);
        }
    }

    #[test]
    fn test_gillespie_different_seeds_differ() {
        let a = gillespie_birth_death(50, 0.5, 0.1, 5.0, 1).unwrap();
        let b = gillespie_birth_death(50, 0.5, 0.1, 5.0, 999).unwrap();
        // Very unlikely to produce identical trajectories
        let same = a.len() == b.len()
            && a.iter()
                .zip(b.iter())
                .all(|(ea, eb)| ea.population == eb.population);
        assert!(
            !same,
            "different seeds should produce different trajectories"
        );
    }

    #[test]
    fn test_gillespie_zero_population() {
        // Start at 0: nothing happens
        let events = gillespie_birth_death(0, 1.0, 0.5, 1.0, 42).unwrap();
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].population, 0);
    }

    #[test]
    fn test_gillespie_time_monotonic() {
        let events = gillespie_birth_death(100, 0.5, 0.1, 5.0, 42).unwrap();
        for i in 1..events.len() {
            assert!(events[i].time >= events[i - 1].time);
        }
    }

    #[test]
    fn test_gillespie_event_serde_roundtrip() {
        let ev = GillespieEvent {
            time: 1.5,
            population: 42,
        };
        let json = serde_json::to_string(&ev).unwrap();
        let back: GillespieEvent = serde_json::from_str(&json).unwrap();
        assert!((ev.time - back.time).abs() < 1e-10);
        assert_eq!(ev.population, back.population);
    }

    // --- Baranyi-Roberts / Gompertz tests ---

    #[test]
    fn test_baranyi_roberts_lag_phase() {
        // During lag, population stays at N0
        let n = baranyi_roberts(100.0, 1e9, 1.0, 2.0, 1.0).unwrap();
        assert!((n - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_baranyi_roberts_after_lag() {
        let n = baranyi_roberts(100.0, 1e9, 1.0, 2.0, 5.0).unwrap();
        assert!(n > 100.0, "should grow after lag");
    }

    #[test]
    fn test_baranyi_roberts_approaches_nmax() {
        let n = baranyi_roberts(100.0, 1e6, 1.0, 1.0, 100.0).unwrap();
        assert!((n - 1e6).abs() < 1.0);
    }

    #[test]
    fn test_gompertz_at_t0() {
        // With long lag, Gompertz at t=0 ≈ N0
        let n = gompertz(100.0, 1e9, 0.5, 10.0, 0.0).unwrap();
        assert!((90.0..=110.0).contains(&n), "n at t=0 = {n}");
    }

    #[test]
    fn test_gompertz_grows() {
        let n0 = gompertz(100.0, 1e6, 1.0, 1.0, 0.0).unwrap();
        let n10 = gompertz(100.0, 1e6, 1.0, 1.0, 10.0).unwrap();
        assert!(n10 > n0, "should grow over time");
    }

    #[test]
    fn test_gompertz_approaches_nmax() {
        let n = gompertz(100.0, 1e6, 1.0, 1.0, 100.0).unwrap();
        assert!((n - 1e6).abs() / 1e6 < 0.01);
    }

    // --- Temperature / pH tests ---

    #[test]
    fn test_ratkowsky_below_tmin() {
        let mu = ratkowsky_temperature(0.04, 0.0, 5.0).unwrap();
        assert!((mu - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_ratkowsky_above_tmin() {
        let mu = ratkowsky_temperature(0.04, 37.0, 5.0).unwrap();
        // b=0.04, T-Tmin=32 → μ = 0.0016 * 1024 = 1.6384
        assert!((mu - 0.04_f64.powi(2) * 32.0_f64.powi(2)).abs() < 1e-10);
    }

    #[test]
    fn test_cardinal_temperature_at_optimum() {
        let ctm = cardinal_temperature(37.0, 5.0, 37.0, 45.0).unwrap();
        assert!(
            (ctm - 1.0).abs() < 1e-6,
            "should be 1.0 at T_opt, got {ctm}"
        );
    }

    #[test]
    fn test_cardinal_temperature_at_extremes() {
        assert!((cardinal_temperature(5.0, 5.0, 37.0, 45.0).unwrap() - 0.0).abs() < 1e-10);
        assert!((cardinal_temperature(45.0, 5.0, 37.0, 45.0).unwrap() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_cardinal_temperature_bad_order() {
        assert!(cardinal_temperature(37.0, 37.0, 5.0, 45.0).is_err());
    }

    #[test]
    fn test_cardinal_ph_at_optimum() {
        let cpm = cardinal_ph(7.0, 4.0, 7.0, 9.0).unwrap();
        assert!(
            (cpm - 1.0).abs() < 1e-6,
            "should be 1.0 at pH_opt, got {cpm}"
        );
    }

    #[test]
    fn test_cardinal_ph_at_extremes() {
        assert!((cardinal_ph(4.0, 4.0, 7.0, 9.0).unwrap() - 0.0).abs() < 1e-10);
        assert!((cardinal_ph(9.0, 4.0, 7.0, 9.0).unwrap() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_cardinal_ph_bad_order() {
        assert!(cardinal_ph(7.0, 7.0, 4.0, 9.0).is_err());
    }

    // --- Maintenance energy tests ---

    #[test]
    fn test_maintenance_substrate_rate() {
        // μ=0.5, Y_true=0.5, m_s=0.1 → q_s = 0.5/0.5 + 0.1 = 1.1
        let q_s = maintenance_substrate_rate(0.5, 0.5, 0.1).unwrap();
        assert!((q_s - 1.1).abs() < 1e-10);
    }

    #[test]
    fn test_maintenance_zero_growth() {
        // At μ=0, only maintenance: q_s = m_s
        let q_s = maintenance_substrate_rate(0.0, 0.5, 0.1).unwrap();
        assert!((q_s - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_apparent_yield_no_maintenance() {
        // m_s=0 → Y_app = Y_true
        let y_app = apparent_yield(0.5, 0.4, 0.0).unwrap();
        assert!((y_app - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_apparent_yield_with_maintenance() {
        // With maintenance, Y_app < Y_true
        let y_app = apparent_yield(0.5, 0.4, 0.1).unwrap();
        assert!(y_app < 0.4);
        assert!(y_app > 0.0);
    }

    // --- Tilman R* / ecology tests ---

    #[test]
    fn test_tilman_r_star() {
        // Same formula as chemostat steady-state substrate
        let r = tilman_r_star(1.0, 0.5, 0.2).unwrap();
        assert!((r - 0.5 * 0.2 / (1.0 - 0.2)).abs() < 1e-10);
    }

    #[test]
    fn test_tilman_r_star_washout() {
        assert!(tilman_r_star(0.5, 0.5, 0.6).is_err());
    }

    #[test]
    fn test_r_star_winner() {
        // Species 1 has lower K_s → lower R* → wins
        let wins = r_star_competition_winner(1.0, 0.1, 1.0, 0.5, 0.2).unwrap();
        assert!(wins, "species 1 (lower K_s) should win");
    }

    #[test]
    fn test_cross_feeding_no_producer() {
        // N_a=0 → no byproduct → no growth for B
        let mu_b = cross_feeding_growth(0.0, 0.5, 0.1, 1.0, 0.5).unwrap();
        assert!((mu_b - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_cross_feeding_with_producer() {
        let mu_b = cross_feeding_growth(1000.0, 0.5, 0.1, 1.0, 0.5).unwrap();
        assert!(mu_b > 0.0, "B should grow on A's byproduct");
    }

    #[test]
    fn test_multi_species_chemostat() {
        // 3 species, species 2 has lowest K_s → wins
        let (winner, r_star) =
            multi_species_chemostat(&[1.0, 1.0, 1.0], &[0.5, 0.1, 0.3], 0.2).unwrap();
        assert_eq!(winner, 1, "species with lowest K_s should win");
        let expected = tilman_r_star(1.0, 0.1, 0.2).unwrap();
        assert!((r_star - expected).abs() < 1e-10);
    }

    #[test]
    fn test_multi_species_chemostat_empty() {
        assert!(multi_species_chemostat(&[], &[], 0.2).is_err());
    }
}
