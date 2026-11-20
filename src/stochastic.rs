//! Stochastic simulation — Gillespie SSA and tau-leaping for reaction networks.
//!
//! Both simulators are seeded and reproducible: the same inputs and seed
//! always produce the same trajectory.

use serde::{Deserialize, Serialize};

use crate::error::{JivanuError, Result, validate_positive};
use crate::rng::Lcg;

/// A stochastic reaction: stoichiometry changes and propensity function index.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct StochasticReaction {
    /// Reaction identifier.
    pub id: String,
    /// State change vector: `(species_index, change)` pairs.
    pub state_change: Vec<(usize, i64)>,
}

/// Snapshot of a stochastic simulation at a point in time.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StochasticState {
    /// Current simulation time.
    pub time: f64,
    /// Species populations.
    pub populations: Vec<u64>,
}

/// Validate the reaction list against the number of species.
fn validate_reactions(reactions: &[StochasticReaction], n_species: usize) -> Result<()> {
    if reactions.is_empty() {
        return Err(JivanuError::ComputationError(
            "at least one reaction required".into(),
        ));
    }
    for rxn in reactions {
        if let Some(&(species, _)) = rxn.state_change.iter().find(|(s, _)| *s >= n_species) {
            return Err(JivanuError::ComputationError(format!(
                "reaction '{}' references species {species}, but only {n_species} exist",
                rxn.id
            )));
        }
    }
    Ok(())
}

/// Evaluate propensities and validate them; returns their sum.
fn checked_propensities(
    propensity: &impl Fn(&[u64]) -> Vec<f64>,
    state: &[u64],
    n_reactions: usize,
) -> Result<(Vec<f64>, f64)> {
    let props = propensity(state);
    if props.len() != n_reactions {
        return Err(JivanuError::ComputationError(format!(
            "propensity function returned {} values for {n_reactions} reactions",
            props.len()
        )));
    }
    // Validate and sum in one pass (`p >= 0.0` is false for NaN).
    let mut total = 0.0;
    for &p in &props {
        if !(p >= 0.0 && p.is_finite()) {
            return Err(JivanuError::ComputationError(format!(
                "propensities must be finite and non-negative, got {p}"
            )));
        }
        total += p;
    }
    Ok((props, total))
}

/// Apply a reaction's state change `firings` times (saturating at 0 and `u64::MAX`).
#[inline]
fn apply_reaction(state: &mut [u64], reaction: &StochasticReaction, firings: u64) {
    for &(species, change) in &reaction.state_change {
        if let Some(pop) = state.get_mut(species) {
            let delta = change.unsigned_abs().saturating_mul(firings);
            *pop = if change >= 0 {
                pop.saturating_add(delta)
            } else {
                pop.saturating_sub(delta)
            };
        }
    }
}

/// Run the Gillespie direct method (SSA) for an arbitrary reaction network.
///
/// # Arguments
///
/// - `initial` — initial species populations
/// - `reactions` — list of reactions with state change vectors
/// - `propensity` — function computing propensity for each reaction given current state.
///   Returns a Vec of propensities (one per reaction).
/// - `t_end` — simulation end time
/// - `max_events` — maximum number of events (prevents runaway)
/// - `seed` — RNG seed
///
/// # Errors
///
/// Returns error if `t_end` is non-positive, there are no reactions, a
/// reaction references a species index outside `initial`, or the propensity
/// function returns the wrong number of values or a negative / non-finite
/// propensity.
#[must_use = "returns the stochastic trajectory without side effects"]
pub fn gillespie_ssa(
    initial: &[u64],
    reactions: &[StochasticReaction],
    propensity: impl Fn(&[u64]) -> Vec<f64>,
    t_end: f64,
    max_events: usize,
    seed: u64,
) -> Result<Vec<StochasticState>> {
    validate_positive(t_end, "t_end")?;
    validate_reactions(reactions, initial.len())?;

    let mut state: Vec<u64> = initial.to_vec();
    let mut t = 0.0;
    let mut rng = Lcg::new(seed);
    let mut trajectory = vec![StochasticState {
        time: t,
        populations: state.clone(),
    }];

    for _ in 0..max_events {
        let (props, total) = checked_propensities(&propensity, &state, reactions.len())?;
        if total <= 0.0 || t >= t_end {
            break;
        }

        // Time to next event
        let u1 = rng.next_f64().max(1e-15);
        t += -u1.ln() / total;
        if t > t_end {
            break;
        }

        // Select reaction: first index whose cumulative propensity exceeds
        // u2 × total. The comparison is strict, so a zero-propensity
        // reaction (which leaves the sum unchanged) is never selected, even
        // when u2 = 0. Rounding fallback: the last reaction that can fire.
        let threshold = rng.next_f64() * total;
        let mut cumulative = 0.0;
        let selected = props
            .iter()
            .position(|&p| {
                cumulative += p;
                cumulative > threshold
            })
            .or_else(|| props.iter().rposition(|&p| p > 0.0));
        if let Some(rxn) = selected.and_then(|i| reactions.get(i)) {
            apply_reaction(&mut state, rxn, 1);
        }

        trajectory.push(StochasticState {
            time: t,
            populations: state.clone(),
        });
    }

    Ok(trajectory)
}

/// Tau-leaping approximate stochastic simulation.
///
/// Instead of simulating one reaction at a time, fires Poisson-distributed
/// numbers of each reaction in a fixed time step `tau`. Much faster than
/// exact SSA for systems with large populations.
///
/// Firings of all reactions within a leap are applied together; populations
/// saturate at zero rather than going negative.
///
/// # Arguments
///
/// - `initial` — initial species populations
/// - `reactions` — reaction list
/// - `propensity` — propensity function
/// - `tau` — leap size (time step)
/// - `t_end` — simulation end time
/// - `seed` — RNG seed
///
/// # Errors
///
/// Returns error if `tau` or `t_end` is non-positive, there are no reactions,
/// a reaction references a species index outside `initial`, or the propensity
/// function returns the wrong number of values or a negative / non-finite
/// propensity.
#[must_use = "returns the stochastic trajectory without side effects"]
pub fn tau_leaping(
    initial: &[u64],
    reactions: &[StochasticReaction],
    propensity: impl Fn(&[u64]) -> Vec<f64>,
    tau: f64,
    t_end: f64,
    seed: u64,
) -> Result<Vec<StochasticState>> {
    validate_positive(tau, "tau")?;
    validate_positive(t_end, "t_end")?;
    validate_reactions(reactions, initial.len())?;

    let mut state: Vec<u64> = initial.to_vec();
    let mut t = 0.0;
    let mut rng = Lcg::new(seed);
    let mut trajectory = vec![StochasticState {
        time: t,
        populations: state.clone(),
    }];
    let mut firings = vec![0u64; reactions.len()];

    while t < t_end {
        let (props, _) = checked_propensities(&propensity, &state, reactions.len())?;
        let step = tau.min(t_end - t);

        // Draw Poisson(propensity × step) firings for every reaction from the
        // same pre-leap state, then apply them together.
        for (count, &prop) in firings.iter_mut().zip(&props) {
            let lambda = prop * step;
            *count = if lambda > 0.0 { rng.poisson(lambda) } else { 0 };
        }
        for (rxn, &count) in reactions.iter().zip(&firings) {
            if count > 0 {
                apply_reaction(&mut state, rxn, count);
            }
        }

        t += step;
        trajectory.push(StochasticState {
            time: t,
            populations: state.clone(),
        });
    }

    Ok(trajectory)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn birth_death_reactions() -> Vec<StochasticReaction> {
        vec![
            StochasticReaction {
                id: "birth".into(),
                state_change: vec![(0, 1)],
            },
            StochasticReaction {
                id: "death".into(),
                state_change: vec![(0, -1)],
            },
        ]
    }

    fn birth_death_propensity(birth_rate: f64, death_rate: f64) -> impl Fn(&[u64]) -> Vec<f64> {
        move |state: &[u64]| vec![birth_rate * state[0] as f64, death_rate * state[0] as f64]
    }

    #[test]
    fn test_gillespie_ssa_pure_birth() {
        let rxns = birth_death_reactions();
        let traj = gillespie_ssa(
            &[10],
            &rxns,
            birth_death_propensity(1.0, 0.0),
            1.0,
            10000,
            42,
        )
        .unwrap();
        assert!(traj.last().unwrap().populations[0] > 10);
    }

    #[test]
    fn test_gillespie_ssa_pure_death() {
        let rxns = birth_death_reactions();
        let traj = gillespie_ssa(
            &[100],
            &rxns,
            birth_death_propensity(0.0, 1.0),
            10.0,
            100000,
            42,
        )
        .unwrap();
        assert!(traj.last().unwrap().populations[0] < 100);
    }

    #[test]
    fn test_gillespie_ssa_reproducible() {
        let rxns = birth_death_reactions();
        let a = gillespie_ssa(
            &[50],
            &rxns,
            birth_death_propensity(0.5, 0.1),
            2.0,
            10000,
            123,
        )
        .unwrap();
        let b = gillespie_ssa(
            &[50],
            &rxns,
            birth_death_propensity(0.5, 0.1),
            2.0,
            10000,
            123,
        )
        .unwrap();
        assert_eq!(a.len(), b.len());
    }

    #[test]
    fn test_gillespie_ssa_two_species() {
        // Predator-prey: prey birth, prey eaten, predator death
        let rxns = vec![
            StochasticReaction {
                id: "prey_birth".into(),
                state_change: vec![(0, 1)],
            },
            StochasticReaction {
                id: "predation".into(),
                state_change: vec![(0, -1), (1, 1)],
            },
            StochasticReaction {
                id: "predator_death".into(),
                state_change: vec![(1, -1)],
            },
        ];
        let prop = |state: &[u64]| {
            vec![
                0.5 * state[0] as f64,                    // prey birth
                0.01 * state[0] as f64 * state[1] as f64, // predation
                0.3 * state[1] as f64,                    // predator death
            ]
        };
        let traj = gillespie_ssa(&[100, 20], &rxns, prop, 5.0, 100000, 42).unwrap();
        assert!(traj.len() > 1);
        assert_eq!(traj[0].populations.len(), 2);
    }

    #[test]
    fn test_gillespie_ssa_empty_reactions() {
        assert!(gillespie_ssa(&[10], &[], |_| vec![], 1.0, 100, 42).is_err());
    }

    #[test]
    fn test_tau_leaping_growth() {
        let rxns = birth_death_reactions();
        let traj = tau_leaping(
            &[100],
            &rxns,
            birth_death_propensity(0.5, 0.1),
            0.01,
            1.0,
            42,
        )
        .unwrap();
        assert!(traj.last().unwrap().populations[0] > 100);
    }

    #[test]
    fn test_tau_leaping_decline() {
        let rxns = birth_death_reactions();
        let traj = tau_leaping(
            &[100],
            &rxns,
            birth_death_propensity(0.0, 1.0),
            0.01,
            1.0,
            42,
        )
        .unwrap();
        assert!(traj.last().unwrap().populations[0] < 100);
    }

    #[test]
    fn test_tau_leaping_trajectory_length() {
        let rxns = birth_death_reactions();
        let traj = tau_leaping(
            &[100],
            &rxns,
            birth_death_propensity(0.5, 0.1),
            0.1,
            1.0,
            42,
        )
        .unwrap();
        // ~10 steps + initial; floating-point may produce 10 or 11 steps
        assert!(
            (10..=12).contains(&traj.len()),
            "expected ~11, got {}",
            traj.len()
        );
    }

    #[test]
    fn test_stochastic_state_serde_roundtrip() {
        let state = StochasticState {
            time: 1.5,
            populations: vec![42, 17],
        };
        let json = serde_json::to_string(&state).unwrap();
        let back: StochasticState = serde_json::from_str(&json).unwrap();
        assert!((state.time - back.time).abs() < 1e-10);
        assert_eq!(state.populations, back.populations);
    }

    #[test]
    fn test_stochastic_reaction_serde_roundtrip() {
        let rxn = StochasticReaction {
            id: "test".into(),
            state_change: vec![(0, 1), (1, -1)],
        };
        let json = serde_json::to_string(&rxn).unwrap();
        let back: StochasticReaction = serde_json::from_str(&json).unwrap();
        assert_eq!(rxn.id, back.id);
    }
    // --- Hardening regressions ---

    #[test]
    fn test_ssa_propensity_length_mismatch_is_error() {
        // A propensity vector longer than the reaction list used to index
        // out of bounds.
        let rxns = birth_death_reactions();
        let r = gillespie_ssa(&[10], &rxns, |_| vec![1.0, 1.0, 1.0], 1.0, 100, 1);
        assert!(r.is_err());
        let r = tau_leaping(&[10], &rxns, |_| vec![1.0], 0.1, 1.0, 1);
        assert!(r.is_err());
    }

    #[test]
    fn test_invalid_propensities_rejected() {
        let rxns = birth_death_reactions();
        assert!(gillespie_ssa(&[10], &rxns, |_| vec![-1.0, 1.0], 1.0, 100, 1).is_err());
        assert!(gillespie_ssa(&[10], &rxns, |_| vec![f64::NAN, 1.0], 1.0, 100, 1).is_err());
        assert!(tau_leaping(&[10], &rxns, |_| vec![f64::INFINITY, 1.0], 0.1, 1.0, 1).is_err());
    }

    #[test]
    fn test_species_index_out_of_range_rejected() {
        let rxns = vec![StochasticReaction {
            id: "bad".into(),
            state_change: vec![(3, 1)],
        }];
        assert!(gillespie_ssa(&[1], &rxns, |_| vec![1.0], 1.0, 10, 1).is_err());
        assert!(tau_leaping(&[1], &rxns, |_| vec![1.0], 0.1, 1.0, 1).is_err());
    }

    #[test]
    fn test_ssa_never_fires_zero_propensity_reaction() {
        let rxns = birth_death_reactions();
        for seed in 0..50 {
            let traj =
                gillespie_ssa(&[5], &rxns, |s| vec![0.0, s[0] as f64], 10.0, 1000, seed).unwrap();
            for w in traj.windows(2) {
                assert!(w[1].populations[0] <= w[0].populations[0], "birth fired");
            }
        }
    }

    #[test]
    fn test_ssa_pure_death_mean_matches_analytic() {
        // E[N(t)] = N0·e^(-d·t) for linear death.
        let rxns = birth_death_reactions();
        let runs = 400;
        let mut total = 0.0;
        for seed in 0..runs {
            let traj = gillespie_ssa(
                &[100],
                &rxns,
                birth_death_propensity(0.0, 1.0),
                1.0,
                10_000,
                seed,
            )
            .unwrap();
            total += traj.last().unwrap().populations[0] as f64;
        }
        let mean = total / f64::from(runs as u32);
        let expected = 100.0 * (-1.0_f64).exp();
        assert!((mean - expected).abs() < 1.5, "mean {mean} vs {expected}");
    }

    #[test]
    fn test_tau_leaping_large_population_mean() {
        // Pure death with N0 = 1e6: Poisson draws with λ ≫ 745 must not be
        // capped. After t = 1, E[N] ≈ 1e6·e^-1 (tau-leaping bias is O(tau)).
        let rxns = birth_death_reactions();
        let traj = tau_leaping(
            &[1_000_000],
            &rxns,
            birth_death_propensity(0.0, 1.0),
            0.001,
            1.0,
            9,
        )
        .unwrap();
        let n = traj.last().unwrap().populations[0] as f64;
        let expected = 1e6 * (-1.0_f64).exp();
        assert!((n - expected).abs() / expected < 0.01, "n = {n}");
    }

    #[test]
    fn test_tau_leaping_reproducible() {
        let rxns = birth_death_reactions();
        let a = tau_leaping(
            &[500],
            &rxns,
            birth_death_propensity(0.6, 0.4),
            0.05,
            2.0,
            77,
        )
        .unwrap();
        let b = tau_leaping(
            &[500],
            &rxns,
            birth_death_propensity(0.6, 0.4),
            0.05,
            2.0,
            77,
        )
        .unwrap();
        assert_eq!(a, b);
    }
}
