//! Stochastic simulation — Gillespie SSA and tau-leaping for reaction networks.

use serde::{Deserialize, Serialize};

use crate::error::{JivanuError, Result, validate_positive};

/// A stochastic reaction: stoichiometry changes and propensity function index.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StochasticReaction {
    /// Reaction identifier.
    pub id: String,
    /// State change vector: `(species_index, change)` pairs.
    pub state_change: Vec<(usize, i64)>,
}

/// Snapshot of a stochastic simulation at a point in time.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StochasticState {
    /// Current simulation time.
    pub time: f64,
    /// Species populations.
    pub populations: Vec<u64>,
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
/// Returns error if inputs are inconsistent.
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
    if reactions.is_empty() {
        return Err(JivanuError::ComputationError(
            "at least one reaction required".into(),
        ));
    }

    let n_species = initial.len();
    let mut state: Vec<u64> = initial.to_vec();
    let mut t = 0.0;
    let mut rng = seed.wrapping_add(1);
    let mut trajectory = Vec::new();

    trajectory.push(StochasticState {
        time: t,
        populations: state.clone(),
    });

    for _ in 0..max_events {
        let props = propensity(&state);
        let total: f64 = props.iter().sum();
        if total <= 0.0 || t >= t_end {
            break;
        }

        // Time to next event
        rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let u1 = ((rng >> 33) as f64 / (1u64 << 31) as f64).max(1e-15);
        let dt = -u1.ln() / total;
        t += dt;

        if t > t_end {
            break;
        }

        // Select reaction
        rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let u2 = (rng >> 33) as f64 / (1u64 << 31) as f64;
        let threshold = u2 * total;
        let mut cumulative = 0.0;
        let mut selected = reactions.len() - 1;
        for (i, p) in props.iter().enumerate() {
            cumulative += p;
            if cumulative >= threshold {
                selected = i;
                break;
            }
        }

        // Apply state change
        for &(species, change) in &reactions[selected].state_change {
            if species < n_species {
                if change >= 0 {
                    state[species] = state[species].saturating_add(change as u64);
                } else {
                    state[species] = state[species].saturating_sub((-change) as u64);
                }
            }
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
/// Returns error if inputs are inconsistent.
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
    if reactions.is_empty() {
        return Err(JivanuError::ComputationError(
            "at least one reaction required".into(),
        ));
    }

    let n_species = initial.len();
    let mut state: Vec<u64> = initial.to_vec();
    let mut t = 0.0;
    let mut rng = seed.wrapping_add(1);
    let mut trajectory = Vec::new();

    trajectory.push(StochasticState {
        time: t,
        populations: state.clone(),
    });

    while t < t_end {
        let props = propensity(&state);
        let step = tau.min(t_end - t);

        // For each reaction, draw Poisson(propensity * tau) firings
        for (rxn_idx, &prop) in props.iter().enumerate() {
            let lambda = prop * step;
            if lambda <= 0.0 {
                continue;
            }
            // Poisson draw via Knuth's algorithm (good for small lambda)
            let firings = poisson_draw(lambda, &mut rng);
            for _ in 0..firings {
                for &(species, change) in &reactions[rxn_idx].state_change {
                    if species < n_species {
                        if change >= 0 {
                            state[species] = state[species].saturating_add(change as u64);
                        } else {
                            state[species] = state[species].saturating_sub((-change) as u64);
                        }
                    }
                }
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

/// Poisson random variate via Knuth's algorithm.
fn poisson_draw(lambda: f64, rng: &mut u64) -> u64 {
    let l = (-lambda).exp();
    let mut k = 0u64;
    let mut p = 1.0;
    loop {
        k += 1;
        *rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let u = (*rng >> 33) as f64 / (1u64 << 31) as f64;
        p *= u;
        if p <= l {
            break;
        }
    }
    k.saturating_sub(1)
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
}
