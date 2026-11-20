//! Metabolism — enzyme kinetics, metabolic pathways, ATP yields.

use serde::{Deserialize, Serialize};

use crate::error::{
    JivanuError, Result, validate_finite, validate_non_negative, validate_positive,
};

/// Michaelis-Menten enzyme kinetics.
///
/// `v = V_max * [S] / (K_m + [S])`
///
/// At `[S] = K_m`, `v = V_max / 2`.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the reaction rate without side effects"]
pub fn michaelis_menten(substrate: f64, v_max: f64, k_m: f64) -> Result<f64> {
    validate_non_negative(substrate, "substrate")?;
    validate_positive(v_max, "v_max")?;
    validate_positive(k_m, "k_m")?;
    Ok(v_max * substrate / (k_m + substrate))
}

/// Lineweaver-Burk double-reciprocal coordinates of one observation.
///
/// Transforms a measured rate `v` at substrate concentration `[S]` into the
/// point `(1/[S], 1/v)` of a Lineweaver-Burk plot, on which Michaelis-Menten
/// data fall on the line `1/v = (K_m / V_max) · (1/[S]) + 1/V_max`
/// (see [`lineweaver_burk_line`]).
///
/// Returns `(x = 1/[S], y = 1/v)`.
///
/// # Errors
///
/// Returns error if v or s is non-positive.
#[inline]
#[must_use = "returns the plot coordinates without side effects"]
pub fn lineweaver_burk(v: f64, s: f64) -> Result<(f64, f64)> {
    validate_positive(v, "v")?;
    validate_positive(s, "s")?;
    Ok((1.0 / s, 1.0 / v))
}

/// Lineweaver-Burk line for given Michaelis-Menten parameters.
///
/// `1/v = (K_m / V_max) · (1/[S]) + 1/V_max`
///
/// Returns `(slope = K_m/V_max, intercept = 1/V_max)`. The x-intercept is
/// `-1/K_m`.
///
/// # Errors
///
/// Returns error if `v_max` or `k_m` is non-positive.
#[inline]
#[must_use = "returns the line parameters without side effects"]
pub fn lineweaver_burk_line(v_max: f64, k_m: f64) -> Result<(f64, f64)> {
    validate_positive(v_max, "v_max")?;
    validate_positive(k_m, "k_m")?;
    Ok((k_m / v_max, 1.0 / v_max))
}

/// ATP yield from glycolysis.
///
/// Net yield: 2 ATP per glucose molecule (substrate-level phosphorylation).
#[inline]
#[must_use]
pub const fn glycolysis_atp() -> u32 {
    2
}

/// ATP yield from oxidative phosphorylation.
///
/// Approximately 34 ATP per glucose (via electron transport chain + chemiosmosis).
#[inline]
#[must_use]
pub const fn oxidative_phosphorylation_atp() -> u32 {
    34
}

/// Total aerobic respiration ATP yield per glucose.
///
/// Glycolysis (2) + Krebs cycle (2) + Oxidative phosphorylation (34) = ~38.
#[inline]
#[must_use]
pub const fn total_aerobic_atp() -> u32 {
    38
}

/// Enzyme inhibition types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum InhibitionType {
    /// Inhibitor competes with substrate for active site.
    Competitive,
    /// Inhibitor binds to enzyme-substrate complex only.
    Uncompetitive,
    /// Inhibitor binds to free enzyme or enzyme-substrate complex.
    Noncompetitive,
}

/// Michaelis-Menten with competitive inhibition.
///
/// `v = V_max * [S] / (K_m * (1 + [I]/K_i) + [S])`
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the inhibited rate without side effects"]
pub fn competitive_inhibition(
    substrate: f64,
    v_max: f64,
    k_m: f64,
    inhibitor: f64,
    k_i: f64,
) -> Result<f64> {
    validate_non_negative(substrate, "substrate")?;
    validate_positive(v_max, "v_max")?;
    validate_positive(k_m, "k_m")?;
    validate_non_negative(inhibitor, "inhibitor")?;
    validate_positive(k_i, "k_i")?;
    let apparent_km = k_m * (1.0 + inhibitor / k_i);
    Ok(v_max * substrate / (apparent_km + substrate))
}

/// Michaelis-Menten with uncompetitive inhibition.
///
/// `v = V_max * [S] / (K_m + [S] * (1 + [I]/K_i))`
///
/// The inhibitor binds only the enzyme-substrate complex, reducing both
/// the apparent V_max and K_m by the same factor.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the inhibited rate without side effects"]
pub fn uncompetitive_inhibition(
    substrate: f64,
    v_max: f64,
    k_m: f64,
    inhibitor: f64,
    k_i: f64,
) -> Result<f64> {
    validate_non_negative(substrate, "substrate")?;
    validate_positive(v_max, "v_max")?;
    validate_positive(k_m, "k_m")?;
    validate_non_negative(inhibitor, "inhibitor")?;
    validate_positive(k_i, "k_i")?;
    let alpha_prime = 1.0 + inhibitor / k_i;
    Ok(v_max * substrate / (k_m + substrate * alpha_prime))
}

/// Michaelis-Menten with noncompetitive inhibition.
///
/// `v = V_max * [S] / ((K_m + [S]) * (1 + [I]/K_i))`
///
/// The inhibitor binds both free enzyme and enzyme-substrate complex
/// equally, reducing V_max without affecting K_m.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the inhibited rate without side effects"]
pub fn noncompetitive_inhibition(
    substrate: f64,
    v_max: f64,
    k_m: f64,
    inhibitor: f64,
    k_i: f64,
) -> Result<f64> {
    validate_non_negative(substrate, "substrate")?;
    validate_positive(v_max, "v_max")?;
    validate_positive(k_m, "k_m")?;
    validate_non_negative(inhibitor, "inhibitor")?;
    validate_positive(k_i, "k_i")?;
    let alpha = 1.0 + inhibitor / k_i;
    Ok(v_max * substrate / ((k_m + substrate) * alpha))
}

/// Fermentation type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum FermentationType {
    /// Ethanol + CO2 (yeast, some bacteria).
    Alcoholic,
    /// Lactic acid (muscle cells, Lactobacillus).
    Lactic,
}

impl FermentationType {
    /// Net ATP yield from fermentation (per glucose).
    #[inline]
    #[must_use]
    pub const fn atp_yield(self) -> u32 {
        // Both types yield 2 ATP (from glycolysis)
        2
    }
}

/// A metabolic reaction with stoichiometric coefficients.
///
/// Negative coefficients are substrates (consumed), positive are products.
/// Each entry is a `(metabolite_name, coefficient)` pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Reaction {
    /// Reaction identifier.
    pub id: String,
    /// Stoichiometric coefficients: metabolite name → coefficient.
    /// Negative = consumed, positive = produced.
    pub stoichiometry: Vec<(String, f64)>,
    /// Whether this reaction is reversible.
    pub reversible: bool,
}

/// A metabolic network: a set of reactions over a set of metabolites.
///
/// Internally constructs the stoichiometric matrix S where `S[i][j]` is the
/// coefficient of metabolite `i` in reaction `j`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetabolicNetwork {
    /// Ordered list of metabolite names (row indices of S).
    pub metabolites: Vec<String>,
    /// Ordered list of reactions (column indices of S).
    pub reactions: Vec<Reaction>,
    /// Stoichiometric matrix (metabolites × reactions), row-major.
    pub s_matrix: Vec<Vec<f64>>,
}

impl MetabolicNetwork {
    /// Build a metabolic network from a set of reactions.
    ///
    /// Automatically discovers all metabolites (in order of first appearance)
    /// and constructs the stoichiometric matrix. Repeated entries for the same
    /// metabolite within a reaction are summed.
    #[must_use]
    pub fn from_reactions(reactions: Vec<Reaction>) -> Self {
        let n_rxn = reactions.len();
        let mut metabolites: Vec<String> = Vec::new();
        let mut s_matrix: Vec<Vec<f64>> = Vec::new();
        let mut index: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
        for (j, rxn) in reactions.iter().enumerate() {
            for (met, coeff) in &rxn.stoichiometry {
                let i = *index.entry(met.as_str()).or_insert_with(|| {
                    metabolites.push(met.clone());
                    s_matrix.push(vec![0.0; n_rxn]);
                    metabolites.len() - 1
                });
                if let Some(cell) = s_matrix.get_mut(i).and_then(|row| row.get_mut(j)) {
                    *cell += *coeff;
                }
            }
        }

        Self {
            metabolites,
            reactions,
            s_matrix,
        }
    }

    /// Number of metabolites (rows of S).
    #[inline]
    #[must_use]
    pub fn n_metabolites(&self) -> usize {
        self.metabolites.len()
    }

    /// Number of reactions (columns of S).
    #[inline]
    #[must_use]
    pub fn n_reactions(&self) -> usize {
        self.reactions.len()
    }

    /// Check that `s_matrix` is `n_metabolites × n_reactions`.
    ///
    /// The fields are public, so a network can be edited into an
    /// inconsistent shape; every computation checks this first.
    fn validate_shape(&self) -> Result<()> {
        let n_rxn = self.n_reactions();
        if self.s_matrix.len() != self.n_metabolites()
            || self.s_matrix.iter().any(|row| row.len() != n_rxn)
        {
            return Err(JivanuError::ComputationError(format!(
                "stoichiometric matrix must be {} × {n_rxn}",
                self.n_metabolites()
            )));
        }
        Ok(())
    }

    fn check_fluxes(&self, fluxes: &[f64]) -> Result<()> {
        if fluxes.len() != self.n_reactions() {
            return Err(JivanuError::ComputationError(format!(
                "flux vector length {} != {} reactions",
                fluxes.len(),
                self.n_reactions()
            )));
        }
        self.validate_shape()
    }

    /// Compute the net production rate for each metabolite: `S · v`.
    ///
    /// A positive value means net production, negative means net consumption.
    ///
    /// # Errors
    ///
    /// Returns error if the flux vector length doesn't match the number of
    /// reactions, or the stoichiometric matrix has an inconsistent shape.
    #[must_use = "returns the net production rates without side effects"]
    pub fn net_production(&self, fluxes: &[f64]) -> Result<Vec<f64>> {
        self.check_fluxes(fluxes)?;
        Ok(self.s_matrix.iter().map(|row| dot(row, fluxes)).collect())
    }

    /// Check whether a flux vector satisfies steady-state mass balance.
    ///
    /// Returns `true` if `|S · v|_∞ < tolerance` (all metabolite
    /// accumulation rates are within tolerance of zero).
    ///
    /// # Errors
    ///
    /// Returns error if the flux vector length doesn't match or `tolerance`
    /// is negative or non-finite.
    #[must_use = "returns whether steady state holds without side effects"]
    pub fn is_steady_state(&self, fluxes: &[f64], tolerance: f64) -> Result<bool> {
        validate_non_negative(tolerance, "tolerance")?;
        self.check_fluxes(fluxes)?;
        Ok(self
            .s_matrix
            .iter()
            .all(|row| dot(row, fluxes).abs() < tolerance))
    }

    /// Compute net ATP yield for a flux distribution.
    ///
    /// Sums `flux[j] * atp_coeff[j]` for each reaction that produces or
    /// consumes ATP. The `atp_metabolite` parameter specifies the name of
    /// the ATP metabolite in the network (e.g., "ATP" or "atp").
    ///
    /// # Errors
    ///
    /// Returns error if the flux vector length doesn't match.
    #[must_use = "returns the net ATP yield without side effects"]
    pub fn net_atp(&self, fluxes: &[f64], atp_metabolite: &str) -> Result<f64> {
        self.check_fluxes(fluxes)?;
        let atp_row = self
            .metabolites
            .iter()
            .position(|m| m == atp_metabolite)
            .and_then(|i| self.s_matrix.get(i));
        // ATP not in network → 0 yield
        Ok(atp_row.map_or(0.0, |row| dot(row, fluxes)))
    }
}

/// Result of flux balance analysis.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FbaResult {
    /// Optimal flux vector (one entry per reaction).
    pub fluxes: Vec<f64>,
    /// Objective function value at optimum.
    pub objective_value: f64,
    /// Whether a feasible solution was found.
    pub feasible: bool,
}

/// Validate FBA inputs and build the bounded LP over the network.
fn fba_problem<'a>(
    network: &'a MetabolicNetwork,
    objective: &[f64],
    lower_bounds: &'a [f64],
    upper_bounds: &'a [f64],
    zeros: &'a [f64],
) -> Result<crate::lp::BoundedLp<'a>> {
    let n_rxn = network.n_reactions();
    if objective.len() != n_rxn || lower_bounds.len() != n_rxn || upper_bounds.len() != n_rxn {
        return Err(JivanuError::ComputationError(
            "objective, lower_bounds, upper_bounds must match reaction count".into(),
        ));
    }
    network.validate_shape()?;
    Ok(crate::lp::BoundedLp {
        rows: &network.s_matrix,
        rhs: zeros,
        lower: lower_bounds,
        upper: upper_bounds,
    })
}

/// Flux balance analysis: maximize a linear objective subject to
/// steady-state mass balance and flux bounds.
///
/// Solves the linear program `max c'v` subject to `Sv = 0`, `lb ≤ v ≤ ub`
/// exactly, using a dense two-phase simplex method (Dantzig pricing with
/// Bland's anti-cycling rule). Suitable for small to medium networks (up to a
/// few hundred reactions).
///
/// Every metabolite in the network is treated as internal (balanced). Model
/// uptake and secretion with explicit exchange reactions (e.g. `→ A` with
/// stoichiometry `[("A", 1.0)]`), as in standard constraint-based modeling.
///
/// Bounds may be infinite (`f64::NEG_INFINITY` / `f64::INFINITY`).
///
/// # Arguments
///
/// - `network` — the metabolic network (provides S matrix)
/// - `objective` — coefficient for each reaction in the objective (length = n_reactions)
/// - `lower_bounds` — minimum flux for each reaction
/// - `upper_bounds` — maximum flux for each reaction
///
/// # Errors
///
/// Returns error if dimensions are inconsistent, a bound is NaN or a lower
/// bound exceeds its upper bound, or the objective is unbounded. An
/// infeasible problem is reported as `Ok` with `feasible == false`.
#[must_use = "returns the FBA result without side effects"]
pub fn flux_balance_analysis(
    network: &MetabolicNetwork,
    objective: &[f64],
    lower_bounds: &[f64],
    upper_bounds: &[f64],
) -> Result<FbaResult> {
    let zeros = vec![0.0; network.n_metabolites()];
    let lp = fba_problem(network, objective, lower_bounds, upper_bounds, &zeros)?;
    match lp.maximize(objective, None)? {
        crate::lp::LpOutcome::Optimal(fluxes) => {
            let objective_value = dot(&fluxes, objective);
            Ok(FbaResult {
                fluxes,
                objective_value,
                feasible: true,
            })
        }
        crate::lp::LpOutcome::Infeasible => {
            let fluxes: Vec<f64> = lower_bounds
                .iter()
                .zip(upper_bounds)
                .map(|(&lo, &hi)| 0.0_f64.clamp(lo, hi))
                .collect();
            let objective_value = dot(&fluxes, objective);
            Ok(FbaResult {
                fluxes,
                objective_value,
                feasible: false,
            })
        }
        crate::lp::LpOutcome::Unbounded => Err(JivanuError::SimulationFailed(
            "FBA objective is unbounded; add finite flux bounds".into(),
        )),
    }
}

#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Flux variability analysis: find the min/max of each flux at (near-)optimal objective.
///
/// First solves FBA for the optimum `Z`, then for each reaction `j` solves two
/// linear programs, minimizing and maximizing `v_j` subject to `Sv = 0`, the
/// flux bounds, and `c'v ≥ Z − (1 − fraction)·|Z|`. With `fraction = 1.0`
/// the ranges are those compatible with the exact optimum; with `0.0` they
/// span the whole feasible flux space.
///
/// Returns `(min_flux, max_flux)` for each reaction. A direction that is
/// unbounded (possible only with infinite bounds) is reported as `±inf`.
///
/// # Errors
///
/// Returns error if FBA fails, the problem is infeasible, dimensions are
/// inconsistent, or `fraction` is outside `[0, 1]`.
#[must_use = "returns the flux ranges without side effects"]
pub fn flux_variability_analysis(
    network: &MetabolicNetwork,
    objective: &[f64],
    lower_bounds: &[f64],
    upper_bounds: &[f64],
    fraction: f64,
) -> Result<Vec<(f64, f64)>> {
    validate_finite(fraction, "fraction")?;
    if !(0.0..=1.0).contains(&fraction) {
        return Err(JivanuError::ComputationError(format!(
            "fraction must be in [0, 1], got {fraction}"
        )));
    }
    let fba = flux_balance_analysis(network, objective, lower_bounds, upper_bounds)?;
    if !fba.feasible {
        return Err(JivanuError::SimulationFailed(
            "FVA requires a feasible FBA problem".into(),
        ));
    }
    let z = fba.objective_value;
    // Small slack keeps the optimal face feasible despite rounding.
    let threshold = z - (1.0 - fraction) * z.abs() - 1e-9 * (1.0 + z.abs());

    let zeros = vec![0.0; network.n_metabolites()];
    let lp = fba_problem(network, objective, lower_bounds, upper_bounds, &zeros)?;
    let n_rxn = network.n_reactions();
    let mut unit = vec![0.0; n_rxn];
    let mut ranges = Vec::with_capacity(n_rxn);
    for j in 0..n_rxn {
        unit[j] = 1.0;
        let max_v = match lp.maximize(&unit, Some((objective, threshold)))? {
            crate::lp::LpOutcome::Optimal(v) => v[j],
            crate::lp::LpOutcome::Unbounded => f64::INFINITY,
            crate::lp::LpOutcome::Infeasible => fba.fluxes[j],
        };
        unit[j] = -1.0;
        let min_v = match lp.maximize(&unit, Some((objective, threshold)))? {
            crate::lp::LpOutcome::Optimal(v) => v[j],
            crate::lp::LpOutcome::Unbounded => f64::NEG_INFINITY,
            crate::lp::LpOutcome::Infeasible => fba.fluxes[j],
        };
        unit[j] = 0.0;
        ranges.push((min_v.min(max_v), max_v.max(min_v)));
    }
    Ok(ranges)
}

/// Hill equation — general cooperative binding / sigmoidal response.
///
/// `response = V_max × [S]^n / (K^n + [S]^n)`
///
/// At `n = 1`, reduces to Michaelis-Menten. `n > 1` gives positive
/// cooperativity (steeper sigmoid). `n < 1` gives negative cooperativity.
///
/// Used throughout microbiology: oxygen binding, gene regulation,
/// drug dose-response (Emax model).
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the response without side effects"]
pub fn hill_equation(substrate: f64, v_max: f64, k: f64, n: f64) -> Result<f64> {
    validate_non_negative(substrate, "substrate")?;
    validate_positive(v_max, "v_max")?;
    validate_positive(k, "k")?;
    validate_positive(n, "n")?;
    Ok(v_max * hill_fraction(substrate, k, n))
}

/// Fractional saturation `S^n / (K^n + S^n)`, evaluated as
/// `1 / (1 + (K/S)^n)` so that large `S` or `n` cannot overflow to `inf/inf`.
#[inline]
fn hill_fraction(s: f64, k: f64, n: f64) -> f64 {
    if s <= 0.0 {
        0.0
    } else {
        1.0 / (1.0 + (k / s).powf(n))
    }
}

/// Emax pharmacodynamic model — sigmoidal drug effect.
///
/// `E = E_max × C^n / (EC50^n + C^n)`
///
/// Standard model for concentration-dependent drug effect.
/// - `E_max` — maximum effect (e.g., maximum kill rate)
/// - `EC50` — concentration producing 50% of E_max
/// - `n` — Hill coefficient (steepness; n=1 hyperbolic, n>1 sigmoidal)
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the drug effect without side effects"]
pub fn emax_model(concentration: f64, e_max: f64, ec50: f64, n: f64) -> Result<f64> {
    validate_non_negative(concentration, "concentration")?;
    validate_positive(e_max, "e_max")?;
    validate_positive(ec50, "ec50")?;
    validate_positive(n, "n")?;
    Ok(e_max * hill_fraction(concentration, ec50, n))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_michaelis_menten_at_km() {
        // At [S] = K_m, v = V_max / 2
        let v = michaelis_menten(1.0, 10.0, 1.0).unwrap();
        assert!((v - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_michaelis_menten_high_substrate() {
        let v = michaelis_menten(1000.0, 10.0, 1.0).unwrap();
        assert!((v - 10.0).abs() < 0.1);
    }

    #[test]
    fn test_michaelis_menten_zero_substrate() {
        let v = michaelis_menten(0.0, 10.0, 1.0).unwrap();
        assert!((v - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_lineweaver_burk() {
        let (inv_s, inv_v) = lineweaver_burk(5.0, 2.0).unwrap();
        assert!((inv_s - 0.5).abs() < 1e-10);
        assert!((inv_v - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_glycolysis_atp() {
        assert_eq!(glycolysis_atp(), 2);
    }

    #[test]
    fn test_oxidative_phosphorylation_atp() {
        assert_eq!(oxidative_phosphorylation_atp(), 34);
    }

    #[test]
    fn test_total_aerobic_atp() {
        assert_eq!(total_aerobic_atp(), 38);
    }

    #[test]
    fn test_competitive_inhibition() {
        // With inhibitor, rate should be less than without
        let v_no_inhib = michaelis_menten(1.0, 10.0, 1.0).unwrap();
        let v_inhib = competitive_inhibition(1.0, 10.0, 1.0, 1.0, 1.0).unwrap();
        assert!(v_inhib < v_no_inhib);
    }

    #[test]
    fn test_competitive_inhibition_zero_inhibitor() {
        // No inhibitor: same as Michaelis-Menten
        let v = competitive_inhibition(1.0, 10.0, 1.0, 0.0, 1.0).unwrap();
        let v_mm = michaelis_menten(1.0, 10.0, 1.0).unwrap();
        assert!((v - v_mm).abs() < 1e-10);
    }

    #[test]
    fn test_uncompetitive_inhibition_reduces_rate() {
        let v_no = michaelis_menten(1.0, 10.0, 1.0).unwrap();
        let v_inh = uncompetitive_inhibition(1.0, 10.0, 1.0, 1.0, 1.0).unwrap();
        assert!(v_inh < v_no);
    }

    #[test]
    fn test_uncompetitive_inhibition_zero_inhibitor() {
        let v = uncompetitive_inhibition(1.0, 10.0, 1.0, 0.0, 1.0).unwrap();
        let v_mm = michaelis_menten(1.0, 10.0, 1.0).unwrap();
        assert!((v - v_mm).abs() < 1e-10);
    }

    #[test]
    fn test_noncompetitive_inhibition_reduces_rate() {
        let v_no = michaelis_menten(1.0, 10.0, 1.0).unwrap();
        let v_inh = noncompetitive_inhibition(1.0, 10.0, 1.0, 1.0, 1.0).unwrap();
        assert!(v_inh < v_no);
    }

    #[test]
    fn test_noncompetitive_inhibition_zero_inhibitor() {
        let v = noncompetitive_inhibition(1.0, 10.0, 1.0, 0.0, 1.0).unwrap();
        let v_mm = michaelis_menten(1.0, 10.0, 1.0).unwrap();
        assert!((v - v_mm).abs() < 1e-10);
    }

    #[test]
    fn test_noncompetitive_does_not_change_km() {
        // At [S] = K_m, noncompetitive should give V_max/(2*alpha)
        let v = noncompetitive_inhibition(1.0, 10.0, 1.0, 1.0, 1.0).unwrap();
        // alpha = 1 + 1/1 = 2, so v = 10*1 / ((1+1)*2) = 10/4 = 2.5
        assert!((v - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_inhibition_types_differ() {
        // All three inhibition types should give different rates
        let v_comp = competitive_inhibition(2.0, 10.0, 1.0, 1.0, 1.0).unwrap();
        let v_uncomp = uncompetitive_inhibition(2.0, 10.0, 1.0, 1.0, 1.0).unwrap();
        let v_noncomp = noncompetitive_inhibition(2.0, 10.0, 1.0, 1.0, 1.0).unwrap();
        // All should be less than uninhibited
        let v_mm = michaelis_menten(2.0, 10.0, 1.0).unwrap();
        assert!(v_comp < v_mm);
        assert!(v_uncomp < v_mm);
        assert!(v_noncomp < v_mm);
        // They should differ from each other
        assert!((v_comp - v_uncomp).abs() > 0.01);
        assert!((v_comp - v_noncomp).abs() > 0.01);
    }

    #[test]
    fn test_fermentation_atp() {
        assert_eq!(FermentationType::Alcoholic.atp_yield(), 2);
        assert_eq!(FermentationType::Lactic.atp_yield(), 2);
    }

    #[test]
    fn test_inhibition_type_serde_roundtrip() {
        let it = InhibitionType::Competitive;
        let json = serde_json::to_string(&it).unwrap();
        let back: InhibitionType = serde_json::from_str(&json).unwrap();
        assert_eq!(it, back);
    }

    #[test]
    fn test_fermentation_serde_roundtrip() {
        let ft = FermentationType::Alcoholic;
        let json = serde_json::to_string(&ft).unwrap();
        let back: FermentationType = serde_json::from_str(&json).unwrap();
        assert_eq!(ft, back);
    }

    // --- Metabolic network tests ---

    /// Build a toy glycolysis-like network:
    /// R1: Glucose → 2 Pyruvate + 2 ATP
    /// R2: Pyruvate → Acetyl-CoA + CO2
    fn toy_network() -> MetabolicNetwork {
        MetabolicNetwork::from_reactions(vec![
            Reaction {
                id: "glycolysis".into(),
                stoichiometry: vec![
                    ("Glucose".into(), -1.0),
                    ("Pyruvate".into(), 2.0),
                    ("ATP".into(), 2.0),
                ],
                reversible: false,
            },
            Reaction {
                id: "pyruvate_decarboxylation".into(),
                stoichiometry: vec![
                    ("Pyruvate".into(), -1.0),
                    ("AcetylCoA".into(), 1.0),
                    ("CO2".into(), 1.0),
                ],
                reversible: false,
            },
        ])
    }

    #[test]
    fn test_network_dimensions() {
        let net = toy_network();
        assert_eq!(net.n_metabolites(), 5); // Glucose, Pyruvate, ATP, AcetylCoA, CO2
        assert_eq!(net.n_reactions(), 2);
        assert_eq!(net.s_matrix.len(), 5);
        assert_eq!(net.s_matrix[0].len(), 2);
    }

    #[test]
    fn test_stoichiometric_matrix_values() {
        let net = toy_network();
        let glucose_idx = net.metabolites.iter().position(|m| m == "Glucose").unwrap();
        let pyruvate_idx = net
            .metabolites
            .iter()
            .position(|m| m == "Pyruvate")
            .unwrap();
        let atp_idx = net.metabolites.iter().position(|m| m == "ATP").unwrap();

        // R1 (col 0): Glucose = -1, Pyruvate = +2, ATP = +2
        assert!((net.s_matrix[glucose_idx][0] - (-1.0)).abs() < 1e-10);
        assert!((net.s_matrix[pyruvate_idx][0] - 2.0).abs() < 1e-10);
        assert!((net.s_matrix[atp_idx][0] - 2.0).abs() < 1e-10);

        // R2 (col 1): Pyruvate = -1
        assert!((net.s_matrix[pyruvate_idx][1] - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_net_production() {
        let net = toy_network();
        // v = [1.0, 2.0]: 1 glycolysis + 2 pyruvate decarboxylations
        let sv = net.net_production(&[1.0, 2.0]).unwrap();
        let glucose_idx = net.metabolites.iter().position(|m| m == "Glucose").unwrap();
        let pyruvate_idx = net
            .metabolites
            .iter()
            .position(|m| m == "Pyruvate")
            .unwrap();
        let atp_idx = net.metabolites.iter().position(|m| m == "ATP").unwrap();

        // Glucose: -1*1 = -1
        assert!((sv[glucose_idx] - (-1.0)).abs() < 1e-10);
        // Pyruvate: 2*1 + (-1)*2 = 0 (balanced!)
        assert!((sv[pyruvate_idx] - 0.0).abs() < 1e-10);
        // ATP: 2*1 + 0*2 = 2
        assert!((sv[atp_idx] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_steady_state_balanced() {
        let net = toy_network();
        // v = [1.0, 2.0] balances pyruvate but not glucose/ATP
        assert!(!net.is_steady_state(&[1.0, 2.0], 1e-10).unwrap());
    }

    #[test]
    fn test_steady_state_zero_flux() {
        let net = toy_network();
        // Zero flux → trivially at steady state
        assert!(net.is_steady_state(&[0.0, 0.0], 1e-10).unwrap());
    }

    #[test]
    fn test_net_atp() {
        let net = toy_network();
        // v = [1.0, 2.0]: ATP yield from glycolysis = 2*1 = 2
        let atp = net.net_atp(&[1.0, 2.0], "ATP").unwrap();
        assert!((atp - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_net_atp_missing_metabolite() {
        let net = toy_network();
        // Nonexistent metabolite → 0
        let atp = net.net_atp(&[1.0, 2.0], "NADH").unwrap();
        assert!((atp - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_net_production_wrong_length() {
        let net = toy_network();
        assert!(net.net_production(&[1.0]).is_err());
    }

    #[test]
    fn test_network_serde_roundtrip() {
        let net = toy_network();
        let json = serde_json::to_string(&net).unwrap();
        let back: MetabolicNetwork = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n_metabolites(), net.n_metabolites());
        assert_eq!(back.n_reactions(), net.n_reactions());
        assert!((back.s_matrix[0][0] - net.s_matrix[0][0]).abs() < 1e-10);
    }

    // --- Hill / Emax tests ---

    #[test]
    fn test_hill_n1_equals_michaelis_menten() {
        let hill = hill_equation(1.0, 10.0, 1.0, 1.0).unwrap();
        let mm = michaelis_menten(1.0, 10.0, 1.0).unwrap();
        assert!((hill - mm).abs() < 1e-10);
    }

    #[test]
    fn test_hill_at_k() {
        // At [S] = K, response = V_max / 2 (for any n)
        let v = hill_equation(5.0, 10.0, 5.0, 3.0).unwrap();
        assert!((v - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_hill_cooperativity() {
        // Higher n → steeper curve → lower response below K, higher above K
        let v_n1 = hill_equation(0.5, 10.0, 1.0, 1.0).unwrap();
        let v_n4 = hill_equation(0.5, 10.0, 1.0, 4.0).unwrap();
        assert!(v_n4 < v_n1, "higher n should give lower response below K");
    }

    #[test]
    fn test_emax_at_ec50() {
        // At C = EC50, E = E_max / 2
        let e = emax_model(5.0, 100.0, 5.0, 2.0).unwrap();
        assert!((e - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_emax_high_concentration() {
        // At very high C, E → E_max
        let e = emax_model(1e6, 100.0, 5.0, 2.0).unwrap();
        assert!((e - 100.0).abs() < 0.01);
    }

    #[test]
    fn test_emax_zero_concentration() {
        let e = emax_model(0.0, 100.0, 5.0, 2.0).unwrap();
        assert!((e - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_duplicate_metabolite_in_reaction() {
        // If a metabolite appears twice, coefficients should accumulate
        let net = MetabolicNetwork::from_reactions(vec![Reaction {
            id: "test".into(),
            stoichiometry: vec![("ATP".into(), 2.0), ("ATP".into(), 1.0)],
            reversible: false,
        }]);
        let atp_idx = net.metabolites.iter().position(|m| m == "ATP").unwrap();
        assert!((net.s_matrix[atp_idx][0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_reaction_serde_roundtrip() {
        let rxn = Reaction {
            id: "test".into(),
            stoichiometry: vec![("A".into(), -1.0), ("B".into(), 1.0)],
            reversible: true,
        };
        let json = serde_json::to_string(&rxn).unwrap();
        let back: Reaction = serde_json::from_str(&json).unwrap();
        assert_eq!(rxn.id, back.id);
        assert_eq!(rxn.reversible, back.reversible);
    }

    // --- FBA tests ---

    fn rxn(id: &str, stoich: &[(&str, f64)]) -> Reaction {
        Reaction {
            id: id.into(),
            stoichiometry: stoich.iter().map(|&(m, c)| (m.into(), c)).collect(),
            reversible: false,
        }
    }

    /// Uptake → A → B → C → secretion.
    fn linear_pathway() -> MetabolicNetwork {
        MetabolicNetwork::from_reactions(vec![
            rxn("EX_A", &[("A", 1.0)]),
            rxn("r1", &[("A", -1.0), ("B", 1.0)]),
            rxn("r2", &[("B", -1.0), ("C", 1.0)]),
            rxn("EX_C", &[("C", -1.0)]),
        ])
    }

    #[test]
    fn test_fba_simple_linear() {
        // Uptake capped at 10 → every flux in the chain is 10.
        let net = linear_pathway();
        let result = flux_balance_analysis(
            &net,
            &[0.0, 0.0, 1.0, 0.0],
            &[0.0; 4],
            &[10.0, 1000.0, 1000.0, 1000.0],
        )
        .unwrap();
        assert!(result.feasible);
        assert!((result.objective_value - 10.0).abs() < 1e-9);
        for v in &result.fluxes {
            assert!((v - 10.0).abs() < 1e-9, "fluxes = {:?}", result.fluxes);
        }
        assert!(net.is_steady_state(&result.fluxes, 1e-9).unwrap());
    }

    #[test]
    fn test_fba_closed_network_has_zero_flux() {
        // A → B → C with no exchange reactions: mass balance forces v = 0.
        let net = MetabolicNetwork::from_reactions(vec![
            rxn("r1", &[("A", -1.0), ("B", 1.0)]),
            rxn("r2", &[("B", -1.0), ("C", 1.0)]),
        ]);
        let result = flux_balance_analysis(&net, &[0.0, 1.0], &[0.0, 0.0], &[10.0, 10.0]).unwrap();
        assert!(result.feasible);
        assert!(result.objective_value.abs() < 1e-12);
    }

    #[test]
    fn test_fba_branch_choice() {
        // A can become B (yield 1) or C (yield 2 ATP); maximize ATP.
        let net = MetabolicNetwork::from_reactions(vec![
            rxn("EX_A", &[("A", 1.0)]),
            rxn("to_B", &[("A", -1.0), ("B", 1.0), ("ATP", 1.0)]),
            rxn("to_C", &[("A", -1.0), ("C", 1.0), ("ATP", 2.0)]),
            rxn("EX_B", &[("B", -1.0)]),
            rxn("EX_C", &[("C", -1.0)]),
            rxn("ATPM", &[("ATP", -1.0)]),
        ]);
        let lb = [0.0; 6];
        let ub = [5.0, 100.0, 100.0, 100.0, 100.0, 100.0];
        let result =
            flux_balance_analysis(&net, &[0.0, 0.0, 0.0, 0.0, 0.0, 1.0], &lb, &ub).unwrap();
        assert!((result.objective_value - 10.0).abs() < 1e-9);
        assert!((result.fluxes[2] - 5.0).abs() < 1e-9, "all flux via to_C");
        assert!(result.fluxes[1].abs() < 1e-9);
    }

    #[test]
    fn test_fba_infeasible() {
        // Secretion forced ≥ 5 but uptake capped at 1.
        let net = linear_pathway();
        let result = flux_balance_analysis(
            &net,
            &[0.0, 0.0, 0.0, 1.0],
            &[0.0, 0.0, 0.0, 5.0],
            &[1.0, 10.0, 10.0, 10.0],
        )
        .unwrap();
        assert!(!result.feasible);
    }

    #[test]
    fn test_fba_unbounded_is_error() {
        let net = linear_pathway();
        let inf = f64::INFINITY;
        assert!(flux_balance_analysis(&net, &[0.0, 0.0, 1.0, 0.0], &[0.0; 4], &[inf; 4]).is_err());
    }

    #[test]
    fn test_fba_rejects_inverted_or_nan_bounds() {
        // Previously panicked inside f64::clamp.
        let net = linear_pathway();
        assert!(flux_balance_analysis(&net, &[1.0; 4], &[0.0, 5.0, 0.0, 0.0], &[1.0; 4]).is_err());
        assert!(
            flux_balance_analysis(&net, &[1.0; 4], &[f64::NAN, 0.0, 0.0, 0.0], &[1.0; 4]).is_err()
        );
    }

    #[test]
    fn test_fba_reversible_reaction() {
        // EX_A: A ⇌ (secretion when positive, uptake when negative).
        let net = MetabolicNetwork::from_reactions(vec![
            rxn("EX_A", &[("A", -1.0)]),
            rxn("r1", &[("A", -1.0), ("B", 1.0)]),
            rxn("EX_B", &[("B", -1.0)]),
        ]);
        let result = flux_balance_analysis(
            &net,
            &[0.0, 1.0, 0.0],
            &[-8.0, 0.0, 0.0],
            &[1000.0, 1000.0, 1000.0],
        )
        .unwrap();
        assert!((result.objective_value - 8.0).abs() < 1e-9);
        assert!((result.fluxes[0] + 8.0).abs() < 1e-9);
    }

    #[test]
    fn test_inconsistent_s_matrix_is_error_not_panic() {
        let mut net = linear_pathway();
        net.s_matrix[0].pop();
        assert!(net.net_production(&[1.0; 4]).is_err());
        assert!(flux_balance_analysis(&net, &[1.0; 4], &[0.0; 4], &[1.0; 4]).is_err());
    }

    #[test]
    fn test_fba_dimension_mismatch() {
        let net = MetabolicNetwork::from_reactions(vec![Reaction {
            id: "r1".into(),
            stoichiometry: vec![("A".into(), -1.0)],
            reversible: false,
        }]);
        assert!(flux_balance_analysis(&net, &[1.0, 2.0], &[0.0], &[10.0]).is_err());
    }

    #[test]
    fn test_fba_result_serde_roundtrip() {
        let result = FbaResult {
            fluxes: vec![1.0, 2.0],
            objective_value: 5.0,
            feasible: true,
        };
        let json = serde_json::to_string(&result).unwrap();
        let back: FbaResult = serde_json::from_str(&json).unwrap();
        assert!((result.objective_value - back.objective_value).abs() < 1e-10);
        assert_eq!(result.feasible, back.feasible);
    }

    #[test]
    fn test_fva_returns_ranges() {
        let net = MetabolicNetwork::from_reactions(vec![
            Reaction {
                id: "r1".into(),
                stoichiometry: vec![("A".into(), -1.0), ("B".into(), 1.0)],
                reversible: false,
            },
            Reaction {
                id: "r2".into(),
                stoichiometry: vec![("B".into(), -1.0), ("C".into(), 1.0)],
                reversible: false,
            },
        ]);
        let ranges =
            flux_variability_analysis(&net, &[0.0, 1.0], &[0.0, 0.0], &[10.0, 10.0], 1.0).unwrap();
        assert_eq!(ranges.len(), 2);
        for (lo, hi) in &ranges {
            assert!(lo <= hi);
        }
    }
    /// Uptake A, two parallel isozymes A → B, secrete B.
    fn parallel_pathway() -> MetabolicNetwork {
        MetabolicNetwork::from_reactions(vec![
            rxn("EX_A", &[("A", 1.0)]),
            rxn("iso1", &[("A", -1.0), ("B", 1.0)]),
            rxn("iso2", &[("A", -1.0), ("B", 1.0)]),
            rxn("EX_B", &[("B", -1.0)]),
        ])
    }

    #[test]
    fn test_fva_parallel_isozymes() {
        let net = parallel_pathway();
        let lb = [0.0; 4];
        let ub = [10.0, 100.0, 100.0, 100.0];
        let c = [0.0, 0.0, 0.0, 1.0];
        let ranges = flux_variability_analysis(&net, &c, &lb, &ub, 1.0).unwrap();
        let expected = [(10.0, 10.0), (0.0, 10.0), (0.0, 10.0), (10.0, 10.0)];
        for (j, ((lo, hi), (elo, ehi))) in ranges.iter().zip(expected).enumerate() {
            assert!((lo - elo).abs() < 1e-7, "reaction {j}: lo {lo} != {elo}");
            assert!((hi - ehi).abs() < 1e-7, "reaction {j}: hi {hi} != {ehi}");
        }
    }

    #[test]
    fn test_fva_fraction_relaxes_ranges() {
        let net = parallel_pathway();
        let lb = [0.0; 4];
        let ub = [10.0, 100.0, 100.0, 100.0];
        let c = [0.0, 0.0, 0.0, 1.0];
        // 50% of optimum: secretion may drop to 5.
        let half = flux_variability_analysis(&net, &c, &lb, &ub, 0.5).unwrap();
        assert!((half[3].0 - 5.0).abs() < 1e-7, "{:?}", half[3]);
        assert!((half[3].1 - 10.0).abs() < 1e-7);
        // 0%: anything feasible.
        let free = flux_variability_analysis(&net, &c, &lb, &ub, 0.0).unwrap();
        assert!(free[0].0.abs() < 1e-7);
    }

    #[test]
    fn test_fva_invalid_fraction() {
        let net = parallel_pathway();
        let lb = [0.0; 4];
        let ub = [10.0; 4];
        let c = [0.0, 0.0, 0.0, 1.0];
        assert!(flux_variability_analysis(&net, &c, &lb, &ub, 1.5).is_err());
        assert!(flux_variability_analysis(&net, &c, &lb, &ub, f64::NAN).is_err());
    }

    #[test]
    fn test_lineweaver_burk_line() {
        // K_m = 2, V_max = 10 → slope 0.2, intercept 0.1
        let (slope, intercept) = lineweaver_burk_line(10.0, 2.0).unwrap();
        assert!((slope - 0.2).abs() < 1e-12);
        assert!((intercept - 0.1).abs() < 1e-12);
        // A point generated by Michaelis-Menten lies on the line.
        let v = michaelis_menten(5.0, 10.0, 2.0).unwrap();
        let (x, y) = lineweaver_burk(v, 5.0).unwrap();
        assert!((y - (slope * x + intercept)).abs() < 1e-12);
    }

    #[test]
    fn test_hill_extreme_inputs_finite() {
        // S^n overflowed to inf/inf = NaN with the naive form.
        let r = hill_equation(1e200, 1.0, 1.0, 4.0).unwrap();
        assert!((r - 1.0).abs() < 1e-12);
        let e = emax_model(1e300, 2.0, 1.0, 3.0).unwrap();
        assert!((e - 2.0).abs() < 1e-12);
        assert_eq!(hill_equation(0.0, 1.0, 1.0, 0.5).unwrap(), 0.0);
    }

    #[test]
    fn test_metabolic_network_serde_roundtrip() {
        let net = linear_pathway();
        let json = serde_json::to_string(&net).unwrap();
        let back: MetabolicNetwork = serde_json::from_str(&json).unwrap();
        assert_eq!(net.metabolites, back.metabolites);
        assert_eq!(net.s_matrix, back.s_matrix);
    }
}
