//! Metabolism — enzyme kinetics, metabolic pathways, ATP yields.

use serde::{Deserialize, Serialize};

use crate::error::{Result, validate_non_negative, validate_positive};

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

/// Lineweaver-Burk double reciprocal plot parameters.
///
/// `1/v = (K_m / V_max) * (1/[S]) + 1/V_max`
///
/// Returns `(slope = K_m/V_max, intercept = 1/V_max)`.
///
/// # Errors
///
/// Returns error if v or s is non-positive.
#[inline]
#[must_use = "returns the plot parameters without side effects"]
pub fn lineweaver_burk(v: f64, s: f64) -> Result<(f64, f64)> {
    validate_positive(v, "v")?;
    validate_positive(s, "s")?;
    Ok((1.0 / s, 1.0 / v))
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
/// The coefficient map keys are metabolite names.
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
    /// Automatically discovers all metabolites and constructs the
    /// stoichiometric matrix.
    #[must_use]
    pub fn from_reactions(reactions: Vec<Reaction>) -> Self {
        // Collect unique metabolite names in order of first appearance
        let mut metabolites: Vec<String> = Vec::new();
        let mut met_index = std::collections::HashMap::new();
        for rxn in &reactions {
            for (met, _) in &rxn.stoichiometry {
                if !met_index.contains_key(met) {
                    met_index.insert(met.clone(), metabolites.len());
                    metabolites.push(met.clone());
                }
            }
        }

        let n_met = metabolites.len();
        let n_rxn = reactions.len();
        let mut s_matrix = vec![vec![0.0; n_rxn]; n_met];

        for (j, rxn) in reactions.iter().enumerate() {
            for (met, coeff) in &rxn.stoichiometry {
                if let Some(&i) = met_index.get(met) {
                    s_matrix[i][j] = *coeff;
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

    /// Compute the net production rate for each metabolite: `S · v`.
    ///
    /// A positive value means net production, negative means net consumption.
    ///
    /// # Errors
    ///
    /// Returns error if the flux vector length doesn't match the number of reactions.
    #[must_use = "returns the net production rates without side effects"]
    pub fn net_production(&self, fluxes: &[f64]) -> Result<Vec<f64>> {
        if fluxes.len() != self.n_reactions() {
            return Err(crate::error::JivanuError::ComputationError(format!(
                "flux vector length {} != {} reactions",
                fluxes.len(),
                self.n_reactions()
            )));
        }
        let mut result = vec![0.0; self.n_metabolites()];
        for (i, row) in self.s_matrix.iter().enumerate() {
            for (j, &coeff) in row.iter().enumerate() {
                result[i] += coeff * fluxes[j];
            }
        }
        Ok(result)
    }

    /// Check whether a flux vector satisfies steady-state mass balance.
    ///
    /// Returns `true` if `|S · v|_∞ < tolerance` (all metabolite
    /// accumulation rates are within tolerance of zero).
    ///
    /// # Errors
    ///
    /// Returns error if the flux vector length doesn't match.
    #[must_use = "returns whether steady state holds without side effects"]
    pub fn is_steady_state(&self, fluxes: &[f64], tolerance: f64) -> Result<bool> {
        let sv = self.net_production(fluxes)?;
        Ok(sv.iter().all(|&v| v.abs() < tolerance))
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
        if fluxes.len() != self.n_reactions() {
            return Err(crate::error::JivanuError::ComputationError(format!(
                "flux vector length {} != {} reactions",
                fluxes.len(),
                self.n_reactions()
            )));
        }
        let atp_row = self.metabolites.iter().position(|m| m == atp_metabolite);
        match atp_row {
            Some(i) => {
                let mut total = 0.0;
                for (j, &coeff) in self.s_matrix[i].iter().enumerate() {
                    total += coeff * fluxes[j];
                }
                Ok(total)
            }
            None => Ok(0.0), // ATP not in network → 0 yield
        }
    }
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
}
