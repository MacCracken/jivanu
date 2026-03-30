//! Metabolism — enzyme kinetics, metabolic pathways, ATP yields.

use serde::{Deserialize, Serialize};

use crate::error::{validate_finite, validate_non_negative, validate_positive, Result};

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
}
