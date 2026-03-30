//! Antibiotic resistance — MIC, kill curves, resistance transfer.

use serde::{Deserialize, Serialize};

use crate::error::{Result, validate_non_negative, validate_positive};

/// Antibiotic classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum AntibioticClass {
    /// Penicillins, cephalosporins, carbapenems.
    BetaLactam,
    /// Gentamicin, tobramycin.
    Aminoglycoside,
    /// Ciprofloxacin, levofloxacin.
    Fluoroquinolone,
    /// Erythromycin, azithromycin.
    Macrolide,
    /// Tetracycline, doxycycline.
    Tetracycline,
    /// Vancomycin.
    Glycopeptide,
}

/// Kill curve: survival fraction as a function of antibiotic concentration.
///
/// `survival = e^(-k * (concentration / MIC - 1))`
///
/// When `concentration = MIC`, survival = 1.0 (boundary definition).
/// Above MIC, survival decreases exponentially.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the survival fraction without side effects"]
pub fn kill_curve(concentration: f64, mic: f64, kill_rate: f64) -> Result<f64> {
    validate_non_negative(concentration, "concentration")?;
    validate_positive(mic, "mic")?;
    validate_positive(kill_rate, "kill_rate")?;

    if concentration <= mic {
        return Ok(1.0); // below MIC: no killing
    }

    let excess = concentration / mic - 1.0;
    Ok((-kill_rate * excess).exp())
}

/// Resistance transfer rate between donor and recipient populations.
///
/// `rate = donor_freq * contact_rate * transfer_efficiency`
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the transfer rate without side effects"]
pub fn resistance_transfer_rate(
    donor_freq: f64,
    contact_rate: f64,
    transfer_efficiency: f64,
) -> Result<f64> {
    validate_non_negative(donor_freq, "donor_freq")?;
    validate_non_negative(contact_rate, "contact_rate")?;
    validate_non_negative(transfer_efficiency, "transfer_efficiency")?;
    Ok(donor_freq * contact_rate * transfer_efficiency)
}

/// Interpretation of a Fractional Inhibitory Concentration (FIC) index.
///
/// Classification follows EUCAST/CLSI consensus thresholds:
/// - Synergy: FIC ≤ 0.5
/// - Additive: 0.5 < FIC ≤ 1.0
/// - Indifferent: 1.0 < FIC ≤ 4.0
/// - Antagonism: FIC > 4.0
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum DrugInteraction {
    /// Combined effect greater than sum of individual effects (FIC ≤ 0.5).
    Synergy,
    /// Combined effect equals sum of individual effects (0.5 < FIC ≤ 1.0).
    Additive,
    /// No meaningful interaction (1.0 < FIC ≤ 4.0).
    Indifferent,
    /// Combined effect less than individual effects (FIC > 4.0).
    Antagonism,
}

/// Fractional Inhibitory Concentration (FIC) index for a two-drug combination.
///
/// `FIC = (MIC_A_combo / MIC_A_alone) + (MIC_B_combo / MIC_B_alone)`
///
/// The FIC index quantifies whether two antimicrobial agents interact
/// synergistically, additively, or antagonistically.
///
/// # Arguments
///
/// - `mic_a_combo` — MIC of drug A in the presence of drug B
/// - `mic_a_alone` — MIC of drug A alone
/// - `mic_b_combo` — MIC of drug B in the presence of drug A
/// - `mic_b_alone` — MIC of drug B alone
///
/// # Errors
///
/// Returns error if any MIC value is non-positive.
#[inline]
#[must_use = "returns the FIC index without side effects"]
pub fn fic_index(
    mic_a_combo: f64,
    mic_a_alone: f64,
    mic_b_combo: f64,
    mic_b_alone: f64,
) -> Result<f64> {
    validate_positive(mic_a_combo, "mic_a_combo")?;
    validate_positive(mic_a_alone, "mic_a_alone")?;
    validate_positive(mic_b_combo, "mic_b_combo")?;
    validate_positive(mic_b_alone, "mic_b_alone")?;
    Ok(mic_a_combo / mic_a_alone + mic_b_combo / mic_b_alone)
}

/// Classify a two-drug interaction from its FIC index.
///
/// Uses EUCAST/CLSI consensus thresholds.
///
/// # Errors
///
/// Returns error if the FIC index is non-positive.
#[inline]
#[must_use = "returns the interaction classification without side effects"]
pub fn classify_interaction(fic: f64) -> Result<DrugInteraction> {
    validate_positive(fic, "fic")?;
    if fic <= 0.5 {
        Ok(DrugInteraction::Synergy)
    } else if fic <= 1.0 {
        Ok(DrugInteraction::Additive)
    } else if fic <= 4.0 {
        Ok(DrugInteraction::Indifferent)
    } else {
        Ok(DrugInteraction::Antagonism)
    }
}

/// Compute the FIC index and its interpretation in one call.
///
/// # Errors
///
/// Returns error if any MIC value is non-positive.
#[must_use = "returns the FIC index and interaction without side effects"]
pub fn fic_interaction(
    mic_a_combo: f64,
    mic_a_alone: f64,
    mic_b_combo: f64,
    mic_b_alone: f64,
) -> Result<(f64, DrugInteraction)> {
    let fic = fic_index(mic_a_combo, mic_a_alone, mic_b_combo, mic_b_alone)?;
    let interaction = classify_interaction(fic)?;
    Ok((fic, interaction))
}

/// Kill curve for a two-drug combination.
///
/// Models the combined bactericidal effect using the Bliss independence model:
///
/// `survival_combo = survival_A × survival_B`
///
/// Each drug's individual survival is computed via [`kill_curve`]. This model
/// assumes independent mechanisms of action — appropriate for drugs from
/// different antibiotic classes.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the combination survival fraction without side effects"]
pub fn combination_kill_curve(
    conc_a: f64,
    mic_a: f64,
    kill_rate_a: f64,
    conc_b: f64,
    mic_b: f64,
    kill_rate_b: f64,
) -> Result<f64> {
    let surv_a = kill_curve(conc_a, mic_a, kill_rate_a)?;
    let surv_b = kill_curve(conc_b, mic_b, kill_rate_b)?;
    Ok(surv_a * surv_b)
}

/// Result of a checkerboard assay for a two-drug combination.
///
/// Contains the FIC index at each grid point (drug A concentrations × drug B
/// concentrations) and the minimum FIC observed (the most synergistic point).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CheckerboardResult {
    /// FIC index at each grid point, row-major: `fic_grid[i][j]` is the FIC
    /// at `conc_a[i]` and `conc_b[j]`.
    pub fic_grid: Vec<Vec<f64>>,
    /// Minimum FIC index observed across the grid.
    pub min_fic: f64,
    /// Classification of the minimum FIC.
    pub interaction: DrugInteraction,
}

/// Run a checkerboard assay over two concentration series.
///
/// For each pair `(conc_a[i], conc_b[j])`, the effective MIC of each drug
/// in combination is modeled using the Loewe additivity null hypothesis:
///
/// `FIC_ij = conc_a[i] / mic_a + conc_b[j] / mic_b`
///
/// This produces the standard checkerboard isobologram grid used in
/// antimicrobial susceptibility testing.
///
/// # Arguments
///
/// - `conc_a` — concentration series for drug A (e.g., twofold dilutions)
/// - `conc_b` — concentration series for drug B
/// - `mic_a` — MIC of drug A alone
/// - `mic_b` — MIC of drug B alone
///
/// # Errors
///
/// Returns error if MICs are non-positive or concentration arrays are empty.
#[must_use = "returns the checkerboard result without side effects"]
pub fn checkerboard(
    conc_a: &[f64],
    conc_b: &[f64],
    mic_a: f64,
    mic_b: f64,
) -> Result<CheckerboardResult> {
    validate_positive(mic_a, "mic_a")?;
    validate_positive(mic_b, "mic_b")?;
    if conc_a.is_empty() || conc_b.is_empty() {
        return Err(crate::error::JivanuError::ComputationError(
            "concentration arrays must not be empty".into(),
        ));
    }

    let mut min_fic = f64::MAX;
    let mut fic_grid = Vec::with_capacity(conc_a.len());

    for &ca in conc_a {
        validate_non_negative(ca, "conc_a element")?;
        let mut row = Vec::with_capacity(conc_b.len());
        for &cb in conc_b {
            validate_non_negative(cb, "conc_b element")?;
            let fic = ca / mic_a + cb / mic_b;
            if fic < min_fic {
                min_fic = fic;
            }
            row.push(fic);
        }
        fic_grid.push(row);
    }

    // min_fic could be 0 if both concentrations are 0; clamp for classification
    let interaction = if min_fic > 0.0 {
        classify_interaction(min_fic)?
    } else {
        DrugInteraction::Synergy
    };

    Ok(CheckerboardResult {
        fic_grid,
        min_fic,
        interaction,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kill_curve_below_mic() {
        let survival = kill_curve(0.5, 1.0, 2.0).unwrap();
        assert!((survival - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_kill_curve_at_mic() {
        let survival = kill_curve(1.0, 1.0, 2.0).unwrap();
        assert!((survival - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_kill_curve_above_mic() {
        let survival = kill_curve(2.0, 1.0, 2.0).unwrap();
        assert!(survival < 1.0);
        assert!(survival > 0.0);
    }

    #[test]
    fn test_kill_curve_high_concentration() {
        let survival = kill_curve(10.0, 1.0, 1.0).unwrap();
        assert!(survival < 0.001);
    }

    #[test]
    fn test_resistance_transfer_rate() {
        let rate = resistance_transfer_rate(0.1, 0.5, 0.01).unwrap();
        assert!((rate - 0.0005).abs() < 1e-10);
    }

    #[test]
    fn test_antibiotic_class_serde_roundtrip() {
        let cls = AntibioticClass::BetaLactam;
        let json = serde_json::to_string(&cls).unwrap();
        let back: AntibioticClass = serde_json::from_str(&json).unwrap();
        assert_eq!(cls, back);
    }

    #[test]
    fn test_fic_index_synergy() {
        // Both MICs drop to 1/4 of alone value → FIC = 0.25 + 0.25 = 0.5
        let fic = fic_index(0.25, 1.0, 0.25, 1.0).unwrap();
        assert!((fic - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_fic_index_additive() {
        // Each at half MIC → FIC = 0.5 + 0.5 = 1.0
        let fic = fic_index(0.5, 1.0, 0.5, 1.0).unwrap();
        assert!((fic - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fic_index_antagonism() {
        // Combo MICs much higher → FIC > 4
        let fic = fic_index(4.0, 1.0, 2.0, 1.0).unwrap();
        assert!((fic - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_fic_index_invalid() {
        assert!(fic_index(0.0, 1.0, 0.5, 1.0).is_err());
        assert!(fic_index(0.5, 0.0, 0.5, 1.0).is_err());
    }

    #[test]
    fn test_classify_interaction_boundaries() {
        assert_eq!(
            classify_interaction(0.25).unwrap(),
            DrugInteraction::Synergy
        );
        assert_eq!(classify_interaction(0.5).unwrap(), DrugInteraction::Synergy);
        assert_eq!(
            classify_interaction(0.75).unwrap(),
            DrugInteraction::Additive
        );
        assert_eq!(
            classify_interaction(1.0).unwrap(),
            DrugInteraction::Additive
        );
        assert_eq!(
            classify_interaction(2.0).unwrap(),
            DrugInteraction::Indifferent
        );
        assert_eq!(
            classify_interaction(4.0).unwrap(),
            DrugInteraction::Indifferent
        );
        assert_eq!(
            classify_interaction(4.1).unwrap(),
            DrugInteraction::Antagonism
        );
    }

    #[test]
    fn test_fic_interaction_combined() {
        let (fic, interaction) = fic_interaction(0.125, 1.0, 0.125, 1.0).unwrap();
        assert!((fic - 0.25).abs() < 1e-10);
        assert_eq!(interaction, DrugInteraction::Synergy);
    }

    #[test]
    fn test_drug_interaction_serde_roundtrip() {
        let di = DrugInteraction::Synergy;
        let json = serde_json::to_string(&di).unwrap();
        let back: DrugInteraction = serde_json::from_str(&json).unwrap();
        assert_eq!(di, back);
    }

    #[test]
    fn test_combination_kill_curve_both_below_mic() {
        // Both below MIC → survival = 1.0 × 1.0
        let surv = combination_kill_curve(0.5, 1.0, 2.0, 0.5, 1.0, 2.0).unwrap();
        assert!((surv - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_combination_kill_curve_one_above() {
        // A above MIC, B below → survival = surv_A × 1.0
        let surv_combo = combination_kill_curve(2.0, 1.0, 2.0, 0.5, 1.0, 2.0).unwrap();
        let surv_a = kill_curve(2.0, 1.0, 2.0).unwrap();
        assert!((surv_combo - surv_a).abs() < 1e-10);
    }

    #[test]
    fn test_combination_kill_curve_both_above() {
        // Both above MIC → survival < either alone
        let surv_a = kill_curve(2.0, 1.0, 1.0).unwrap();
        let surv_b = kill_curve(2.0, 1.0, 1.0).unwrap();
        let surv_combo = combination_kill_curve(2.0, 1.0, 1.0, 2.0, 1.0, 1.0).unwrap();
        assert!(surv_combo < surv_a);
        assert!((surv_combo - surv_a * surv_b).abs() < 1e-10);
    }

    #[test]
    fn test_checkerboard_basic() {
        // Twofold dilutions: 0, 0.25, 0.5, 1.0 × MIC
        let conc_a = [0.0, 0.25, 0.5, 1.0];
        let conc_b = [0.0, 0.25, 0.5, 1.0];
        let result = checkerboard(&conc_a, &conc_b, 1.0, 1.0).unwrap();
        assert_eq!(result.fic_grid.len(), 4);
        assert_eq!(result.fic_grid[0].len(), 4);
        // Corner (0,0) → FIC = 0
        assert!((result.fic_grid[0][0] - 0.0).abs() < 1e-10);
        // (1.0, 1.0) → FIC = 2.0
        assert!((result.fic_grid[3][3] - 2.0).abs() < 1e-10);
        // (0.25, 0.25) → FIC = 0.5
        assert!((result.fic_grid[1][1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_checkerboard_min_fic() {
        let conc_a = [0.125, 0.25, 0.5];
        let conc_b = [0.125, 0.25, 0.5];
        let result = checkerboard(&conc_a, &conc_b, 1.0, 1.0).unwrap();
        // Min is at (0.125, 0.125) → FIC = 0.25
        assert!((result.min_fic - 0.25).abs() < 1e-10);
        assert_eq!(result.interaction, DrugInteraction::Synergy);
    }

    #[test]
    fn test_checkerboard_empty_error() {
        assert!(checkerboard(&[], &[0.5], 1.0, 1.0).is_err());
        assert!(checkerboard(&[0.5], &[], 1.0, 1.0).is_err());
    }

    #[test]
    fn test_checkerboard_invalid_mic() {
        assert!(checkerboard(&[0.5], &[0.5], 0.0, 1.0).is_err());
    }

    #[test]
    fn test_checkerboard_serde_roundtrip() {
        let result = checkerboard(&[0.25, 0.5], &[0.25, 0.5], 1.0, 1.0).unwrap();
        let json = serde_json::to_string(&result).unwrap();
        let back: CheckerboardResult = serde_json::from_str(&json).unwrap();
        assert!((result.min_fic - back.min_fic).abs() < 1e-10);
        assert_eq!(result.interaction, back.interaction);
    }
}
