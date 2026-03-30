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
}
