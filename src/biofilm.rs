//! Biofilm — formation stages, quorum sensing, nutrient diffusion.

use serde::{Deserialize, Serialize};

use crate::error::{validate_finite, validate_non_negative, validate_positive, Result};

/// Stages of biofilm development.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[non_exhaustive]
pub enum BiofilmStage {
    /// Initial reversible attachment to surface.
    Attachment,
    /// Formation of small clusters.
    Microcolony,
    /// Three-dimensional structure with channels.
    Maturation,
    /// Cells detach and colonize new surfaces.
    Dispersal,
}

/// Quorum sensing: does signal concentration exceed threshold?
///
/// Quorum sensing activates gene expression when the autoinducer
/// concentration reaches a critical threshold.
#[inline]
#[must_use = "returns whether quorum is reached without side effects"]
pub fn quorum_sensing(signal_concentration: f64, threshold: f64) -> Result<bool> {
    validate_non_negative(signal_concentration, "signal_concentration")?;
    validate_positive(threshold, "threshold")?;
    Ok(signal_concentration >= threshold)
}

/// Nutrient diffusion through biofilm matrix (Fick's first law application).
///
/// `J = -D * dC/dx ≈ D * C_surface / thickness`
///
/// Returns the flux (concentration/time).
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the diffusion flux without side effects"]
pub fn diffusion_through_matrix(
    nutrient: f64,
    thickness: f64,
    diffusivity: f64,
) -> Result<f64> {
    validate_non_negative(nutrient, "nutrient")?;
    validate_positive(thickness, "thickness")?;
    validate_positive(diffusivity, "diffusivity")?;
    Ok(diffusivity * nutrient / thickness)
}

/// Attachment rate as a function of surface energy and flow rate.
///
/// Higher surface energy promotes attachment; higher flow opposes it.
///
/// # Errors
///
/// Returns error if parameters are non-finite.
#[inline]
#[must_use = "returns the attachment rate without side effects"]
pub fn attachment_rate(surface_energy: f64, flow_rate: f64) -> Result<f64> {
    validate_finite(surface_energy, "surface_energy")?;
    validate_non_negative(flow_rate, "flow_rate")?;
    let se = surface_energy.clamp(0.0, 1.0);
    // Simple model: attachment increases with surface energy, decreases with flow
    Ok(se / (1.0 + flow_rate))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quorum_sensing_above() {
        assert!(quorum_sensing(10.0, 5.0).unwrap());
    }

    #[test]
    fn test_quorum_sensing_below() {
        assert!(!quorum_sensing(3.0, 5.0).unwrap());
    }

    #[test]
    fn test_quorum_sensing_at_threshold() {
        assert!(quorum_sensing(5.0, 5.0).unwrap());
    }

    #[test]
    fn test_diffusion() {
        let flux = diffusion_through_matrix(10.0, 2.0, 0.5).unwrap();
        assert!((flux - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_attachment_rate_high_energy() {
        let rate = attachment_rate(1.0, 0.0).unwrap();
        assert!((rate - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_attachment_rate_high_flow() {
        let low = attachment_rate(0.5, 10.0).unwrap();
        let high = attachment_rate(0.5, 0.0).unwrap();
        assert!(low < high);
    }

    #[test]
    fn test_biofilm_stage_ordering() {
        assert!(BiofilmStage::Attachment < BiofilmStage::Microcolony);
        assert!(BiofilmStage::Microcolony < BiofilmStage::Maturation);
    }

    #[test]
    fn test_biofilm_stage_serde_roundtrip() {
        let stage = BiofilmStage::Maturation;
        let json = serde_json::to_string(&stage).unwrap();
        let back: BiofilmStage = serde_json::from_str(&json).unwrap();
        assert_eq!(stage, back);
    }
}
