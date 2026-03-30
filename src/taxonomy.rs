//! Taxonomy — microbial classification, morphology, metabolic types.

use serde::{Deserialize, Serialize};

/// Three domains of life.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Domain {
    /// True bacteria.
    Bacteria,
    /// Archaea (extremophiles, methanogens).
    Archaea,
    /// Eukaryotes (fungi, protists, plants, animals).
    Eukarya,
}

/// Gram stain classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum GramStain {
    /// Thick peptidoglycan layer, retains crystal violet.
    Positive,
    /// Thin peptidoglycan, outer membrane, pink safranin.
    Negative,
}

/// Cell shape (morphology).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum CellShape {
    /// Spherical.
    Coccus,
    /// Rod-shaped.
    Bacillus,
    /// Rigid spiral.
    Spirillum,
    /// Comma-shaped.
    Vibrio,
    /// Flexible spiral with axial filament.
    Spirochete,
}

/// Oxygen requirement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum OxygenRequirement {
    /// Requires oxygen for growth.
    ObligateAerobe,
    /// Cannot tolerate oxygen.
    ObligateAnaerobe,
    /// Can grow with or without oxygen.
    Facultative,
    /// Requires low oxygen levels.
    Microaerophilic,
}

/// Basic microbial organism profile.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MicrobialProfile {
    /// Domain of life.
    pub domain: Domain,
    /// Gram stain result.
    pub gram_stain: GramStain,
    /// Cell morphology.
    pub shape: CellShape,
    /// Oxygen requirement.
    pub oxygen: OxygenRequirement,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_domain_serde_roundtrip() {
        let d = Domain::Bacteria;
        let json = serde_json::to_string(&d).unwrap();
        let back: Domain = serde_json::from_str(&json).unwrap();
        assert_eq!(d, back);
    }

    #[test]
    fn test_gram_stain_serde_roundtrip() {
        let gs = GramStain::Negative;
        let json = serde_json::to_string(&gs).unwrap();
        let back: GramStain = serde_json::from_str(&json).unwrap();
        assert_eq!(gs, back);
    }

    #[test]
    fn test_cell_shape_serde_roundtrip() {
        let shape = CellShape::Spirochete;
        let json = serde_json::to_string(&shape).unwrap();
        let back: CellShape = serde_json::from_str(&json).unwrap();
        assert_eq!(shape, back);
    }

    #[test]
    fn test_oxygen_requirement_serde_roundtrip() {
        let oxy = OxygenRequirement::Facultative;
        let json = serde_json::to_string(&oxy).unwrap();
        let back: OxygenRequirement = serde_json::from_str(&json).unwrap();
        assert_eq!(oxy, back);
    }

    #[test]
    fn test_microbial_profile_serde_roundtrip() {
        let profile = MicrobialProfile {
            domain: Domain::Bacteria,
            gram_stain: GramStain::Negative,
            shape: CellShape::Bacillus,
            oxygen: OxygenRequirement::Facultative,
        };
        let json = serde_json::to_string(&profile).unwrap();
        let back: MicrobialProfile = serde_json::from_str(&json).unwrap();
        assert_eq!(profile.domain, back.domain);
        assert_eq!(profile.shape, back.shape);
    }
}
