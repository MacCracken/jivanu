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

impl Domain {
    /// Whether organisms in this domain have a nucleus.
    #[inline]
    #[must_use]
    pub const fn is_eukaryotic(self) -> bool {
        matches!(self, Self::Eukarya)
    }

    /// Whether organisms in this domain are prokaryotic.
    #[inline]
    #[must_use]
    pub const fn is_prokaryotic(self) -> bool {
        matches!(self, Self::Bacteria | Self::Archaea)
    }
}

impl GramStain {
    /// Whether the cell wall contains an outer membrane.
    ///
    /// Gram-negative bacteria have an outer membrane with lipopolysaccharide (LPS).
    #[inline]
    #[must_use]
    pub const fn has_outer_membrane(self) -> bool {
        matches!(self, Self::Negative)
    }
}

impl CellShape {
    /// Whether the morphology is elongated (rod-like or spiral).
    #[inline]
    #[must_use]
    pub const fn is_elongated(self) -> bool {
        matches!(
            self,
            Self::Bacillus | Self::Spirillum | Self::Vibrio | Self::Spirochete
        )
    }
}

impl OxygenRequirement {
    /// Whether the organism can grow in the presence of oxygen.
    #[inline]
    #[must_use]
    pub const fn tolerates_oxygen(self) -> bool {
        matches!(
            self,
            Self::ObligateAerobe | Self::Facultative | Self::Microaerophilic
        )
    }

    /// Whether the organism can grow without oxygen.
    #[inline]
    #[must_use]
    pub const fn tolerates_anaerobic(self) -> bool {
        matches!(self, Self::ObligateAnaerobe | Self::Facultative)
    }
}

impl MicrobialProfile {
    /// Whether this organism is likely susceptible to beta-lactam antibiotics.
    ///
    /// Beta-lactams target peptidoglycan synthesis. Gram-positive bacteria
    /// (thick peptidoglycan) are generally more susceptible. Archaea and
    /// eukaryotes lack peptidoglycan entirely.
    #[inline]
    #[must_use]
    pub fn beta_lactam_susceptible(&self) -> bool {
        matches!(self.domain, Domain::Bacteria) && matches!(self.gram_stain, GramStain::Positive)
    }

    /// Whether this organism can form endospores.
    ///
    /// Endospore formation is restricted to Gram-positive, rod-shaped
    /// bacteria (Bacillus, Clostridium). This is a simplified heuristic.
    #[inline]
    #[must_use]
    pub fn can_form_endospores(&self) -> bool {
        matches!(self.domain, Domain::Bacteria)
            && matches!(self.gram_stain, GramStain::Positive)
            && matches!(self.shape, CellShape::Bacillus)
    }

    /// Whether this organism is likely motile.
    ///
    /// Spiral/vibrio morphologies strongly correlate with flagellar motility.
    /// Cocci are typically nonmotile. This is a simplified heuristic.
    #[inline]
    #[must_use]
    pub fn likely_motile(&self) -> bool {
        matches!(
            self.shape,
            CellShape::Spirillum | CellShape::Vibrio | CellShape::Spirochete
        )
    }
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
    fn test_domain_eukaryotic() {
        assert!(Domain::Eukarya.is_eukaryotic());
        assert!(!Domain::Bacteria.is_eukaryotic());
        assert!(!Domain::Archaea.is_eukaryotic());
    }

    #[test]
    fn test_domain_prokaryotic() {
        assert!(Domain::Bacteria.is_prokaryotic());
        assert!(Domain::Archaea.is_prokaryotic());
        assert!(!Domain::Eukarya.is_prokaryotic());
    }

    #[test]
    fn test_gram_stain_outer_membrane() {
        assert!(GramStain::Negative.has_outer_membrane());
        assert!(!GramStain::Positive.has_outer_membrane());
    }

    #[test]
    fn test_cell_shape_elongated() {
        assert!(CellShape::Bacillus.is_elongated());
        assert!(CellShape::Spirillum.is_elongated());
        assert!(!CellShape::Coccus.is_elongated());
    }

    #[test]
    fn test_oxygen_tolerates() {
        assert!(OxygenRequirement::Facultative.tolerates_oxygen());
        assert!(OxygenRequirement::Facultative.tolerates_anaerobic());
        assert!(!OxygenRequirement::ObligateAerobe.tolerates_anaerobic());
        assert!(!OxygenRequirement::ObligateAnaerobe.tolerates_oxygen());
    }

    #[test]
    fn test_beta_lactam_susceptible() {
        let gram_pos = MicrobialProfile {
            domain: Domain::Bacteria,
            gram_stain: GramStain::Positive,
            shape: CellShape::Coccus,
            oxygen: OxygenRequirement::Facultative,
        };
        assert!(gram_pos.beta_lactam_susceptible());

        let gram_neg = MicrobialProfile {
            domain: Domain::Bacteria,
            gram_stain: GramStain::Negative,
            shape: CellShape::Bacillus,
            oxygen: OxygenRequirement::Facultative,
        };
        assert!(!gram_neg.beta_lactam_susceptible());
    }

    #[test]
    fn test_can_form_endospores() {
        let bacillus = MicrobialProfile {
            domain: Domain::Bacteria,
            gram_stain: GramStain::Positive,
            shape: CellShape::Bacillus,
            oxygen: OxygenRequirement::Facultative,
        };
        assert!(bacillus.can_form_endospores());

        let coccus = MicrobialProfile {
            domain: Domain::Bacteria,
            gram_stain: GramStain::Positive,
            shape: CellShape::Coccus,
            oxygen: OxygenRequirement::Facultative,
        };
        assert!(!coccus.can_form_endospores());
    }

    #[test]
    fn test_likely_motile() {
        let spirillum = MicrobialProfile {
            domain: Domain::Bacteria,
            gram_stain: GramStain::Negative,
            shape: CellShape::Spirillum,
            oxygen: OxygenRequirement::Microaerophilic,
        };
        assert!(spirillum.likely_motile());

        let coccus = MicrobialProfile {
            domain: Domain::Bacteria,
            gram_stain: GramStain::Positive,
            shape: CellShape::Coccus,
            oxygen: OxygenRequirement::Facultative,
        };
        assert!(!coccus.likely_motile());
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
