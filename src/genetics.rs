//! Microbial genetics — mutation rates, Hardy-Weinberg, GC content, gene transfer.

use serde::{Deserialize, Serialize};

use crate::error::{JivanuError, Result, validate_finite, validate_positive};

/// Mutation rate per base per generation.
///
/// `rate = mutations / (bases * generations)`
///
/// # Errors
///
/// Returns error if bases or generations is zero.
#[inline]
#[must_use = "returns the mutation rate without side effects"]
pub fn mutation_rate(mutations: u64, bases: u64, generations: u64) -> Result<f64> {
    if bases == 0 {
        return Err(JivanuError::ComputationError("bases must be > 0".into()));
    }
    if generations == 0 {
        return Err(JivanuError::ComputationError(
            "generations must be > 0".into(),
        ));
    }
    Ok(mutations as f64 / (bases as f64 * generations as f64))
}

/// Hardy-Weinberg equilibrium genotype frequencies.
///
/// Given allele frequency `p` (dominant), computes:
/// - `p^2` (homozygous dominant)
/// - `2pq` (heterozygous)
/// - `q^2` (homozygous recessive)
///
/// where `q = 1 - p`.
///
/// # Errors
///
/// Returns error if p is outside [0, 1].
#[must_use = "returns genotype frequencies (p2, 2pq, q2) without side effects"]
pub fn hardy_weinberg(p: f64) -> Result<(f64, f64, f64)> {
    validate_finite(p, "p")?;
    if !(0.0..=1.0).contains(&p) {
        return Err(JivanuError::ComputationError(
            "allele frequency p must be in [0, 1]".into(),
        ));
    }
    let q = 1.0 - p;
    Ok((p * p, 2.0 * p * q, q * q))
}

/// GC content: fraction of guanine and cytosine bases in a DNA sequence.
///
/// # Errors
///
/// Returns error if the sequence is empty or contains non-DNA characters.
#[must_use = "returns the GC content fraction without side effects"]
pub fn gc_content(dna: &str) -> Result<f64> {
    if dna.is_empty() {
        return Err(JivanuError::ComputationError(
            "DNA sequence must not be empty".into(),
        ));
    }
    let mut gc = 0usize;
    let mut total = 0usize;
    for c in dna.chars() {
        match c.to_ascii_uppercase() {
            'G' | 'C' => {
                gc += 1;
                total += 1;
            }
            'A' | 'T' => {
                total += 1;
            }
            _ => {
                return Err(JivanuError::ComputationError(format!(
                    "invalid DNA character: {c}"
                )));
            }
        }
    }
    Ok(gc as f64 / total as f64)
}

/// Selection coefficient: relative fitness difference.
///
/// `s = (w_mutant - w_wildtype) / w_wildtype`
///
/// s > 0: mutant is fitter. s < 0: mutant is less fit.
///
/// # Errors
///
/// Returns error if wildtype fitness is non-positive.
#[inline]
#[must_use = "returns the selection coefficient without side effects"]
pub fn selection_coefficient(fitness_mutant: f64, fitness_wildtype: f64) -> Result<f64> {
    validate_finite(fitness_mutant, "fitness_mutant")?;
    validate_positive(fitness_wildtype, "fitness_wildtype")?;
    Ok((fitness_mutant - fitness_wildtype) / fitness_wildtype)
}

/// Horizontal gene transfer mechanisms.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum GeneTransferMechanism {
    /// Direct cell-to-cell DNA transfer via pili.
    Conjugation,
    /// Phage-mediated DNA transfer.
    Transduction,
    /// Uptake of free DNA from the environment.
    Transformation,
}

/// Simple codon table for standard genetic code.
///
/// Returns the single-letter amino acid code for a DNA codon (3 bases).
///
/// # Errors
///
/// Returns error if the codon is invalid.
#[must_use = "returns the amino acid without side effects"]
pub fn translate_codon(codon: &str) -> Result<char> {
    if codon.len() != 3 {
        return Err(JivanuError::ComputationError(
            "codon must be exactly 3 bases".into(),
        ));
    }
    let codon_upper: String = codon.chars().map(|c| c.to_ascii_uppercase()).collect();
    match codon_upper.as_str() {
        // Phenylalanine
        "TTT" | "TTC" => Ok('F'),
        // Leucine
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => Ok('L'),
        // Isoleucine
        "ATT" | "ATC" | "ATA" => Ok('I'),
        // Methionine (start)
        "ATG" => Ok('M'),
        // Valine
        "GTT" | "GTC" | "GTA" | "GTG" => Ok('V'),
        // Serine
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => Ok('S'),
        // Proline
        "CCT" | "CCC" | "CCA" | "CCG" => Ok('P'),
        // Threonine
        "ACT" | "ACC" | "ACA" | "ACG" => Ok('T'),
        // Alanine
        "GCT" | "GCC" | "GCA" | "GCG" => Ok('A'),
        // Tyrosine
        "TAT" | "TAC" => Ok('Y'),
        // Stop codons
        "TAA" | "TAG" | "TGA" => Ok('*'),
        // Histidine
        "CAT" | "CAC" => Ok('H'),
        // Glutamine
        "CAA" | "CAG" => Ok('Q'),
        // Asparagine
        "AAT" | "AAC" => Ok('N'),
        // Lysine
        "AAA" | "AAG" => Ok('K'),
        // Aspartic acid
        "GAT" | "GAC" => Ok('D'),
        // Glutamic acid
        "GAA" | "GAG" => Ok('E'),
        // Cysteine
        "TGT" | "TGC" => Ok('C'),
        // Tryptophan
        "TGG" => Ok('W'),
        // Arginine
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => Ok('R'),
        // Glycine
        "GGT" | "GGC" | "GGA" | "GGG" => Ok('G'),
        _ => Err(JivanuError::ComputationError(format!(
            "unknown codon: {codon_upper}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mutation_rate() {
        let rate = mutation_rate(10, 1000, 100).unwrap();
        assert!((rate - 0.0001).abs() < 1e-10);
    }

    #[test]
    fn test_mutation_rate_zero_bases() {
        assert!(mutation_rate(1, 0, 1).is_err());
    }

    #[test]
    fn test_hardy_weinberg_p06() {
        let (p2, pq2, q2) = hardy_weinberg(0.6).unwrap();
        assert!((p2 - 0.36).abs() < 1e-10);
        assert!((pq2 - 0.48).abs() < 1e-10);
        assert!((q2 - 0.16).abs() < 1e-10);
    }

    #[test]
    fn test_hardy_weinberg_sum_to_one() {
        let (p2, pq2, q2) = hardy_weinberg(0.3).unwrap();
        assert!((p2 + pq2 + q2 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hardy_weinberg_invalid_p() {
        assert!(hardy_weinberg(1.5).is_err());
        assert!(hardy_weinberg(-0.1).is_err());
    }

    #[test]
    fn test_gc_content_balanced() {
        let gc = gc_content("ATGC").unwrap();
        assert!((gc - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_gc_content_all_gc() {
        let gc = gc_content("GGCC").unwrap();
        assert!((gc - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gc_content_all_at() {
        let gc = gc_content("AATT").unwrap();
        assert!((gc - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_gc_content_empty() {
        assert!(gc_content("").is_err());
    }

    #[test]
    fn test_gc_content_invalid() {
        assert!(gc_content("ATXG").is_err());
    }

    #[test]
    fn test_selection_coefficient_neutral() {
        let s = selection_coefficient(1.0, 1.0).unwrap();
        assert!((s - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_selection_coefficient_beneficial() {
        let s = selection_coefficient(1.1, 1.0).unwrap();
        assert!((s - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_translate_codon_met() {
        assert_eq!(translate_codon("ATG").unwrap(), 'M');
    }

    #[test]
    fn test_translate_codon_stop() {
        assert_eq!(translate_codon("TAA").unwrap(), '*');
    }

    #[test]
    fn test_translate_codon_invalid() {
        assert!(translate_codon("XX").is_err());
    }

    #[test]
    fn test_gene_transfer_serde_roundtrip() {
        let mech = GeneTransferMechanism::Conjugation;
        let json = serde_json::to_string(&mech).unwrap();
        let back: GeneTransferMechanism = serde_json::from_str(&json).unwrap();
        assert_eq!(mech, back);
    }
}
