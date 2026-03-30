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

/// Standard amino acids encoded by the genetic code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[non_exhaustive]
pub enum AminoAcid {
    /// Alanine (A / Ala).
    Alanine,
    /// Arginine (R / Arg).
    Arginine,
    /// Asparagine (N / Asn).
    Asparagine,
    /// Aspartic acid (D / Asp).
    AsparticAcid,
    /// Cysteine (C / Cys).
    Cysteine,
    /// Glutamic acid (E / Glu).
    GlutamicAcid,
    /// Glutamine (Q / Gln).
    Glutamine,
    /// Glycine (G / Gly).
    Glycine,
    /// Histidine (H / His).
    Histidine,
    /// Isoleucine (I / Ile).
    Isoleucine,
    /// Leucine (L / Leu).
    Leucine,
    /// Lysine (K / Lys).
    Lysine,
    /// Methionine (M / Met) — also the standard start codon.
    Methionine,
    /// Phenylalanine (F / Phe).
    Phenylalanine,
    /// Proline (P / Pro).
    Proline,
    /// Serine (S / Ser).
    Serine,
    /// Threonine (T / Thr).
    Threonine,
    /// Tryptophan (W / Trp).
    Tryptophan,
    /// Tyrosine (Y / Tyr).
    Tyrosine,
    /// Valine (V / Val).
    Valine,
    /// Stop codon (translation termination).
    Stop,
}

impl AminoAcid {
    /// Single-letter IUPAC code.
    #[inline]
    #[must_use]
    pub const fn one_letter(self) -> char {
        match self {
            Self::Alanine => 'A',
            Self::Arginine => 'R',
            Self::Asparagine => 'N',
            Self::AsparticAcid => 'D',
            Self::Cysteine => 'C',
            Self::GlutamicAcid => 'E',
            Self::Glutamine => 'Q',
            Self::Glycine => 'G',
            Self::Histidine => 'H',
            Self::Isoleucine => 'I',
            Self::Leucine => 'L',
            Self::Lysine => 'K',
            Self::Methionine => 'M',
            Self::Phenylalanine => 'F',
            Self::Proline => 'P',
            Self::Serine => 'S',
            Self::Threonine => 'T',
            Self::Tryptophan => 'W',
            Self::Tyrosine => 'Y',
            Self::Valine => 'V',
            Self::Stop => '*',
        }
    }

    /// Three-letter abbreviation.
    #[inline]
    #[must_use]
    pub const fn three_letter(self) -> &'static str {
        match self {
            Self::Alanine => "Ala",
            Self::Arginine => "Arg",
            Self::Asparagine => "Asn",
            Self::AsparticAcid => "Asp",
            Self::Cysteine => "Cys",
            Self::GlutamicAcid => "Glu",
            Self::Glutamine => "Gln",
            Self::Glycine => "Gly",
            Self::Histidine => "His",
            Self::Isoleucine => "Ile",
            Self::Leucine => "Leu",
            Self::Lysine => "Lys",
            Self::Methionine => "Met",
            Self::Phenylalanine => "Phe",
            Self::Proline => "Pro",
            Self::Serine => "Ser",
            Self::Threonine => "Thr",
            Self::Tryptophan => "Trp",
            Self::Tyrosine => "Tyr",
            Self::Valine => "Val",
            Self::Stop => "Ter",
        }
    }

    /// Full amino acid name.
    #[inline]
    #[must_use]
    pub const fn full_name(self) -> &'static str {
        match self {
            Self::Alanine => "Alanine",
            Self::Arginine => "Arginine",
            Self::Asparagine => "Asparagine",
            Self::AsparticAcid => "Aspartic acid",
            Self::Cysteine => "Cysteine",
            Self::GlutamicAcid => "Glutamic acid",
            Self::Glutamine => "Glutamine",
            Self::Glycine => "Glycine",
            Self::Histidine => "Histidine",
            Self::Isoleucine => "Isoleucine",
            Self::Leucine => "Leucine",
            Self::Lysine => "Lysine",
            Self::Methionine => "Methionine",
            Self::Phenylalanine => "Phenylalanine",
            Self::Proline => "Proline",
            Self::Serine => "Serine",
            Self::Threonine => "Threonine",
            Self::Tryptophan => "Tryptophan",
            Self::Tyrosine => "Tyrosine",
            Self::Valine => "Valine",
            Self::Stop => "Stop",
        }
    }

    /// Average molecular weight of the free amino acid in daltons (Da).
    ///
    /// Values from standard biochemistry references (NIST / UniProt).
    /// For peptide mass calculation, subtract 18.015 Da (water) per
    /// peptide bond formed.
    /// Stop returns 0.0.
    #[inline]
    #[must_use]
    pub const fn molecular_weight(self) -> f64 {
        match self {
            Self::Alanine => 89.094,
            Self::Arginine => 174.203,
            Self::Asparagine => 132.119,
            Self::AsparticAcid => 133.104,
            Self::Cysteine => 121.159,
            Self::GlutamicAcid => 147.130,
            Self::Glutamine => 146.146,
            Self::Glycine => 75.032,
            Self::Histidine => 155.156,
            Self::Isoleucine => 131.175,
            Self::Leucine => 131.175,
            Self::Lysine => 146.189,
            Self::Methionine => 149.208,
            Self::Phenylalanine => 165.192,
            Self::Proline => 115.132,
            Self::Serine => 105.093,
            Self::Threonine => 119.119,
            Self::Tryptophan => 204.228,
            Self::Tyrosine => 181.191,
            Self::Valine => 117.148,
            Self::Stop => 0.0,
        }
    }

    /// Look up an amino acid from its single-letter IUPAC code.
    ///
    /// # Errors
    ///
    /// Returns error if the character is not a valid amino acid code.
    #[must_use = "returns the amino acid without side effects"]
    pub fn from_one_letter(code: char) -> Result<Self> {
        match code.to_ascii_uppercase() {
            'A' => Ok(Self::Alanine),
            'R' => Ok(Self::Arginine),
            'N' => Ok(Self::Asparagine),
            'D' => Ok(Self::AsparticAcid),
            'C' => Ok(Self::Cysteine),
            'E' => Ok(Self::GlutamicAcid),
            'Q' => Ok(Self::Glutamine),
            'G' => Ok(Self::Glycine),
            'H' => Ok(Self::Histidine),
            'I' => Ok(Self::Isoleucine),
            'L' => Ok(Self::Leucine),
            'K' => Ok(Self::Lysine),
            'M' => Ok(Self::Methionine),
            'F' => Ok(Self::Phenylalanine),
            'P' => Ok(Self::Proline),
            'S' => Ok(Self::Serine),
            'T' => Ok(Self::Threonine),
            'W' => Ok(Self::Tryptophan),
            'Y' => Ok(Self::Tyrosine),
            'V' => Ok(Self::Valine),
            '*' => Ok(Self::Stop),
            _ => Err(JivanuError::ComputationError(format!(
                "unknown amino acid code: {code}"
            ))),
        }
    }
}

/// Amino acid side-chain charge classification at physiological pH (~7.4).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum ChargeClass {
    /// Positively charged at pH 7.4 (Arg, His, Lys).
    Positive,
    /// Negatively charged at pH 7.4 (Asp, Glu).
    Negative,
    /// Uncharged polar (Asn, Cys, Gln, Ser, Thr, Tyr).
    Polar,
    /// Nonpolar / hydrophobic (Ala, Gly, Ile, Leu, Met, Phe, Pro, Trp, Val).
    Nonpolar,
}

impl AminoAcid {
    /// Side-chain charge class at physiological pH (~7.4).
    ///
    /// Histidine (pKa ~6.0) is classified as positive because it is
    /// partially protonated at pH 7.4 and its charge state is
    /// biologically significant for enzyme catalysis.
    #[inline]
    #[must_use]
    pub const fn charge_class(self) -> ChargeClass {
        match self {
            Self::Arginine | Self::Histidine | Self::Lysine => ChargeClass::Positive,
            Self::AsparticAcid | Self::GlutamicAcid => ChargeClass::Negative,
            Self::Asparagine
            | Self::Cysteine
            | Self::Glutamine
            | Self::Serine
            | Self::Threonine
            | Self::Tyrosine => ChargeClass::Polar,
            Self::Alanine
            | Self::Glycine
            | Self::Isoleucine
            | Self::Leucine
            | Self::Methionine
            | Self::Phenylalanine
            | Self::Proline
            | Self::Tryptophan
            | Self::Valine => ChargeClass::Nonpolar,
            Self::Stop => ChargeClass::Nonpolar,
        }
    }

    /// Kyte-Doolittle hydrophobicity index.
    ///
    /// Scale ranges from −4.5 (most hydrophilic, Arg) to +4.5 (most
    /// hydrophobic, Ile). Reference: Kyte & Doolittle (1982) J. Mol. Biol.
    /// 157:105-132.
    ///
    /// Stop returns 0.0.
    #[inline]
    #[must_use]
    pub const fn hydrophobicity(self) -> f64 {
        match self {
            Self::Isoleucine => 4.5,
            Self::Valine => 4.2,
            Self::Leucine => 3.8,
            Self::Phenylalanine => 2.8,
            Self::Cysteine => 2.5,
            Self::Methionine => 1.9,
            Self::Alanine => 1.8,
            Self::Glycine => -0.4,
            Self::Threonine => -0.7,
            Self::Serine => -0.8,
            Self::Tryptophan => -0.9,
            Self::Tyrosine => -1.3,
            Self::Proline => -1.6,
            Self::Histidine => -3.2,
            Self::GlutamicAcid => -3.5,
            Self::Glutamine => -3.5,
            Self::AsparticAcid => -3.5,
            Self::Asparagine => -3.5,
            Self::Lysine => -3.9,
            Self::Arginine => -4.5,
            Self::Stop => 0.0,
        }
    }

    /// Isoelectric point (pI) — pH at which the amino acid carries no net charge.
    ///
    /// Values from standard biochemistry references (Lehninger, Stryer).
    ///
    /// Stop returns 0.0.
    #[inline]
    #[must_use]
    pub const fn isoelectric_point(self) -> f64 {
        match self {
            Self::Alanine => 6.00,
            Self::Arginine => 10.76,
            Self::Asparagine => 5.41,
            Self::AsparticAcid => 2.77,
            Self::Cysteine => 5.07,
            Self::GlutamicAcid => 3.22,
            Self::Glutamine => 5.65,
            Self::Glycine => 5.97,
            Self::Histidine => 7.59,
            Self::Isoleucine => 6.02,
            Self::Leucine => 5.98,
            Self::Lysine => 9.74,
            Self::Methionine => 5.74,
            Self::Phenylalanine => 5.48,
            Self::Proline => 6.30,
            Self::Serine => 5.68,
            Self::Threonine => 5.60,
            Self::Tryptophan => 5.89,
            Self::Tyrosine => 5.66,
            Self::Valine => 5.96,
            Self::Stop => 0.0,
        }
    }
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

/// Translate a DNA codon to its amino acid using the standard genetic code.
///
/// Returns the [`AminoAcid`] encoded by a 3-base DNA codon. All 64 codons
/// of the standard genetic code are covered.
///
/// # Errors
///
/// Returns error if the codon is not exactly 3 valid DNA bases.
#[inline]
#[must_use = "returns the amino acid without side effects"]
pub fn translate_codon_to_aa(codon: &str) -> Result<AminoAcid> {
    if codon.len() != 3 {
        return Err(JivanuError::ComputationError(
            "codon must be exactly 3 bases".into(),
        ));
    }
    let mut buf = [0u8; 3];
    for (i, c) in codon.bytes().enumerate() {
        buf[i] = c.to_ascii_uppercase();
    }
    match &buf {
        // Phenylalanine
        b"TTT" | b"TTC" => Ok(AminoAcid::Phenylalanine),
        // Leucine
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => Ok(AminoAcid::Leucine),
        // Isoleucine
        b"ATT" | b"ATC" | b"ATA" => Ok(AminoAcid::Isoleucine),
        // Methionine (start)
        b"ATG" => Ok(AminoAcid::Methionine),
        // Valine
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => Ok(AminoAcid::Valine),
        // Serine
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => Ok(AminoAcid::Serine),
        // Proline
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => Ok(AminoAcid::Proline),
        // Threonine
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => Ok(AminoAcid::Threonine),
        // Alanine
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => Ok(AminoAcid::Alanine),
        // Tyrosine
        b"TAT" | b"TAC" => Ok(AminoAcid::Tyrosine),
        // Stop codons
        b"TAA" | b"TAG" | b"TGA" => Ok(AminoAcid::Stop),
        // Histidine
        b"CAT" | b"CAC" => Ok(AminoAcid::Histidine),
        // Glutamine
        b"CAA" | b"CAG" => Ok(AminoAcid::Glutamine),
        // Asparagine
        b"AAT" | b"AAC" => Ok(AminoAcid::Asparagine),
        // Lysine
        b"AAA" | b"AAG" => Ok(AminoAcid::Lysine),
        // Aspartic acid
        b"GAT" | b"GAC" => Ok(AminoAcid::AsparticAcid),
        // Glutamic acid
        b"GAA" | b"GAG" => Ok(AminoAcid::GlutamicAcid),
        // Cysteine
        b"TGT" | b"TGC" => Ok(AminoAcid::Cysteine),
        // Tryptophan
        b"TGG" => Ok(AminoAcid::Tryptophan),
        // Arginine
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => Ok(AminoAcid::Arginine),
        // Glycine
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => Ok(AminoAcid::Glycine),
        _ => Err(JivanuError::ComputationError(format!(
            "unknown codon: {}",
            core::str::from_utf8(&buf).unwrap_or("???")
        ))),
    }
}

/// Translate a DNA codon to its single-letter amino acid code.
///
/// Convenience wrapper around [`translate_codon_to_aa`] that returns the
/// IUPAC one-letter code as a `char`.
///
/// # Errors
///
/// Returns error if the codon is not exactly 3 valid DNA bases.
#[inline]
#[must_use = "returns the amino acid without side effects"]
pub fn translate_codon(codon: &str) -> Result<char> {
    translate_codon_to_aa(codon).map(|aa| aa.one_letter())
}

/// Number of codons that encode a given amino acid (codon degeneracy).
///
/// In the standard genetic code, degeneracy ranges from 1 (Met, Trp) to 6
/// (Arg, Leu, Ser).
#[inline]
#[must_use]
pub const fn codon_degeneracy(aa: AminoAcid) -> u8 {
    match aa {
        AminoAcid::Methionine | AminoAcid::Tryptophan => 1,
        AminoAcid::Phenylalanine
        | AminoAcid::Tyrosine
        | AminoAcid::Histidine
        | AminoAcid::Glutamine
        | AminoAcid::Asparagine
        | AminoAcid::Lysine
        | AminoAcid::AsparticAcid
        | AminoAcid::GlutamicAcid
        | AminoAcid::Cysteine => 2,
        AminoAcid::Isoleucine | AminoAcid::Stop => 3,
        AminoAcid::Valine
        | AminoAcid::Proline
        | AminoAcid::Threonine
        | AminoAcid::Alanine
        | AminoAcid::Glycine => 4,
        AminoAcid::Leucine | AminoAcid::Arginine | AminoAcid::Serine => 6,
    }
}

/// Reverse complement of a DNA sequence.
///
/// `A↔T`, `G↔C`. The result is reversed and complemented.
///
/// # Errors
///
/// Returns error if the sequence contains non-DNA characters.
#[must_use = "returns the reverse complement without side effects"]
pub fn reverse_complement(dna: &str) -> Result<String> {
    let mut result = Vec::with_capacity(dna.len());
    for b in dna.bytes().rev() {
        let comp = match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            _ => {
                return Err(JivanuError::ComputationError(format!(
                    "invalid DNA character: {}",
                    b as char
                )));
            }
        };
        result.push(comp);
    }
    // All output bytes are ASCII (A, T, G, C), so this cannot fail.
    String::from_utf8(result)
        .map_err(|e| JivanuError::ComputationError(format!("internal error: invalid UTF-8: {e}")))
}

/// Translate a DNA coding sequence (open reading frame) to a protein sequence.
///
/// Reads codons in-frame from the start of the sequence. Stops at the first
/// stop codon or end of sequence. Partial trailing codons are ignored.
///
/// Returns the amino acid sequence as single-letter IUPAC codes.
///
/// # Errors
///
/// Returns error if any codon contains invalid DNA characters.
#[must_use = "returns the protein sequence without side effects"]
pub fn translate_orf(dna: &str) -> Result<String> {
    let bytes = dna.as_bytes();
    let mut protein = String::with_capacity(bytes.len() / 3);
    let mut i = 0;
    while i + 3 <= bytes.len() {
        let codon = &dna[i..i + 3];
        let aa = translate_codon_to_aa(codon)?;
        if aa == AminoAcid::Stop {
            break;
        }
        protein.push(aa.one_letter());
        i += 3;
    }
    Ok(protein)
}

/// Estimated molecular weight of a protein from its amino acid sequence.
///
/// Sums individual amino acid weights and subtracts water (18.015 Da)
/// for each peptide bond formed.
///
/// Accepts single-letter IUPAC codes. Stop codons (`*`) are ignored.
///
/// # Errors
///
/// Returns error if the sequence is empty or contains invalid characters.
#[must_use = "returns the molecular weight without side effects"]
pub fn protein_molecular_weight(sequence: &str) -> Result<f64> {
    if sequence.is_empty() {
        return Err(JivanuError::ComputationError(
            "protein sequence must not be empty".into(),
        ));
    }
    let mut total = 0.0;
    let mut count = 0u64;
    for c in sequence.chars() {
        if c == '*' {
            continue;
        }
        let aa = AminoAcid::from_one_letter(c)?;
        total += aa.molecular_weight();
        count += 1;
    }
    if count == 0 {
        return Err(JivanuError::ComputationError(
            "protein sequence contains no amino acids".into(),
        ));
    }
    // Subtract water for each peptide bond (count - 1 bonds)
    if count > 1 {
        total -= (count - 1) as f64 * 18.015;
    }
    Ok(total)
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

    #[test]
    fn test_amino_acid_serde_roundtrip_all() {
        let all = [
            AminoAcid::Alanine,
            AminoAcid::Arginine,
            AminoAcid::Asparagine,
            AminoAcid::AsparticAcid,
            AminoAcid::Cysteine,
            AminoAcid::GlutamicAcid,
            AminoAcid::Glutamine,
            AminoAcid::Glycine,
            AminoAcid::Histidine,
            AminoAcid::Isoleucine,
            AminoAcid::Leucine,
            AminoAcid::Lysine,
            AminoAcid::Methionine,
            AminoAcid::Phenylalanine,
            AminoAcid::Proline,
            AminoAcid::Serine,
            AminoAcid::Threonine,
            AminoAcid::Tryptophan,
            AminoAcid::Tyrosine,
            AminoAcid::Valine,
            AminoAcid::Stop,
        ];
        for aa in &all {
            let json = serde_json::to_string(aa).unwrap();
            let back: AminoAcid = serde_json::from_str(&json).unwrap();
            assert_eq!(*aa, back);
        }
    }

    #[test]
    fn test_amino_acid_one_letter_codes() {
        assert_eq!(AminoAcid::Alanine.one_letter(), 'A');
        assert_eq!(AminoAcid::Tryptophan.one_letter(), 'W');
        assert_eq!(AminoAcid::Stop.one_letter(), '*');
    }

    #[test]
    fn test_amino_acid_three_letter_codes() {
        assert_eq!(AminoAcid::Alanine.three_letter(), "Ala");
        assert_eq!(AminoAcid::Tryptophan.three_letter(), "Trp");
        assert_eq!(AminoAcid::Stop.three_letter(), "Ter");
    }

    #[test]
    fn test_amino_acid_full_name() {
        assert_eq!(AminoAcid::AsparticAcid.full_name(), "Aspartic acid");
        assert_eq!(AminoAcid::GlutamicAcid.full_name(), "Glutamic acid");
    }

    #[test]
    fn test_amino_acid_molecular_weight_known_values() {
        // Glycine: lightest amino acid, ~75 Da
        assert!((AminoAcid::Glycine.molecular_weight() - 75.032).abs() < 0.001);
        // Tryptophan: heaviest standard amino acid, ~204 Da
        assert!((AminoAcid::Tryptophan.molecular_weight() - 204.228).abs() < 0.001);
        // Leucine and Isoleucine are isomers: same MW
        assert!(
            (AminoAcid::Leucine.molecular_weight() - AminoAcid::Isoleucine.molecular_weight())
                .abs()
                < 1e-10
        );
        // Stop has no weight
        assert!((AminoAcid::Stop.molecular_weight() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_amino_acid_from_one_letter_roundtrip() {
        let all = [
            'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T',
            'W', 'Y', 'V', '*',
        ];
        for code in &all {
            let aa = AminoAcid::from_one_letter(*code).unwrap();
            assert_eq!(aa.one_letter(), *code);
        }
    }

    #[test]
    fn test_amino_acid_from_one_letter_case_insensitive() {
        assert_eq!(AminoAcid::from_one_letter('a').unwrap(), AminoAcid::Alanine);
        assert_eq!(
            AminoAcid::from_one_letter('w').unwrap(),
            AminoAcid::Tryptophan
        );
    }

    #[test]
    fn test_amino_acid_from_one_letter_invalid() {
        assert!(AminoAcid::from_one_letter('X').is_err());
        assert!(AminoAcid::from_one_letter('Z').is_err());
    }

    #[test]
    fn test_translate_codon_to_aa() {
        assert_eq!(translate_codon_to_aa("ATG").unwrap(), AminoAcid::Methionine);
        assert_eq!(translate_codon_to_aa("TAA").unwrap(), AminoAcid::Stop);
        assert_eq!(translate_codon_to_aa("TGG").unwrap(), AminoAcid::Tryptophan);
        assert_eq!(
            translate_codon_to_aa("TTT").unwrap(),
            AminoAcid::Phenylalanine
        );
    }

    #[test]
    fn test_translate_codon_to_aa_case_insensitive() {
        assert_eq!(translate_codon_to_aa("atg").unwrap(), AminoAcid::Methionine);
        assert_eq!(translate_codon_to_aa("Atg").unwrap(), AminoAcid::Methionine);
    }

    #[test]
    fn test_translate_codon_to_aa_invalid() {
        assert!(translate_codon_to_aa("XX").is_err());
        assert!(translate_codon_to_aa("ATGC").is_err());
    }

    #[test]
    fn test_translate_codon_backwards_compat() {
        // Old char API still works, consistent with new AA API
        assert_eq!(translate_codon("ATG").unwrap(), 'M');
        assert_eq!(
            translate_codon("ATG").unwrap(),
            translate_codon_to_aa("ATG").unwrap().one_letter()
        );
    }

    #[test]
    fn test_all_64_codons_covered() {
        let bases = ['A', 'T', 'G', 'C'];
        let mut count = 0;
        for &b1 in &bases {
            for &b2 in &bases {
                for &b3 in &bases {
                    let codon = format!("{b1}{b2}{b3}");
                    let result = translate_codon_to_aa(&codon);
                    assert!(
                        result.is_ok(),
                        "codon {codon} not in table: {:?}",
                        result.err()
                    );
                    count += 1;
                }
            }
        }
        assert_eq!(count, 64);
    }

    #[test]
    fn test_codon_degeneracy_known_values() {
        // Met and Trp have exactly 1 codon each
        assert_eq!(codon_degeneracy(AminoAcid::Methionine), 1);
        assert_eq!(codon_degeneracy(AminoAcid::Tryptophan), 1);
        // Leu, Arg, Ser have 6 codons each
        assert_eq!(codon_degeneracy(AminoAcid::Leucine), 6);
        assert_eq!(codon_degeneracy(AminoAcid::Arginine), 6);
        assert_eq!(codon_degeneracy(AminoAcid::Serine), 6);
        // Stop has 3 codons (TAA, TAG, TGA)
        assert_eq!(codon_degeneracy(AminoAcid::Stop), 3);
    }

    #[test]
    fn test_codon_degeneracy_sums_to_64() {
        let all = [
            AminoAcid::Alanine,
            AminoAcid::Arginine,
            AminoAcid::Asparagine,
            AminoAcid::AsparticAcid,
            AminoAcid::Cysteine,
            AminoAcid::GlutamicAcid,
            AminoAcid::Glutamine,
            AminoAcid::Glycine,
            AminoAcid::Histidine,
            AminoAcid::Isoleucine,
            AminoAcid::Leucine,
            AminoAcid::Lysine,
            AminoAcid::Methionine,
            AminoAcid::Phenylalanine,
            AminoAcid::Proline,
            AminoAcid::Serine,
            AminoAcid::Threonine,
            AminoAcid::Tryptophan,
            AminoAcid::Tyrosine,
            AminoAcid::Valine,
            AminoAcid::Stop,
        ];
        let total: u32 = all.iter().map(|aa| codon_degeneracy(*aa) as u32).sum();
        assert_eq!(total, 64, "total codon count must equal 64");
    }

    #[test]
    fn test_amino_acid_charge_classes() {
        assert_eq!(AminoAcid::Arginine.charge_class(), ChargeClass::Positive);
        assert_eq!(AminoAcid::Histidine.charge_class(), ChargeClass::Positive);
        assert_eq!(AminoAcid::Lysine.charge_class(), ChargeClass::Positive);
        assert_eq!(
            AminoAcid::AsparticAcid.charge_class(),
            ChargeClass::Negative
        );
        assert_eq!(
            AminoAcid::GlutamicAcid.charge_class(),
            ChargeClass::Negative
        );
        assert_eq!(AminoAcid::Serine.charge_class(), ChargeClass::Polar);
        assert_eq!(AminoAcid::Alanine.charge_class(), ChargeClass::Nonpolar);
        assert_eq!(AminoAcid::Isoleucine.charge_class(), ChargeClass::Nonpolar);
    }

    #[test]
    fn test_hydrophobicity_extremes() {
        // Isoleucine is the most hydrophobic
        assert!((AminoAcid::Isoleucine.hydrophobicity() - 4.5).abs() < 1e-10);
        // Arginine is the most hydrophilic
        assert!((AminoAcid::Arginine.hydrophobicity() - (-4.5)).abs() < 1e-10);
    }

    #[test]
    fn test_hydrophobicity_scale_range() {
        let all = [
            AminoAcid::Alanine,
            AminoAcid::Arginine,
            AminoAcid::Asparagine,
            AminoAcid::AsparticAcid,
            AminoAcid::Cysteine,
            AminoAcid::GlutamicAcid,
            AminoAcid::Glutamine,
            AminoAcid::Glycine,
            AminoAcid::Histidine,
            AminoAcid::Isoleucine,
            AminoAcid::Leucine,
            AminoAcid::Lysine,
            AminoAcid::Methionine,
            AminoAcid::Phenylalanine,
            AminoAcid::Proline,
            AminoAcid::Serine,
            AminoAcid::Threonine,
            AminoAcid::Tryptophan,
            AminoAcid::Tyrosine,
            AminoAcid::Valine,
        ];
        for aa in &all {
            let h = aa.hydrophobicity();
            assert!(
                (-4.5..=4.5).contains(&h),
                "{:?} hydrophobicity {h} out of Kyte-Doolittle range",
                aa
            );
        }
    }

    #[test]
    fn test_isoelectric_point_known_values() {
        // Asp and Glu: acidic, pI < 4
        assert!(AminoAcid::AsparticAcid.isoelectric_point() < 4.0);
        assert!(AminoAcid::GlutamicAcid.isoelectric_point() < 4.0);
        // Arg and Lys: basic, pI > 9
        assert!(AminoAcid::Arginine.isoelectric_point() > 9.0);
        assert!(AminoAcid::Lysine.isoelectric_point() > 9.0);
        // His: weakly basic, pI ~7.6
        assert!((AminoAcid::Histidine.isoelectric_point() - 7.59).abs() < 0.01);
    }

    #[test]
    fn test_charge_class_serde_roundtrip() {
        let cc = ChargeClass::Positive;
        let json = serde_json::to_string(&cc).unwrap();
        let back: ChargeClass = serde_json::from_str(&json).unwrap();
        assert_eq!(cc, back);
    }

    // --- Sequence utility tests ---

    #[test]
    fn test_reverse_complement_simple() {
        assert_eq!(reverse_complement("ATGC").unwrap(), "GCAT");
    }

    #[test]
    fn test_reverse_complement_palindrome() {
        // ATAT → reverse = TATA → complement = ATAT
        assert_eq!(reverse_complement("ATAT").unwrap(), "ATAT");
    }

    #[test]
    fn test_reverse_complement_case_insensitive() {
        assert_eq!(reverse_complement("atgc").unwrap(), "GCAT");
    }

    #[test]
    fn test_reverse_complement_invalid() {
        assert!(reverse_complement("ATXG").is_err());
    }

    #[test]
    fn test_translate_orf_met_only() {
        assert_eq!(translate_orf("ATG").unwrap(), "M");
    }

    #[test]
    fn test_translate_orf_with_stop() {
        // ATG GCT TAA → M A (stop)
        assert_eq!(translate_orf("ATGGCTTAA").unwrap(), "MA");
    }

    #[test]
    fn test_translate_orf_no_stop() {
        // ATG GCT GCT → M A A (no stop, reads to end)
        assert_eq!(translate_orf("ATGGCTGCT").unwrap(), "MAA");
    }

    #[test]
    fn test_translate_orf_partial_codon_ignored() {
        // ATG GC → M (trailing 2 bases ignored)
        assert_eq!(translate_orf("ATGGC").unwrap(), "M");
    }

    #[test]
    fn test_translate_orf_empty() {
        assert_eq!(translate_orf("").unwrap(), "");
    }

    #[test]
    fn test_translate_orf_invalid() {
        assert!(translate_orf("XYZ").is_err());
    }

    #[test]
    fn test_protein_molecular_weight_single_aa() {
        // Single glycine: no peptide bonds
        let mw = protein_molecular_weight("G").unwrap();
        assert!((mw - 75.032).abs() < 0.001);
    }

    #[test]
    fn test_protein_molecular_weight_dipeptide() {
        // G-A dipeptide: 75.032 + 89.094 - 18.015 = 146.111
        let mw = protein_molecular_weight("GA").unwrap();
        assert!((mw - 146.111).abs() < 0.001);
    }

    #[test]
    fn test_protein_molecular_weight_ignores_stop() {
        let mw_with = protein_molecular_weight("GA*").unwrap();
        let mw_without = protein_molecular_weight("GA").unwrap();
        assert!((mw_with - mw_without).abs() < 1e-10);
    }

    #[test]
    fn test_protein_molecular_weight_empty() {
        assert!(protein_molecular_weight("").is_err());
    }

    #[test]
    fn test_protein_molecular_weight_invalid() {
        assert!(protein_molecular_weight("GXA").is_err());
    }

    #[test]
    fn test_translate_orf_to_protein_mw_pipeline() {
        // Full pipeline: DNA → protein → MW
        let protein = translate_orf("ATGGCTTAA").unwrap(); // MA
        let mw = protein_molecular_weight(&protein).unwrap();
        // M (149.208) + A (89.094) - water (18.015) = 220.287
        assert!((mw - 220.287).abs() < 0.001);
    }

    #[test]
    fn test_amino_acid_ordering() {
        // Alphabetical by variant name
        assert!(AminoAcid::Alanine < AminoAcid::Arginine);
        assert!(AminoAcid::Valine < AminoAcid::Stop);
    }

    #[test]
    fn test_amino_acid_count() {
        // 20 standard amino acids + Stop = 21 variants
        let all = [
            AminoAcid::Alanine,
            AminoAcid::Arginine,
            AminoAcid::Asparagine,
            AminoAcid::AsparticAcid,
            AminoAcid::Cysteine,
            AminoAcid::GlutamicAcid,
            AminoAcid::Glutamine,
            AminoAcid::Glycine,
            AminoAcid::Histidine,
            AminoAcid::Isoleucine,
            AminoAcid::Leucine,
            AminoAcid::Lysine,
            AminoAcid::Methionine,
            AminoAcid::Phenylalanine,
            AminoAcid::Proline,
            AminoAcid::Serine,
            AminoAcid::Threonine,
            AminoAcid::Tryptophan,
            AminoAcid::Tyrosine,
            AminoAcid::Valine,
            AminoAcid::Stop,
        ];
        assert_eq!(all.len(), 21);
        // All unique
        let mut set = std::collections::HashSet::new();
        for aa in &all {
            assert!(set.insert(aa));
        }
    }
}
