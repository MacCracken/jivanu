//! Integration tests for jivanu.

use jivanu::biofilm;
use jivanu::epidemiology;
use jivanu::genetics;
use jivanu::growth;
use jivanu::metabolism;
use jivanu::pharmacokinetics;
use jivanu::resistance;
use jivanu::taxonomy;

// --- Growth ---

#[test]
fn test_monod_at_high_substrate() {
    let mu = growth::monod_kinetics(1000.0, 0.5, 1.0).unwrap();
    assert!((mu - 0.5).abs() < 0.001);
}

#[test]
fn test_doubling_time_1hr() {
    let td = growth::doubling_time(core::f64::consts::LN_2).unwrap();
    assert!((td - 1.0).abs() < 1e-10);
}

#[test]
fn test_exponential_growth_doubles() {
    let n = growth::exponential_growth(100.0, core::f64::consts::LN_2, 1.0).unwrap();
    assert!((n - 200.0).abs() < 1e-6);
}

#[test]
fn test_competition_outcome_coexistence() {
    let params = growth::CompetitionParams {
        r1: 0.5,
        r2: 0.4,
        k1: 1000.0,
        k2: 800.0,
        alpha12: 0.5,
        alpha21: 0.5,
    };
    assert_eq!(
        growth::competition_outcome(&params).unwrap(),
        growth::CompetitionOutcome::Coexistence
    );
}

// --- Metabolism ---

#[test]
fn test_michaelis_menten_at_km() {
    let v = metabolism::michaelis_menten(1.0, 10.0, 1.0).unwrap();
    assert!((v - 5.0).abs() < 1e-10);
}

#[test]
fn test_glycolysis_2_atp() {
    assert_eq!(metabolism::glycolysis_atp(), 2);
}

#[test]
fn test_inhibition_reduces_rate() {
    let v_mm = metabolism::michaelis_menten(1.0, 10.0, 1.0).unwrap();
    let v_comp = metabolism::competitive_inhibition(1.0, 10.0, 1.0, 1.0, 1.0).unwrap();
    let v_uncomp = metabolism::uncompetitive_inhibition(1.0, 10.0, 1.0, 1.0, 1.0).unwrap();
    let v_noncomp = metabolism::noncompetitive_inhibition(1.0, 10.0, 1.0, 1.0, 1.0).unwrap();
    assert!(v_comp < v_mm);
    assert!(v_uncomp < v_mm);
    assert!(v_noncomp < v_mm);
}

#[test]
fn test_metabolic_network_steady_state() {
    // A → B → C with equal flux: B is at steady state
    let net = metabolism::MetabolicNetwork::from_reactions(vec![
        metabolism::Reaction {
            id: "r1".into(),
            stoichiometry: vec![("A".into(), -1.0), ("B".into(), 1.0)],
            reversible: false,
        },
        metabolism::Reaction {
            id: "r2".into(),
            stoichiometry: vec![("B".into(), -1.0), ("C".into(), 1.0)],
            reversible: false,
        },
    ]);
    let sv = net.net_production(&[1.0, 1.0]).unwrap();
    let b_idx = net.metabolites.iter().position(|m| m == "B").unwrap();
    assert!((sv[b_idx] - 0.0).abs() < 1e-10);
}

// --- Genetics ---

#[test]
fn test_hardy_weinberg_p06() {
    let (p2, pq2, q2) = genetics::hardy_weinberg(0.6).unwrap();
    assert!((p2 - 0.36).abs() < 1e-10);
    assert!((pq2 - 0.48).abs() < 1e-10);
    assert!((q2 - 0.16).abs() < 1e-10);
}

#[test]
fn test_gc_content_50_percent() {
    let gc = genetics::gc_content("ATGC").unwrap();
    assert!((gc - 0.5).abs() < 1e-10);
}

#[test]
fn test_dna_to_protein_pipeline() {
    let protein = genetics::translate_orf("ATGGCTTGGTAA").unwrap();
    assert_eq!(protein, "MAW");
    let mw = genetics::protein_molecular_weight(&protein).unwrap();
    assert!(mw > 0.0);
}

#[test]
fn test_reverse_complement_involution() {
    let dna = "ATGCGATCGA";
    let rc = genetics::reverse_complement(dna).unwrap();
    let rc2 = genetics::reverse_complement(&rc).unwrap();
    assert_eq!(rc2, dna.to_ascii_uppercase());
}

// --- Epidemiology ---

#[test]
fn test_r_naught_2_5() {
    let r0 = epidemiology::r_naught(0.5, 0.2).unwrap();
    assert!((r0 - 2.5).abs() < 1e-10);
}

#[test]
fn test_herd_immunity_60_percent() {
    let h = epidemiology::herd_immunity_threshold(2.5).unwrap();
    assert!((h - 0.6).abs() < 1e-10);
}

#[test]
fn test_vaccination_suppresses_r_eff() {
    let r_eff = epidemiology::effective_r(3.0, 0.7).unwrap();
    assert!(r_eff < 1.0, "70% coverage should suppress R0=3 below 1");
}

// --- Pharmacokinetics ---

#[test]
fn test_iv_bolus_halves_at_half_life() {
    let k_e = 0.1;
    let t_half = pharmacokinetics::half_life(k_e).unwrap();
    let c0 = pharmacokinetics::iv_bolus_concentration(500.0, 50.0, k_e, 0.0).unwrap();
    let c_half = pharmacokinetics::iv_bolus_concentration(500.0, 50.0, k_e, t_half).unwrap();
    assert!((c_half - c0 / 2.0).abs() < 1e-6);
}

#[test]
fn test_oral_cmax_less_than_iv_c0() {
    let cmax = pharmacokinetics::oral_cmax(500.0, 1.0, 50.0, 1.0, 0.1).unwrap();
    let c0_iv = 500.0 / 50.0;
    assert!(cmax < c0_iv);
}

// --- Resistance ---

#[test]
fn test_fic_synergy() {
    let (fic, interaction) = resistance::fic_interaction(0.125, 1.0, 0.125, 1.0).unwrap();
    assert!(fic <= 0.5);
    assert_eq!(interaction, resistance::DrugInteraction::Synergy);
}

#[test]
fn test_kill_curve_above_mic() {
    let survival = resistance::kill_curve(2.0, 1.0, 2.0).unwrap();
    assert!(survival < 1.0);
    assert!(survival > 0.0);
}

// --- Biofilm ---

#[test]
fn test_quorum_sensing_threshold() {
    assert!(biofilm::quorum_sensing(10.0, 5.0).unwrap());
    assert!(!biofilm::quorum_sensing(3.0, 5.0).unwrap());
}

#[test]
fn test_diffusion_proportional_to_nutrient() {
    let low = biofilm::diffusion_through_matrix(1.0, 1.0, 0.5).unwrap();
    let high = biofilm::diffusion_through_matrix(10.0, 1.0, 0.5).unwrap();
    assert!(high > low);
}

// --- Taxonomy ---

#[test]
fn test_taxonomy_serde_roundtrip() {
    let profile = taxonomy::MicrobialProfile {
        domain: taxonomy::Domain::Bacteria,
        gram_stain: taxonomy::GramStain::Negative,
        shape: taxonomy::CellShape::Bacillus,
        oxygen: taxonomy::OxygenRequirement::Facultative,
    };
    let json = serde_json::to_string(&profile).unwrap();
    let back: taxonomy::MicrobialProfile = serde_json::from_str(&json).unwrap();
    assert_eq!(profile.domain, back.domain);
    assert_eq!(profile.gram_stain, back.gram_stain);
    assert_eq!(profile.shape, back.shape);
    assert_eq!(profile.oxygen, back.oxygen);
}
