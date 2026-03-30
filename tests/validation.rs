//! Validation tests against published microbiology reference values.
//!
//! Each test cites its source. These ensure the models produce outputs
//! consistent with established literature.

// --- Growth kinetics ---

/// E. coli doubling time ~20 min (0.333 hr) in rich media at 37°C.
/// Growth rate μ = ln(2) / 0.333 ≈ 2.08 /hr.
/// Source: Neidhardt, Ingraham & Schaechter, "Physiology of the Bacterial Cell" (1990).
#[test]
fn validate_ecoli_doubling_time() {
    let mu = core::f64::consts::LN_2 / (20.0 / 60.0); // 20 min in hours
    let td = jivanu::growth::doubling_time(mu).unwrap();
    assert!((td - 20.0 / 60.0).abs() < 0.001, "td = {td}");
}

/// Monod kinetics: E. coli on glucose, K_s ≈ 2-4 mg/L.
/// At S = K_s, μ = μ_max / 2.
/// Source: Monod (1949) "The Growth of Bacterial Cultures", Annual Review of Microbiology.
#[test]
fn validate_monod_half_saturation() {
    let mu_max = 2.0;
    let k_s = 3.0; // mg/L, midpoint of published range
    let mu = jivanu::growth::monod_kinetics(k_s, mu_max, k_s).unwrap();
    assert!((mu - mu_max / 2.0).abs() < 1e-10);
}

// --- Enzyme kinetics ---

/// Michaelis-Menten: at [S] = K_m, v = V_max/2.
/// Universal property of MM kinetics.
/// Source: Michaelis & Menten (1913) Biochemische Zeitschrift.
#[test]
fn validate_michaelis_menten_identity() {
    let v_max = 100.0;
    let k_m = 5.0;
    let v = jivanu::metabolism::michaelis_menten(k_m, v_max, k_m).unwrap();
    assert!((v - v_max / 2.0).abs() < 1e-10);
}

/// Competitive inhibition increases apparent K_m by factor (1 + [I]/K_i)
/// but does not change V_max.
/// Source: Lehninger "Principles of Biochemistry" (8th ed).
#[test]
fn validate_competitive_inhibition_vmax_unchanged() {
    // At very high [S], competitive inhibition is overcome → v → V_max
    let v_max = 100.0;
    let v = jivanu::metabolism::competitive_inhibition(1e6, v_max, 1.0, 10.0, 1.0).unwrap();
    assert!((v - v_max).abs() < 0.1, "v = {v}");
}

/// Noncompetitive inhibition reduces V_max but does not change K_m.
/// At [S] >> K_m: v → V_max / (1 + [I]/K_i).
/// Source: Lehninger "Principles of Biochemistry" (8th ed).
#[test]
fn validate_noncompetitive_reduces_vmax() {
    let v_max = 100.0;
    let inhibitor = 5.0;
    let k_i = 5.0;
    // Expected apparent V_max = 100 / (1 + 5/5) = 50
    let v = jivanu::metabolism::noncompetitive_inhibition(1e6, v_max, 1.0, inhibitor, k_i).unwrap();
    assert!((v - 50.0).abs() < 0.1, "v = {v}");
}

/// Total aerobic ATP yield: ~38 per glucose.
/// Source: Stryer "Biochemistry" (9th ed), Table 18.4.
#[test]
fn validate_aerobic_atp_yield() {
    assert_eq!(jivanu::metabolism::total_aerobic_atp(), 38);
    assert_eq!(jivanu::metabolism::glycolysis_atp(), 2);
    assert_eq!(jivanu::metabolism::oxidative_phosphorylation_atp(), 34);
}

// --- Epidemiology ---

/// Measles: R0 ≈ 12-18, herd immunity threshold ≈ 92-95%.
/// Source: Fine, Eames & Heymann (2011) Clin Infect Dis 52(7):911-916.
#[test]
fn validate_measles_herd_immunity() {
    let h = jivanu::epidemiology::herd_immunity_threshold(15.0).unwrap();
    assert!(h > 0.90 && h < 0.97, "measles HIT = {h}");
}

/// COVID-19 original strain: R0 ≈ 2.5, herd immunity ≈ 60%.
/// Source: Liu et al. (2020) J Travel Med 27(2).
#[test]
fn validate_covid_herd_immunity() {
    let h = jivanu::epidemiology::herd_immunity_threshold(2.5).unwrap();
    assert!((h - 0.6).abs() < 1e-10);
}

/// SIR conservation: S + I + R must sum to 1.0 over trajectory.
/// Fundamental property of closed-population compartmental models.
#[test]
fn validate_sir_population_conservation() {
    let traj =
        jivanu::epidemiology::sir_trajectory(0.999, 0.001, 0.0, 0.5, 0.1, 0.01, 1000).unwrap();
    for state in &traj {
        let total = state.s + state.i + state.r;
        assert!((total - 1.0).abs() < 0.01, "S+I+R = {total} at some step");
    }
}

/// Vaccination: R_eff < 1 when coverage exceeds herd immunity threshold.
/// Source: Anderson & May "Infectious Diseases of Humans" (1991).
#[test]
fn validate_vaccination_suppresses_transmission() {
    let r0 = 5.0;
    let vc = jivanu::epidemiology::critical_vaccination_coverage(r0).unwrap();
    let r_eff = jivanu::epidemiology::effective_r(r0, vc + 0.01).unwrap();
    assert!(r_eff < 1.0, "R_eff = {r_eff} should be < 1 above Vc");
}

// --- Pharmacokinetics ---

/// Half-life / k_e relationship: t½ = ln(2)/k_e.
/// Universal PK identity.
#[test]
fn validate_half_life_identity() {
    for k_e in [0.01, 0.1, 0.5, 1.0, 5.0] {
        let t_half = jivanu::pharmacokinetics::half_life(k_e).unwrap();
        let k_back = jivanu::pharmacokinetics::elimination_rate(t_half).unwrap();
        assert!((k_e - k_back).abs() < 1e-10, "k_e={k_e}, k_back={k_back}");
    }
}

/// Gentamicin PK: V_d ≈ 0.25 L/kg, t½ ≈ 2 hr, dose = 5 mg/kg.
/// For 70 kg patient: V_d = 17.5 L, dose = 350 mg, k_e ≈ 0.347.
/// C_0 = 350/17.5 = 20 mg/L. At t = 2 hr, C ≈ 10 mg/L.
/// Source: Bauer "Clinical Pharmacokinetics Handbook" (2008).
#[test]
fn validate_gentamicin_pk() {
    let dose = 350.0;
    let v_d = 17.5;
    let k_e = core::f64::consts::LN_2 / 2.0; // t½ = 2 hr
    let c0 = jivanu::pharmacokinetics::iv_bolus_concentration(dose, v_d, k_e, 0.0).unwrap();
    assert!((c0 - 20.0).abs() < 0.1, "C0 = {c0}");
    let c2 = jivanu::pharmacokinetics::iv_bolus_concentration(dose, v_d, k_e, 2.0).unwrap();
    assert!((c2 - 10.0).abs() < 0.1, "C(2hr) = {c2}");
}

// --- Genetics ---

/// E. coli GC content ≈ 50.8%.
/// Source: Blattner et al. (1997) Science 277(5331):1453-1462 (E. coli K-12 genome).
#[test]
fn validate_ecoli_gc_range() {
    // Test with a representative E. coli-like sequence (balanced GC)
    let gc = jivanu::genetics::gc_content("ATGCATGCATGC").unwrap();
    assert!((gc - 0.5).abs() < 0.01);
}

/// Standard genetic code: ATG = Met (start), TAA/TAG/TGA = Stop.
/// Universal across nearly all organisms.
/// Source: Crick (1968) J Mol Biol 38:367-379.
#[test]
fn validate_standard_genetic_code() {
    assert_eq!(
        jivanu::genetics::translate_codon_to_aa("ATG").unwrap(),
        jivanu::genetics::AminoAcid::Methionine
    );
    assert_eq!(
        jivanu::genetics::translate_codon_to_aa("TAA").unwrap(),
        jivanu::genetics::AminoAcid::Stop
    );
    assert_eq!(
        jivanu::genetics::translate_codon_to_aa("TAG").unwrap(),
        jivanu::genetics::AminoAcid::Stop
    );
    assert_eq!(
        jivanu::genetics::translate_codon_to_aa("TGA").unwrap(),
        jivanu::genetics::AminoAcid::Stop
    );
}

/// Glycine: lightest amino acid, MW ≈ 75.03 Da.
/// Tryptophan: heaviest, MW ≈ 204.23 Da.
/// Source: NIST Chemistry WebBook.
#[test]
fn validate_amino_acid_mw_extremes() {
    let gly = jivanu::genetics::AminoAcid::Glycine.molecular_weight();
    let trp = jivanu::genetics::AminoAcid::Tryptophan.molecular_weight();
    assert!((gly - 75.032).abs() < 0.01);
    assert!((trp - 204.228).abs() < 0.01);
    // Gly is lightest, Trp is heaviest
    assert!(gly < trp);
}

/// Kyte-Doolittle: Ile is most hydrophobic (+4.5), Arg most hydrophilic (−4.5).
/// Source: Kyte & Doolittle (1982) J Mol Biol 157:105-132.
#[test]
fn validate_kyte_doolittle_extremes() {
    assert!((jivanu::genetics::AminoAcid::Isoleucine.hydrophobicity() - 4.5).abs() < 1e-10);
    assert!((jivanu::genetics::AminoAcid::Arginine.hydrophobicity() - (-4.5)).abs() < 1e-10);
}

// --- Resistance ---

/// FIC index: synergy defined as ≤ 0.5 per EUCAST/CLSI.
/// Source: Odds (2003) J Antimicrob Chemother 52:1.
#[test]
fn validate_fic_synergy_threshold() {
    let interaction = jivanu::resistance::classify_interaction(0.5).unwrap();
    assert_eq!(interaction, jivanu::resistance::DrugInteraction::Synergy);
    let interaction = jivanu::resistance::classify_interaction(0.51).unwrap();
    assert_eq!(interaction, jivanu::resistance::DrugInteraction::Additive);
}

// --- PK/PD bridge ---

/// Aminoglycosides: target Cmax/MIC ≥ 8-10 for efficacy.
/// Source: Moore, Lietman & Smith (1987) J Infect Dis 155:93-99.
#[test]
fn validate_aminoglycoside_pk_pd_target() {
    // Gentamicin: dose=350mg, Vd=17.5L, MIC=2 → Cmax/MIC = 10
    let ratio = jivanu::bridge::cmax_mic_ratio(350.0, 17.5, 2.0).unwrap();
    assert!(ratio >= 8.0, "Cmax/MIC = {ratio}, target ≥ 8");
}
