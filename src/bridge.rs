//! Cross-module bridges — composing models from different domains.
//!
//! Functions that connect pharmacokinetics, resistance, growth, and other
//! modules into integrated simulations.

use serde::{Deserialize, Serialize};

use crate::error::{Result, validate_non_negative, validate_positive};
use crate::{biofilm, growth, pharmacokinetics, resistance};

/// A single time point in a PK-driven time-kill simulation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TimeKillPoint {
    /// Time after drug administration (hr).
    pub time: f64,
    /// Plasma drug concentration (mg/L).
    pub concentration: f64,
    /// Surviving fraction of bacteria (0–1).
    pub survival: f64,
}

/// Simulate bacterial survival over time given IV bolus drug administration.
///
/// At each time step, computes the plasma concentration from the
/// one-compartment IV bolus model, then applies the kill curve to
/// determine bacterial survival.
///
/// This bridges [`pharmacokinetics::iv_bolus_concentration`] with
/// [`resistance::kill_curve`] to produce a full time-kill trajectory.
///
/// # Arguments
///
/// - `dose` — administered dose (mg)
/// - `v_d` — volume of distribution (L)
/// - `k_e` — elimination rate constant (1/hr)
/// - `mic` — minimum inhibitory concentration (mg/L)
/// - `kill_rate` — bacterial kill rate constant
/// - `dt` — time step (hr)
/// - `steps` — number of time steps
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the time-kill trajectory without side effects"]
pub fn iv_time_kill(
    dose: f64,
    v_d: f64,
    k_e: f64,
    mic: f64,
    kill_rate: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<TimeKillPoint>> {
    validate_positive(dose, "dose")?;
    validate_positive(v_d, "v_d")?;
    validate_positive(k_e, "k_e")?;
    validate_positive(mic, "mic")?;
    validate_positive(kill_rate, "kill_rate")?;
    validate_positive(dt, "dt")?;

    let mut trajectory = Vec::with_capacity(steps + 1);
    let mut cumulative_survival = 1.0;

    for step in 0..=steps {
        let t = dt * step as f64;
        let conc = pharmacokinetics::iv_bolus_concentration(dose, v_d, k_e, t)?;
        // Instantaneous survival fraction at this concentration
        let instant_survival = resistance::kill_curve(conc, mic, kill_rate)?;
        // Cumulative: survival decreases over time, never recovers
        cumulative_survival *= instant_survival.powf(dt);

        trajectory.push(TimeKillPoint {
            time: t,
            concentration: conc,
            survival: cumulative_survival,
        });
    }

    Ok(trajectory)
}

/// Parameters for an oral PK time-kill simulation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct OralTimeKillParams {
    /// Administered dose (mg).
    pub dose: f64,
    /// Fraction absorbed (0–1).
    pub bioavailability: f64,
    /// Volume of distribution (L).
    pub v_d: f64,
    /// Absorption rate constant (1/hr).
    pub k_a: f64,
    /// Elimination rate constant (1/hr).
    pub k_e: f64,
    /// Minimum inhibitory concentration (mg/L).
    pub mic: f64,
    /// Bacterial kill rate constant.
    pub kill_rate: f64,
}

/// Simulate bacterial survival over time given oral drug administration.
///
/// Bridges [`pharmacokinetics::oral_concentration`] with
/// [`resistance::kill_curve`].
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the time-kill trajectory without side effects"]
pub fn oral_time_kill(
    params: &OralTimeKillParams,
    dt: f64,
    steps: usize,
) -> Result<Vec<TimeKillPoint>> {
    validate_positive(params.dose, "dose")?;
    validate_positive(params.bioavailability, "bioavailability")?;
    validate_positive(params.v_d, "v_d")?;
    validate_positive(params.k_a, "k_a")?;
    validate_positive(params.k_e, "k_e")?;
    validate_positive(params.mic, "mic")?;
    validate_positive(params.kill_rate, "kill_rate")?;
    validate_positive(dt, "dt")?;

    let mut trajectory = Vec::with_capacity(steps + 1);
    let mut cumulative_survival = 1.0;

    for step in 0..=steps {
        let t = dt * step as f64;
        let conc = pharmacokinetics::oral_concentration(
            params.dose,
            params.bioavailability,
            params.v_d,
            params.k_a,
            params.k_e,
            t,
        )?;
        let instant_survival = resistance::kill_curve(conc, params.mic, params.kill_rate)?;
        cumulative_survival *= instant_survival.powf(dt);

        trajectory.push(TimeKillPoint {
            time: t,
            concentration: conc,
            survival: cumulative_survival,
        });
    }

    Ok(trajectory)
}

/// Time above MIC for an IV bolus.
///
/// `T>MIC = ln(C_0 / MIC) / k_e`
///
/// The duration for which drug concentration exceeds the MIC — a key
/// PK/PD index for time-dependent antibiotics (beta-lactams, macrolides).
///
/// # Errors
///
/// Returns error if parameters are invalid or C_0 ≤ MIC.
#[inline]
#[must_use = "returns the time above MIC without side effects"]
pub fn time_above_mic(dose: f64, v_d: f64, k_e: f64, mic: f64) -> Result<f64> {
    validate_positive(dose, "dose")?;
    validate_positive(v_d, "v_d")?;
    validate_positive(k_e, "k_e")?;
    validate_positive(mic, "mic")?;

    let c0 = dose / v_d;
    if c0 <= mic {
        return Ok(0.0); // Never exceeds MIC
    }
    Ok((c0 / mic).ln() / k_e)
}

/// Peak concentration to MIC ratio (Cmax/MIC).
///
/// Key PK/PD index for concentration-dependent antibiotics
/// (aminoglycosides, fluoroquinolones). Higher ratios predict
/// better bactericidal activity.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the Cmax/MIC ratio without side effects"]
pub fn cmax_mic_ratio(dose: f64, v_d: f64, mic: f64) -> Result<f64> {
    validate_positive(dose, "dose")?;
    validate_positive(v_d, "v_d")?;
    validate_positive(mic, "mic")?;
    Ok(dose / (v_d * mic))
}

/// AUC/MIC ratio for IV bolus.
///
/// Key PK/PD index that integrates both time and concentration exposure.
/// `AUC/MIC = dose / (V_d × k_e × MIC)`
///
/// Target values vary by antibiotic: fluoroquinolones typically require
/// AUC/MIC > 125 for Gram-negative infections.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the AUC/MIC ratio without side effects"]
pub fn auc_mic_ratio(dose: f64, v_d: f64, k_e: f64, mic: f64) -> Result<f64> {
    validate_positive(mic, "mic")?;
    let auc = pharmacokinetics::auc_iv_bolus(dose, v_d, k_e)?;
    Ok(auc / mic)
}

/// Determine whether regrowth occurs after IV bolus dosing.
///
/// Combines growth (Monod kinetics) with PK-driven killing. Returns the
/// time at which the drug concentration drops below MIC, after which
/// surviving bacteria can resume exponential growth.
///
/// Returns `None` if C_0 never exceeds MIC (no killing phase).
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the regrowth time without side effects"]
pub fn regrowth_time(dose: f64, v_d: f64, k_e: f64, mic: f64) -> Result<Option<f64>> {
    validate_positive(dose, "dose")?;
    validate_positive(v_d, "v_d")?;
    validate_positive(k_e, "k_e")?;
    validate_positive(mic, "mic")?;

    let c0 = dose / v_d;
    if c0 <= mic {
        return Ok(None);
    }
    Ok(Some((c0 / mic).ln() / k_e))
}

/// Minimum dose to achieve a target Cmax/MIC ratio (IV bolus).
///
/// `dose = target_ratio × V_d × MIC`
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the minimum dose without side effects"]
pub fn dose_for_cmax_mic(target_ratio: f64, v_d: f64, mic: f64) -> Result<f64> {
    validate_positive(target_ratio, "target_ratio")?;
    validate_positive(v_d, "v_d")?;
    validate_positive(mic, "mic")?;
    Ok(target_ratio * v_d * mic)
}

/// Multi-dose IV bolus concentration at time `t`.
///
/// Superposition of `n_doses` given at regular intervals of `tau` hours.
///
/// `C(t) = Σ_{k=0}^{n-1} (dose/V_d) × e^(-k_e × (t - k×τ))` for `t ≥ k×τ`
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the concentration without side effects"]
pub fn multi_dose_concentration(
    dose: f64,
    v_d: f64,
    k_e: f64,
    tau: f64,
    n_doses: usize,
    t: f64,
) -> Result<f64> {
    validate_positive(dose, "dose")?;
    validate_positive(v_d, "v_d")?;
    validate_positive(k_e, "k_e")?;
    validate_positive(tau, "tau")?;
    validate_non_negative(t, "t")?;

    let c0 = dose / v_d;
    let mut total = 0.0;
    for k in 0..n_doses {
        let t_dose = tau * k as f64;
        if t >= t_dose {
            total += c0 * (-k_e * (t - t_dose)).exp();
        }
    }
    Ok(total)
}

/// Steady-state trough concentration for repeated IV bolus dosing.
///
/// `C_trough = (dose/V_d) × e^(-k_e×τ) / (1 - e^(-k_e×τ))`
///
/// The minimum concentration just before the next dose at steady state.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the trough concentration without side effects"]
pub fn steady_state_trough(dose: f64, v_d: f64, k_e: f64, tau: f64) -> Result<f64> {
    validate_positive(dose, "dose")?;
    validate_positive(v_d, "v_d")?;
    validate_positive(k_e, "k_e")?;
    validate_positive(tau, "tau")?;
    let c0 = dose / v_d;
    let decay = (-k_e * tau).exp();
    Ok(c0 * decay / (1.0 - decay))
}

/// Growth rate modifier from biofilm nutrient limitation.
///
/// Combines Monod kinetics with biofilm diffusion to compute an effective
/// growth rate. Nutrient availability inside the biofilm is reduced by
/// diffusion limitation through the matrix.
///
/// `μ_eff = μ_max × S_eff / (K_s + S_eff)`
///
/// where `S_eff = diffusivity × S_bulk / thickness` (flux-based approximation).
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the effective growth rate without side effects"]
pub fn biofilm_limited_growth(
    mu_max: f64,
    k_s: f64,
    substrate_bulk: f64,
    biofilm_thickness: f64,
    diffusivity: f64,
) -> Result<f64> {
    validate_positive(mu_max, "mu_max")?;
    validate_positive(k_s, "k_s")?;
    validate_non_negative(substrate_bulk, "substrate_bulk")?;
    validate_positive(biofilm_thickness, "biofilm_thickness")?;
    validate_positive(diffusivity, "diffusivity")?;

    let flux = biofilm::diffusion_through_matrix(substrate_bulk, biofilm_thickness, diffusivity)?;
    growth::monod_kinetics(flux, mu_max, k_s)
}

/// Growth rate multiplier based on biofilm maturation stage.
///
/// Biofilm stage affects growth rate:
/// - Attachment: reduced (establishing, 0.5×)
/// - Microcolony: near-normal (0.9×)
/// - Maturation: reduced by diffusion limitation (0.6×)
/// - Dispersal: maximal planktonic growth (1.0×)
///
/// Based on general biofilm physiology (Stewart & Franklin, 2008).
#[inline]
#[must_use]
pub const fn biofilm_stage_growth_modifier(stage: biofilm::BiofilmStage) -> f64 {
    match stage {
        biofilm::BiofilmStage::Attachment => 0.5,
        biofilm::BiofilmStage::Microcolony => 0.9,
        biofilm::BiofilmStage::Maturation => 0.6,
        biofilm::BiofilmStage::Dispersal => 1.0,
    }
}

/// Antibiotic efficacy modifier within a biofilm.
///
/// Biofilms reduce antibiotic penetration, increasing the effective MIC.
/// The MIC multiplier depends on maturation stage:
/// - Attachment: 2× (early protection)
/// - Microcolony: 10× (developing matrix)
/// - Maturation: 100–1000× (full EPS matrix)
/// - Dispersal: 1× (planktonic, full susceptibility)
///
/// Based on Mah & O'Toole (2001) Trends in Microbiology.
#[inline]
#[must_use]
pub const fn biofilm_mic_multiplier(stage: biofilm::BiofilmStage) -> f64 {
    match stage {
        biofilm::BiofilmStage::Attachment => 2.0,
        biofilm::BiofilmStage::Microcolony => 10.0,
        biofilm::BiofilmStage::Maturation => 100.0,
        biofilm::BiofilmStage::Dispersal => 1.0,
    }
}

/// Kill curve adjusted for biofilm protection.
///
/// Applies the biofilm MIC multiplier to compute survival within a
/// biofilm at a given maturation stage.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the biofilm-adjusted survival without side effects"]
pub fn biofilm_kill_curve(
    concentration: f64,
    planktonic_mic: f64,
    kill_rate: f64,
    stage: biofilm::BiofilmStage,
) -> Result<f64> {
    let effective_mic = planktonic_mic * biofilm_mic_multiplier(stage);
    resistance::kill_curve(concentration, effective_mic, kill_rate)
}

/// Post-antibiotic effect (PAE): delayed regrowth time after drug removal.
///
/// After antibiotic exposure ceases, surviving bacteria remain growth-
/// suppressed for a duration proportional to the peak exposure.
///
/// `PAE = C × ln(Cmax / MIC)` hours, where C is a drug-class constant.
///
/// Typical C values: aminoglycosides ~1.5, fluoroquinolones ~1.0,
/// beta-lactams ~0.5.
///
/// Returns the PAE duration in hours. Returns 0 if Cmax ≤ MIC.
///
/// Reference: Craig (1993) Clin Infect Dis 17:S235-S243.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[inline]
#[must_use = "returns the PAE duration without side effects"]
pub fn post_antibiotic_effect(cmax: f64, mic: f64, pae_constant: f64) -> Result<f64> {
    validate_non_negative(cmax, "cmax")?;
    validate_positive(mic, "mic")?;
    validate_non_negative(pae_constant, "pae_constant")?;
    if cmax <= mic {
        return Ok(0.0);
    }
    Ok(pae_constant * (cmax / mic).ln())
}

/// Time-kill ODE: bacterial population dynamics under antibiotic exposure.
///
/// `dN/dt = k_growth × N × (1 - N/K) - k_kill × Emax(C) × N`
///
/// Combines logistic growth with Emax-driven killing. Uses the Hill/Emax
/// model for concentration-dependent bactericidal effect.
///
/// Returns the population after one time step.
///
/// # Errors
///
/// Returns error if parameters are invalid.
#[must_use = "returns the population without side effects"]
#[allow(clippy::too_many_arguments)]
pub fn time_kill_ode_step(
    population: f64,
    capacity: f64,
    growth_rate: f64,
    concentration: f64,
    e_max: f64,
    ec50: f64,
    hill_n: f64,
    dt: f64,
) -> Result<f64> {
    validate_non_negative(population, "population")?;
    validate_positive(capacity, "capacity")?;
    validate_non_negative(growth_rate, "growth_rate")?;
    validate_non_negative(concentration, "concentration")?;
    validate_non_negative(e_max, "e_max")?;
    validate_positive(ec50, "ec50")?;
    validate_positive(hill_n, "hill_n")?;
    validate_positive(dt, "dt")?;

    let kill_effect = crate::metabolism::emax_model(concentration, e_max, ec50, hill_n)?;
    let growth = growth_rate * population * (1.0 - population / capacity);
    let killing = kill_effect * population;
    let dn = growth - killing;
    Ok((population + dn * dt).max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iv_time_kill_trajectory() {
        let traj = iv_time_kill(500.0, 50.0, 0.1, 5.0, 1.0, 0.5, 10).unwrap();
        assert_eq!(traj.len(), 11);
        // First point: t=0, conc = 10 mg/L (above MIC=5)
        assert!((traj[0].concentration - 10.0).abs() < 1e-6);
        assert!((traj[0].time - 0.0).abs() < 1e-10);
        // Survival should decrease over time while conc > MIC
        assert!(traj[5].survival < 1.0);
    }

    #[test]
    fn test_iv_time_kill_below_mic() {
        // Dose so low that C0 < MIC → no killing
        let traj = iv_time_kill(100.0, 50.0, 0.1, 5.0, 1.0, 0.5, 5).unwrap();
        // C0 = 2 < MIC=5, survival stays 1.0
        for pt in &traj {
            assert!((pt.survival - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_oral_time_kill_trajectory() {
        let params = OralTimeKillParams {
            dose: 500.0,
            bioavailability: 0.8,
            v_d: 50.0,
            k_a: 1.0,
            k_e: 0.1,
            mic: 5.0,
            kill_rate: 1.0,
        };
        let traj = oral_time_kill(&params, 0.5, 20).unwrap();
        assert_eq!(traj.len(), 21);
        // At t=0, oral concentration is 0
        assert!((traj[0].concentration - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_time_above_mic() {
        // C0 = 500/50 = 10, MIC = 5, k_e = 0.1
        // T>MIC = ln(10/5) / 0.1 = ln(2) / 0.1 ≈ 6.93 hr
        let t = time_above_mic(500.0, 50.0, 0.1, 5.0).unwrap();
        assert!((t - core::f64::consts::LN_2 / 0.1).abs() < 1e-6);
    }

    #[test]
    fn test_time_above_mic_below_threshold() {
        // C0 = 2 < MIC = 5
        let t = time_above_mic(100.0, 50.0, 0.1, 5.0).unwrap();
        assert!((t - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_cmax_mic_ratio() {
        // C0 = 500/50 = 10, MIC = 2 → ratio = 5
        let ratio = cmax_mic_ratio(500.0, 50.0, 2.0).unwrap();
        assert!((ratio - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_auc_mic_ratio() {
        // AUC = 500/(50*0.1) = 100, MIC = 2 → ratio = 50
        let ratio = auc_mic_ratio(500.0, 50.0, 0.1, 2.0).unwrap();
        assert!((ratio - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_regrowth_time() {
        let t = regrowth_time(500.0, 50.0, 0.1, 5.0).unwrap();
        assert!(t.is_some());
        // Same as time_above_mic
        let t_above = time_above_mic(500.0, 50.0, 0.1, 5.0).unwrap();
        assert!((t.unwrap() - t_above).abs() < 1e-10);
    }

    #[test]
    fn test_regrowth_time_never_exceeds() {
        let t = regrowth_time(100.0, 50.0, 0.1, 5.0).unwrap();
        assert!(t.is_none());
    }

    #[test]
    fn test_dose_for_cmax_mic() {
        // Want Cmax/MIC = 10, Vd = 50, MIC = 2 → dose = 1000
        let dose = dose_for_cmax_mic(10.0, 50.0, 2.0).unwrap();
        assert!((dose - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_multi_dose_first_dose() {
        // Single dose at t=0, check at t=0
        let c = multi_dose_concentration(500.0, 50.0, 0.1, 8.0, 1, 0.0).unwrap();
        assert!((c - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_multi_dose_accumulation() {
        // Second dose adds to remaining first dose
        let c1 = multi_dose_concentration(500.0, 50.0, 0.1, 8.0, 1, 8.0).unwrap();
        let c2 = multi_dose_concentration(500.0, 50.0, 0.1, 8.0, 2, 8.0).unwrap();
        // c2 should be c1 (residual from first) + C0 (new dose)
        assert!(c2 > c1);
        assert!((c2 - c1 - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_multi_dose_before_dose() {
        // At t=4, only first dose given (tau=8), second not yet
        let c = multi_dose_concentration(500.0, 50.0, 0.1, 8.0, 2, 4.0).unwrap();
        let c_single = multi_dose_concentration(500.0, 50.0, 0.1, 8.0, 1, 4.0).unwrap();
        assert!((c - c_single).abs() < 1e-10);
    }

    #[test]
    fn test_steady_state_trough() {
        let c = steady_state_trough(500.0, 50.0, 0.1, 8.0).unwrap();
        assert!(c > 0.0);
        // Trough should be less than C0
        assert!(c < 10.0);
    }

    #[test]
    fn test_steady_state_trough_short_interval() {
        // Shorter interval → higher trough (less time to eliminate)
        let c_short = steady_state_trough(500.0, 50.0, 0.1, 4.0).unwrap();
        let c_long = steady_state_trough(500.0, 50.0, 0.1, 12.0).unwrap();
        assert!(c_short > c_long);
    }

    // --- Biofilm-competition bridge tests ---

    #[test]
    fn test_biofilm_limited_growth() {
        // Thick biofilm reduces effective substrate → lower growth rate
        let thin = biofilm_limited_growth(1.0, 1.0, 10.0, 0.1, 0.5).unwrap();
        let thick = biofilm_limited_growth(1.0, 1.0, 10.0, 1.0, 0.5).unwrap();
        assert!(thick < thin, "thicker biofilm should limit growth more");
    }

    #[test]
    fn test_biofilm_stage_growth_modifier() {
        assert!(
            (biofilm_stage_growth_modifier(biofilm::BiofilmStage::Dispersal) - 1.0).abs() < 1e-10
        );
        assert!(biofilm_stage_growth_modifier(biofilm::BiofilmStage::Maturation) < 1.0);
        assert!(
            biofilm_stage_growth_modifier(biofilm::BiofilmStage::Attachment)
                < biofilm_stage_growth_modifier(biofilm::BiofilmStage::Microcolony)
        );
    }

    #[test]
    fn test_biofilm_mic_multiplier() {
        assert!((biofilm_mic_multiplier(biofilm::BiofilmStage::Dispersal) - 1.0).abs() < 1e-10);
        assert!(
            biofilm_mic_multiplier(biofilm::BiofilmStage::Maturation)
                > biofilm_mic_multiplier(biofilm::BiofilmStage::Microcolony)
        );
    }

    #[test]
    fn test_biofilm_kill_curve_protection() {
        // Same concentration: mature biofilm survives better than planktonic
        let planktonic = resistance::kill_curve(5.0, 1.0, 1.0).unwrap();
        let biofilm_surv =
            biofilm_kill_curve(5.0, 1.0, 1.0, biofilm::BiofilmStage::Maturation).unwrap();
        assert!(biofilm_surv > planktonic, "biofilm should protect bacteria");
    }

    #[test]
    fn test_biofilm_kill_curve_dispersal_equals_planktonic() {
        let planktonic = resistance::kill_curve(5.0, 1.0, 1.0).unwrap();
        let dispersal =
            biofilm_kill_curve(5.0, 1.0, 1.0, biofilm::BiofilmStage::Dispersal).unwrap();
        assert!((planktonic - dispersal).abs() < 1e-10);
    }

    // --- PAE + time-kill ODE tests ---

    #[test]
    fn test_pae_above_mic() {
        let pae = post_antibiotic_effect(10.0, 2.0, 1.0).unwrap();
        // PAE = 1.0 * ln(10/2) = ln(5) ≈ 1.609
        assert!((pae - 5.0_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_pae_below_mic() {
        let pae = post_antibiotic_effect(1.0, 2.0, 1.0).unwrap();
        assert!((pae - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_time_kill_ode_growth_no_drug() {
        // No drug → pure logistic growth
        let n = time_kill_ode_step(100.0, 1000.0, 0.5, 0.0, 1.0, 5.0, 1.0, 0.1).unwrap();
        assert!(n > 100.0, "should grow without drug");
    }

    #[test]
    fn test_time_kill_ode_high_drug() {
        // High drug concentration → population declines
        let n = time_kill_ode_step(1000.0, 10000.0, 0.5, 100.0, 5.0, 5.0, 2.0, 0.1).unwrap();
        assert!(n < 1000.0, "should decline under heavy drug");
    }

    #[test]
    fn test_time_kill_ode_no_growth() {
        // No growth, no drug → no change
        let n = time_kill_ode_step(100.0, 1000.0, 0.0, 0.0, 1.0, 5.0, 1.0, 0.1).unwrap();
        assert!((n - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_oral_time_kill_params_serde_roundtrip() {
        let params = OralTimeKillParams {
            dose: 500.0,
            bioavailability: 0.8,
            v_d: 50.0,
            k_a: 1.0,
            k_e: 0.1,
            mic: 5.0,
            kill_rate: 1.0,
        };
        let json = serde_json::to_string(&params).unwrap();
        let back: OralTimeKillParams = serde_json::from_str(&json).unwrap();
        assert!((params.dose - back.dose).abs() < 1e-10);
        assert!((params.bioavailability - back.bioavailability).abs() < 1e-10);
    }

    #[test]
    fn test_time_kill_point_serde_roundtrip() {
        let pt = TimeKillPoint {
            time: 1.0,
            concentration: 5.0,
            survival: 0.8,
        };
        let json = serde_json::to_string(&pt).unwrap();
        let back: TimeKillPoint = serde_json::from_str(&json).unwrap();
        assert!((pt.survival - back.survival).abs() < 1e-10);
    }
}
