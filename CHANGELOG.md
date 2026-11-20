# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.1.0] - 2026-11-19

Toolchain and dependency refresh plus a P(-1) audit: correctness repairs,
panic hardening, and optimizations. Most changes fix wrong results; see the
**Changed** and **Fixed** sections for differences in numeric output.

### Changed

- **MSRV raised to Rust 1.100** (from 1.89); edition 2024. Lets the resolver
  pick the latest dependency releases.
- Dependencies: `hisab` 1.4, `kimiya` 1.2, `serde` 1.0.229, `thiserror` 2.0.21;
  dev: `criterion` 0.8.2, `serde_json` 1.0.154.
- Removed unused dependencies `tracing`, `tracing-subscriber`, and `pramana`
  (no code used them). The `logging` and `pramana` features remain as
  **deprecated no-ops** so existing feature lists keep compiling.
- `metabolism::flux_balance_analysis` now solves the linear program exactly
  with a dense two-phase simplex (Dantzig pricing, Bland anti-cycling). The
  previous gradient/projection heuristic could return infeasible or
  suboptimal fluxes. Infinite bounds are supported; an unbounded objective
  returns an error. As in standard FBA, every metabolite is balanced, so
  model uptake/secretion with exchange reactions.
- `metabolism::flux_variability_analysis` now solves a min and a max LP per
  reaction subject to `c·v ≥ Z − (1 − fraction)·|Z|`. `fraction` was
  previously ignored and the ranges were not computed from the LP.
- `pharmacokinetics::oral_concentration`, `oral_tmax`, `oral_cmax` accept
  `k_a == k_e` (exact limit) and use cancellation-free forms that stay
  accurate when the rates are close and for flip-flop kinetics at long times.
- `bridge::iv_time_kill` / `oral_time_kill`: survival is exactly 1 at `t = 0`
  and accumulates by the trapezoidal rule on the log-survival rate.
  Previously a full step of killing was applied at `t = 0`.
- `metabolism::hill_equation` / `emax_model` use the overflow-safe form
  `1 / (1 + (K/S)^n)` (one `powf` instead of two; no `inf/inf = NaN`).
- `stochastic::tau_leaping` Poisson draws use Hörmann's PTRS for `λ ≥ 30`
  (exact, constant time). Knuth's method saturated near 745 for large `λ`
  and cost O(λ). Seeded trajectories with `λ < 30` are unchanged; those with
  larger `λ` differ from 1.0. Firings within a leap are applied in one batch.
- `stochastic::gillespie_ssa` uses a strict cumulative comparison, so a
  zero-propensity reaction can never fire.
- All seeded simulators share one internal RNG; its stream is pinned by a test
  against the original generator.
- `growth::n_strain_competition_step` no longer clones its parameters every
  step; `MetabolicNetwork::from_reactions` builds the matrix in one pass.
- `genetics::gc_content`, `translate_orf`, and `reverse_complement` work on
  bytes, and their errors name the offending character.

### Performance

Criterion medians, 1.0.0 code vs 1.1.0, on the same machine (`benchmarks.md`):

| Benchmark | 1.0.0 | 1.1.0 | Change |
|-----------|-------|-------|--------|
| `metabolism/fba_core_17rxn` | 3.93 ms | 5.15 µs | ~760× faster (and now exact) |
| `metabolism/fva_core_17rxn` | 3.94 ms | 195 µs | ~20× faster (and now exact) |
| `stochastic/tau_leaping_1e6_pop` | 34.1 µs | 818 ns | ~42× faster (and no longer capped) |
| `metabolism/hill_equation` | 32.5 ns | 18.6 ns | −43% |
| `growth/n_strain_step_5` | 286 ns | 195 ns | −32% |
| `growth/gillespie_birth_death_100` | 1.19 µs | 813 ns | −32% |
| `metabolism/net_production` | 25.9 ns | 19.3 ns | −25% |
| `bridge/multi_dose_concentration_30` | 60.1 ns | 47.5 ns | −21% |
| `genetics/gc_content_32bp` | 41.7 ns | 34.0 ns | −18% |
| `pharmacokinetics/oral_concentration` | 21.6 ns | 18.8 ns | −13% |
| `stochastic/gillespie_ssa_100` | 2.39 µs | 2.51 µs | +5% (per-event propensity validation) |

Other benchmarks are within ±5%.

### Fixed

- `genetics::AminoAcid::Glycine.molecular_weight()` returned the monoisotopic
  mass (75.032 Da) while every other residue uses average mass; now 75.067 Da.
  Affects `protein_molecular_weight` for glycine-containing sequences.
- `genetics::translate_orf` panicked on non-ASCII input (byte-index slicing of
  `&str`); it now returns an error.
- `metabolism::flux_balance_analysis` panicked (inside `f64::clamp`) when a
  lower bound exceeded its upper bound or a bound was NaN; it now returns an error.
- `stochastic::gillespie_ssa` / `tau_leaping` panicked when the propensity
  function returned more values than reactions. They now return an error for
  wrong-length, negative, or non-finite propensities, and for reactions
  that reference a species index outside the state (previously silently ignored).
- `epidemiology::sir_trajectory` / `seir_trajectory` with `steps = 0` errored
  with `hisab` but not without; both now return the initial state.
- `growth::gompertz` produced NaN when `n_max == n0`; it now errors when `n_max <= n0`.
- `epidemiology::case_fatality_rate` rejects `deaths > cases`.
- `pharmacokinetics::auc_trapezoidal` rejects non-finite times and
  negative/non-finite concentrations (a NaN time silently yielded NaN).
- Trajectory functions return an error instead of overflowing on absurd step counts.
- `MetabolicNetwork` methods return an error instead of panicking when the
  public `s_matrix` field has an inconsistent shape.
- Documentation: `lineweaver_burk` claimed to return slope/intercept but
  returns plot coordinates (doc corrected; see `lineweaver_burk_line`);
  `regrowth_time`, `baranyi_roberts` (sharp-lag form), `checkerboard`
  (model grid), and `herd_immunity_threshold` (`R0 ≤ 1`) clarified.

### Added

- `metabolism::lineweaver_burk_line(v_max, k_m)` → `(slope, intercept)`.
- `PartialEq` on state, parameter, and result structs (`SirState`, `SeirState`,
  `SeirParams`, `SirsState`, `SirsParams`, `TwoStrainState`, `CompetitionParams`,
  `GillespieEvent`, `TwoCompartmentState`, `TwoCompartmentParams`,
  `TimeKillPoint`, `OralTimeKillParams`, `FbaResult`, `StochasticState`);
  `PartialEq + Eq` on `StochasticReaction`; `Copy + Eq + Hash` on
  `MicrobialProfile`; `Hash` on the taxonomy enums.
- Manifest lints: `unsafe_code = "forbid"`, `missing_docs = "warn"`, and
  `clippy::{unwrap_used, expect_used, panic} = "deny"`, which enforce
  panic-free library code (allowed in tests via `clippy.toml`).
- CI also runs the test suite with default and no-default features (the
  forward-Euler paths were previously untested in CI). GitHub Actions bumped
  (`checkout` v7, `codecov-action` v7, `action-gh-release` v3).
- 11 new benchmarks (FBA/FVA, tau-leaping at large populations, ORF
  translation, N-strain competition, SIR trajectory, Tmax, time-kill,
  multi-dose, Hill, SSA).
- 59 new tests (334 → 393), including a randomized cross-check of the simplex solver
  against brute-force vertex enumeration and analytic checks (SSA mean decay,
  SIR final-size relation, time-kill integral, Bateman agreement).

## [1.0.0] - 2026-03-31

### Added

#### Pharmacokinetics (new module)
- One-compartment IV bolus model (`iv_bolus_concentration`)
- One-compartment oral model with Bateman equation (`oral_concentration`)
- Tmax, Cmax calculation for oral dosing (`oral_tmax`, `oral_cmax`)
- AUC via trapezoidal rule and analytical IV bolus formula
- Two-compartment IV bolus model with RK4 integration (`two_compartment_step`)
- Half-life / elimination rate interconversion
- Plasma concentration helper

#### Epidemiology
- SIRS model with waning immunity and vaccination (`sirs_step`)
- Effective reproduction number with vaccination (`effective_r`)
- Critical vaccination coverage calculation
- SEIR trajectory function (parity with existing SIR trajectory)
- RK4 integration (hisab feature) for SIR, SEIR, and SIRS models
- Conservation tests verifying S+I+R drift < 1e-10 with RK4

#### Genetics
- `AminoAcid` enum — 20 standard amino acids + Stop
- Amino acid properties: molecular weight, one/three-letter codes, full name
- Biochemical properties: charge class, Kyte-Doolittle hydrophobicity, isoelectric point
- `translate_codon_to_aa` — returns `AminoAcid` (zero-alloc bytes-based)
- `codon_degeneracy` — verified sum to 64
- `ChargeClass` enum for side-chain classification
- `reverse_complement` — DNA reverse complement
- `translate_orf` — translate open reading frame to protein sequence
- `protein_molecular_weight` — estimate protein MW from sequence

#### Growth
- Two-strain Lotka-Volterra competition (`competition_step`, `competition_outcome`)
- N-strain generalized competition with interaction matrix (`n_strain_competition_step`)
- `CompetitionOutcome` enum (Coexistence, Strain1Wins, Strain2Wins, Unstable)
- RK4 integration (hisab feature) for competition ODE systems

#### Resistance
- Fractional Inhibitory Concentration index (`fic_index`, `classify_interaction`)
- Drug interaction classification (Synergy/Additive/Indifferent/Antagonism)
- Checkerboard assay model with isobologram grid
- Two-drug combination kill curve (Bliss independence model)

#### Metabolism
- Stoichiometric matrix representation (`Reaction`, `MetabolicNetwork`)
- Steady-state mass balance check (`is_steady_state`)
- Net metabolite production (`net_production` — S·v computation)
- Net ATP yield from flux distribution
- Uncompetitive enzyme inhibition
- Noncompetitive enzyme inhibition
- Hill equation / sigmoidal cooperative binding
- Emax pharmacodynamic model
- Flux balance analysis with iterative projection solver
- Flux variability analysis (FVA)

#### Growth (expanded)
- Baranyi-Roberts lag-phase growth model
- Gompertz growth model
- Ratkowsky square-root temperature model
- Cardinal temperature model (Rosso)
- Cardinal pH model (Rosso)
- Herbert-Pirt maintenance energy model and apparent yield
- Tilman R* resource competition and winner prediction
- Cross-feeding / syntrophy growth model
- Multi-species chemostat dynamics with competitive exclusion

#### Stochastic (new module)
- General Gillespie SSA for arbitrary reaction networks
- Tau-leaping approximate stochastic simulation
- Poisson-distributed event firing

#### Bridge (new module)
- PK-driven IV time-kill simulation (`iv_time_kill`)
- PK-driven oral time-kill simulation (`oral_time_kill`)
- PK/PD indices: T>MIC, Cmax/MIC, AUC/MIC
- Multi-dose IV bolus concentration with superposition
- Steady-state trough concentration
- Regrowth time prediction
- Dose calculation for target Cmax/MIC ratio
- Post-antibiotic effect (PAE) duration model
- Time-kill ODE (logistic growth + Emax killing)
- Biofilm-limited growth rate (diffusion + Monod)
- Biofilm stage growth modifiers and MIC multipliers
- Biofilm-adjusted kill curve
- Kimiya bridges (feature-gated): Arrhenius→growth modifier, pH→growth modifier,
  temperature-adjusted Michaelis-Menten, chemical half-life→PK

#### Taxonomy
- Domain trait queries (`is_eukaryotic`, `is_prokaryotic`)
- Gram stain properties (`has_outer_membrane`)
- Cell shape classification (`is_elongated`)
- Oxygen tolerance queries (`tolerates_oxygen`, `tolerates_anaerobic`)
- Microbial profile heuristics (`beta_lactam_susceptible`, `can_form_endospores`, `likely_motile`)

#### Validation
- 18 literature-validated tests with cited sources (Monod, Michaelis-Menten, measles/COVID herd immunity, gentamicin PK, Kyte-Doolittle, EUCAST FIC)

### Changed

- Epidemiology ODE integration refactored to use `hisab::num::rk4` when the
  `hisab` feature is enabled (4th-order accuracy, Neumaier-compensated
  accumulation). Forward Euler retained as fallback.
- `translate_codon` reimplemented as wrapper around `translate_codon_to_aa`
  (bytes-based, zero heap allocation)

## [0.1.0] - 2026-03-29

### Added

- Initial scaffold with all domain modules
- Full test suite with known-good reference values
- Criterion benchmarks
- CI/CD workflows
