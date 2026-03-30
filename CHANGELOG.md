# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

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

#### Bridge (new module)
- PK-driven IV time-kill simulation (`iv_time_kill`)
- PK-driven oral time-kill simulation (`oral_time_kill`)
- PK/PD indices: T>MIC, Cmax/MIC, AUC/MIC
- Multi-dose IV bolus concentration with superposition
- Steady-state trough concentration
- Regrowth time prediction
- Dose calculation for target Cmax/MIC ratio

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
