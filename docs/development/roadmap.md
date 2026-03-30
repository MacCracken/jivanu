# Roadmap

## v0.1.0 (scaffold)

- [x] Growth kinetics (exponential, logistic, Monod, chemostat)
- [x] Metabolism (Michaelis-Menten, ATP, competitive inhibition)
- [x] Genetics (mutation rate, Hardy-Weinberg, GC content, codons)
- [x] Epidemiology (SIR, SEIR, R0, herd immunity)
- [x] Biofilm (stages, quorum sensing, diffusion)
- [x] Resistance (kill curves, MIC, transfer)
- [x] Taxonomy (classification enums)

## Unreleased

### Done

- [x] RK4 integration via hisab for all ODE-based models
- [x] Full codon table with amino acid properties (MW, charge, hydrophobicity, pI)
- [x] Multi-strain competition (2-strain + N-strain Lotka-Volterra)
- [x] Antibiotic combination effects (FIC, checkerboard, Bliss independence)
- [x] Metabolic flux analysis basics (stoichiometric matrix, steady-state, net ATP)
- [x] Pharmacokinetics module (IV bolus, oral PK, two-compartment, AUC)
- [x] SIRS model with waning immunity and vaccination
- [x] Sequence utilities (reverse complement, ORF translation, protein MW)
- [x] Uncompetitive and noncompetitive enzyme inhibition
- [x] SEIR trajectory function

- [x] Cross-module bridge: PK → resistance time-kill simulations, PK/PD indices
- [x] Taxonomy functions (domain/gram/shape/oxygen queries, profile heuristics)
- [x] Model validation test suite (18 tests with literature citations)
- [x] Integration tests expanded to cover all 10 modules
- [x] Gillespie birth-death stochastic simulation
- [x] Competition → biofilm coupling (nutrient limitation, stage modifiers, biofilm MIC)
- [x] 16 benchmarks covering all modules

### Remaining for v1.0

#### Critical
- [ ] Hill equation / Emax PD model (standard sigmoidal dose-response)
- [ ] Baranyi-Roberts lag-phase growth model (predictive food microbiology standard)
- [ ] Temperature/pH growth rate modifiers (cardinal models)

#### High
- [ ] General Gillespie SSA for arbitrary reaction networks
- [ ] Tau-leaping approximate stochastic simulation
- [ ] Tilman R* resource competition model
- [ ] Cross-feeding / syntrophy models
- [ ] Gompertz growth model
- [ ] Post-antibiotic effect (PAE) in bridge module
- [ ] Multi-species chemostat dynamics
- [ ] Herbert-Pirt maintenance energy model

#### Medium
- [ ] Flux bounds + objective function on MetabolicNetwork
- [ ] LP-based flux balance analysis solver (pure Rust simplex)
- [ ] Flux variability analysis (FVA)
- [ ] Time-kill ODE model (dN/dt with growth + Emax killing)
- [ ] kimiya bridge module

## v1.0 Criteria

- 80%+ test coverage
- All models validated against published microbiology data
- Stable public API (#[non_exhaustive] on all public structs — design decision pending)
- Complete integration test suite
- Benchmarks for all modules

## Post-v1.0

### v1.1–v1.2: High Value Extensions
- Genome-scale FBA (1000+ reactions, sparse matrices)
- Dynamic FBA (dFBA) — coupled FBA + ODE
- Gene-reaction rules (GPR) and gene essentiality screens
- SBML Level 3 import/export (feature-gated)
- Diauxic growth / substrate switching

### v1.3+: Advanced / Specialized
- Spatial reaction-diffusion (lattice-based, for biofilm modeling)
- Agent-based individual cell modeling
- Gene regulatory networks (Hill-based transcription)
- Population genetics (Wright-Fisher, Moran, coalescent)
- HGT dynamics (conjugation/transduction kinetics coupled to growth)
- Metabolic control analysis (MCA)
- Community FBA (SteadyCom)
- Antimicrobial resistance evolution (within-host PK/PD → between-host epi)
- Chemical Langevin equation
- Whole-cell modeling components (macromolecular budgets, growth laws)
