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

### Remaining

- [ ] Stochastic growth models (Gillespie algorithm via pramana)
- [ ] Competition → biofilm coupling
- [ ] Benchmarks for all public functions
- [ ] kimiya bridge module

## v1.0 Criteria

- 80%+ test coverage
- All models validated against published microbiology data
- Stable public API
- Complete integration test suite
- Benchmarks for all public functions
