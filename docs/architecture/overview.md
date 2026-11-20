# Architecture Overview

## Module Map

```
jivanu/
├── growth.rs             — Exponential/logistic growth, Monod, chemostat, competition
├── metabolism.rs          — Michaelis-Menten, ATP, inhibition, stoichiometric networks
├── genetics.rs            — Mutation, Hardy-Weinberg, GC, codons, amino acids, sequences
├── epidemiology.rs        — SIR/SEIR/SIRS, vaccination, R0, herd immunity
├── pharmacokinetics.rs    — IV/oral PK, two-compartment, AUC, half-life
├── biofilm.rs             — Formation stages, quorum sensing, diffusion
├── resistance.rs          — Kill curves, MIC, FIC, checkerboard, combinations
├── stochastic.rs          — Gillespie SSA, tau-leaping
├── bridge.rs              — PK→resistance time-kill, PK/PD indices, multi-dose, kimiya
├── taxonomy.rs            — Domain, Gram stain, morphology, oxygen, trait queries
├── error.rs               — JivanuError enum, validation helpers
├── lp.rs        (private) — dense two-phase simplex backing FBA/FVA
└── rng.rs       (private) — seeded LCG + Poisson (Knuth / PTRS) for stochastic models
```

## Optional Dependencies

- `hisab` — RK4 ODE integration for epidemiology, competition, and PK models
- `kimiya` — chemistry bridges (Arrhenius, pH, temperature-adjusted kinetics, half-life)
- `pramana`, `logging` — deprecated no-op features (dependencies removed in 1.1.0)

## Invariants

- Library code is panic-free: no `unwrap`/`expect`/`panic!` (enforced by
  manifest lints), invalid input returns `JivanuError`.
- Seeded stochastic simulations are reproducible for a given seed.

## Consumers

- sangha: epidemiological models feed social dynamics
- kimiya: biochemistry overlap
- kiran/joshua: ecosystem simulation, disease mechanics
- medical/health applications: PK/PD, resistance modeling
