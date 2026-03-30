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
├── taxonomy.rs            — Domain, Gram stain, morphology, oxygen
└── error.rs               — JivanuError enum
```

## Optional Dependencies

- `hisab` — RK4 ODE integration for epidemiology, competition, and PK models
- `pramana` — statistics/probability (planned: stochastic models)
- `kimiya` — chemistry bridge (planned)

## Consumers

- sangha: epidemiological models feed social dynamics
- kimiya: biochemistry overlap
- kiran/joshua: ecosystem simulation, disease mechanics
- medical/health applications: PK/PD, resistance modeling
