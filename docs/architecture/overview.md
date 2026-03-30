# Architecture Overview

## Module Map

```
jivanu/
├── growth.rs        — Exponential/logistic growth, Monod, doubling time
├── metabolism.rs    — Michaelis-Menten, ATP yields, inhibition
├── genetics.rs      — Mutation rates, Hardy-Weinberg, GC content, codons
├── epidemiology.rs  — SIR/SEIR models, R0, herd immunity
├── biofilm.rs       — Formation stages, quorum sensing, diffusion
├── resistance.rs    — Kill curves, MIC, resistance transfer
├── taxonomy.rs      — Domain, Gram stain, morphology, oxygen
└── error.rs         — JivanuError enum
```

## Consumers

- sangha: epidemiological models feed social dynamics
- kimiya: biochemistry overlap
- kiran/joshua: ecosystem simulation, disease mechanics
