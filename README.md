# Jivanu

**जीवाणु** (Hindi: microbe, bacterium) — Microbiology engine for growth kinetics, metabolism, genetics, and epidemiology.

Part of the [AGNOS](https://github.com/MacCracken/agnosticos) science crate ecosystem.

## Key Capabilities

- **Growth Kinetics**: Exponential/logistic/Gompertz/Baranyi growth, Monod kinetics, chemostat,
  Lotka-Volterra competition, cardinal temperature/pH models
- **Metabolism**: Michaelis-Menten, enzyme inhibition, Hill/Emax, ATP yields, stoichiometric
  networks, exact flux balance analysis (FBA) and flux variability analysis (FVA)
- **Genetics**: Mutation rates, Hardy-Weinberg equilibrium, GC content, codon translation,
  amino acid properties
- **Epidemiology**: SIR/SEIR/SIRS models, R0, herd immunity, vaccination, case fatality rate
- **Pharmacokinetics**: IV/oral one-compartment, two-compartment, AUC, half-life
- **Stochastic**: Gillespie SSA and tau-leaping for arbitrary reaction networks
- **Biofilm**: Formation stages, quorum sensing, nutrient diffusion
- **Resistance**: Kill curves, MIC, FIC index, checkerboard, drug combinations
- **Bridges**: PK-driven time-kill, PK/PD indices (T>MIC, Cmax/MIC, AUC/MIC), biofilm effects
- **Taxonomy**: Domain classification, Gram stain, morphology, oxygen requirements

## Quick Start

```rust
use jivanu::{growth, genetics, epidemiology};

// Monod kinetics: growth rate at half-saturation
let mu = growth::monod_kinetics(1.0, 0.5, 1.0).unwrap();
assert!((mu - 0.25).abs() < 1e-10); // mu_max/2 at S=K_s

// GC content
let gc = genetics::gc_content("ATGC").unwrap();
assert!((gc - 0.5).abs() < 1e-10);

// Herd immunity for R0 = 2.5
let h = epidemiology::herd_immunity_threshold(2.5).unwrap();
assert!((h - 0.6).abs() < 1e-10); // 60%
```

## Feature Flags

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | Yes | Standard library support |
| `hisab` | No | RK4 ODE integration via hisab (forward Euler otherwise) |
| `kimiya` | No | Chemistry bridges via kimiya (Arrhenius, pH, half-life) |
| `full` | No | All features |
| `pramana`, `logging` | No | Deprecated no-ops (kept for compatibility) |

Minimum supported Rust version: **1.100**.

## License

GPL-3.0-only
