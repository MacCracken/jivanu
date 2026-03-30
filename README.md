# Jivanu

**जीवाणु** (Hindi: microbe, bacterium) — Microbiology engine for growth kinetics, metabolism, genetics, and epidemiology.

Part of the [AGNOS](https://github.com/MacCracken/agnosticos) science crate ecosystem.

## Key Capabilities

- **Growth Kinetics**: Exponential/logistic growth, Monod kinetics, doubling time, chemostat
- **Metabolism**: Michaelis-Menten enzyme kinetics, ATP yields, enzyme inhibition
- **Genetics**: Mutation rates, Hardy-Weinberg equilibrium, GC content, codon translation
- **Epidemiology**: SIR/SEIR models, R0, herd immunity threshold, case fatality rate
- **Biofilm**: Formation stages, quorum sensing, nutrient diffusion
- **Resistance**: Kill curves, MIC, resistance transfer rates
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
| `hisab` | No | Advanced math via hisab |
| `pramana` | No | Statistics via pramana |
| `kimiya` | No | Chemistry via kimiya |
| `logging` | No | Tracing subscriber |
| `full` | No | All features |

## License

GPL-3.0-only
