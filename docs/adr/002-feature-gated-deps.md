# ADR 002: Feature-Gated Optional Dependencies

## Status

Accepted

## Context

Jivanu is consumed by multiple downstream crates (sangha, kimiya, kiran/joshua)
with different dependency budgets. Not all consumers need ODE integration (hisab),
chemistry bridges (kimiya), or statistics (pramana).

## Decision

All external AGNOS crate dependencies are optional and feature-gated:

| Feature | Dependency | Provides |
|---------|-----------|----------|
| `hisab` | hisab 1.x | RK4 ODE integration for epi/competition/PK |
| `kimiya` | kimiya 1.x | Arrhenius, pH, temperature-adjusted kinetics |
| `pramana` | pramana 1.x | Statistics/probability (planned) |
| `std` | serde/std, thiserror/std | Standard library support (default) |
| `full` | all of the above | Everything |

Core functionality (growth kinetics, enzyme kinetics, genetics, taxonomy)
works with zero optional dependencies.

## Consequences

- Downstream crates opt in to exactly the features they need
- Feature-gated code requires `#[cfg(feature = "...")]` on functions and tests
- CI must test all feature combinations (default, hisab, all-features)
- Bridge module functions are split between always-available and feature-gated
