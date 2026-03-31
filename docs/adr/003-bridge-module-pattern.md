# ADR 003: Cross-Module Bridge Pattern

## Status

Accepted

## Context

Jivanu's domain modules (pharmacokinetics, resistance, growth, biofilm,
metabolism) model different aspects of microbiology. Real-world scenarios
require composing these models: drug concentration over time drives bacterial
killing, biofilm stage modifies antibiotic efficacy, nutrient diffusion
limits growth rate within biofilms.

Placing cross-module logic inside either source module creates coupling.

## Decision

A dedicated `bridge` module contains all cross-module composition functions.
Bridge functions take primitive values or types from multiple modules and
combine their computations without modifying either source module.

Categories:
- **PK-resistance bridges**: time-kill simulations, PK/PD indices
- **Biofilm-growth bridges**: diffusion-limited growth, stage modifiers
- **Biofilm-resistance bridges**: MIC multipliers, adjusted kill curves
- **Kimiya bridges** (feature-gated): temperature/pH chemistry to growth

Each bridge function documents which source functions it composes.

## Consequences

- Source modules remain independent and testable in isolation
- New cross-module integrations go in one place
- Feature-gated chemistry bridges coexist with always-available internal bridges
- Bridge module grows as the crate matures (post-v1.0 additions land here)
