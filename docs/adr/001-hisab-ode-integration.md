# ADR 001: Hisab RK4 ODE Integration

## Status

Accepted

## Context

Jivanu's epidemiology, competition, and pharmacokinetics models require
numerical integration of ordinary differential equations (SIR, SEIR, SIRS,
Lotka-Volterra, two-compartment PK). The initial scaffold used forward Euler
(first-order accuracy), which produces measurable conservation drift at
coarse time steps.

The AGNOS ecosystem already provides `hisab`, a mathematics library with a
fourth-order Runge-Kutta (RK4) integrator featuring Neumaier-compensated
accumulation.

## Decision

Use `hisab::num::rk4` and `hisab::num::rk4_trajectory` for all ODE-based
models when the `hisab` feature is enabled. Retain forward Euler as a fallback
behind `#[cfg(not(feature = "hisab"))]` for builds without the dependency.

Each ODE function follows the pattern:
- Public function validates inputs, delegates to `*_inner`
- `#[cfg(feature = "hisab")] fn *_inner(...)` uses RK4
- `#[cfg(not(feature = "hisab"))] fn *_inner(...)` uses Euler

## Consequences

- S+I+R conservation drift drops from ~1e-3 (Euler) to ~1e-10 (RK4) at dt=0.1
- No mandatory dependency on hisab for downstream consumers
- Two code paths to maintain per ODE model
- Consistent pattern makes adding new ODE models straightforward
