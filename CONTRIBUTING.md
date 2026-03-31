# Contributing to Jivanu

Thank you for your interest in contributing to Jivanu.

## Development Workflow

1. Fork and clone the repository
2. Create a feature branch from `main`
3. Make your changes
4. Run `make check` to validate
5. Open a pull request

## Prerequisites

- Rust stable (MSRV 1.89)
- Components: `rustfmt`, `clippy`
- Optional: `cargo-audit`, `cargo-deny`, `cargo-llvm-cov`, `cargo-tarpaulin`

## Makefile Targets

| Command | Description |
|---------|-------------|
| `make check` | fmt + clippy + test + audit |
| `make fmt` | Check formatting |
| `make clippy` | Lint with `-D warnings` |
| `make test` | Run test suite (all features + no-default) |
| `make audit` | Security audit |
| `make deny` | Supply chain checks |
| `make bench` | Run benchmarks with history tracking |
| `make coverage` | Generate HTML coverage report |
| `make coverage-check` | Verify coverage >= 70% threshold |
| `make doc` | Build documentation (strict) |

## Adding a Module

1. Create `src/module_name.rs` with module doc comment
2. Add `pub mod module_name;` to `src/lib.rs` (alphabetical order)
3. Add unit tests in the module (`#[cfg(test)] mod tests`)
4. Add serde roundtrip tests for all Serialize+Deserialize types
5. Add integration tests in `tests/integration.rs`
6. Add benchmarks in `benches/benchmarks.rs` for hot-path functions
7. Update README feature list
8. Update `docs/architecture/overview.md` module map

If the module requires an external dependency, gate it behind a feature flag
in `Cargo.toml` and use `#[cfg(feature = "...")]` on the module and its tests.

## Feature Flags

| Feature | Description |
|---------|-------------|
| `std` | Standard library support (default) |
| `hisab` | RK4 ODE integration for epidemiology, competition, PK models |
| `kimiya` | Chemistry bridges (Arrhenius, pH, temperature-adjusted kinetics) |
| `pramana` | Statistics/probability (planned: stochastic models) |
| `full` | Enables all optional features |

## Code Style

- `cargo fmt` — mandatory
- `cargo clippy -- -D warnings` — zero warnings
- Doc comments on all public items with `# Errors` section for fallible functions
- `#[non_exhaustive]` on public enums
- `#[must_use]` on all pure functions
- `#[inline]` on hot-path functions
- No `unwrap()` or `panic!()` in library code
- Every Serialize+Deserialize type must have a serde roundtrip test
- Real microbiology formulas with known-good reference values

## Testing

- Unit tests colocated in modules (`#[cfg(test)] mod tests`)
- Integration tests in `tests/integration.rs`
- Validation tests against published literature in `tests/validation.rs`
- Feature-gated tests with `#[cfg(feature = "...")]`
- Target: 80%+ line coverage (CI gate at 70%)
- Use `< 1e-10` epsilon for float comparisons

## Benchmarks

- All benchmarks use Criterion in `benches/benchmarks.rs`
- Run `make bench` to record results to `bench-history.csv` and generate `benchmarks.md`
- Performance claims in CHANGELOG must include benchmark numbers
- Never skip benchmarks before claiming improvements

## Commits

- Use conventional-style messages
- One logical change per commit

## License

By contributing, you agree that your contributions will be licensed under GPL-3.0-only (see [LICENSE](LICENSE)).
