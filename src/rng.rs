//! Internal pseudo-random number generation for the stochastic simulators.
//!
//! A 64-bit linear congruential generator (Knuth MMIX multiplier) whose high
//! 31 bits form a uniform variate. Seeded simulations are reproducible: the
//! same seed always yields the same sequence.

/// Knuth MMIX LCG multiplier.
const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;

/// 2^31, the resolution of [`Lcg::next_f64`].
const TWO_POW_31: f64 = 2_147_483_648.0;

/// Below this mean, Poisson variates use Knuth's multiplication method;
/// at or above it, Hörmann's PTRS transformed-rejection method.
const PTRS_THRESHOLD: f64 = 30.0;

/// Minimal seeded LCG used by the Gillespie and tau-leaping simulators.
#[derive(Debug, Clone)]
pub(crate) struct Lcg {
    state: u64,
}

impl Lcg {
    /// Create a generator from a seed (offset by one so a zero seed is usable).
    #[inline]
    #[must_use]
    pub(crate) const fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }

    /// Uniform variate in `[0, 1)` with 31 bits of resolution.
    #[inline]
    pub(crate) fn next_f64(&mut self) -> f64 {
        self.state = self.state.wrapping_mul(LCG_MULTIPLIER).wrapping_add(1);
        (self.state >> 33) as f64 / TWO_POW_31
    }

    /// Poisson variate with mean `lambda` (must be finite and positive).
    ///
    /// Exact for every `lambda`: Knuth's method below [`PTRS_THRESHOLD`],
    /// Hörmann's PTRS (1993) above it, which runs in constant expected time.
    #[inline]
    pub(crate) fn poisson(&mut self, lambda: f64) -> u64 {
        if lambda < PTRS_THRESHOLD {
            self.poisson_knuth(lambda)
        } else {
            self.poisson_ptrs(lambda)
        }
    }

    fn poisson_knuth(&mut self, lambda: f64) -> u64 {
        let l = (-lambda).exp();
        let mut k = 0u64;
        let mut p = 1.0;
        loop {
            k += 1;
            p *= self.next_f64();
            if p <= l {
                break;
            }
        }
        k.saturating_sub(1)
    }

    /// Hörmann, W. (1993). The transformed rejection method for generating
    /// Poisson random variables. Insurance: Mathematics and Economics 12:39-45.
    fn poisson_ptrs(&mut self, lambda: f64) -> u64 {
        let slam = lambda.sqrt();
        let loglam = lambda.ln();
        let b = 0.931 + 2.53 * slam;
        let a = -0.059 + 0.02483 * b;
        let ln_invalpha = (1.1239 + 1.1328 / (b - 3.4)).ln();
        let vr = 0.9277 - 3.6224 / (b - 2.0);
        loop {
            let u = self.next_f64() - 0.5;
            let v = self.next_f64();
            let us = 0.5 - u.abs();
            let k = ((2.0 * a / us + b) * u + lambda + 0.43).floor();
            if us >= 0.07 && v <= vr {
                return k as u64;
            }
            if k < 0.0 || (us < 0.013 && v > us) {
                continue;
            }
            let lhs = v.ln() + ln_invalpha - (a / (us * us) + b).ln();
            if lhs <= -lambda + k * loglam - ln_factorial(k) {
                return k as u64;
            }
        }
    }
}

/// `ln(k!)` for a non-negative integer-valued `k`.
fn ln_factorial(k: f64) -> f64 {
    if k < 20.0 {
        let mut product = 1.0;
        let mut i = 2.0;
        while i <= k {
            product *= i;
            i += 1.0;
        }
        product.ln()
    } else {
        // Stirling series for ln Γ(x) at x = k + 1; error < 1e-14 for x ≥ 21.
        let x = k + 1.0;
        let x2 = x * x;
        (x - 0.5) * x.ln() - x + 0.5 * (2.0 * core::f64::consts::PI).ln() + 1.0 / (12.0 * x)
            - 1.0 / (360.0 * x * x2)
            + 1.0 / (1260.0 * x * x2 * x2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniform_range() {
        let mut rng = Lcg::new(7);
        for _ in 0..10_000 {
            let u = rng.next_f64();
            assert!((0.0..1.0).contains(&u));
        }
    }

    #[test]
    fn test_reproducible() {
        let mut a = Lcg::new(42);
        let mut b = Lcg::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_f64().to_bits(), b.next_f64().to_bits());
        }
    }

    #[test]
    fn test_matches_original_inline_generator() {
        // The simulators previously inlined this generator; seeded results
        // must not change.
        for seed in [0u64, 1, 42, u64::MAX] {
            let mut state = seed.wrapping_add(1);
            let mut rng = Lcg::new(seed);
            for _ in 0..1000 {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                let expected = (state >> 33) as f64 / (1u64 << 31) as f64;
                assert_eq!(rng.next_f64().to_bits(), expected.to_bits());
            }
        }
    }

    #[test]
    fn test_ln_factorial_matches_direct_sum() {
        for k in 0..60u32 {
            let direct: f64 = (2..=k).map(|i| f64::from(i).ln()).sum();
            let approx = ln_factorial(f64::from(k));
            assert!(
                (direct - approx).abs() < 1e-10,
                "k={k}: {direct} vs {approx}"
            );
        }
    }

    fn sample_moments(lambda: f64, n: usize) -> (f64, f64) {
        let mut rng = Lcg::new(12345);
        let draws: Vec<f64> = (0..n).map(|_| rng.poisson(lambda) as f64).collect();
        let mean = draws.iter().sum::<f64>() / n as f64;
        let var = draws.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
        (mean, var)
    }

    #[test]
    fn test_poisson_small_lambda_moments() {
        let (mean, var) = sample_moments(4.0, 50_000);
        assert!((mean - 4.0).abs() < 0.05, "mean = {mean}");
        assert!((var - 4.0).abs() < 0.15, "var = {var}");
    }

    #[test]
    fn test_poisson_ptrs_moments() {
        let (mean, var) = sample_moments(100.0, 50_000);
        assert!((mean - 100.0).abs() < 0.3, "mean = {mean}");
        assert!((var - 100.0).abs() < 3.0, "var = {var}");
    }

    #[test]
    fn test_poisson_large_lambda_not_capped() {
        // Knuth's method alone saturates near 745 because e^-λ underflows.
        let (mean, _) = sample_moments(10_000.0, 5_000);
        assert!((mean - 10_000.0).abs() < 10.0, "mean = {mean}");
    }
}
