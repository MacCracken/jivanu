//! Dense two-phase simplex solver for small linear programs.
//!
//! Internal to the crate; backs flux balance and flux variability analysis.
//! Problems are stated over bounded variables,
//!
//! ```text
//! maximize c·v   subject to   A v = b,   lb ≤ v ≤ ub,   [g·v ≥ h]
//! ```
//!
//! and converted to standard form (`min f·x, M x = d, x ≥ 0`) by shifting,
//! reflecting, or splitting each variable according to which of its bounds
//! are finite. Finite upper bounds become explicit slack rows. The tableau is
//! dense, so the solver suits networks up to a few hundred reactions.

use crate::error::{JivanuError, Result};

/// Pivot / reduced-cost tolerance.
const EPS: f64 = 1e-9;

/// After this many consecutive degenerate pivots, switch from Dantzig's rule
/// to Bland's rule, which cannot cycle.
const DEGENERATE_SWITCH: usize = 50;

/// Result of a linear program.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum LpOutcome {
    /// Optimal primal solution in the original variables.
    Optimal(Vec<f64>),
    /// No point satisfies the constraints.
    Infeasible,
    /// The objective is unbounded over the feasible set.
    Unbounded,
}

/// A linear program over bounded variables with equality constraints.
pub(crate) struct BoundedLp<'a> {
    /// Equality constraint rows of `A` (each of length `n`).
    pub rows: &'a [Vec<f64>],
    /// Right-hand side `b` (one per row).
    pub rhs: &'a [f64],
    /// Lower bounds (may be `-inf`).
    pub lower: &'a [f64],
    /// Upper bounds (may be `+inf`).
    pub upper: &'a [f64],
}

/// How an original variable `v_j` maps onto standard-form columns.
#[derive(Debug, Clone, Copy)]
struct VarMap {
    /// `v_j = offset + sign * x[col] (- x[neg])`.
    offset: f64,
    col: usize,
    sign: f64,
    neg: Option<usize>,
}

impl BoundedLp<'_> {
    /// Validate dimensions and bounds.
    fn validate(&self) -> Result<usize> {
        let n = self.lower.len();
        if self.upper.len() != n {
            return Err(JivanuError::ComputationError(
                "lower and upper bounds must have the same length".into(),
            ));
        }
        if self.rhs.len() != self.rows.len() {
            return Err(JivanuError::ComputationError(
                "constraint rows and right-hand side must have the same length".into(),
            ));
        }
        if let Some(i) = self.rows.iter().position(|r| r.len() != n) {
            return Err(JivanuError::ComputationError(format!(
                "constraint row {i} has length {}, expected {n}",
                self.rows[i].len()
            )));
        }
        for (j, (&lo, &hi)) in self.lower.iter().zip(self.upper).enumerate() {
            if lo.is_nan() || hi.is_nan() || lo == f64::INFINITY || hi == f64::NEG_INFINITY {
                return Err(JivanuError::ComputationError(format!(
                    "bounds of variable {j} are invalid: [{lo}, {hi}]"
                )));
            }
            if lo > hi {
                return Err(JivanuError::ComputationError(format!(
                    "lower bound exceeds upper bound for variable {j}: {lo} > {hi}"
                )));
            }
        }
        if self
            .rows
            .iter()
            .flatten()
            .chain(self.rhs)
            .any(|v| !v.is_finite())
        {
            return Err(JivanuError::ComputationError(
                "constraint coefficients must be finite".into(),
            ));
        }
        Ok(n)
    }

    /// Maximize `c·v`, optionally subject to an extra constraint `g·v ≥ h`.
    pub(crate) fn maximize(&self, c: &[f64], extra_ge: Option<(&[f64], f64)>) -> Result<LpOutcome> {
        let n = self.validate()?;
        if c.len() != n || c.iter().any(|v| !v.is_finite()) {
            return Err(JivanuError::ComputationError(format!(
                "objective must have {n} finite coefficients"
            )));
        }
        if let Some((g, h)) = extra_ge
            && (g.len() != n || g.iter().any(|v| !v.is_finite()) || !h.is_finite())
        {
            return Err(JivanuError::ComputationError(format!(
                "extra constraint must have {n} finite coefficients and a finite bound"
            )));
        }

        // --- Map variables to non-negative standard-form columns ---
        let mut maps = Vec::with_capacity(n);
        let mut n_x = 0usize;
        let mut ub_rows: Vec<(usize, f64)> = Vec::new();
        for (&lo, &hi) in self.lower.iter().zip(self.upper) {
            let map = match (lo.is_finite(), hi.is_finite()) {
                (true, finite_hi) => {
                    if finite_hi {
                        ub_rows.push((n_x, hi - lo));
                    }
                    VarMap {
                        offset: lo,
                        col: n_x,
                        sign: 1.0,
                        neg: None,
                    }
                }
                (false, true) => VarMap {
                    offset: hi,
                    col: n_x,
                    sign: -1.0,
                    neg: None,
                },
                (false, false) => {
                    n_x += 1;
                    VarMap {
                        offset: 0.0,
                        col: n_x - 1,
                        sign: 1.0,
                        neg: Some(n_x),
                    }
                }
            };
            n_x += 1;
            maps.push(map);
        }

        let n_extra = usize::from(extra_ge.is_some());
        let n_cols = n_x + ub_rows.len() + n_extra;
        let n_rows = self.rows.len() + ub_rows.len() + n_extra;
        let mut a = vec![0.0; n_rows * n_cols];
        let mut d = vec![0.0; n_rows];

        // Express a linear form over v as coefficients over x plus a constant.
        let write_form = |row: usize, coeffs: &[f64], a: &mut [f64]| -> f64 {
            let mut constant = 0.0;
            for (map, &coef) in maps.iter().zip(coeffs) {
                if coef == 0.0 {
                    continue;
                }
                constant += coef * map.offset;
                a[row * n_cols + map.col] += coef * map.sign;
                if let Some(neg) = map.neg {
                    a[row * n_cols + neg] -= coef;
                }
            }
            constant
        };

        for (i, (row, &b)) in self.rows.iter().zip(self.rhs).enumerate() {
            let constant = write_form(i, row, &mut a);
            d[i] = b - constant;
        }
        let base = self.rows.len();
        for (k, &(col, width)) in ub_rows.iter().enumerate() {
            let r = base + k;
            a[r * n_cols + col] = 1.0;
            a[r * n_cols + n_x + k] = 1.0;
            d[r] = width;
        }
        if let Some((g, h)) = extra_ge {
            let r = n_rows - 1;
            let constant = write_form(r, g, &mut a);
            a[r * n_cols + n_cols - 1] = -1.0;
            d[r] = h - constant;
        }

        // Objective: maximize c·v  ⇔  minimize -(c mapped)·x.
        let mut cost = vec![0.0; n_cols];
        for (map, &cj) in maps.iter().zip(c) {
            cost[map.col] -= cj * map.sign;
            if let Some(neg) = map.neg {
                cost[neg] += cj;
            }
        }

        let outcome = solve_standard(a, d, &cost, n_rows, n_cols)?;
        Ok(match outcome {
            StdOutcome::Optimal(x) => {
                let v = maps
                    .iter()
                    .zip(self.lower.iter().zip(self.upper))
                    .map(|(map, (&lo, &hi))| {
                        let mut value = map.offset + map.sign * x[map.col];
                        if let Some(neg) = map.neg {
                            value -= x[neg];
                        }
                        value.clamp(lo, hi)
                    })
                    .collect();
                LpOutcome::Optimal(v)
            }
            StdOutcome::Infeasible => LpOutcome::Infeasible,
            StdOutcome::Unbounded => LpOutcome::Unbounded,
        })
    }
}

enum StdOutcome {
    Optimal(Vec<f64>),
    Infeasible,
    Unbounded,
}

enum SimplexEnd {
    Optimal,
    Unbounded,
}

/// Dense simplex tableau with an explicit basis.
struct Tableau {
    /// `m × width` row-major; last column is the right-hand side.
    t: Vec<f64>,
    /// Reduced-cost row (`width` entries; last is `-objective`).
    z: Vec<f64>,
    basis: Vec<usize>,
    m: usize,
    width: usize,
}

impl Tableau {
    #[inline]
    fn rhs(&self, r: usize) -> f64 {
        self.t[r * self.width + self.width - 1]
    }

    fn pivot(&mut self, r: usize, k: usize) {
        let w = self.width;
        let p = self.t[r * w + k];
        let inv = 1.0 / p;
        for v in &mut self.t[r * w..(r + 1) * w] {
            *v *= inv;
        }
        self.t[r * w + k] = 1.0;
        let (before, rest) = self.t.split_at_mut(r * w);
        let (pivot_row, after) = rest.split_at_mut(w);
        for row in before.chunks_exact_mut(w).chain(after.chunks_exact_mut(w)) {
            let factor = row[k];
            if factor != 0.0 {
                for (x, &pv) in row.iter_mut().zip(pivot_row.iter()) {
                    *x -= factor * pv;
                }
                row[k] = 0.0;
                if let Some(rhs) = row.last_mut()
                    && *rhs < 0.0
                    && *rhs > -EPS
                {
                    *rhs = 0.0;
                }
            }
        }
        let factor = self.z[k];
        if factor != 0.0 {
            for (x, &pv) in self.z.iter_mut().zip(pivot_row.iter()) {
                *x -= factor * pv;
            }
            self.z[k] = 0.0;
        }
        self.basis[r] = k;
    }

    /// Run simplex iterations, letting only columns `< allowed` enter.
    fn run(&mut self, allowed: usize) -> Result<SimplexEnd> {
        let w = self.width;
        let max_iter = 50 * (self.m + w).max(20);
        let mut degenerate = 0usize;
        for _ in 0..max_iter {
            let bland = degenerate >= DEGENERATE_SWITCH;
            let mut entering = None;
            let mut best = -EPS;
            for (k, &zk) in self.z[..allowed].iter().enumerate() {
                if zk < best {
                    entering = Some(k);
                    if bland {
                        break;
                    }
                    best = zk;
                }
            }
            let Some(e) = entering else {
                return Ok(SimplexEnd::Optimal);
            };

            let mut leave: Option<(usize, f64)> = None;
            for r in 0..self.m {
                let coef = self.t[r * w + e];
                if coef > EPS {
                    let ratio = self.rhs(r) / coef;
                    leave = match leave {
                        Some((lr, lratio))
                            if ratio > lratio + EPS
                                || ((ratio - lratio).abs() <= EPS
                                    && self.basis[r] >= self.basis[lr]) =>
                        {
                            Some((lr, lratio))
                        }
                        _ => Some((r, ratio)),
                    };
                }
            }
            let Some((r, ratio)) = leave else {
                return Ok(SimplexEnd::Unbounded);
            };
            if ratio <= EPS {
                degenerate += 1;
            } else {
                degenerate = 0;
            }
            self.pivot(r, e);
        }
        Err(JivanuError::SimulationFailed(
            "simplex iteration limit reached".into(),
        ))
    }
}

/// Solve `min cost·x  s.t.  a x = d, x ≥ 0` (two-phase simplex).
fn solve_standard(
    mut a: Vec<f64>,
    mut d: Vec<f64>,
    cost: &[f64],
    m: usize,
    n: usize,
) -> Result<StdOutcome> {
    // Make every right-hand side non-negative.
    for (r, dr) in d.iter_mut().enumerate() {
        if *dr < 0.0 {
            *dr = -*dr;
            for v in &mut a[r * n..(r + 1) * n] {
                *v = -*v;
            }
        }
    }

    // Reuse unit columns (slacks) as the starting basis where possible.
    let mut basis: Vec<Option<usize>> = vec![None; m];
    for k in 0..n {
        let mut unit_row = None;
        let mut is_unit = true;
        for r in 0..m {
            let v = a[r * n + k];
            if v == 0.0 {
                continue;
            }
            if v == 1.0 && unit_row.is_none() {
                unit_row = Some(r);
            } else {
                is_unit = false;
                break;
            }
        }
        if is_unit
            && let Some(r) = unit_row
            && basis[r].is_none()
        {
            basis[r] = Some(k);
        }
    }

    let artificial_rows: Vec<usize> = (0..m).filter(|&r| basis[r].is_none()).collect();
    let n_art = artificial_rows.len();
    let width = n + n_art + 1;
    let mut t = vec![0.0; m * width];
    for r in 0..m {
        t[r * width..r * width + n].copy_from_slice(&a[r * n..(r + 1) * n]);
        t[r * width + width - 1] = d[r];
    }
    let mut full_basis = vec![0usize; m];
    for (r, b) in basis.iter().enumerate() {
        if let Some(k) = b {
            full_basis[r] = *k;
        }
    }
    for (idx, &r) in artificial_rows.iter().enumerate() {
        t[r * width + n + idx] = 1.0;
        full_basis[r] = n + idx;
    }

    // Phase I: minimize the sum of artificials.
    let mut z = vec![0.0; width];
    for &r in &artificial_rows {
        for (zk, &tk) in z.iter_mut().zip(&t[r * width..(r + 1) * width]) {
            *zk -= tk;
        }
    }
    for zk in &mut z[n..n + n_art] {
        *zk = 0.0;
    }
    let mut tab = Tableau {
        t,
        z,
        basis: full_basis,
        m,
        width,
    };

    if n_art > 0 {
        if let SimplexEnd::Unbounded = tab.run(n + n_art)? {
            return Err(JivanuError::SimulationFailed(
                "phase I simplex unexpectedly unbounded".into(),
            ));
        }
        let infeasibility = -tab.z[width - 1];
        let scale = 1.0 + d.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
        if infeasibility > 1e-7 * scale {
            return Ok(StdOutcome::Infeasible);
        }
        // Drive remaining (zero-valued) artificials out of the basis.
        for r in 0..m {
            if tab.basis[r] >= n
                && let Some(k) = (0..n).find(|&k| tab.t[r * width + k].abs() > EPS)
            {
                tab.pivot(r, k);
            }
        }
    }

    // Phase II: reduced costs for the true objective.
    let mut z = vec![0.0; width];
    z[..n].copy_from_slice(cost);
    for r in 0..m {
        let b = tab.basis[r];
        let cb = if b < n { cost[b] } else { 0.0 };
        if cb != 0.0 {
            for (zk, &tk) in z.iter_mut().zip(&tab.t[r * width..(r + 1) * width]) {
                *zk -= cb * tk;
            }
        }
    }
    for zk in &mut z[n..n + n_art] {
        *zk = 0.0;
    }
    tab.z = z;

    if let SimplexEnd::Unbounded = tab.run(n)? {
        return Ok(StdOutcome::Unbounded);
    }

    let mut x = vec![0.0; n];
    for r in 0..m {
        let b = tab.basis[r];
        if b < n {
            x[b] = tab.rhs(r).max(0.0);
        }
    }
    Ok(StdOutcome::Optimal(x))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn optimal(outcome: LpOutcome) -> Vec<f64> {
        match outcome {
            LpOutcome::Optimal(v) => v,
            other => panic!("expected optimal, got {other:?}"),
        }
    }

    #[test]
    fn test_textbook_lp() {
        // max 3x + 5y  s.t. x ≤ 4, 2y ≤ 12, 3x + 2y ≤ 18  →  (2, 6), z = 36
        // Third constraint via slack s: 3x + 2y + s = 18, s ≥ 0.
        let rows = vec![vec![3.0, 2.0, 1.0]];
        let lp = BoundedLp {
            rows: &rows,
            rhs: &[18.0],
            lower: &[0.0, 0.0, 0.0],
            upper: &[4.0, 6.0, f64::INFINITY],
        };
        let v = optimal(lp.maximize(&[3.0, 5.0, 0.0], None).unwrap());
        assert!((v[0] - 2.0).abs() < 1e-9);
        assert!((v[1] - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_free_and_reflected_variables() {
        // max -|x| style: max y s.t. y - x = 0, x ∈ (-inf, 3], y free → 3
        let rows = vec![vec![-1.0, 1.0]];
        let lp = BoundedLp {
            rows: &rows,
            rhs: &[0.0],
            lower: &[f64::NEG_INFINITY, f64::NEG_INFINITY],
            upper: &[3.0, f64::INFINITY],
        };
        let v = optimal(lp.maximize(&[0.0, 1.0], None).unwrap());
        assert!((v[1] - 3.0).abs() < 1e-9);
        // Minimizing y is unbounded below.
        assert_eq!(
            lp.maximize(&[0.0, -1.0], None).unwrap(),
            LpOutcome::Unbounded
        );
    }

    #[test]
    fn test_infeasible() {
        // x + y = 10 with x, y ∈ [0, 2]
        let rows = vec![vec![1.0, 1.0]];
        let lp = BoundedLp {
            rows: &rows,
            rhs: &[10.0],
            lower: &[0.0, 0.0],
            upper: &[2.0, 2.0],
        };
        assert_eq!(
            lp.maximize(&[1.0, 0.0], None).unwrap(),
            LpOutcome::Infeasible
        );
    }

    #[test]
    fn test_extra_ge_constraint() {
        // max -x s.t. x + y = 4 (y ∈ [0, 4]), with extra x ≥ 1.5 → x = 1.5
        let rows = vec![vec![1.0, 1.0]];
        let lp = BoundedLp {
            rows: &rows,
            rhs: &[4.0],
            lower: &[0.0, 0.0],
            upper: &[10.0, 4.0],
        };
        let v = optimal(lp.maximize(&[-1.0, 0.0], Some((&[1.0, 0.0], 1.5))).unwrap());
        assert!((v[0] - 1.5).abs() < 1e-9);
    }

    #[test]
    fn test_redundant_rows() {
        // Duplicate equality rows must not break phase I.
        let rows = vec![vec![1.0, -1.0], vec![2.0, -2.0]];
        let lp = BoundedLp {
            rows: &rows,
            rhs: &[0.0, 0.0],
            lower: &[0.0, 0.0],
            upper: &[5.0, 7.0],
        };
        let v = optimal(lp.maximize(&[1.0, 1.0], None).unwrap());
        assert!((v[0] - 5.0).abs() < 1e-9);
        assert!((v[1] - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_invalid_bounds_rejected() {
        let rows: Vec<Vec<f64>> = vec![];
        let lp = BoundedLp {
            rows: &rows,
            rhs: &[],
            lower: &[1.0],
            upper: &[0.0],
        };
        assert!(lp.maximize(&[1.0], None).is_err());
        let lp = BoundedLp {
            rows: &rows,
            rhs: &[],
            lower: &[f64::NAN],
            upper: &[0.0],
        };
        assert!(lp.maximize(&[1.0], None).is_err());
    }
    /// Solve a small dense system by Gaussian elimination (test helper).
    fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
        let n = b.len();
        for col in 0..n {
            let piv = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
            if a[piv][col].abs() < 1e-12 {
                return None;
            }
            a.swap(col, piv);
            b.swap(col, piv);
            for row in 0..n {
                if row != col {
                    let f = a[row][col] / a[col][col];
                    for k in col..n {
                        a[row][k] -= f * a[col][k];
                    }
                    b[row] -= f * b[col];
                }
            }
        }
        Some((0..n).map(|i| b[i] / a[i][i]).collect())
    }

    fn rank(rows: &[Vec<f64>]) -> usize {
        let mut a = rows.to_vec();
        let mut r = 0;
        let n = a.first().map_or(0, Vec::len);
        for col in 0..n {
            let Some(piv) = (r..a.len()).find(|&i| a[i][col].abs() > 1e-12) else {
                continue;
            };
            a.swap(r, piv);
            for i in 0..a.len() {
                if i != r {
                    let f = a[i][col] / a[r][col];
                    for k in 0..n {
                        a[i][k] -= f * a[r][k];
                    }
                }
            }
            r += 1;
        }
        r
    }

    /// Brute-force optimum by enumerating every vertex of
    /// `{A v = b, lb ≤ v ≤ ub}` (all bounds finite).
    fn brute_force_max(
        rows: &[Vec<f64>],
        rhs: &[f64],
        lb: &[f64],
        ub: &[f64],
        c: &[f64],
    ) -> Option<f64> {
        let n = lb.len();
        let m = rows.len();
        let mut best: Option<f64> = None;
        for mask in 0u32..(1 << n) {
            if mask.count_ones() as usize != m {
                continue;
            }
            let basic: Vec<usize> = (0..n).filter(|j| mask & (1 << j) != 0).collect();
            let nonbasic: Vec<usize> = (0..n).filter(|j| mask & (1 << j) == 0).collect();
            for choice in 0u32..(1 << nonbasic.len()) {
                let mut v = vec![0.0; n];
                for (k, &j) in nonbasic.iter().enumerate() {
                    v[j] = if choice & (1 << k) != 0 { ub[j] } else { lb[j] };
                }
                let a: Vec<Vec<f64>> = rows
                    .iter()
                    .map(|r| basic.iter().map(|&j| r[j]).collect())
                    .collect();
                let b: Vec<f64> = rows
                    .iter()
                    .zip(rhs)
                    .map(|(r, &bi)| bi - nonbasic.iter().map(|&j| r[j] * v[j]).sum::<f64>())
                    .collect();
                let sol = if m == 0 {
                    Some(vec![])
                } else {
                    solve_dense(a, b)
                };
                let Some(sol) = sol else { continue };
                for (&j, &x) in basic.iter().zip(&sol) {
                    v[j] = x;
                }
                if v.iter()
                    .zip(lb.iter().zip(ub))
                    .all(|(&x, (&l, &u))| x >= l - 1e-7 && x <= u + 1e-7)
                {
                    let obj: f64 = v.iter().zip(c).map(|(x, y)| x * y).sum();
                    best = Some(best.map_or(obj, |b| b.max(obj)));
                }
            }
        }
        best
    }

    #[test]
    fn test_random_lps_match_vertex_enumeration() {
        let mut rng = crate::rng::Lcg::new(2026);
        let mut int = |lo: i32, hi: i32| lo + (rng.next_f64() * f64::from(hi - lo + 1)) as i32;
        let mut checked = 0;
        for _ in 0..400 {
            let n = int(2, 6) as usize;
            let m = int(0, (n as i32 - 1).min(3)) as usize;
            // Small integer coefficients make degenerate vertices common.
            let rows: Vec<Vec<f64>> = (0..m)
                .map(|_| (0..n).map(|_| f64::from(int(-2, 2))).collect())
                .collect();
            let lb: Vec<f64> = (0..n).map(|_| f64::from(int(-3, 1))).collect();
            let ub: Vec<f64> = lb.iter().map(|&l| l + f64::from(int(0, 4))).collect();
            // Pick a right-hand side reachable from a random interior point
            // half the time, so both feasible and infeasible cases occur.
            let rhs: Vec<f64> = if int(0, 1) == 0 {
                let p: Vec<f64> = lb
                    .iter()
                    .zip(&ub)
                    .map(|(&l, &u)| f64::from(int(0, 2)) * 0.5 * (u - l) + l)
                    .collect();
                rows.iter()
                    .map(|r| r.iter().zip(&p).map(|(a, x)| a * x).sum())
                    .collect()
            } else {
                (0..m).map(|_| f64::from(int(-4, 4))).collect()
            };
            let c: Vec<f64> = (0..n).map(|_| f64::from(int(-3, 3))).collect();
            // Vertex enumeration needs full row rank; redundant rows are
            // covered by `test_redundant_rows`.
            if rank(&rows) < m {
                continue;
            }

            let lp = BoundedLp {
                rows: &rows,
                rhs: &rhs,
                lower: &lb,
                upper: &ub,
            };
            let expected = brute_force_max(&rows, &rhs, &lb, &ub, &c);
            match (lp.maximize(&c, None).unwrap(), expected) {
                (LpOutcome::Optimal(v), Some(best)) => {
                    let obj: f64 = v.iter().zip(&c).map(|(x, y)| x * y).sum();
                    assert!(
                        (obj - best).abs() < 1e-6,
                        "rows={rows:?} rhs={rhs:?} lb={lb:?} ub={ub:?} c={c:?}: {obj} vs {best}"
                    );
                    for (r, &bi) in rows.iter().zip(&rhs) {
                        let lhs: f64 = r.iter().zip(&v).map(|(a, x)| a * x).sum();
                        assert!((lhs - bi).abs() < 1e-6);
                    }
                    checked += 1;
                }
                (LpOutcome::Infeasible, None) => {}
                (got, want) => panic!(
                    "mismatch rows={rows:?} rhs={rhs:?} lb={lb:?} ub={ub:?} c={c:?}: got {got:?}, brute force {want:?}"
                ),
            }
        }
        assert!(checked > 100, "only {checked} feasible cases exercised");
    }
}
