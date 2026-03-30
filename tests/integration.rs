//! Integration tests for jivanu.

use jivanu::epidemiology;
use jivanu::genetics;
use jivanu::growth;
use jivanu::metabolism;

#[test]
fn test_monod_at_high_substrate() {
    let mu = growth::monod_kinetics(1000.0, 0.5, 1.0).unwrap();
    assert!((mu - 0.5).abs() < 0.001);
}

#[test]
fn test_doubling_time_1hr() {
    let td = growth::doubling_time(core::f64::consts::LN_2).unwrap();
    assert!((td - 1.0).abs() < 1e-10);
}

#[test]
fn test_michaelis_menten_at_km() {
    let v = metabolism::michaelis_menten(1.0, 10.0, 1.0).unwrap();
    assert!((v - 5.0).abs() < 1e-10);
}

#[test]
fn test_glycolysis_2_atp() {
    assert_eq!(metabolism::glycolysis_atp(), 2);
}

#[test]
fn test_hardy_weinberg_p06() {
    let (p2, pq2, q2) = genetics::hardy_weinberg(0.6).unwrap();
    assert!((p2 - 0.36).abs() < 1e-10);
    assert!((pq2 - 0.48).abs() < 1e-10);
    assert!((q2 - 0.16).abs() < 1e-10);
}

#[test]
fn test_gc_content_50_percent() {
    let gc = genetics::gc_content("ATGC").unwrap();
    assert!((gc - 0.5).abs() < 1e-10);
}

#[test]
fn test_r_naught_2_5() {
    let r0 = epidemiology::r_naught(0.5, 0.2).unwrap();
    assert!((r0 - 2.5).abs() < 1e-10);
}

#[test]
fn test_herd_immunity_60_percent() {
    let h = epidemiology::herd_immunity_threshold(2.5).unwrap();
    assert!((h - 0.6).abs() < 1e-10);
}

#[test]
fn test_exponential_growth_doubles() {
    let n = growth::exponential_growth(100.0, core::f64::consts::LN_2, 1.0).unwrap();
    assert!((n - 200.0).abs() < 1e-6);
}
