use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

fn bench_monod(c: &mut Criterion) {
    c.bench_function("growth/monod_kinetics", |b| {
        b.iter(|| jivanu::growth::monod_kinetics(black_box(5.0), black_box(0.5), black_box(1.0)))
    });
}

fn bench_michaelis_menten(c: &mut Criterion) {
    c.bench_function("metabolism/michaelis_menten", |b| {
        b.iter(|| {
            jivanu::metabolism::michaelis_menten(black_box(5.0), black_box(10.0), black_box(1.0))
        })
    });
}

fn bench_sir_step(c: &mut Criterion) {
    c.bench_function("epidemiology/sir_step", |b| {
        b.iter(|| {
            jivanu::epidemiology::sir_step(
                black_box(0.99),
                black_box(0.01),
                black_box(0.0),
                black_box(0.5),
                black_box(0.1),
                black_box(0.01),
            )
        })
    });
}

fn bench_gc_content(c: &mut Criterion) {
    let dna = "ATGCATGCATGCATGCATGCATGCATGCATGC";
    c.bench_function("genetics/gc_content_32bp", |b| {
        b.iter(|| jivanu::genetics::gc_content(black_box(dna)))
    });
}

fn bench_exponential_growth(c: &mut Criterion) {
    c.bench_function("growth/exponential_growth", |b| {
        b.iter(|| {
            jivanu::growth::exponential_growth(black_box(100.0), black_box(0.5), black_box(10.0))
        })
    });
}

fn bench_translate_codon(c: &mut Criterion) {
    c.bench_function("genetics/translate_codon_to_aa", |b| {
        b.iter(|| jivanu::genetics::translate_codon_to_aa(black_box("ATG")))
    });
}

fn bench_fic_index(c: &mut Criterion) {
    c.bench_function("resistance/fic_index", |b| {
        b.iter(|| {
            jivanu::resistance::fic_index(
                black_box(0.25),
                black_box(1.0),
                black_box(0.25),
                black_box(1.0),
            )
        })
    });
}

fn bench_kill_curve(c: &mut Criterion) {
    c.bench_function("resistance/kill_curve", |b| {
        b.iter(|| jivanu::resistance::kill_curve(black_box(2.0), black_box(1.0), black_box(1.0)))
    });
}

criterion_group!(
    benches,
    bench_monod,
    bench_michaelis_menten,
    bench_sir_step,
    bench_gc_content,
    bench_exponential_growth,
    bench_translate_codon,
    bench_fic_index,
    bench_kill_curve,
);
criterion_main!(benches);
