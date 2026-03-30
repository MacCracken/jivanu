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

fn bench_competition_step(c: &mut Criterion) {
    let state = jivanu::growth::TwoStrainState {
        n1: 100.0,
        n2: 50.0,
    };
    let params = jivanu::growth::CompetitionParams {
        r1: 0.5,
        r2: 0.4,
        k1: 1000.0,
        k2: 800.0,
        alpha12: 0.6,
        alpha21: 0.3,
    };
    c.bench_function("growth/competition_step", |b| {
        b.iter(|| {
            jivanu::growth::competition_step(black_box(&state), black_box(&params), black_box(0.1))
        })
    });
}

fn bench_net_production(c: &mut Criterion) {
    let net = jivanu::metabolism::MetabolicNetwork::from_reactions(vec![
        jivanu::metabolism::Reaction {
            id: "r1".into(),
            stoichiometry: vec![("A".into(), -1.0), ("B".into(), 2.0), ("C".into(), 1.0)],
            reversible: false,
        },
        jivanu::metabolism::Reaction {
            id: "r2".into(),
            stoichiometry: vec![("B".into(), -1.0), ("D".into(), 1.0)],
            reversible: false,
        },
    ]);
    let fluxes = [1.0, 2.0];
    c.bench_function("metabolism/net_production", |b| {
        b.iter(|| net.net_production(black_box(&fluxes)))
    });
}

fn bench_iv_bolus(c: &mut Criterion) {
    c.bench_function("pharmacokinetics/iv_bolus", |b| {
        b.iter(|| {
            jivanu::pharmacokinetics::iv_bolus_concentration(
                black_box(500.0),
                black_box(50.0),
                black_box(0.1),
                black_box(2.0),
            )
        })
    });
}

fn bench_reverse_complement(c: &mut Criterion) {
    let dna = "ATGCATGCATGCATGCATGCATGCATGCATGC";
    c.bench_function("genetics/reverse_complement_32bp", |b| {
        b.iter(|| jivanu::genetics::reverse_complement(black_box(dna)))
    });
}

fn bench_sirs_step(c: &mut Criterion) {
    let state = jivanu::epidemiology::SirsState {
        s: 0.8,
        i: 0.1,
        r: 0.1,
    };
    let params = jivanu::epidemiology::SirsParams {
        beta: 0.5,
        gamma: 0.1,
        delta: 0.05,
        nu: 0.02,
        dt: 0.01,
    };
    c.bench_function("epidemiology/sirs_step", |b| {
        b.iter(|| jivanu::epidemiology::sirs_step(black_box(&state), black_box(&params)))
    });
}

fn bench_oral_concentration(c: &mut Criterion) {
    c.bench_function("pharmacokinetics/oral_concentration", |b| {
        b.iter(|| {
            jivanu::pharmacokinetics::oral_concentration(
                black_box(500.0),
                black_box(0.8),
                black_box(50.0),
                black_box(1.0),
                black_box(0.1),
                black_box(3.0),
            )
        })
    });
}

fn bench_biofilm_kill_curve(c: &mut Criterion) {
    c.bench_function("bridge/biofilm_kill_curve", |b| {
        b.iter(|| {
            jivanu::bridge::biofilm_kill_curve(
                black_box(5.0),
                black_box(1.0),
                black_box(1.0),
                black_box(jivanu::biofilm::BiofilmStage::Maturation),
            )
        })
    });
}

fn bench_gillespie(c: &mut Criterion) {
    c.bench_function("growth/gillespie_birth_death_100", |b| {
        b.iter(|| {
            jivanu::growth::gillespie_birth_death(
                black_box(100),
                black_box(0.5),
                black_box(0.1),
                black_box(1.0),
                black_box(42),
            )
        })
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
    bench_competition_step,
    bench_net_production,
    bench_iv_bolus,
    bench_reverse_complement,
    bench_sirs_step,
    bench_oral_concentration,
    bench_biofilm_kill_curve,
    bench_gillespie,
);
criterion_main!(benches);
