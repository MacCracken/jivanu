//! Criterion benchmarks for jivanu hot paths.

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

fn rxn(id: &str, stoich: &[(&str, f64)]) -> jivanu::metabolism::Reaction {
    jivanu::metabolism::Reaction {
        id: id.into(),
        stoichiometry: stoich.iter().map(|&(m, c)| (m.into(), c)).collect(),
        reversible: false,
    }
}

/// Synthetic 17-reaction core network: uptake, a 10-step chain with ATP
/// yield, two bypass branches, byproduct secretion, biomass, maintenance.
fn core_network() -> jivanu::metabolism::MetabolicNetwork {
    let mut rxns = vec![rxn("EX_glc", &[("glc", 1.0)])];
    let mut prev = "glc".to_string();
    for i in 0..10 {
        let next = format!("m{i}");
        let atp = if i % 3 == 0 { 2.0 } else { 0.0 };
        rxns.push(rxn(
            &format!("r{i}"),
            &[(&prev, -1.0), (&next, 1.0), ("atp", atp)],
        ));
        prev = next;
    }
    rxns.push(rxn("bypass_2_5", &[("m2", -1.0), ("m5", 1.0)]));
    rxns.push(rxn(
        "bypass_4_8",
        &[("m4", -1.0), ("m8", 1.0), ("atp", 1.0)],
    ));
    rxns.push(rxn("EX_m3", &[("m3", -1.0)]));
    rxns.push(rxn("EX_m6", &[("m6", -1.0)]));
    rxns.push(rxn("biomass", &[("m9", -1.0), ("atp", -3.0)]));
    rxns.push(rxn("ATPM", &[("atp", -1.0)]));
    jivanu::metabolism::MetabolicNetwork::from_reactions(rxns)
}

fn fba_inputs(net: &jivanu::metabolism::MetabolicNetwork) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = net.n_reactions();
    let mut objective = vec![0.0; n];
    let mut lb = vec![0.0; n];
    let mut ub = vec![1000.0; n];
    ub[0] = 10.0; // glucose uptake
    lb[n - 1] = 1.0; // ATP maintenance
    objective[n - 2] = 1.0; // biomass
    (objective, lb, ub)
}

fn bench_fba(c: &mut Criterion) {
    let net = core_network();
    let (obj, lb, ub) = fba_inputs(&net);
    c.bench_function("metabolism/fba_core_17rxn", |b| {
        b.iter(|| {
            jivanu::metabolism::flux_balance_analysis(
                black_box(&net),
                black_box(&obj),
                black_box(&lb),
                black_box(&ub),
            )
        })
    });
}

fn bench_fva(c: &mut Criterion) {
    let net = core_network();
    let (obj, lb, ub) = fba_inputs(&net);
    c.bench_function("metabolism/fva_core_17rxn", |b| {
        b.iter(|| {
            jivanu::metabolism::flux_variability_analysis(
                black_box(&net),
                black_box(&obj),
                black_box(&lb),
                black_box(&ub),
                black_box(0.9),
            )
        })
    });
}

fn bench_hill(c: &mut Criterion) {
    c.bench_function("metabolism/hill_equation", |b| {
        b.iter(|| {
            jivanu::metabolism::hill_equation(
                black_box(3.0),
                black_box(10.0),
                black_box(2.0),
                black_box(2.5),
            )
        })
    });
}

fn bench_translate_orf(c: &mut Criterion) {
    let orf = "ATGGCTAAAGGTCTGCGTTCTGAAGATACCCAGTGGATTCACTATTTCAACCCG".repeat(6);
    c.bench_function("genetics/translate_orf_324bp", |b| {
        b.iter(|| jivanu::genetics::translate_orf(black_box(&orf)))
    });
}

fn bench_n_strain(c: &mut Criterion) {
    let pops = [100.0, 80.0, 60.0, 40.0, 20.0];
    let rates = [0.5, 0.45, 0.4, 0.35, 0.3];
    let caps = [1000.0, 900.0, 800.0, 700.0, 600.0];
    let alpha: Vec<Vec<f64>> = (0..5)
        .map(|i| (0..5).map(|j| if i == j { 1.0 } else { 0.3 }).collect())
        .collect();
    c.bench_function("growth/n_strain_step_5", |b| {
        b.iter(|| {
            jivanu::growth::n_strain_competition_step(
                black_box(&pops),
                black_box(&rates),
                black_box(&caps),
                black_box(&alpha),
                black_box(0.1),
            )
        })
    });
}

fn birth_death() -> Vec<jivanu::stochastic::StochasticReaction> {
    vec![
        jivanu::stochastic::StochasticReaction {
            id: "birth".into(),
            state_change: vec![(0, 1)],
        },
        jivanu::stochastic::StochasticReaction {
            id: "death".into(),
            state_change: vec![(0, -1)],
        },
    ]
}

fn bench_gillespie_ssa(c: &mut Criterion) {
    let rxns = birth_death();
    c.bench_function("stochastic/gillespie_ssa_100", |b| {
        b.iter(|| {
            jivanu::stochastic::gillespie_ssa(
                black_box(&[100]),
                black_box(&rxns),
                |s: &[u64]| vec![0.5 * s[0] as f64, 0.1 * s[0] as f64],
                black_box(1.0),
                black_box(100_000),
                black_box(42),
            )
        })
    });
}

fn bench_tau_leaping_large(c: &mut Criterion) {
    // λ per leap ≈ 1e5: exercises the PTRS Poisson sampler.
    let rxns = birth_death();
    c.bench_function("stochastic/tau_leaping_1e6_pop", |b| {
        b.iter(|| {
            jivanu::stochastic::tau_leaping(
                black_box(&[1_000_000]),
                black_box(&rxns),
                |s: &[u64]| vec![0.6 * s[0] as f64, 0.5 * s[0] as f64],
                black_box(0.1),
                black_box(1.0),
                black_box(7),
            )
        })
    });
}

fn bench_sir_trajectory(c: &mut Criterion) {
    c.bench_function("epidemiology/sir_trajectory_1000", |b| {
        b.iter(|| {
            jivanu::epidemiology::sir_trajectory(
                black_box(0.99),
                black_box(0.01),
                black_box(0.0),
                black_box(0.3),
                black_box(0.1),
                black_box(0.1),
                black_box(1000),
            )
        })
    });
}

fn bench_oral_tmax(c: &mut Criterion) {
    c.bench_function("pharmacokinetics/oral_tmax", |b| {
        b.iter(|| jivanu::pharmacokinetics::oral_tmax(black_box(1.2), black_box(0.15)))
    });
}

fn bench_iv_time_kill(c: &mut Criterion) {
    c.bench_function("bridge/iv_time_kill_100", |b| {
        b.iter(|| {
            jivanu::bridge::iv_time_kill(
                black_box(500.0),
                black_box(50.0),
                black_box(0.1),
                black_box(5.0),
                black_box(1.0),
                black_box(0.25),
                black_box(100),
            )
        })
    });
}

fn bench_multi_dose(c: &mut Criterion) {
    c.bench_function("bridge/multi_dose_concentration_30", |b| {
        b.iter(|| {
            jivanu::bridge::multi_dose_concentration(
                black_box(500.0),
                black_box(50.0),
                black_box(0.1),
                black_box(8.0),
                black_box(30),
                black_box(36.0),
            )
        })
    });
}

criterion_group!(
    benches,
    // growth
    bench_monod,
    bench_exponential_growth,
    bench_competition_step,
    bench_n_strain,
    bench_gillespie,
    // metabolism
    bench_michaelis_menten,
    bench_hill,
    bench_net_production,
    bench_fba,
    bench_fva,
    // genetics
    bench_gc_content,
    bench_translate_codon,
    bench_reverse_complement,
    bench_translate_orf,
    // epidemiology
    bench_sir_step,
    bench_sirs_step,
    bench_sir_trajectory,
    // pharmacokinetics
    bench_iv_bolus,
    bench_oral_concentration,
    bench_oral_tmax,
    // resistance
    bench_fic_index,
    bench_kill_curve,
    // stochastic
    bench_gillespie_ssa,
    bench_tau_leaping_large,
    // bridge
    bench_biofilm_kill_curve,
    bench_iv_time_kill,
    bench_multi_dose,
);
criterion_main!(benches);
