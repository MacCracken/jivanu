//! Basic examples of jivanu microbiology formulas.

fn main() -> Result<(), jivanu::JivanuError> {
    // Growth: Monod kinetics
    let mu = jivanu::growth::monod_kinetics(5.0, 0.5, 1.0)?;
    println!("Monod growth rate (S=5, mu_max=0.5, K_s=1): {mu:.3} /hr");

    // Growth: doubling time
    let td = jivanu::growth::doubling_time(mu)?;
    println!("Doubling time: {td:.2} hr");

    // Metabolism: Michaelis-Menten
    let v = jivanu::metabolism::michaelis_menten(2.0, 10.0, 1.0)?;
    println!("Michaelis-Menten rate (S=2, Vmax=10, Km=1): {v:.2}");

    // Genetics: GC content
    let gc = jivanu::genetics::gc_content("ATGCATGCATGC")?;
    println!("GC content of ATGCATGCATGC: {:.1}%", gc * 100.0);

    // Epidemiology: R0 and herd immunity
    let r0 = jivanu::epidemiology::r_naught(0.5, 0.2)?;
    let h = jivanu::epidemiology::herd_immunity_threshold(r0)?;
    println!("R0 = {r0:.1}, herd immunity = {:.1}%", h * 100.0);

    // Hardy-Weinberg
    let (p2, pq2, q2) = jivanu::genetics::hardy_weinberg(0.6)?;
    println!("Hardy-Weinberg (p=0.6): AA={p2:.2}, Aa={pq2:.2}, aa={q2:.2}");

    Ok(())
}
