//! # Jivanu — Microbiology Engine
//!
//! **जीवाणु** (Hindi: microbe, bacterium)
//!
//! A Rust library for computational microbiology: growth kinetics,
//! metabolism, genetics, epidemiology, biofilm, and antibiotic resistance.
//!
//! ## Modules
//!
//! - [`growth`] — Exponential/logistic/Gompertz growth, Monod kinetics, chemostat,
//!   competition, cardinal temperature/pH models
//! - [`metabolism`] — Michaelis-Menten, enzyme inhibition, Hill/Emax, ATP yields,
//!   stoichiometric networks, flux balance and flux variability analysis
//! - [`genetics`] — Mutation rates, Hardy-Weinberg, GC content, codon table,
//!   amino acid properties, translation
//! - [`epidemiology`] — SIR/SEIR/SIRS models, R0, herd immunity, vaccination
//! - [`biofilm`] — Formation stages, quorum sensing, diffusion
//! - [`resistance`] — Kill curves, MIC, FIC index, checkerboard, drug combinations
//! - [`pharmacokinetics`] — IV/oral one-compartment, two-compartment, AUC
//! - [`stochastic`] — Gillespie SSA and tau-leaping for reaction networks
//! - [`bridge`] — Cross-module models (PK-driven time-kill, PK/PD indices, biofilm)
//! - [`taxonomy`] — Domain, Gram stain, morphology, oxygen requirements
//!
//! ## Feature flags
//!
//! - `hisab` — fourth-order Runge-Kutta integration for the ODE models
//!   (forward Euler otherwise)
//! - `kimiya` — chemistry bridges in [`bridge`] (Arrhenius, pH, half-life)
//!
//! ## Example
//!
//! ```
//! use jivanu::growth;
//!
//! // Monod kinetics: growth rate at half-saturation
//! let mu = growth::monod_kinetics(1.0, 0.5, 1.0)?;
//! assert!((mu - 0.25).abs() < 1e-10); // mu_max/2 at S=K_s
//!
//! // Doubling time
//! let td = growth::doubling_time(0.693)?;
//! assert!((td - 1.0).abs() < 0.01); // ~1 hour
//! # Ok::<(), jivanu::JivanuError>(())
//! ```

pub mod biofilm;
pub mod bridge;
pub mod epidemiology;
pub mod error;
pub mod genetics;
pub mod growth;
mod lp;
pub mod metabolism;
pub mod pharmacokinetics;
pub mod resistance;
mod rng;
pub mod stochastic;
pub mod taxonomy;

pub use error::JivanuError;
