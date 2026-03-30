//! # Jivanu — Microbiology Engine
//!
//! **जीवाणु** (Hindi: microbe, bacterium)
//!
//! A Rust library for computational microbiology: growth kinetics,
//! metabolism, genetics, epidemiology, biofilm, and antibiotic resistance.
//!
//! ## Modules
//!
//! - [`growth`] — Exponential/logistic growth, Monod kinetics, doubling time
//! - [`metabolism`] — Michaelis-Menten, ATP yields, enzyme inhibition
//! - [`genetics`] — Mutation rates, Hardy-Weinberg, GC content, codon table
//! - [`epidemiology`] — SIR/SEIR models, R0, herd immunity
//! - [`biofilm`] — Formation stages, quorum sensing, diffusion
//! - [`resistance`] — Kill curves, MIC, resistance transfer, drug combinations
//! - [`pharmacokinetics`] — Drug concentration models, PK parameters, AUC
//! - [`taxonomy`] — Domain, Gram stain, morphology, oxygen requirements
//!
//! ## Example
//!
//! ```
//! use jivanu::growth;
//!
//! // Monod kinetics: growth rate at half-saturation
//! let mu = growth::monod_kinetics(1.0, 0.5, 1.0).unwrap();
//! assert!((mu - 0.25).abs() < 1e-10); // mu_max/2 at S=K_s
//!
//! // Doubling time
//! let td = growth::doubling_time(0.693).unwrap();
//! assert!((td - 1.0).abs() < 0.01); // ~1 hour
//! ```

#![warn(missing_docs)]

pub mod biofilm;
pub mod epidemiology;
pub mod error;
pub mod genetics;
pub mod growth;
pub mod metabolism;
pub mod pharmacokinetics;
pub mod resistance;
pub mod taxonomy;

pub use error::JivanuError;
