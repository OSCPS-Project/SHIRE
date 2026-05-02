//! #SingleComponent
//!
//! Contains lib file for property solvers for a single component system
use crate::thermodynamics::ideal::BaseEOSModel;
use uom::si::f64::*;

/// #SingleComponentSystemSolver
pub struct SingleComponentSystemSolver {
    /// Equation of state model
    pub model: Option<Box<dyn BaseEOSModel>>,
    /// volume
    pub volume: Volume,
    /// Temperature
    pub temperature: ThermodynamicTemperature,
    /// moles
    pub moles: AmountOfSubstance

}

pub trait SingleComponentBaseMethods {

}
