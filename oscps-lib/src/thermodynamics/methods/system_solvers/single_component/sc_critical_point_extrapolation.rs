//! # SCCriticalPointExtrapolation

// internal imports
use crate::thermodynamics::ideal::BaseEOSModel;

// external imports
use uom::si::f64::*;

/// # SCCriticalPointExtrapolation
/// Extrapolates the critical point temperature and pressure
pub struct SCCriticalPointExtrapolation {
    /// Equation of state model
    pub model: Option<Box<dyn BaseEOSModel>>,
    /// liquid phase volume
    pub v_L: Volume,
    /// vapor phase volume
    pub v_V: Volume,
    /// Temperature
    pub temperature: ThermodynamicTemperature,
    /// moles
    pub moles: AmountOfSubstance
}

impl SCCriticalPointExtrapolation {
    fn solver_setup() {

    }

    fn tsat_computation() {

    }

    fn psat_computation() {

    }
}


