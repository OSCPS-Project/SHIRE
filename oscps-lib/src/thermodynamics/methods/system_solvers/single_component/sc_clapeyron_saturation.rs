//! # SCClapeyronSaturation

// internal imports
use crate::thermodynamics::ideal::BaseEOSModel;

// external imports
use uom::si::f64::*;

/// # SCClapeyronSaturation
/// Calculates Saturated temperature using the Clapeyron equation
pub struct SCClapeyronSaturation {
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

impl SCClapeyronSaturation {
    pub fn tsat_computation() {

    }

    fn solver_setup() {

    }
}

