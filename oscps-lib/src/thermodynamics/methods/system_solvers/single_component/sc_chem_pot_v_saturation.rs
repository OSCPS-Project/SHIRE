//! # SCChemPotVSaturation

// internal imports
use crate::thermodynamics::ideal::BaseEOSModel;

// external imports
use uom::si::f64::*;

/// Uses equality of Chemical Potentials with a volume basis to compute sat pressure. 
/// Default for sat pressure calculations.
pub struct SCChemPotVSaturation {
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

impl SCChemPotVSaturation {
    /// saturated pressure calculations
    fn psat_computation() {

    }

    fn solver_setup() {

    }
}
