//! # SCAntoineSaturation

// internal imports
use crate::thermodynamics::ideal::BaseEOSModel;

// external imports
use uom::si::f64::*;


/// # SCAntoineSaturation
/// Antoine's equation for saturation temperature calcs.
/// Default for sat temp calcs.
pub struct SCAntoineSaturation {
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

impl SCAntoineSaturation {
    /// saturated temp calcalation
    fn tsat_computation() {
        // will call upon the 'pub fn sc_compute_equilibrium_temperature()' function to solve for
        // the Tsat 
    }
    /// setting up solver with correct parameters prior to solving for Tsat
    fn solver_setup() {

    }
}

