//! #SingleComponent
//!
//! Contains lib file for property solvers for a single component system

// internal imports
use crate::thermodynamics::ideal::BaseEOSModel;

// modules that are part of the single component system solver
pub mod saturation;

// external imports
use uom::si::f64::*;

/// #SingleComponentSystemSolver
pub struct SingleComponentSystemSolver {
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

impl SingleComponentSystemSolver {
    /// checks if values are positive
    fn _check_positive_values() {

    }
    
    ///thermodynamic derivative routine that returns how molar volume changes 
    ///with pressure at fixed temperature (and composition). Used in the validity
    ///checks function
    fn _dv_dp_calculation() {

    }

    /// validity checks for mechanical stability of volumes and that liquid and vapor phase
    /// volumes are different by a threshold
    pub fn sc_validity_checks() {
        
    }

    /// base function equilibrium pressure calculation
    pub fn sc_compute_equilibrium_pressure() {

    }

    /// base function equilibrium temperature calculation
    pub fn sc_compute_equilibrium_temperature() {

    }
}


