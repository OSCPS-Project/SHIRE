//! # Thermodynamics
//!
//! This module will hold all the functions related to calculating 
//! themrodynamic properties for the blocks and chemical species.
//!
//! Inspired by: https://github.com/ClapeyronThermo/Clapeyron.jl

///Importing EOSModels
pub mod ideal;
pub mod cubic;

///Importing Supporting Thermodynamic Methods

/// Importing Chemical Properties Used by Thermo Packages
use crate::properties::Chemical;
use crate::stream::ComponentData;

///Importing External Packages
use std::sync::Arc;
use uom::si::f64::*;
use uom::si::heat_capacity;
use uom::si::mass;
use uom::si::pressure;
use uom::si::thermodynamic_temperature;
use uom::si::energy;
use uom::si::amount_of_substance;
use nalgebra::DMatrix;
use std::collections::HashMap;

#[allow(dead_code)]
///# ThermodynamicConstants
///
/// Struct for storing physical constants for thermodynamics.
/// TODO: Reimplement the use of uom for dimensional analysis.
pub enum ThermodynamicConstants {
    /// The Universal gas constant in J/(mol*K)
    UniversalGasConstant, // J/(mol*K)
    /// Standard temperature in K
    StandardTemperature,  // T_0
    /// Standard pressure in Pa
    StandardPressure,     // P_0
    /// Avogadro's number in mol^-1
    AvogadroNumber,       // N_A
    /// Boltzmann Constant
    BoltzmannConstant     // k_B
}

#[allow(dead_code)] 
/// Implements values of thermodynamic constants.
impl ThermodynamicConstants {
    /// Returns the value of the thermodynamic constant with its appropriate type.
    pub fn value(&self) -> Box<dyn std::any::Any> {
        match self {
            ThermodynamicConstants::UniversalGasConstant => {
                let r = 8.314462618;
                let constant = Energy::new::<energy::joule>(r) / (ThermodynamicTemperature::new::<thermodynamic_temperature::kelvin>(1.0)* AmountOfSubstance::new::<amount_of_substance::mole>(1.0));
                Box::new(constant)
            },
            ThermodynamicConstants::StandardTemperature => {
                Box::new(ThermodynamicTemperature::new::<thermodynamic_temperature::kelvin>(273.15))
            }
            ThermodynamicConstants::StandardPressure => {
                Box::new(Pressure::new::<pressure::pascal>(101325.0))
            },
            ThermodynamicConstants::AvogadroNumber => Box::new(6.02214076e23), //Units: particles/mole
            ThermodynamicConstants::BoltzmannConstant => Box::new(HeatCapacity::new::<heat_capacity::joule_per_kelvin>(1.380_649e-23))
        }
    }
}

///# EOSParams
///
/// Will hold in the parameters information from the database that the various EOS models will use
#[derive(Clone)]
pub enum EOSParams {
    ///Will hold like parameters for a single chemical/groups (e.g acentric factor)
    SingleParameterType(SingleParameter),

    ///Will hold like and unlike parameters for a pair of chemical species/groups
    BinaryParameterType(BinaryParameter),

    ///Parameters for associating interactions between two sites on 2 species/groups(e.g
    ///associating energy or bonding volume)
    AssociatingParameterType(AssociatingParameter),

    ///Parameters associated to groups 
    GroupContributionParameterType(EOSGroupContributionParameter),

    ///Will hold parameters from the reference model
    ReferenceStateParameterType(ReferenceStateParameter)
}

///# SingleParameter
///
///Pure component thermodynamic properties
#[derive(Clone)]
pub struct SingleParameter {
    
}

///# BinaryParameter
///
///Thermodynamic properties for binary interactions
#[derive(Clone)]
pub struct BinaryParameter {

}

///# AssociatingParameter
///
///Thermodynamic properties for associating interactions
#[derive(Clone)]
pub struct AssociatingParameter {

}

///# ReferenceState
///
/// Enumeration that will contain the types of reference states that will be used by the different
/// equation of state models.
#[derive(Clone)]
pub enum ReferenceStateType{}

///# ReferenceStateParams
/// 
/// Will contain the parameters that will be part of each enum member in the ``ReferenceState``
/// enumeration.
#[derive(Clone)]
pub struct ReferenceStateParameter{}

///# EOSGroupContributionTypes
///
/// Enumeration to hold the type of groups used within ``EOSGroupContributionParameters``
#[derive(Clone)]
pub enum EOSGroupContributionTypes{}

///# EOSGroupContributionParameters
///
/// Contain struct definition for Group Contributions. Specifically, this struct will perform
/// calculations to study properties of groups (such as CH2, CH3, etc within a long-chain
/// hydrocarbon) and the interactions between groups.
///
/// This can be used to estimate thermodynamic properties using a molecule's functional groups.
///
/// Derived from ClapeyronThermo (GroupParams.jl)
///
#[derive(Clone)]
pub struct EOSGroupContributionParameter {
    ///Type of group contribution
    pub group_type : Arc<EOSGroupContributionTypes>, 
    /// list of the components
    pub components : Arc<Vec<ComponentData>>, 
    /// A list of all the connections between groups
    pub n_intragroups : Arc<Vec<DMatrix<i64>>>, 
    /// A list of all unique groups
    pub flattened_groups : Arc<Vec<String>>, 
    ///multiplicitiy of each unique group for each chemical species 
    pub n_flattened_groups : Arc<Vec<Vec<i64>>>, 
    /// references
    pub sourcecsvs : Arc<Vec<String>>
}

impl EOSGroupContributionParameter {
    ///Constructor for ``EOSGroupsContributionParameters``
    pub fn new(
        group_type : Arc<EOSGroupContributionTypes>,
        components : Arc<Vec<ComponentData>>,
        sourcecsvs : Arc<Vec<String>>,
    ) ->Self {
        let flattened_groups : Vec<String> = Vec::new();
        let n_flattened_groups : Vec<Vec<i64>> = Vec::new();
        let empty_intergroup: DMatrix<i64> = DMatrix::zeros(0, 0);
        let n_intragroups: Vec<DMatrix<i64>> = vec![empty_intergroup.clone(); components.len()];

        return EOSGroupContributionParameter { 
            group_type: group_type, 
            components: components, 
            n_intragroups: Arc::new(n_intragroups), 
            flattened_groups: Arc::new(flattened_groups), 
            n_flattened_groups: Arc::new(n_flattened_groups),
            sourcecsvs: sourcecsvs 
        };
    }
    
    pub fn build_gc_groups(self) {
        
    }
}


#[cfg(test)]
mod tests_thermo_module {}
