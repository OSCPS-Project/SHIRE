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
use uom::si::action;

use nalgebra::DMatrix;
use std::collections::HashMap;

#[allow(dead_code)]
///# ThermodynamicConstants
///
/// Struct for storing physical constants for thermodynamics.
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
    BoltzmannConstant,    // k_B
    ///Planck's Constant
    PlancksConstant       //h
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
            ThermodynamicConstants::StandardTemperature => {Box::new(ThermodynamicTemperature::new::<thermodynamic_temperature::kelvin>(273.15))} // K
            ThermodynamicConstants::StandardPressure => {Box::new(Pressure::new::<pressure::pascal>(101325.0))}, //Pa
            ThermodynamicConstants::AvogadroNumber => Box::new(6.02214076e23), //Units: particles/mole
            ThermodynamicConstants::BoltzmannConstant => Box::new(HeatCapacity::new::<heat_capacity::joule_per_kelvin>(1.380_649e-23)),
            ThermodynamicConstants::PlancksConstant => Box::new(Action::new::<action::joule_second>(6.62607015e-34)) //J*s
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
    value : f64    
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
pub struct ReferenceStateParameter{

}

///# EOSGroupContributionTypes
///
/// Enumeration to hold the type of groups used within ``EOSGroupContributionParameters``
#[derive(Clone)]
pub enum EOSGroupContributionTypes{
    ///Group contribution for WalkerIdeal EOS
    WalkerIdeal
}

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
    pub n_intragroups : Arc<Vec<Vec<Vec<f64>>>>, 
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
        let n_intragroups : Vec<Vec<Vec<f64>>> = Vec::new();

        return EOSGroupContributionParameter { 
            group_type: group_type, 
            components: components, 
            n_intragroups: Arc::new(n_intragroups), 
            flattened_groups: Arc::new(flattened_groups), 
            n_flattened_groups: Arc::new(n_flattened_groups),
            sourcecsvs: sourcecsvs 
        };
    }
    /// Builds the intragroups for the Group Param
    pub fn build_intragroups(
        &mut self, 
        db_group_contribution_intragroups : &Vec<Vec<((String, String), f64)>>
    ) -> Vec<Vec<Vec<f64>>> {
       let group_names = self.flattened_groups.as_ref().clone();
       let n_groups = group_names.len();
       let n_components = self.components.as_ref().len();

       let mut n_intragroups: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_components);

       for i in 0..n_components {
           // Create an n_groups x n_groups matrix filled with 0.0
           let mut matrix: Vec<Vec<f64>> = vec![vec![0.0; n_groups]; n_groups];
           let gc_pair_i = &db_group_contribution_intragroups[i];
           assert!(!gc_pair_i.is_empty(), "Intragroup information was requested, but is missing from component {}", i);
           for pair_ik in gc_pair_i.iter() {
               let ((k1, k2), val) = pair_ik; // destructure the tuple
               let n1 = group_names
                   .iter()
                   .position(|x| *x == *k1)
                   .expect(&format!("group {} not found", k1));
               let n2 = group_names
                   .iter()
                   .position(|x| *x == *k2)
                   .expect(&format!("group {} not found", k2));

               matrix[n1][n2] = *val;
               matrix[n2][n1] = *val;

           }
           n_intragroups.push(matrix);
       }
       self.n_intragroups = Arc::new(n_intragroups.clone());
       n_intragroups
    }
}


#[cfg(test)]
mod tests_thermo_module {}
