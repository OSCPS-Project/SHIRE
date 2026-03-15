//! # Walker Ideal Model
//!
//! Contains implementation of the ``WalkerModel`` a type of ideal Equation of State Model that
//! also account for the vibrational and rotational modes of molecules in thermo calculations.

//External Imports
use std::sync::Arc;
use uom::si::f64::*;
use uom::si::molar_energy;
use uom::si::molar_heat_capacity;
use uom::si::pressure;
use uom::si::thermodynamic_temperature;
use uom::si::energy;
use uom::si::amount_of_substance;
use uom::si::volume;
use uom::si::ratio;
use uom::si::molar_mass;

//Internal Imports
use crate::thermodynamics::EOSParams;
use crate::thermodynamics::ReferenceStateParameter;
use crate::thermodynamics::ideal::BaseEOSModel;
use crate::thermodynamics::EOSGroupContributionParameter;
use crate::stream::ComponentData;

///# WalkerModel
///
///Another type of ideal equation of state model that will also account for vibrational and
///rotational modes of molecules in thermodynamic calculations.
pub struct WalkerModel {
    ///List of components (coming from 'Stream' struct)
    pub components: Arc<Vec<ComponentData>>,
    /// Overall Molecular weight (SingleParameter)
    pub molecular_weight : Arc<EOSParams>,
    /// Nrot - param for Walker model(SingleParameter)
    pub n_rot : Arc<EOSParams>, 
    /// θ1 - param for Walker model(SingleParameter)
    pub theta_1 : Arc<EOSParams>, 
    /// θ2 - param for Walker model(SingleParameter)
    pub theta_2 : Arc<EOSParams>,
    /// θ3 - param for Walker model(SingleParameter)
    pub theta_3 : Arc<EOSParams>, 
    /// θ4 - param for Walker model(SingleParameter)
    pub theta_4 : Arc<EOSParams>,
    /// deg1 - param for Walker model(SingleParameter)
    pub deg_1 : Arc<EOSParams>,
    /// deg2 - param for Walker model(SingleParameter)
    pub deg_2 : Arc<EOSParams>,
    /// deg3 - param for Walker model(SingleParameter)
    pub deg_3 : Arc<EOSParams>,
    /// deg4 - param for Walker model(SingleParameter)
    pub deg_4 : Arc<EOSParams>,
    /// reference state for the EOS model
    pub reference_state : Arc<ReferenceStateParameter>,
    /// group contributions
    pub eos_groups : Arc<EOSGroupContributionParameter>

}

//Implementing the `BaseEOSModel` trait for the BasIdeal EOS Package. This will the default method
//for calculating the ideal helmholtz free energy.
impl BaseEOSModel for WalkerModel {
    fn components(&self) -> Arc<Vec<ComponentData>> {
        return Arc::clone(&self.components);
    }
    // Overriding the default function for the Walker Ideal Model
    fn ideal_helmholtz(&self, V: Volume, T: ThermodynamicTemperature, z: Vec<AmountOfSubstance>) -> Energy {
        let a_ideal = 0.0;
        let rotational_modes = vec![self.theta_1.as_ref(), self.theta_2.as_ref(), self.theta_3.as_ref(), self.theta_4.as_ref()];
        let vibrational_modes = vec![self.deg_1.as_ref(), self.deg_2.as_ref(), self.deg_3.as_ref(), self.deg_4.as_ref()];
        let n_groups = &self.eos_groups.n_flattened_groups.clone();
        let residual = 0.0;
        let sum_moles = self.total_moles();
        // loop over components
        for (ni, &zi) in n_groups.iter().zip(z.iter()) {
            // Molecular weight for this component
            let molecular_weight_i = ni.iter()
                .map(|&n| n as f64 * sum_moles)
                .collect();
            // Rotational contribution
            //FIXME:Need to figure out how to extract value from self.n_rot
            let n_rotational_mode_i = ni.iter().map(|&n| n as f64 * self.n_rot).collect()/(ni.iter().sum()); 
            // Thermal wavelength
            let lambda = H / (K_B * T.value() * mwi / N_A).sqrt();

            // Ideal translational term
            a_ideal += xlogx(zi.get::<mole>(), N_A / V.value() * lambda.powi(3));

            // Rotational term
            a_ideal += zi.get::<mole>() * (-nroti / 2.0 * T.value().ln());

            // Vibrational term
            let mut vib_sum = 0.0;

            for (k, &ni_k) in ni.iter().enumerate() {
                let vib_inner: f64 = (0..4)
                    .map(|v| g_vib[v] * walker_fi(theta_vib[v], T.value()))
                    .sum();
                vib_sum += ni_k as f64 * vib_inner;
            }

            a_ideal += zi.get::<mole>() * vib_sum;
        }

        // normalize by total moles if desired
        a_ideal /= sum_moles;

        return Energy::new::<energy::joule>(a_ideal);
    }
}

impl WalkerModel {
    ///Constructor for the ``WalkerModel`` struct
    ///
    /// # Arguments
    /// * species - The list of species
    /// * molec_weight - The total molecular weight
    /// * The n_rot, theta, and deg are empirical values that are required for the thermodynamic
    /// calculations
    ///
    /// # Returns
    /// Instance of the ``WalkerModel`` struct
    pub fn new(
        species : Arc<Vec<ComponentData>>, 
        molec_weight : Arc<EOSParams>, 
        n_rot : Arc<EOSParams>, 
        theta_values : Arc<Vec<EOSParams>>, 
        deg_values : Arc<Vec<EOSParams>>, 
        reference_state : Arc<ReferenceStateParameter>,
        eos_groups : Arc<EOSGroupContributionParameter>)
        -> Self {
            return WalkerModel { 
                components: species, 
                molecular_weight: molec_weight, 
                n_rot: n_rot, 
                theta_1: Arc::new(theta_values[0].clone()), 
                theta_2: Arc::new(theta_values[1].clone()), 
                theta_3: Arc::new(theta_values[2].clone()), 
                theta_4: Arc::new(theta_values[3].clone()), 
                deg_1: Arc::new(deg_values[0].clone()), 
                deg_2: Arc::new(deg_values[1].clone()), 
                deg_3: Arc::new(deg_values[2].clone()), 
                deg_4: Arc::new(deg_values[3].clone()), 
                reference_state: reference_state,
                eos_groups: eos_groups };
    }
}


