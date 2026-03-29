//! # Walker Ideal Model
//!
//! Contains implementation of the ``WalkerModel`` a type of ideal Equation of State Model that
//! also account for the vibrational and rotational modes of molecules in thermo calculations.
//! 
//! Ref
//! ---
//! Walker, P. J., & Haslam, A. J. (2020). A new predictive group-contribution 
//! ideal-heat-capacity model and its influence on second-derivative properties 
//! calculated using a free-energy equation of state. 
//! Journal of Chemical and Engineering Data, 65(12), 5809–5829. 
//! [doi:10.1021/acs.jced.0c00723](https://doi.org/10.1021/acs.jced.0c00723)

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
use uom::si::action;
use uom::si::heat_capacity;

//Internal Imports
use crate::thermodynamics::EOSParams;
use crate::thermodynamics::ReferenceStateParameter;
use crate::thermodynamics::ideal::BaseEOSModel;
use crate::thermodynamics::EOSGroupContributionParameter;
use crate::thermodynamics::ThermodynamicConstants::*;
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
    // FIXME: Check whether implemented correctly!
    fn ideal_helmholtz(&self, V: Volume, T: ThermodynamicTemperature, z: Vec<AmountOfSubstance>) -> Energy {
        let k_b : f64 = BoltzmannConstant.value().downcast_ref::<HeatCapacity>().unwrap().get::<heat_capacity::joule_per_kelvin>();
        let h : f64 = PlancksConstant.value().downcast_ref::<Action>().unwrap().get::<action::joule_second>();
        let n_a : f64 = *AvogadroNumber.value().downcast_ref::<f64>().unwrap();
        let vol : f64 = V.get::<volume::cubic_meter>();
        let temp : f64 = T.get::<thermodynamic_temperature::kelvin>();
        
        let mut a_ideal = 0.0;
        let rotational_modes = vec![self.theta_1.as_ref(), self.theta_2.as_ref(), self.theta_3.as_ref(), self.theta_4.as_ref()];
        let vibrational_modes = vec![self.deg_1.as_ref(), self.deg_2.as_ref(), self.deg_3.as_ref(), self.deg_4.as_ref()];
        let n_groups = &self.eos_groups.n_flattened_groups.clone();
        let sum_moles = self.total_moles();
        // loop over components
        for (ni, &zi) in n_groups.iter().zip(z.iter()) {
            // Molecular weight for this component
            let molecular_weight_i : f64 = ni.iter()
                .map(|&n| (sum_moles * n as f64).get::<amount_of_substance::mole>()) // f64
                .sum();
            // Rotational contribution
            let n_rot_value = match &*self.n_rot {
                EOSParams::SingleParameterType(param) => param.value.unwrap(), // adjust field name
                _ => panic!("Expected SingleParameterType for n_rot"),
            };
            let n_rotational_mode_i: f64 = ni.iter()
                .map(|&n| n as f64 * n_rot_value)
                .sum::<f64>()   // sum the weighted values
                / ni.iter().map(|&n| n as f64).sum::<f64>(); // sum of n
            // Thermal wavelength
            let lambda : f64 = h / (k_b * temp * molecular_weight_i / n_a).sqrt();

            // Ideal translational term
            a_ideal += self.xlogx(zi.get::<amount_of_substance::mole>(), n_a / vol * lambda.powi(3));

            // Rotational term
            a_ideal += zi.get::<amount_of_substance::mole>() * (-n_rotational_mode_i / 2.0 * temp.ln());

            // Vibrational term
            let mut vib_sum = 0.0;

            for (k, &ni_k) in ni.iter().enumerate() {
                let vib_inner: f64 = (0..4)
                    .map(|v| {
                        let vi: f64 = match &*vibrational_modes[v] {
                            EOSParams::SingleParameterType(param) => param.value, // adjust field name
                            _ => panic!("Expected SingleParameterType for n_rot"),
                        };
                        vi * self.walker_fi(rotational_modes[v], temp) // <-- no extra ')'
                    })
                .sum();
                vib_sum += ni_k as f64 * vib_inner;
            }

            a_ideal += zi.get::<amount_of_substance::mole>() * vib_sum;
        }

        // normalize by total moles if desired
        a_ideal /= sum_moles.get::<amount_of_substance::mole>();

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
    fn xlogx(&self, x: f64, y: f64) -> f64 {
        if x <= 0.0 || y <= 0.0 {
            0.0 // convention: 0 * ln(anything) = 0
        } else {
            x * y.ln()
        }
    }

    fn walker_fi(&self, theta : &EOSParams , T : f64) -> f64 {
        let theta_val : f64 = match &*theta {
                EOSParams::SingleParameterType(param) => param.value.unwrap(), // adjust field name
                _ => panic!("Expected SingleParameterType for n_rot"),
            };
        if(theta_val != 0.0) {
            return (1.0 - (-theta_val / T).exp()).ln() + theta_val / (2.0 * T);
        }
        else {
            return 0.0;
        }
    }
}


