//! # Properties
//!
//! Contains chemical properties for species in the simulation.


extern crate uom;
extern crate pubchem;
use anyhow::Result;
use uom::si::f64::*;
use std::{thread,time::Duration};
use serde::{Serialize, Deserialize};

#[derive(Clone)]
#[allow(dead_code)]
/// Used by the "Chemical" struct to create the pubchem::Compound obj based on
/// either the chemical name or the pubchem id of the chemical
pub enum ChemicalIdentifier {
    /// The PubChem ID of the component.
    PubchemID(u32),
    /// The actual name of the component.
    CompoundName(String),
}

#[derive(Clone)]
#[allow(dead_code)]
/// A struct to store information regarding the chemical properties of a 
/// particular functional group.
pub struct FunctionalGroup {
    /// The (PubChem)[<https://pubchem.ncbi.nlm.nih.gov/>] CID of a compound.
    pub pubchem_obj: pubchem::Compound,
    /// Physical properties of a compound.
    pub properties: ChemicalProperties,
}

#[derive(Clone)]
#[allow(dead_code)]
/// A struct to store information regarding the chemical properties of a 
/// particular substance. The "Chemical" struct is a wrapper for the 
/// pubchem::Compound object
pub struct Chemical {
    /// The (PubChem)[<https://pubchem.ncbi.nlm.nih.gov/>] CID of a compound.
    pub pubchem_obj: pubchem::Compound,
    /// Physical properties of a compound.
    pub properties: ChemicalProperties,
    /// functional groups present for this chemical (group name, db index, multiplicity)
    pub groups: Vec<(FunctionalGroup, i32, i64)>
}

#[allow(dead_code)]
/// Implementation of the chemical of interest.
impl Chemical {
    /// Constructs a new chemical.
    pub fn new(identifier: ChemicalIdentifier) -> Result<Self> {
        let pubchem_chemical_object = match identifier {
            ChemicalIdentifier::PubchemID(id) => pubchem::Compound::new(id),
            ChemicalIdentifier::CompoundName(name) => pubchem::Compound::with_name(name.as_str()),
        };
        let mut request_counter = 0;
        let mut cid_vec = None;
        while request_counter <= 10 {
            match pubchem_chemical_object.cids(){
                Ok(cid_list) => {
                    cid_vec = Some(cid_list);
                    break;
                },
                _ => {
                    request_counter += 1;
                    thread::sleep(Duration::from_secs(10));
                }
            };
        }

        // let cid_vec = pubchem_chemical_object.cids().unwrap();
        let cid: i32 = cid_vec.unwrap()[0];
        let prop = ChemicalProperties::new(cid);
        let groups : Vec<(FunctionalGroup, i32, i64)> = Vec::new();
        Ok(Chemical {
            pubchem_obj: pubchem_chemical_object,
            properties: prop,
            groups : groups
        })
    }
    /// Returns the pubchem object for the compound.
    pub fn get_pubchem_obj(&self) -> &pubchem::Compound {
        &self.pubchem_obj
    }

    /// Returns the "ChemicalProperties" object for the "Chemical" object.
    pub fn get_properties(&self) -> &ChemicalProperties {
        &self.properties
    }
}

#[derive(Clone)]
#[allow(dead_code)]
/// Struct containing physical properties of a chemical species
pub struct ChemicalProperties {}

impl ChemicalProperties{
    /// constructor for the ``ChemicalProperties`` struct
    pub fn new(_cid: i32) -> ChemicalProperties {
        return ChemicalProperties {};
    }
}


#[cfg(test)]
mod chemical_species_tests {
}
