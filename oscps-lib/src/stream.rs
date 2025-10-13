//! # Stream

use std::fmt;

// NOTE: Temporarily disabled until the thermodynamics crate is thread-safe.
// use crate::thermodynamics::ThermoState;
use crate::simulation::BlockReference;

// HACK: Streams must be connected to something. They cannot be floating.
/// # Stream
///
/// Struct to hold stream information
pub struct Stream {
    /// Instance of ThermoState struct that holds thermodynamic information.
    // pub thermo: Option<ThermoState>, // HACK: Temporarily disable to enable thread-safety.
    /// ID of source block
    pub from: BlockReference,
    /// ID of destination block
    pub to: BlockReference,
}

impl fmt::Debug for Stream {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Stream")
            .field("from", &"BlockReference") // Don't recursively debug
            .field("to", &"BlockReference")
            .finish()
    }
}

impl Stream {
    /// Constructor for 'Stream' struct
    pub fn new(from: BlockReference, to: BlockReference) -> Stream {
        Stream { from, to }
    }
}
