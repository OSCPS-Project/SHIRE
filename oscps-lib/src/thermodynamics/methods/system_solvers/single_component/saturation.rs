//! #SingleComponentSaturationMethods



/// Antoine's equation for saturation temperature calcs (default for sat temp calcs)
pub trait SCAntoineSat {
    /// saturated temp calcalation
    fn tsat_computation() {

    }
}

/// Uses equality of Chemical Potentials with a volume basis to compute sat pressure. Default for
/// sat pressure calculations
pub trait SCChemPotSat {
    /// saturated pressure calculations
    fn psat_computation() {

    }
}

