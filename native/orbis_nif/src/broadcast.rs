//! Rustler boundary for the `astrodynamics-gnss` broadcast-navigation product.
//!
//! Pure glue: `broadcast_parse/1` decodes RINEX navigation text, calls
//! [`BroadcastEphemeris::from_nav`], and returns a resource handle holding the
//! parsed records. `broadcast_position/4` evaluates one satellite's orbit and
//! clock at an instant via the crate's [`EphemerisSource`] contract; the
//! single-point-positioning solve consumes the same handle via
//! `spp::spp_solve_broadcast/15`. No parsing grammar or orbit math lives here.

use astrodynamics_gnss::ephemeris::{BroadcastEphemeris, EphemerisSource};
use astrodynamics_gnss::{GnssSatelliteId, GnssSystem};
use rustler::{Encoder, Env, Error, NifResult, ResourceArc, Term};

/// Resource handle holding a parsed broadcast-navigation product across calls.
pub struct BroadcastResource {
    pub store: BroadcastEphemeris,
}

#[rustler::resource_impl]
impl rustler::Resource for BroadcastResource {}

/// Map a GNSS single-letter system identifier (e.g. `"G"`) onto the crate's
/// [`GnssSystem`]. Pure identifier translation; mirrors `sp3::system_from_letter`.
fn system_from_letter(letter: &str) -> NifResult<GnssSystem> {
    let c = letter
        .chars()
        .next()
        .ok_or_else(|| Error::Term(Box::new("empty GNSS system letter")))?;
    GnssSystem::from_letter(c)
        .ok_or_else(|| Error::Term(Box::new(format!("unknown GNSS system letter {letter:?}"))))
}

/// Parse RINEX 3.x/4.xx navigation text into a broadcast-ephemeris resource handle.
///
/// Dirty-CPU: parsing a full daily multi-GNSS file is unbounded relative to the
/// 1 ms NIF budget. On a malformed file returns the parser's error as a term.
#[rustler::nif(schedule = "DirtyCpu")]
fn broadcast_parse(text: String) -> NifResult<ResourceArc<BroadcastResource>> {
    let store =
        BroadcastEphemeris::from_nav(&text).map_err(|e| Error::Term(Box::new(e.to_string())))?;
    Ok(ResourceArc::new(BroadcastResource { store }))
}

/// Evaluate `sat`'s broadcast orbit and clock at `t_j2000_s` against a loaded
/// handle.
///
/// `t_j2000_s` is the query instant as a continuous second-of-J2000 in the
/// GPST-aligned scale the crate's [`EphemerisSource`] contract expects (it maps
/// that onto each system's own time — BDT for BeiDou, UTC-referenced for GLONASS
/// — internally). Returns `{x_m, y_m, z_m, clock_s}` — ECEF meters and the
/// satellite clock offset in seconds — or the atom `nil` when the product has no
/// usable ephemeris for that satellite at that instant (the crate returns
/// `None`). The miss is encoded as an atom rather than a tuple of NaNs, which the
/// BEAM cannot represent. Pure glue over
/// [`EphemerisSource::position_clock_at_j2000_s`]; no orbit math or file I/O
/// lives here.
#[rustler::nif]
fn broadcast_position<'a>(
    env: Env<'a>,
    handle: ResourceArc<BroadcastResource>,
    system_letter: String,
    prn: u8,
    t_j2000_s: f64,
) -> NifResult<Term<'a>> {
    let system = system_from_letter(&system_letter)?;
    let sat = GnssSatelliteId::new(system, prn);

    match handle.store.position_clock_at_j2000_s(sat, t_j2000_s) {
        Some(([x_m, y_m, z_m], clock_s)) => Ok((x_m, y_m, z_m, clock_s).encode(env)),
        None => Ok(rustler::types::atom::nil().encode(env)),
    }
}
