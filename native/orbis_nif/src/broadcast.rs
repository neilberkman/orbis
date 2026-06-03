//! Rustler boundary for the `astrodynamics-gnss` broadcast-navigation product.
//!
//! Pure glue: `broadcast_parse/1` decodes RINEX navigation text, calls
//! [`BroadcastStore::from_nav`], and returns a resource handle holding the
//! parsed records. The single-point-positioning solve consumes that handle via
//! `spp::spp_solve_broadcast/15`. No parsing grammar or orbit math lives here.

use astrodynamics_gnss::BroadcastStore;
use rustler::{Error, NifResult, ResourceArc};

/// Resource handle holding a parsed broadcast-navigation product across calls.
pub struct BroadcastResource {
    pub store: BroadcastStore,
}

#[rustler::resource_impl]
impl rustler::Resource for BroadcastResource {}

/// Parse RINEX 3 navigation text into a broadcast-ephemeris resource handle.
///
/// Dirty-CPU: parsing a full daily multi-GNSS file is unbounded relative to the
/// 1 ms NIF budget. On a malformed file returns the parser's error as a term.
#[rustler::nif(schedule = "DirtyCpu")]
fn broadcast_parse(text: String) -> NifResult<ResourceArc<BroadcastResource>> {
    let store =
        BroadcastStore::from_nav(&text).map_err(|e| Error::Term(Box::new(e.to_string())))?;
    Ok(ResourceArc::new(BroadcastResource { store }))
}
