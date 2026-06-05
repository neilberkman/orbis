//! Rustler boundary for the `astrodynamics-gnss` SP3 precise-ephemeris product.
//!
//! This module is **pure glue**: it decodes Erlang terms, calls
//! the `astrodynamics_gnss::ephemeris` public APIs, manages the parsed product as a
//! Rustler resource handle, and encodes results back. No SP3 grammar, no unit
//! conversion, and no interpolation numerics live here — those are the crate's
//! responsibility. In particular:
//!
//! - `sp3_parse/1` decodes a byte buffer, calls [`Sp3::parse`], and returns a
//!   [`ResourceArc`] wrapping the parsed product. The bytes are parsed exactly
//!   once; nothing stores a path to re-open per call.
//! - `sp3_position/6` operates on that handle plus a decoded epoch; it never
//!   touches the filesystem.
//! - `sp3_satellite_ids/1` exposes only the parsed header satellite tokens, so
//!   Elixir validation code can compare product identity without re-reading the
//!   file or probing interpolation.

use astrodynamics::time::model::{Instant, JulianDateSplit, TimeScale};
use astrodynamics_gnss::ephemeris::Sp3;
use astrodynamics_gnss::{GnssSatelliteId, GnssSystem};
use rustler::{Encoder, Env, Error, NifResult, ResourceArc, Term};

/// Resource handle holding a parsed SP3 product across NIF calls.
///
/// The parsed [`Sp3`] is read-only after construction, so the handle is shared
/// (`ResourceArc`) and evaluation borrows it immutably. The BEAM GC drops it
/// when the last Elixir reference is collected.
pub struct Sp3Resource {
    pub sp3: Sp3,
}

#[rustler::resource_impl]
impl rustler::Resource for Sp3Resource {}

/// Map a GNSS single-letter system identifier (as the Elixir side passes it,
/// e.g. `"G"`) onto the crate's [`GnssSystem`]. Pure identifier translation.
fn system_from_letter(letter: &str) -> NifResult<GnssSystem> {
    let c = letter
        .chars()
        .next()
        .ok_or_else(|| Error::Term(Box::new("empty GNSS system letter")))?;
    GnssSystem::from_letter(c)
        .ok_or_else(|| Error::Term(Box::new(format!("unknown GNSS system letter {letter:?}"))))
}

/// Map a time-scale abbreviation onto the core [`TimeScale`]. Pure translation;
/// used so an Elixir caller can name the epoch's scale explicitly when it is not
/// the file's own header scale.
fn time_scale_from_abbrev(abbrev: &str) -> NifResult<TimeScale> {
    Ok(match abbrev {
        "UTC" => TimeScale::Utc,
        "TAI" => TimeScale::Tai,
        "TT" => TimeScale::Tt,
        "TDB" => TimeScale::Tdb,
        "GPST" => TimeScale::Gpst,
        "GST" => TimeScale::Gst,
        "BDT" => TimeScale::Bdt,
        other => {
            return Err(Error::Term(Box::new(format!(
                "unknown time scale {other:?}"
            ))))
        }
    })
}

/// Parse an SP3-c / SP3-d byte buffer into a resource handle.
///
/// Dirty-CPU: parsing a full IGS day file is unbounded relative to the 1 ms NIF
/// budget. On success returns the [`Sp3Resource`] handle; on a
/// malformed buffer returns the crate's parse-error reason as an Erlang term.
#[rustler::nif(schedule = "DirtyCpu")]
fn sp3_parse(bytes: rustler::Binary) -> NifResult<ResourceArc<Sp3Resource>> {
    let sp3 = Sp3::parse(bytes.as_slice()).map_err(|e| Error::Term(Box::new(e.to_string())))?;
    Ok(ResourceArc::new(Sp3Resource { sp3 }))
}

/// The file's own header time-scale abbreviation (e.g. `"GPST"`), so the Elixir
/// wrapper can tag a query epoch in the product's native scale.
#[rustler::nif]
fn sp3_time_scale(handle: ResourceArc<Sp3Resource>) -> NifResult<String> {
    Ok(handle.sp3.header.time_scale.abbrev().to_string())
}

/// The SP3/RINEX satellite tokens declared in the product header, e.g. `"G01"`.
#[rustler::nif]
fn sp3_satellite_ids(handle: ResourceArc<Sp3Resource>) -> NifResult<Vec<String>> {
    Ok(handle
        .sp3
        .satellites()
        .iter()
        .map(|sat| sat.to_string())
        .collect())
}

/// Evaluate `sat`'s interpolated state at `epoch` against a loaded handle.
///
/// The epoch is a split Julian date `(jd_whole, jd_fraction)` in the named
/// `scale`. Returns `{x_m, y_m, z_m, clock}` where `clock` is the satellite
/// clock offset in seconds, or the atom `nil` when the satellite has no clock
/// estimate at the epoch (the crate returns `None`). The clock is encoded as a
/// term rather than a float so a missing clock is not forced through `NaN`,
/// which the BEAM cannot represent.
///
/// Operates only on the resource handle — no file I/O.
#[rustler::nif]
fn sp3_position<'a>(
    env: Env<'a>,
    handle: ResourceArc<Sp3Resource>,
    system_letter: String,
    prn: u8,
    scale: String,
    jd_whole: f64,
    jd_fraction: f64,
) -> NifResult<Term<'a>> {
    let system = system_from_letter(&system_letter)?;
    let sat = GnssSatelliteId::new(system, prn);
    let scale = time_scale_from_abbrev(&scale)?;
    let epoch = Instant::from_julian_date(scale, JulianDateSplit::new(jd_whole, jd_fraction));

    let state = handle
        .sp3
        .position(sat, epoch)
        .map_err(|e| Error::Term(Box::new(e.to_string())))?;

    // Encode clock as `nil` when absent so a fixed-arity tuple never carries a
    // NaN float (unrepresentable on the BEAM); the Elixir wrapper maps `nil`
    // straight through to `clock_s: nil`.
    let clock_term: Term<'a> = match state.clock_s {
        Some(c) => c.encode(env),
        None => rustler::types::atom::nil().encode(env),
    };

    Ok((
        state.position.x_m,
        state.position.y_m,
        state.position.z_m,
        clock_term,
    )
        .encode(env))
}
