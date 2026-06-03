//! Rustler boundary for the `astrodynamics-gnss` ionospheric delay models.
//!
//! This module is **pure glue**: it decodes Erlang terms, calls the
//! `astrodynamics_gnss::ionex` public APIs, manages the parsed IONEX product as
//! a Rustler resource handle, and encodes the results back. No Klobuchar
//! polynomial, no single-layer-model geometry, and no grid interpolation lives
//! here — those are the crate's responsibility.
//!
//! - `klobuchar_delay/7` evaluates the GPS broadcast Klobuchar L1 model scaled
//!   to the requested carrier, taking radians at the boundary (the public Elixir
//!   wrapper converts its degree inputs).
//! - `ionex_parse/1` decodes a byte buffer, calls [`Ionex::parse`], and returns
//!   a [`ResourceArc`] wrapping the parsed grid; the bytes are parsed once.
//! - `ionex_slant_delay/7` operates on that handle plus an integer J2000-second
//!   epoch; it never touches the filesystem.

use astrodynamics_gnss::ionex::Ionex;
use astrodynamics_gnss::Wgs84Geodetic;
use astrodynamics_gnss::{ionex_slant_delay, klobuchar_native, KlobucharParams};
use rustler::{Error, NifResult, ResourceArc};

/// Resource handle holding a parsed IONEX product across NIF calls.
///
/// The parsed [`Ionex`] grid is read-only after construction, so the handle is
/// shared (`ResourceArc`) and evaluation borrows it immutably. The BEAM GC
/// drops it when the last Elixir reference is collected.
pub struct IonexResource {
    pub ionex: Ionex,
}

#[rustler::resource_impl]
impl rustler::Resource for IonexResource {}

/// GPS broadcast Klobuchar L1 ionospheric group delay (positive meters).
///
/// All inputs arrive in the model's native boundary units: receiver
/// latitude/longitude and satellite azimuth/elevation in **degrees**, and the
/// GPS **second-of-day** in `[0, 86400)`. The Elixir wrapper supplies these
/// directly (it has the degree inputs and forms the second-of-day from the
/// epoch's integer clock fields), so no angle or time conversion happens at this
/// boundary and the delay is bit-exact to the model reference. `frequency_hz` is
/// the carrier on which the delay is reported (the model is dispersive).
#[rustler::nif]
#[allow(clippy::too_many_arguments)]
fn klobuchar_delay(
    lat_deg: f64,
    lon_deg: f64,
    azimuth_deg: f64,
    elevation_deg: f64,
    t_gps_s: f64,
    frequency_hz: f64,
    alpha: (f64, f64, f64, f64),
    beta: (f64, f64, f64, f64),
) -> NifResult<f64> {
    let params = KlobucharParams {
        alpha: [alpha.0, alpha.1, alpha.2, alpha.3],
        beta: [beta.0, beta.1, beta.2, beta.3],
    };
    Ok(klobuchar_native(
        &params,
        lat_deg,
        lon_deg,
        azimuth_deg,
        elevation_deg,
        t_gps_s,
        frequency_hz,
    ))
}

/// Parse an IONEX byte buffer into a resource handle.
///
/// Dirty-CPU: a full daily IONEX map set is unbounded relative to the 1 ms NIF
/// budget. On success returns the [`IonexResource`] handle; on a malformed
/// buffer returns the crate's parse-error reason as an Erlang term.
#[rustler::nif(schedule = "DirtyCpu")]
fn ionex_parse(bytes: rustler::Binary) -> NifResult<ResourceArc<IonexResource>> {
    let ionex = Ionex::parse(bytes.as_slice()).map_err(|e| Error::Term(Box::new(e.to_string())))?;
    Ok(ResourceArc::new(IonexResource { ionex }))
}

/// IONEX vertical-TEC-grid slant ionospheric group delay (positive meters).
///
/// Operates on the parsed handle plus the receiver geodetic latitude/longitude
/// and the satellite azimuth/elevation in radians. `epoch_j2000_s` is integer
/// seconds since the J2000 epoch so it lands exactly on the product's own epoch
/// axis with no float-rounded time entering the temporal bracket. `frequency_hz`
/// is the carrier on which the delay is reported. No file I/O.
#[rustler::nif]
#[allow(clippy::too_many_arguments)]
fn ionex_slant(
    handle: ResourceArc<IonexResource>,
    lat_rad: f64,
    lon_rad: f64,
    elevation_rad: f64,
    azimuth_rad: f64,
    epoch_j2000_s: i64,
    frequency_hz: f64,
) -> NifResult<f64> {
    let receiver = Wgs84Geodetic::new(lat_rad, lon_rad, 0.0);
    Ok(ionex_slant_delay(
        &handle.ionex,
        receiver,
        elevation_rad,
        azimuth_rad,
        epoch_j2000_s,
        frequency_hz,
    ))
}
