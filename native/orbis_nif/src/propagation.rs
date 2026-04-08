//! SGP4 orbit propagation via the in-house `astrodynamics::sgp4` module
//! (pure-Rust port of Vallado SGP4, bit-exact to non-FMA Vallado at 0 ULP).

use astrodynamics::sgp4::{ElementSet, JulianDate, OpsMode, Satellite};
use rustler::{Encoder, Env, NifResult, Term};

type DateTuple = (i32, i32, i32);
type TimeTuple = (i32, i32, i32, i32);
struct DatetimeComponents {
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: i32,
    microsecond: i32,
}

fn parse_datetime_tuple(term: Term) -> NifResult<DatetimeComponents> {
    let ((y, m, d), (h, min, s, us)): (DateTuple, TimeTuple) = term.decode()?;
    Ok(DatetimeComponents {
        year: y,
        month: m,
        day: d,
        hour: h,
        minute: min,
        second: s,
        microsecond: us,
    })
}

fn get_map_val<'a, T: rustler::Decoder<'a>>(
    env: Env<'a>,
    map: Term<'a>,
    key: &str,
) -> NifResult<T> {
    let atom = rustler::types::atom::Atom::from_str(env, key)?;
    let val = map.map_get(atom.to_term(env))?;
    val.decode()
}

pub(crate) fn propagate_with_elements_impl<'a>(
    env: Env<'a>,
    tle_map: Term<'a>,
    datetime_tuple: Term<'a>,
) -> NifResult<Term<'a>> {
    let ok = rustler::types::atom::Atom::from_str(env, "ok")?;
    let error = rustler::types::atom::Atom::from_str(env, "error")?;

    let elements = ElementSet {
        epoch_year_two_digit: get_map_val(env, tle_map, "epochyr")?,
        epoch_days: get_map_val(env, tle_map, "epochdays")?,
        bstar: get_map_val(env, tle_map, "bstar")?,
        mean_motion_dot: get_map_val(env, tle_map, "mean_motion_dot")?,
        mean_motion_double_dot: get_map_val(env, tle_map, "mean_motion_double_dot")?,
        eccentricity: get_map_val(env, tle_map, "eccentricity")?,
        argument_of_perigee_deg: get_map_val(env, tle_map, "arg_perigee_deg")?,
        inclination_deg: get_map_val(env, tle_map, "inclination_deg")?,
        mean_anomaly_deg: get_map_val(env, tle_map, "mean_anomaly_deg")?,
        mean_motion_rev_per_day: get_map_val(env, tle_map, "mean_motion")?,
        right_ascension_deg: get_map_val(env, tle_map, "raan_deg")?,
        catalog_number: 0,
    };

    let dt = parse_datetime_tuple(datetime_tuple)?;

    // AFSPC opsmode: matches the historical orbis behavior (which used the
    // third-party sgp4 crate's `_afspc_compatibility_mode` functions). The
    // Skyfield reference values stored in oracle tests were calibrated to
    // AFSPC, so we preserve that mode for compatibility.
    let satellite = match Satellite::from_elements_with_opsmode(&elements, OpsMode::Afspc) {
        Ok(s) => s,
        Err(e) => return Ok((error, format!("SGP4 init: {e}")).encode(env)),
    };

    // Compute the target Julian Date from the supplied UTC components and let
    // Satellite::propagate_jd subtract the satrec's *cached* exact epoch
    // internally. This avoids any precision drift that would otherwise come
    // from computing the epoch JD a second time on this side.
    let target_jd = match utc_components_to_jd_split(&dt) {
        Some(jd) => jd,
        None => return Ok((error, "invalid datetime").encode(env)),
    };

    match satellite.propagate_jd(target_jd) {
        Ok(pred) => {
            let pos = (pred.position[0], pred.position[1], pred.position[2]);
            let vel = (pred.velocity[0], pred.velocity[1], pred.velocity[2]);
            Ok((ok, (pos, vel)).encode(env))
        }
        Err(e) => Ok((error, format!("SGP4 propagate: {e}")).encode(env)),
    }
}

/// Build a split-form Julian date `(jd_whole, jd_fraction)` from a UTC
/// calendar tuple. Returns `None` for an invalid date.
///
/// Splits at midnight: `jd_whole` is the integer JD for the calendar day at
/// noon UTC, and `jd_fraction` is the fractional part (negative for AM,
/// positive for PM-of-the-noon-day).
fn utc_components_to_jd_split(dt: &DatetimeComponents) -> Option<JulianDate> {
    if dt.month < 1 || dt.month > 12 || dt.day < 1 || dt.day > 31 {
        return None;
    }
    let a = (14 - dt.month) / 12;
    let y = dt.year + 4800 - a;
    let m = dt.month + 12 * a - 3;
    let jdn = dt.day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    // jdn is the Julian Day Number at *noon* of the calendar date.
    // Convert to midnight-anchored JD by subtracting 0.5.
    let jd_midnight = jdn as f64 - 0.5;
    let frac = (dt.hour as f64) / 24.0
        + (dt.minute as f64) / 1440.0
        + (dt.second as f64) / 86400.0
        + (dt.microsecond as f64) / 86_400_000_000.0;
    Some(JulianDate(jd_midnight, frac))
}
