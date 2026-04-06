//! SGP4 orbit propagation via the published `sgp4` crate.

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

    let epochyr: i32 = get_map_val(env, tle_map, "epochyr")?;
    let epochdays: f64 = get_map_val(env, tle_map, "epochdays")?;
    let bstar: f64 = get_map_val(env, tle_map, "bstar")?;
    let ndot: f64 = get_map_val(env, tle_map, "mean_motion_dot")?;
    let nddot: f64 = get_map_val(env, tle_map, "mean_motion_double_dot")?;
    let ecco: f64 = get_map_val(env, tle_map, "eccentricity")?;
    let argpo: f64 = get_map_val(env, tle_map, "arg_perigee_deg")?;
    let inclo: f64 = get_map_val(env, tle_map, "inclination_deg")?;
    let mo: f64 = get_map_val(env, tle_map, "mean_anomaly_deg")?;
    let no_kozai: f64 = get_map_val(env, tle_map, "mean_motion")?;
    let nodeo: f64 = get_map_val(env, tle_map, "raan_deg")?;

    let dt = parse_datetime_tuple(datetime_tuple)?;

    // Build epoch datetime
    let year_full = if epochyr < 57 { 2000 + epochyr } else { 1900 + epochyr };
    let epoch_dt = epoch_to_naive_datetime(year_full, epochdays);

    let elements = sgp4::Elements {
        object_name: None,
        international_designator: None,
        norad_id: 0,
        classification: sgp4::Classification::Unclassified,
        datetime: epoch_dt,
        mean_motion_dot: ndot,
        mean_motion_ddot: nddot,
        drag_term: bstar,
        element_set_number: 0,
        inclination: inclo,
        right_ascension: nodeo,
        eccentricity: ecco,
        argument_of_perigee: argpo,
        mean_anomaly: mo,
        mean_motion: no_kozai,
        revolution_number: 0,
        ephemeris_type: 0,
    };

    let constants = match sgp4::Constants::from_elements_afspc_compatibility_mode(&elements) {
        Ok(c) => c,
        Err(e) => return Ok((error, format!("SGP4 init: {e}")).encode(env)),
    };

    // Build target NaiveDateTime
    let sec_total = dt.second as f64 + dt.microsecond as f64 / 1_000_000.0;
    let whole_sec = sec_total.floor() as u32;
    let nanos = ((sec_total - whole_sec as f64) * 1e9) as u32;
    let target_dt = chrono::NaiveDate::from_ymd_opt(dt.year, dt.month as u32, dt.day as u32)
        .and_then(|d| d.and_hms_nano_opt(dt.hour as u32, dt.minute as u32, whole_sec, nanos))
        .ok_or_else(|| rustler::Error::RaiseTerm(Box::new("invalid datetime")))?;

    // Use the crate's own epoch calculation for consistent tsince
    let tsince = match elements.datetime_to_minutes_since_epoch(&target_dt) {
        Ok(t) => t,
        Err(e) => return Ok((error, format!("tsince: {e}")).encode(env)),
    };

    match constants.propagate_afspc_compatibility_mode(tsince) {
        Ok(pred) => {
            let pos = (pred.position[0], pred.position[1], pred.position[2]);
            let vel = (pred.velocity[0], pred.velocity[1], pred.velocity[2]);
            Ok((ok, (pos, vel)).encode(env))
        }
        Err(e) => Ok((error, format!("SGP4 propagate: {e}")).encode(env)),
    }
}

fn epoch_to_naive_datetime(year: i32, epochdays: f64) -> chrono::NaiveDateTime {
    let jan1 = chrono::NaiveDate::from_ymd_opt(year, 1, 1)
        .unwrap()
        .and_hms_opt(0, 0, 0)
        .unwrap();
    let days_from_jan1 = epochdays - 1.0;
    let whole_days = days_from_jan1.floor() as i64;
    let frac_day = days_from_jan1 - whole_days as f64;
    let frac_nanos = (frac_day * 86400.0 * 1e9) as i64;
    jan1 + chrono::TimeDelta::days(whole_days) + chrono::TimeDelta::nanoseconds(frac_nanos)
}
