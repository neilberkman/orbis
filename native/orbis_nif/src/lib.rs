mod atmosphere;
mod conjunction;
mod coordinates;
mod doppler;
mod ephemeris;
mod gauss;
mod iau2000a_data;
mod iers_data;
mod iod;
mod lambert;
mod matrix;
mod nutation;
mod precession;
mod propagation;
mod time_scales;

use rustler::{Env, NifResult, Term};

#[rustler::nif]
fn propagate_with_elements<'a>(
    env: Env<'a>,
    tle_map: Term<'a>,
    datetime_tuple: Term<'a>,
) -> NifResult<Term<'a>> {
    propagation::propagate_with_elements_impl(env, tle_map, datetime_tuple)
}

type Vec3 = (f64, f64, f64);

#[rustler::nif]
#[allow(clippy::too_many_arguments)]
fn teme_to_gcrs(
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    datetime_tuple: Term,
    skyfield_compat: bool,
) -> NifResult<(Vec3, Vec3)> {
    coordinates::teme_to_gcrs_impl(x, y, z, vx, vy, vz, datetime_tuple, skyfield_compat)
}

#[rustler::nif]
fn gcrs_to_itrs(
    x: f64,
    y: f64,
    z: f64,
    datetime_tuple: Term,
    skyfield_compat: bool,
) -> NifResult<(f64, f64, f64)> {
    coordinates::gcrs_to_itrs_impl(x, y, z, datetime_tuple, skyfield_compat)
}

#[rustler::nif]
fn itrs_to_geodetic(
    x: f64,
    y: f64,
    z: f64,
) -> NifResult<(f64, f64, f64)> {
    coordinates::itrs_to_geodetic_impl(x, y, z)
}

#[rustler::nif]
#[allow(clippy::too_many_arguments)]
fn gcrs_to_topocentric(
    sat_x: f64,
    sat_y: f64,
    sat_z: f64,
    station_lat_deg: f64,
    station_lon_deg: f64,
    station_alt_km: f64,
    datetime_tuple: Term,
    skyfield_compat: bool,
) -> NifResult<(f64, f64, f64)> {
    coordinates::gcrs_to_topocentric_impl(
        sat_x, sat_y, sat_z,
        station_lat_deg, station_lon_deg, station_alt_km,
        datetime_tuple, skyfield_compat,
    )
}

#[rustler::nif]
#[allow(clippy::too_many_arguments)]
fn atmosphere_density(
    lat_deg: f64,
    lon_deg: f64,
    alt_km: f64,
    year: i32,
    doy: i32,
    sec: f64,
    f107: f64,
    f107a: f64,
    ap: f64,
) -> NifResult<(f64, f64)> {
    atmosphere::atmosphere_density_impl(lat_deg, lon_deg, alt_km, year, doy, sec, f107, f107a, ap)
}

#[rustler::nif(schedule = "DirtyCpu")]
fn get_body_position(
    file_path: String,
    target_name: String,
    observer_name: String,
    jd_whole: f64,
    jd_fraction: f64,
    skyfield_compat: bool,
) -> NifResult<(f64, f64, f64)> {
    ephemeris::get_body_position_impl(file_path, target_name, observer_name, jd_whole, jd_fraction, skyfield_compat)
}

#[rustler::nif]
fn utc_to_tdb_jd_split(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
) -> NifResult<(f64, f64)> {
    let ts = time_scales::TimeScales::from_utc(year, month, day, hour, minute, second);
    Ok((ts.jd_whole, ts.tdb_fraction))
}

#[rustler::nif]
fn utc_to_tdb_jd(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
) -> NifResult<f64> {
    let ts = time_scales::TimeScales::from_utc(year, month, day, hour, minute, second);
    Ok(ts.jd_tdb)
}

#[rustler::nif]
#[allow(clippy::too_many_arguments)]
fn doppler_compute(
    sat_x: f64,
    sat_y: f64,
    sat_z: f64,
    sat_vx: f64,
    sat_vy: f64,
    sat_vz: f64,
    station_lat_deg: f64,
    station_lon_deg: f64,
    station_alt_km: f64,
    datetime_tuple: Term,
) -> NifResult<(f64, f64)> {
    doppler::doppler_compute_impl(
        sat_x,
        sat_y,
        sat_z,
        sat_vx,
        sat_vy,
        sat_vz,
        station_lat_deg,
        station_lon_deg,
        station_alt_km,
        datetime_tuple,
    )
}

#[rustler::nif]
fn iod_gibbs(r1: Vec3, r2: Vec3, r3: Vec3) -> NifResult<(Vec3, f64, f64, f64)> {
    iod::gibbs_impl(r1, r2, r3)
}

#[rustler::nif]
#[allow(clippy::too_many_arguments)]
fn iod_hgibbs(
    r1: Vec3,
    r2: Vec3,
    r3: Vec3,
    jd1: f64,
    jd2: f64,
    jd3: f64,
) -> NifResult<(Vec3, f64, f64, f64)> {
    iod::hgibbs_impl(r1, r2, r3, jd1, jd2, jd3)
}

#[rustler::nif]
#[allow(clippy::too_many_arguments)]
fn iod_gauss(
    decl1: f64, decl2: f64, decl3: f64,
    rtasc1: f64, rtasc2: f64, rtasc3: f64,
    jd1: f64, jdf1: f64,
    jd2: f64, jdf2: f64,
    jd3: f64, jdf3: f64,
    rseci1: Vec3, rseci2: Vec3, rseci3: Vec3,
) -> NifResult<(Vec3, Vec3)> {
    gauss::gauss_impl(
        decl1, decl2, decl3, rtasc1, rtasc2, rtasc3,
        jd1, jdf1, jd2, jdf2, jd3, jdf3,
        rseci1, rseci2, rseci3,
    )
}

#[rustler::nif]
#[allow(clippy::too_many_arguments)]
fn lambert_battin(
    r1: Vec3,
    r2: Vec3,
    v1: Vec3,
    dm: i32,
    de: i32,
    nrev: i32,
    dtsec: f64,
) -> NifResult<(Vec3, Vec3)> {
    lambert::lambert_battin_impl(r1, r2, v1, dm, de, nrev, dtsec)
}

#[rustler::nif(schedule = "DirtyCpu")]
#[allow(clippy::too_many_arguments)]
fn find_conjunctions<'a>(
    env: Env<'a>,
    line1_a: String,
    line2_a: String,
    line1_b: String,
    line2_b: String,
    start_min: f64,
    end_min: f64,
    step_min: f64,
    threshold_km: f64,
) -> NifResult<Term<'a>> {
    conjunction::conjunction_impl(
        env, &line1_a, &line2_a, &line1_b, &line2_b,
        start_min, end_min, step_min, threshold_km,
    )
}

rustler::init!("Elixir.Orbis.NIF");
