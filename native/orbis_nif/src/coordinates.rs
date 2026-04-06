//! Coordinate transformation pipeline.
//!
//! TEME → GCRS replicates Skyfield's exact computation path including AU/day
//! unit scaling for bit-exact (0 ULP) parity.
//!
//! Also provides GCRS → ITRS, ITRS → geodetic (WGS84), and topocentric
//! (az/el/range) transformations.

use rustler::{NifResult, Term};

use crate::matrix::{inline_mxmxm, inline_rxr, inline_tr, Mat3};
use crate::nutation::{
    build_skyfield_nutation_matrix, skyfield_equation_of_the_equinoxes_complimentary_terms,
    skyfield_iau2000a_radians, skyfield_mean_obliquity_radians,
};
use crate::precession::{build_icrs_to_j2000, compute_skyfield_precession_matrix};
use crate::time_scales::TimeScales;

const AU_KM: f64 = 149597870.700;
const DAY_S: f64 = 86400.0;
const TAU: f64 = std::f64::consts::TAU;
const T0: f64 = 2451545.0;

// WGS84 ellipsoid constants
const WGS84_A: f64 = 6378.137; // semi-major axis in km
const WGS84_F: f64 = 1.0 / 298.257223563; // flattening
const WGS84_E2: f64 = 2.0 * WGS84_F - WGS84_F * WGS84_F; // first eccentricity squared

/// Final matrix-vector multiply using explicit FMA.
/// This matches numpy's vectorized behavior and is the ONLY place
/// where f64::mul_add() should be used.
fn mat3_vec3_mul_fma(r: &Mat3, p: &[f64; 3]) -> [f64; 3] {
    let mut result = [0.0_f64; 3];
    for i in 0..3 {
        let sum = r[i][0] * p[0];
        let sum = f64::mul_add(r[i][1], p[1], sum);
        let sum = f64::mul_add(r[i][2], p[2], sum);
        result[i] = sum;
    }
    result
}

fn build_rot_z(angle: f64) -> Mat3 {
    let c = angle.cos();
    let s = angle.sin();
    [
        [c, -s, 0.0],
        [s, c, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

fn earth_rotation_angle(jd_whole: f64, ut1_fraction: f64) -> f64 {
    let days_since_j2000 = jd_whole - T0 + ut1_fraction;
    // Force separate rounded operations to match Skyfield/Python's path.
    let spins_since_j2000: f64 = {
        let v = 0.00273781191135448 * days_since_j2000;
        // Use black_box-like pattern to prevent optimization
        let v_stored: f64 = v;
        v_stored
    };
    let th = 0.7790572732640 + spins_since_j2000;
    let mut result = (th % 1.0 + jd_whole % 1.0 + ut1_fraction) % 1.0;
    if result < 0.0 {
        result += 1.0;
    }
    result
}

fn compute_theta_gmst1982(jd_whole: f64, ut1_fraction: f64) -> f64 {
    let t = (jd_whole - T0 + ut1_fraction) / 36525.0;
    let g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t) * t) * t;
    let mut theta = ((jd_whole % 1.0) + ut1_fraction + (g / DAY_S) % 1.0) % 1.0 * TAU;
    if theta < 0.0 {
        theta += TAU;
    }
    theta
}

fn sidereal_time_hours(jd_whole: f64, ut1_fraction: f64, tdb_fraction: f64) -> f64 {
    let theta = earth_rotation_angle(jd_whole, ut1_fraction);
    let t = (jd_whole - T0 + tdb_fraction) / 36525.0;
    let st = 0.014506
        + ((((-0.0000000368 * t - 0.000029956) * t - 0.00000044) * t + 1.3915817) * t
            + 4612.156534)
            * t;
    let mut result = (st / 54000.0 + theta * 24.0) % 24.0;
    if result < 0.0 {
        result += 24.0;
    }
    result
}

fn gast_radians(ts: &TimeScales, dpsi: f64) -> f64 {
    let gmst_hours = sidereal_time_hours(ts.jd_whole, ts.ut1_fraction, ts.tdb_fraction);
    let mean_ob = skyfield_mean_obliquity_radians(ts.jd_tdb);
    let c_terms = skyfield_equation_of_the_equinoxes_complimentary_terms(ts.jd_tt);
    let eq_eq = dpsi * mean_ob.cos() + c_terms;
    let mut gast_hours = (gmst_hours + eq_eq / TAU * 24.0) % 24.0;
    if gast_hours < 0.0 {
        gast_hours += 24.0;
    }
    gast_hours / 24.0 * TAU
}

/// Build the TEME→GCRS rotation matrix T from time scales.
fn build_teme_to_gcrs_matrix(ts: &TimeScales, skyfield_compat: bool) -> Mat3 {
    let (dpsi, deps) = skyfield_iau2000a_radians(ts.jd_tt);
    let mean_ob = skyfield_mean_obliquity_radians(ts.jd_tdb);
    let true_ob = mean_ob + deps;

    let n = build_skyfield_nutation_matrix(mean_ob, true_ob, dpsi);
    let p = compute_skyfield_precession_matrix(ts.jd_tdb);
    let b = build_icrs_to_j2000();

    // Skyfield uses Kahan-compensated triple product (matching numpy einsum).
    // Direct mode uses standard sequential multiply (more precise).
    let m = if skyfield_compat {
        inline_mxmxm(&n, &p, &b)
    } else {
        let np = inline_rxr(&n, &p);
        inline_rxr(&np, &b)
    };

    let gast = gast_radians(ts, dpsi);
    let theta = compute_theta_gmst1982(ts.jd_whole, ts.ut1_fraction);
    let angle = theta - gast;

    let r = build_rot_z(angle);
    let g = inline_rxr(&r, &m);
    inline_tr(&g)
}

/// Standard (non-FMA) matrix-vector multiply.
pub(crate) fn mat3_vec3_mul(r: &Mat3, p: &[f64; 3]) -> [f64; 3] {
    let mut result = [0.0_f64; 3];
    for i in 0..3 {
        let mut sum = 0.0;
        for j in 0..3 {
            sum += r[i][j] * p[j];
        }
        result[i] = sum;
    }
    result
}

/// Core TEME→GCRS transform. Returns ((px,py,pz), (vx,vy,vz)).
#[allow(clippy::too_many_arguments)]
pub(crate) fn teme_to_gcrs_compute(
    x: f64, y: f64, z: f64,
    vx: f64, vy: f64, vz: f64,
    ts: &TimeScales,
    skyfield_compat: bool,
) -> (Vec3, Vec3) {
    let t = build_teme_to_gcrs_matrix(ts, skyfield_compat);

    if skyfield_compat {
        // AU/day scaling + FMA multiply matching Skyfield's _at() path.
        let r_au = [x / AU_KM, y / AU_KM, z / AU_KM];
        let r_gcrs_au = mat3_vec3_mul_fma(&t, &r_au);
        let r_gcrs = (r_gcrs_au[0] * AU_KM, r_gcrs_au[1] * AU_KM, r_gcrs_au[2] * AU_KM);

        let v_au_d = [vx / AU_KM * DAY_S, vy / AU_KM * DAY_S, vz / AU_KM * DAY_S];
        let v_gcrs_au_d = mat3_vec3_mul_fma(&t, &v_au_d);
        let v_gcrs = (
            v_gcrs_au_d[0] * AU_KM / DAY_S,
            v_gcrs_au_d[1] * AU_KM / DAY_S,
            v_gcrs_au_d[2] * AU_KM / DAY_S,
        );
        (r_gcrs, v_gcrs)
    } else {
        // Direct km/s multiply — no AU round-trip, no FMA.
        let r_teme = [x, y, z];
        let r_g = mat3_vec3_mul(&t, &r_teme);
        let v_teme = [vx, vy, vz];
        let v_g = mat3_vec3_mul(&t, &v_teme);
        ((r_g[0], r_g[1], r_g[2]), (v_g[0], v_g[1], v_g[2]))
    }
}

/// NIF entry point for teme_to_gcrs.
#[allow(clippy::too_many_arguments)]
pub(crate) fn teme_to_gcrs_impl(
    x: f64, y: f64, z: f64,
    vx: f64, vy: f64, vz: f64,
    datetime_tuple: Term,
    skyfield_compat: bool,
) -> NifResult<(Vec3, Vec3)> {
    let (year, month, day, hour, minute, second, microsecond) =
        parse_datetime_tuple(datetime_tuple)?;

    let second_with_micro = second as f64 + microsecond as f64 / 1_000_000.0;
    let ts = TimeScales::from_utc(year, month, day, hour, minute, second_with_micro);

    Ok(teme_to_gcrs_compute(x, y, z, vx, vy, vz, &ts, skyfield_compat))
}

type DateTuple = (i32, i32, i32);
type TimeTuple = (i32, i32, i32, i32);
pub(crate) type Vec3 = (f64, f64, f64);

pub(crate) fn parse_datetime_tuple(term: Term) -> NifResult<(i32, i32, i32, i32, i32, i32, i32)> {
    let (date, time): (DateTuple, TimeTuple) = term.decode()?;
    Ok((date.0, date.1, date.2, time.0, time.1, time.2, time.3))
}

// ---------------------------------------------------------------------------
// GCRS → ITRS (Earth-fixed / ECEF)
// ---------------------------------------------------------------------------

/// Build the full GCRS→ITRS rotation matrix for a given time.
///
/// This is the transpose of the ITRS→GCRS matrix, which combines
/// precession, nutation, and Earth rotation.
pub(crate) fn gcrs_to_itrs_matrix(ts: &TimeScales) -> Mat3 {
    let (dpsi, deps) = skyfield_iau2000a_radians(ts.jd_tt);
    let mean_ob = skyfield_mean_obliquity_radians(ts.jd_tdb);
    let true_ob = mean_ob + deps;

    let n = build_skyfield_nutation_matrix(mean_ob, true_ob, dpsi);
    let p = compute_skyfield_precession_matrix(ts.jd_tdb);
    let b = build_icrs_to_j2000();

    // Celestial-to-terrestrial: combine precession, nutation, frame bias
    let m = inline_mxmxm(&n, &p, &b);

    let gast = gast_radians(ts, dpsi);

    // GAST rotation takes us from true-equator-equinox to ITRS
    let r_gast = build_rot_z(-gast);

    // GCRS→ITRS = R_z(-GAST) * (N * P * B)
    inline_rxr(&r_gast, &m)
}

/// Core GCRS→ITRS transform. Returns (x, y, z) in km.
pub(crate) fn gcrs_to_itrs_compute(
    x: f64, y: f64, z: f64,
    ts: &TimeScales,
    skyfield_compat: bool,
) -> (f64, f64, f64) {
    let mat = gcrs_to_itrs_matrix(ts);

    if skyfield_compat {
        // Skyfield: mxv(R, pos_au) in AU, then convert to km.
        // For ITRS, scalar (non-FMA) multiply matches einsum's rounding.
        // (Unlike TEME→GCRS where FMA is needed — the difference is due to
        // the specific matrix/vector values and how rounding interacts.)
        let pos_au = [x / AU_KM, y / AU_KM, z / AU_KM];
        let r = mat3_vec3_mul(&mat, &pos_au);
        (r[0] * AU_KM, r[1] * AU_KM, r[2] * AU_KM)
    } else {
        let pos = [x, y, z];
        let r = mat3_vec3_mul(&mat, &pos);
        (r[0], r[1], r[2])
    }
}

/// NIF entry point for gcrs_to_itrs.
pub(crate) fn gcrs_to_itrs_impl(
    x: f64, y: f64, z: f64,
    datetime_tuple: Term,
    skyfield_compat: bool,
) -> NifResult<(f64, f64, f64)> {
    let (year, month, day, hour, minute, second, microsecond) =
        parse_datetime_tuple(datetime_tuple)?;
    let second_with_micro = second as f64 + microsecond as f64 / 1_000_000.0;
    let ts = TimeScales::from_utc(year, month, day, hour, minute, second_with_micro);
    Ok(gcrs_to_itrs_compute(x, y, z, &ts, skyfield_compat))
}

// ---------------------------------------------------------------------------
// ITRS → Geodetic (WGS84 lat/lon/alt)
// ---------------------------------------------------------------------------

/// Convert ECEF/ITRS (km) to geodetic coordinates.
/// Returns (latitude_deg, longitude_deg, altitude_km).
///
/// Replicates Skyfield's exact algorithm (wgs84.subpoint / _compute_latitude)
/// which works in AU with exactly 3 iterations.
pub(crate) fn itrs_to_geodetic_compute(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    // Convert to AU to match Skyfield's computation path.
    let x_au = x / AU_KM;
    let y_au = y / AU_KM;
    let z_au = z / AU_KM;

    let a_au = WGS84_A / AU_KM; // Earth equatorial radius in AU
    let r_xy = (x_au * x_au + y_au * y_au).sqrt();

    // Longitude: match Skyfield's exact normalization:
    // (arctan2(y, x) - pi) % tau - pi
    // Python's % always returns positive; Rust's can be negative.
    let lon_raw = y_au.atan2(x_au);
    let pi = std::f64::consts::PI;
    let mut lon_shifted = (lon_raw - pi) % TAU;
    if lon_shifted < 0.0 { lon_shifted += TAU; }
    let lon = lon_shifted - pi;

    // Latitude: 3 iterations matching Skyfield exactly
    let mut lat = z_au.atan2(r_xy);
    let mut a_c = 0.0_f64;
    let mut hyp = 0.0_f64;

    for _ in 0..3 {
        let sin_lat = lat.sin();
        let e2_sin_lat = WGS84_E2 * sin_lat;
        a_c = a_au / (1.0 - e2_sin_lat * sin_lat).sqrt();
        hyp = z_au + a_c * e2_sin_lat;
        lat = hyp.atan2(r_xy);
    }

    // Elevation in AU, then convert to km
    let height_au = (hyp * hyp + r_xy * r_xy).sqrt() - a_c;
    let alt = height_au * AU_KM;

    // Skyfield's Angle.degrees uses: radians * 360.0 / tau
    // This gives different rounding than radians * (180.0 / PI).
    (lat * 360.0 / TAU, lon * 360.0 / TAU, alt)
}

/// NIF entry point for itrs_to_geodetic.
pub(crate) fn itrs_to_geodetic_impl(
    x: f64,
    y: f64,
    z: f64,
) -> NifResult<(f64, f64, f64)> {
    Ok(itrs_to_geodetic_compute(x, y, z))
}

// ---------------------------------------------------------------------------
// Topocentric (az/el/range) from ground station to satellite
// ---------------------------------------------------------------------------

/// Convert geodetic (lat_deg, lon_deg, alt_km) to ECEF/ITRS (km).
pub(crate) fn geodetic_to_itrs(lat_deg: f64, lon_deg: f64, alt_km: f64) -> (f64, f64, f64) {
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let sin_lon = lon.sin();
    let cos_lon = lon.cos();

    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();

    let x = (n + alt_km) * cos_lat * cos_lon;
    let y = (n + alt_km) * cos_lat * sin_lon;
    let z = (n * (1.0 - WGS84_E2) + alt_km) * sin_lat;

    (x, y, z)
}

/// Compute station ECEF/ITRS position directly in AU.
/// Matches Skyfield's Geoid.latlon which works in AU from the start,
/// avoiding the km→AU_KM division that introduces 1 ULP rounding.
fn geodetic_to_itrs_au(lat_deg: f64, lon_deg: f64, alt_km: f64) -> [f64; 3] {
    let lat = lat_deg * TAU / 360.0;
    let lon = lon_deg * TAU / 360.0;

    let sinphi = lat.sin();
    let cosphi = lat.cos();

    let radius_au = WGS84_A / AU_KM;
    let elevation_au = alt_km / AU_KM;

    let omf2 = (1.0 - WGS84_F) * (1.0 - WGS84_F);
    let c = 1.0 / (cosphi * cosphi + sinphi * sinphi * omf2).sqrt();
    let s = omf2 * c;

    let radius_xy = radius_au * c;
    let xy = (radius_xy + elevation_au) * cosphi;
    let x = xy * lon.cos();
    let y = xy * lon.sin();

    let radius_z = radius_au * s;
    let z = (radius_z + elevation_au) * sinphi;

    [x, y, z]
}

/// Build the ECEF→ENU rotation matrix for a given geodetic position.
fn ecef_to_enu_matrix(lat_deg: f64, lon_deg: f64) -> Mat3 {
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let sin_lon = lon.sin();
    let cos_lon = lon.cos();

    // ENU rotation matrix:
    // E = [-sin(lon),           cos(lon),          0       ]
    // N = [-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)]
    // U = [ cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)]
    [
        [-sin_lon, cos_lon, 0.0],
        [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
        [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat],
    ]
}

/// Compute topocentric az/el/range from a ground station to a satellite.
///
/// Returns (azimuth_deg, elevation_deg, range_km).
#[allow(clippy::too_many_arguments)]
pub(crate) fn gcrs_to_topocentric_compute(
    sat_x: f64,
    sat_y: f64,
    sat_z: f64,
    station_lat_deg: f64,
    station_lon_deg: f64,
    station_alt_km: f64,
    ts: &TimeScales,
    skyfield_compat: bool,
) -> (f64, f64, f64) {
    if skyfield_compat {
        return gcrs_to_topocentric_skyfield(
            sat_x, sat_y, sat_z,
            station_lat_deg, station_lon_deg, station_alt_km, ts,
        );
    }

    // Standard path: GCRS→ITRS→subtract→ENU
    let (sat_itrs_x, sat_itrs_y, sat_itrs_z) =
        gcrs_to_itrs_compute(sat_x, sat_y, sat_z, ts, false);

    let (stn_x, stn_y, stn_z) =
        geodetic_to_itrs(station_lat_deg, station_lon_deg, station_alt_km);

    let dx = sat_itrs_x - stn_x;
    let dy = sat_itrs_y - stn_y;
    let dz = sat_itrs_z - stn_z;

    let enu_mat = ecef_to_enu_matrix(station_lat_deg, station_lon_deg);
    let enu = mat3_vec3_mul(&enu_mat, &[dx, dy, dz]);
    let east = enu[0];
    let north = enu[1];
    let up = enu[2];

    // Range
    let range = (east * east + north * north + up * up).sqrt();

    // Elevation
    let elevation = (up / range).asin().to_degrees();

    // Azimuth (measured clockwise from north)
    let mut azimuth = east.atan2(north).to_degrees();
    if azimuth < 0.0 {
        azimuth += 360.0;
    }

    (azimuth, elevation, range)
}

/// Skyfield-compatible topocentric: stays in GCRS AU the entire time.
///
/// Replicates Skyfield's altaz computation:
/// 1. R_lat = rot_y(lat)[::-1]  (row-reversed Y rotation)
/// 2. R_latlon = mxm(R_lat, rot_z(-lon))
/// 3. R_full = mxm(R_latlon, itrs_rotation)
/// 4. station_gcrs_au = transpose(itrs_rotation) * station_itrs_au
/// 5. diff_au = sat_gcrs_au - station_gcrs_au
/// 6. enu_au = mxv(R_full, diff_au)
/// 7. to_spherical(enu_au) → (range_au, elevation_rad, azimuth_rad)
fn gcrs_to_topocentric_skyfield(
    sat_x: f64, sat_y: f64, sat_z: f64,
    station_lat_deg: f64, station_lon_deg: f64, station_alt_km: f64,
    ts: &TimeScales,
) -> (f64, f64, f64) {
    let lat_rad = station_lat_deg * TAU / 360.0;
    let lon_rad = station_lon_deg * TAU / 360.0;

    // Build R_lat = rot_y(lat)[::-1]  (rows reversed)
    let cy = lat_rad.cos();
    let sy = lat_rad.sin();
    // rot_y(lat) = [[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]]
    // [::-1] reverses rows: [[-sy, 0, cy], [0, 1, 0], [cy, 0, sy]]
    let r_lat: Mat3 = [
        [-sy, 0.0, cy],
        [0.0, 1.0, 0.0],
        [cy, 0.0, sy],
    ];

    // R_latlon = mxm(R_lat, rot_z(-lon))
    let rz_neg_lon = build_rot_z(-lon_rad);
    let r_latlon = inline_rxr(&r_lat, &rz_neg_lon);

    // R_full = mxm(R_latlon, itrs_rotation)
    let r_itrs = gcrs_to_itrs_matrix(ts);
    let r_full = inline_rxr(&r_latlon, &r_itrs);

    // Station ITRS position directly in AU, matching Skyfield's Geoid.latlon
    // which computes in AU from the start (not km then / AU_KM).
    let stn_itrs_au = geodetic_to_itrs_au(station_lat_deg, station_lon_deg, station_alt_km);

    // Station GCRS AU = transpose(R_itrs) * station_itrs_au
    let r_itrs_t = inline_tr(&r_itrs);
    let stn_gcrs_au = mat3_vec3_mul(&r_itrs_t, &stn_itrs_au);

    // Satellite GCRS in AU
    let sat_au = [sat_x / AU_KM, sat_y / AU_KM, sat_z / AU_KM];

    // Difference vector in GCRS AU
    let diff_au = [
        sat_au[0] - stn_gcrs_au[0],
        sat_au[1] - stn_gcrs_au[1],
        sat_au[2] - stn_gcrs_au[2],
    ];

    // Rotate to ENU-ish frame: mxv(R_full, diff_au)
    let enu_au = mat3_vec3_mul(&r_full, &diff_au);

    // to_spherical: r, theta (elevation), phi (azimuth)
    let ex = enu_au[0];
    let ey = enu_au[1];
    let ez = enu_au[2];

    let r_au = (ex * ex + ey * ey + ez * ez).sqrt();
    let elevation_rad = ez.atan2((ex * ex + ey * ey).sqrt());
    let mut azimuth_rad = ey.atan2(ex) % TAU;
    if azimuth_rad < 0.0 { azimuth_rad += TAU; }

    let range_km = r_au * AU_KM;
    let elevation_deg = elevation_rad * 360.0 / TAU;
    let azimuth_deg = azimuth_rad * 360.0 / TAU;

    (azimuth_deg, elevation_deg, range_km)
}

/// NIF entry point for gcrs_to_topocentric.
#[allow(clippy::too_many_arguments)]
pub(crate) fn gcrs_to_topocentric_impl(
    sat_x: f64, sat_y: f64, sat_z: f64,
    station_lat_deg: f64, station_lon_deg: f64, station_alt_km: f64,
    datetime_tuple: Term,
    skyfield_compat: bool,
) -> NifResult<(f64, f64, f64)> {
    let (year, month, day, hour, minute, second, microsecond) =
        parse_datetime_tuple(datetime_tuple)?;
    let second_with_micro = second as f64 + microsecond as f64 / 1_000_000.0;
    let ts = TimeScales::from_utc(year, month, day, hour, minute, second_with_micro);
    Ok(gcrs_to_topocentric_compute(
        sat_x, sat_y, sat_z,
        station_lat_deg, station_lon_deg, station_alt_km,
        &ts, skyfield_compat,
    ))
}
