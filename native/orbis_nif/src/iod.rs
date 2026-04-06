//! Initial Orbit Determination (IOD) methods.
//!
//! Implements Gibbs and Herrick-Gibbs methods from Vallado's
//! "Fundamentals of Astrodynamics and Applications" (2022).

use rustler::NifResult;

/// Earth gravitational parameter (km^3/s^2), matching Vallado.
const MU: f64 = 398600.4415;

/// Cross product of two 3-vectors.
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product of two 3-vectors.
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Magnitude of a 3-vector.
fn mag(a: &[f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Unit vector.
fn unit(a: &[f64; 3]) -> [f64; 3] {
    let m = mag(a);
    [a[0] / m, a[1] / m, a[2] / m]
}

/// Vector addition: a + b.
fn vadd(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scalar-vector multiply.
fn smul(s: f64, a: &[f64; 3]) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}

/// Gibbs method: determine velocity at r2 from three coplanar position vectors.
///
/// Algorithm 54, Vallado 2022, pp. 460-467.
///
/// Returns (v2, theta12_rad, theta23_rad, copa_rad) or an error string.
pub fn gibbs(
    r1: &[f64; 3],
    r2: &[f64; 3],
    r3: &[f64; 3],
) -> Result<([f64; 3], f64, f64, f64), String> {
    let magr1 = mag(r1);
    let magr2 = mag(r2);
    let magr3 = mag(r3);

    // Cross products
    let p = cross(r2, r3);
    let q = cross(r3, r1);
    let w = cross(r1, r2);

    // Coplanarity angle
    let copa = dot(&unit(&p), &unit(r1)).asin();

    // D = P + Q + W
    let d = vadd(&vadd(&p, &q), &w);
    let magd = mag(&d);

    // N = |r1|*P + |r2|*Q + |r3|*W
    let n = vadd(&vadd(&smul(magr1, &p), &smul(magr2, &q)), &smul(magr3, &w));
    let magn = mag(&n);

    if magd < 1e-6 || magn < 1e-6 {
        return Err("orbit determination not possible".to_string());
    }

    // Angles between position vectors
    let theta12 = (dot(r1, r2) / (magr1 * magr2)).clamp(-1.0, 1.0).acos();
    let theta23 = (dot(r2, r3) / (magr2 * magr3)).clamp(-1.0, 1.0).acos();

    // S vector
    let r1mr2 = magr1 - magr2;
    let r3mr1 = magr3 - magr1;
    let r2mr3 = magr2 - magr3;
    let s = vadd(
        &vadd(&smul(r1mr2, r3), &smul(r3mr1, r2)),
        &smul(r2mr3, r1),
    );

    // B = D × r2
    let b = cross(&d, r2);

    // Scaling factor
    let lg = (MU / (magd * magn)).sqrt();

    // v2 = (lg / |r2|) * B + lg * S
    let v2 = vadd(&smul(lg / magr2, &b), &smul(lg, &s));

    Ok((v2, theta12, theta23, copa))
}

/// Herrick-Gibbs method: determine velocity at r2 from three closely-spaced
/// position vectors with timestamps.
///
/// Algorithm 55, Vallado 2022, pp. 467-472.
///
/// Times are in Julian days (or any consistent unit — only differences matter).
///
/// Returns (v2, theta12_rad, theta23_rad, copa_rad) or an error string.
pub fn hgibbs(
    r1: &[f64; 3],
    r2: &[f64; 3],
    r3: &[f64; 3],
    jd1: f64,
    jd2: f64,
    jd3: f64,
) -> Result<([f64; 3], f64, f64, f64), String> {
    let magr1 = mag(r1);
    let magr2 = mag(r2);
    let magr3 = mag(r3);

    // Time differences (in seconds if jd is in days)
    let dt21 = (jd2 - jd1) * 86400.0;
    let dt31 = (jd3 - jd1) * 86400.0;
    let dt32 = (jd3 - jd2) * 86400.0;

    // Cross products for coplanarity check
    let p = cross(r2, r3);

    // Coplanarity angle
    let copa = dot(&unit(&p), &unit(r1)).asin();

    // Angles between position vectors
    let theta12 = (dot(r1, r2) / (magr1 * magr2)).clamp(-1.0, 1.0).acos();
    let theta23 = (dot(r2, r3) / (magr2 * magr3)).clamp(-1.0, 1.0).acos();

    // Herrick-Gibbs velocity approximation
    let term1 = smul(
        -dt32 * (1.0 / (dt21 * dt31) + MU / (12.0 * magr1.powi(3))),
        r1,
    );
    let term2 = smul(
        (dt32 - dt21) * (1.0 / (dt21 * dt32) + MU / (12.0 * magr2.powi(3))),
        r2,
    );
    let term3 = smul(
        dt21 * (1.0 / (dt32 * dt31) + MU / (12.0 * magr3.powi(3))),
        r3,
    );

    let v2 = vadd(&vadd(&term1, &term2), &term3);

    Ok((v2, theta12, theta23, copa))
}

type Vec3 = (f64, f64, f64);

pub(crate) fn gibbs_impl(
    r1: Vec3,
    r2: Vec3,
    r3: Vec3,
) -> NifResult<(Vec3, f64, f64, f64)> {
    let r1a = [r1.0, r1.1, r1.2];
    let r2a = [r2.0, r2.1, r2.2];
    let r3a = [r3.0, r3.1, r3.2];

    match gibbs(&r1a, &r2a, &r3a) {
        Ok((v2, theta12, theta23, copa)) => {
            Ok(((v2[0], v2[1], v2[2]), theta12, theta23, copa))
        }
        Err(_) => Ok(((0.0, 0.0, 0.0), 0.0, 0.0, 0.0)),
    }
}

pub(crate) fn hgibbs_impl(
    r1: Vec3,
    r2: Vec3,
    r3: Vec3,
    jd1: f64,
    jd2: f64,
    jd3: f64,
) -> NifResult<(Vec3, f64, f64, f64)> {
    let r1a = [r1.0, r1.1, r1.2];
    let r2a = [r2.0, r2.1, r2.2];
    let r3a = [r3.0, r3.1, r3.2];

    match hgibbs(&r1a, &r2a, &r3a, jd1, jd2, jd3) {
        Ok((v2, theta12, theta23, copa)) => {
            Ok(((v2[0], v2[1], v2[2]), theta12, theta23, copa))
        }
        Err(_) => Ok(((0.0, 0.0, 0.0), 0.0, 0.0, 0.0)),
    }
}
