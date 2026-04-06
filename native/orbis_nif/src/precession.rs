//! IAU 2006 precession matrix and ICRS-to-J2000 frame bias, ported from the
//! C++ Skyfield-compatible implementation.
//!
//! All arithmetic uses plain operators (no `f64::mul_add`) so that
//! rounding matches CPython / Skyfield compiled without FMA contraction.

use crate::matrix::Mat3;

const T0: f64 = 2_451_545.0;
#[allow(clippy::excessive_precision)]
const ASEC2RAD: f64 = 4.848_136_811_095_359_935_899_141e-6;

// ---------------------------------------------------------------------------
// Precession matrix (IAU 2006)
// ---------------------------------------------------------------------------

/// Compute the 3x3 precession rotation matrix for the given TDB Julian date,
/// using the IAU 2006 Fukushima-Williams parameterisation.
pub(crate) fn compute_skyfield_precession_matrix(jd_tdb: f64) -> Mat3 {
    const EPS0_ARCSEC: f64 = 84381.406;

    let t = (jd_tdb - T0) / 36525.0;

    let psia = ((((-0.0000000951 * t + 0.000132851) * t - 0.00114045) * t - 1.0790069) * t
        + 5038.481507)
        * t;

    let omegaa = ((((0.0000003337 * t - 0.000000467) * t - 0.00772503) * t + 0.0512623) * t
        - 0.025754)
        * t
        + EPS0_ARCSEC;

    let chia = ((((-0.0000000560 * t + 0.000170663) * t - 0.00121197) * t - 2.3814292) * t
        + 10.556403)
        * t;

    let eps0 = EPS0_ARCSEC * ASEC2RAD;
    let psia_rad = psia * ASEC2RAD;
    let omegaa_rad = omegaa * ASEC2RAD;
    let chia_rad = chia * ASEC2RAD;

    let sa = eps0.sin();
    let ca = eps0.cos();
    let sb = (-psia_rad).sin();
    let cb = (-psia_rad).cos();
    let sc = (-omegaa_rad).sin();
    let cc = (-omegaa_rad).cos();
    let sd = chia_rad.sin();
    let cd = chia_rad.cos();

    [
        [
            cd * cb - sb * sd * cc,
            cd * sb * ca + sd * cc * cb * ca - sa * sd * sc,
            cd * sb * sa + sd * cc * cb * sa + ca * sd * sc,
        ],
        [
            -sd * cb - sb * cd * cc,
            -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc,
            -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc,
        ],
        [
            sb * sc,
            -sc * cb * ca - sa * cc,
            -sc * cb * sa + cc * ca,
        ],
    ]
}

// ---------------------------------------------------------------------------
// ICRS-to-J2000 frame bias matrix
// ---------------------------------------------------------------------------

/// Build the ICRS-to-J2000 frame bias rotation matrix.
///
/// This accounts for the small offset between the ICRS axes and the
/// mean J2000.0 dynamical frame, parameterised by the three bias angles
/// xi_0, eta_0, and da_0 from the IAU 2006 conventions.
pub(crate) fn build_icrs_to_j2000() -> Mat3 {
    let xi0 = -0.0166170 * ASEC2RAD;
    let eta0 = -0.0068192 * ASEC2RAD;
    let da0 = -0.01460 * ASEC2RAD;

    let yx = -da0;
    let zx = xi0;
    let xy = da0;
    let zy = eta0;
    let xz = -xi0;
    let yz = -eta0;

    [
        [1.0 - 0.5 * (yx * yx + zx * zx), xy, xz],
        [yx, 1.0 - 0.5 * (yx * yx + zy * zy), yz],
        [zx, zy, 1.0 - 0.5 * (zy * zy + zx * zx)],
    ]
}
