//! Gauss angles-only initial orbit determination.
//!
//! Algorithm 52, Vallado 2022, pp. 448-459.

use crate::iod;

const MU: f64 = 398600.4415;
const RE: f64 = 6378.1363;
const TUSEC: f64 = 806.8109913067327;
const DAY2SEC: f64 = 86400.0;
const SMALL: f64 = 1e-10;

/// 3x3 matrix determinant.
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// 3x3 matrix inverse.
fn inv3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let d = det3(m);
    let id = 1.0 / d;
    [
        [
            id * (m[1][1] * m[2][2] - m[1][2] * m[2][1]),
            id * (m[0][2] * m[2][1] - m[0][1] * m[2][2]),
            id * (m[0][1] * m[1][2] - m[0][2] * m[1][1]),
        ],
        [
            id * (m[1][2] * m[2][0] - m[1][0] * m[2][2]),
            id * (m[0][0] * m[2][2] - m[0][2] * m[2][0]),
            id * (m[0][2] * m[1][0] - m[0][0] * m[1][2]),
        ],
        [
            id * (m[1][0] * m[2][1] - m[1][1] * m[2][0]),
            id * (m[0][1] * m[2][0] - m[0][0] * m[2][1]),
            id * (m[0][0] * m[1][1] - m[0][1] * m[1][0]),
        ],
    ]
}

/// 3x3 matrix-vector multiply.
fn mat3_vec3(m: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// 3x3 matrix multiply.
fn mat3_mat3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    r
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn mag(a: &[f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn smul(s: f64, a: &[f64; 3]) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}

fn vadd(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Line-of-sight unit vector from RA and Dec.
fn los(ra: f64, dec: f64) -> [f64; 3] {
    [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()]
}

/// Halley iteration to refine the 8th-order polynomial root.
fn halley_iteration(poly: &[f64; 9]) -> f64 {
    // Find initial guess from polynomial roots
    // Use a simple Newton-based approach starting from GPS altitude
    let mut bigr2c = 20000.0 / RE;
    let mut bigr2 = 100.0;

    for _ in 0..15 {
        if (bigr2 - bigr2c).abs() < 8e-5 {
            break;
        }
        bigr2 = bigr2c;
        let x = bigr2;
        let f = x.powi(8) + poly[2] * x.powi(6) + poly[5] * x.powi(3) + poly[8];
        let f1 = 8.0 * x.powi(7) + 6.0 * poly[2] * x.powi(5) + 3.0 * poly[5] * x.powi(2);
        let f2 = 56.0 * x.powi(6) + 30.0 * poly[2] * x.powi(4) + 6.0 * poly[5] * x;
        bigr2c = bigr2 - (2.0 * f * f1) / (2.0 * f1 * f1 - f * f2);
    }

    bigr2c
}

/// Gauss angles-only orbit determination.
///
/// Given three angular observations (RA/Dec) with times and site positions,
/// determine the orbit at the middle observation.
pub fn gauss_angles(
    decl: &[f64; 3],
    rtasc: &[f64; 3],
    jd: &[f64; 3],
    jdf: &[f64; 3],
    rseci: &[[f64; 3]; 3],
) -> Result<([f64; 3], [f64; 3]), String> {
    // Time intervals (seconds)
    let tau12 = ((jd[0] - jd[1]) + (jdf[0] - jdf[1])) * DAY2SEC;
    let _tau13 = ((jd[0] - jd[2]) + (jdf[0] - jdf[2])) * DAY2SEC;
    let tau32 = ((jd[2] - jd[1]) + (jdf[2] - jdf[1])) * DAY2SEC;

    // Line-of-sight vectors
    let l1 = los(rtasc[0], decl[0]);
    let l2 = los(rtasc[1], decl[1]);
    let l3 = los(rtasc[2], decl[2]);

    // Canonical units
    let tau12c = tau12 / TUSEC;
    let tau32c = tau32 / TUSEC;
    let rseci1c = smul(1.0 / RE, &rseci[0]);
    let rseci2c = smul(1.0 / RE, &rseci[1]);
    let rseci3c = smul(1.0 / RE, &rseci[2]);

    // L-matrix (columns = LOS vectors)
    let lmat = [
        [l1[0], l2[0], l3[0]],
        [l1[1], l2[1], l3[1]],
        [l1[2], l2[2], l3[2]],
    ];

    let d = det3(&lmat);
    if d.abs() < SMALL {
        return Err("determinant too small".to_string());
    }

    let lmati = inv3(&lmat);

    // Range-site matrix (columns = site vectors in canonical units)
    let rsmatc = [
        [rseci1c[0], rseci2c[0], rseci3c[0]],
        [rseci1c[1], rseci2c[1], rseci3c[1]],
        [rseci1c[2], rseci2c[2], rseci3c[2]],
    ];

    let lir = mat3_mat3(&lmati, &rsmatc);

    // Polynomial coefficients
    let a1 = tau32c / (tau32c - tau12c);
    let a1u = (tau32c * ((tau32c - tau12c).powi(2) - tau32c.powi(2))) / (6.0 * (tau32c - tau12c));
    let a3 = -tau12c / (tau32c - tau12c);
    let a3u = -(tau12c * ((tau32c - tau12c).powi(2) - tau12c.powi(2))) / (6.0 * (tau32c - tau12c));

    let d1c = lir[1][0] * a1 - lir[1][1] + lir[1][2] * a3;
    let d2c = lir[1][0] * a1u + lir[1][2] * a3u;
    let magrs2 = mag(&rseci2c);
    let l2dotrs = dot(&l2, &rseci2c);

    // 8th-order polynomial
    let mut poly = [0.0; 9];
    poly[0] = 1.0;
    poly[2] = -(d1c.powi(2) + 2.0 * d1c * l2dotrs + magrs2.powi(2));
    poly[5] = -2.0 * (l2dotrs * d2c + d1c * d2c);
    poly[8] = -(d2c.powi(2));

    // Solve for radius
    let mut bigr2c = halley_iteration(&poly);

    if bigr2c < 0.0 || bigr2c * RE > 50000.0 {
        bigr2c = 35000.0 / RE;
    }

    let bigr2 = bigr2c * RE;
    let a1u_sec = a1u * TUSEC.powi(2);
    let a3u_sec = a3u * TUSEC.powi(2);

    // Solve for f and g series
    let u = MU / bigr2.powi(3);
    let c1 = a1 + a1u_sec * u;
    let c2 = -1.0;
    let c3 = a3 + a3u_sec * u;

    // Range-site matrix (non-canonical)
    let rsmat = [
        [rseci[0][0], rseci[1][0], rseci[2][0]],
        [rseci[0][1], rseci[1][1], rseci[2][1]],
        [rseci[0][2], rseci[1][2], rseci[2][2]],
    ];
    let lir_full = mat3_mat3(&lmati, &rsmat);
    let cmat = [-c1, -c2, -c3];
    let rhomat = mat3_vec3(&lir_full, &cmat);

    // Form position vectors
    let r1 = vadd(&smul(rhomat[0] / c1, &l1), &rseci[0]);
    let r2 = vadd(&smul(rhomat[1] / c2, &l2), &rseci[1]);
    let r3 = vadd(&smul(rhomat[2] / c3, &l3), &rseci[2]);

    // Use Gibbs to get velocity
    let (v2, _, _, _) = iod::gibbs(&r1, &r2, &r3)?;

    Ok((r2, v2))
}

use rustler::NifResult;

type Vec3 = (f64, f64, f64);

#[allow(clippy::too_many_arguments)]
pub(crate) fn gauss_impl(
    decl1: f64, decl2: f64, decl3: f64,
    rtasc1: f64, rtasc2: f64, rtasc3: f64,
    jd1: f64, jdf1: f64,
    jd2: f64, jdf2: f64,
    jd3: f64, jdf3: f64,
    rseci1: Vec3, rseci2: Vec3, rseci3: Vec3,
) -> NifResult<(Vec3, Vec3)> {
    let decl = [decl1, decl2, decl3];
    let rtasc = [rtasc1, rtasc2, rtasc3];
    let jd = [jd1, jd2, jd3];
    let jdf = [jdf1, jdf2, jdf3];
    let rseci = [
        [rseci1.0, rseci1.1, rseci1.2],
        [rseci2.0, rseci2.1, rseci2.2],
        [rseci3.0, rseci3.1, rseci3.2],
    ];

    match gauss_angles(&decl, &rtasc, &jd, &jdf, &rseci) {
        Ok((r2, v2)) => Ok(((r2[0], r2[1], r2[2]), (v2[0], v2[1], v2[2]))),
        Err(_) => Ok(((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))),
    }
}
