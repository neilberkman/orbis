//! Lambert problem solver using Battin's method.
//!
//! Algorithm 61, Vallado 2022, pp. 505-510.

use std::f64::consts::PI;

const MU: f64 = 398600.4415;
const TWOPI: f64 = 2.0 * PI;
const SMALL: f64 = 1e-10;

fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
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

/// Continued fraction for Battin's method (Vallado Eq. 7-69).
fn seebatt(v: f64) -> f64 {
    let c = [
        9.0 / 7.0,
        16.0 / 63.0,
        25.0 / 99.0,
        36.0 / 143.0,
        49.0 / 195.0,
        64.0 / 255.0,
        81.0 / 323.0,
        100.0 / 399.0,
        121.0 / 483.0,
        144.0 / 575.0,
        169.0 / 675.0,
        196.0 / 783.0,
        225.0 / 899.0,
        256.0 / 1023.0,
        289.0 / 1155.0,
        324.0 / 1295.0,
        361.0 / 1443.0,
        400.0 / 1599.0,
        441.0 / 1763.0,
        484.0 / 1935.0,
    ];

    let sqrtopv = (1.0 + v).sqrt();
    let eta = v / (1.0 + sqrtopv).powi(2);

    let mut term2 = 1.0 + c[19] * eta;
    for j in (0..19).rev() {
        term2 = 1.0 + c[j] * eta / term2;
    }

    8.0 * (1.0 + sqrtopv) / (3.0 + 1.0 / (5.0 + eta + (9.0 / 7.0) * eta / term2))
}

/// Continued fraction for Battin's method (Vallado Eq. 7-70).
fn kbatt(v: f64) -> f64 {
    let d = [
        1.0 / 3.0,
        4.0 / 27.0,
        8.0 / 27.0,
        2.0 / 9.0,
        22.0 / 81.0,
        208.0 / 891.0,
        340.0 / 1287.0,
        418.0 / 1755.0,
        598.0 / 2295.0,
        700.0 / 2907.0,
        928.0 / 3591.0,
        1054.0 / 4347.0,
        1330.0 / 5175.0,
        1480.0 / 6075.0,
        1804.0 / 7047.0,
        1978.0 / 8091.0,
        2350.0 / 9207.0,
        2548.0 / 10395.0,
        2968.0 / 11655.0,
        3190.0 / 12987.0,
        3658.0 / 14391.0,
    ];

    // Forward pass
    let mut sum1: f64 = d[0];
    let mut delold: f64 = 1.0;
    let mut termold: f64 = d[0];
    let ktr = 21;

    for di in d.iter().take(ktr).skip(1) {
        if termold.abs() <= 1e-8 {
            break;
        }
        let del = 1.0 / (1.0 + di * v * delold);
        let term = termold * (del - 1.0);
        sum1 += term;
        delold = del;
        termold = term;
    }
    let _ = sum1; // forward pass result not used in final — backward pass is

    // Backward pass
    let mut term2 = 1.0 + d[ktr - 1] * v;
    for i in 0..(ktr - 2) {
        let sum2 = d[ktr - i - 2] * v / term2;
        term2 = 1.0 + sum2;
    }

    d[0] / term2
}

/// Hodograph velocity transfer (Thompson 2013/2018).
fn hodograph(
    r1: &[f64; 3],
    r2: &[f64; 3],
    v1: &[f64; 3],
    p: f64,
    ecc: f64,
    dnu: f64,
    dtsec: f64,
) -> ([f64; 3], [f64; 3]) {
    let magr1 = mag(r1);
    let magr2 = mag(r2);

    let a = MU * (1.0 / magr1 - 1.0 / p);
    let b = (MU * ecc / p).powi(2) - a * a;
    let x1_abs = if b <= 0.0 { 0.0 } else { b.sqrt() };
    let mut x1 = -x1_abs;

    let nvec;
    if dnu.sin().abs() < SMALL {
        // 180-degree transfer
        let cp = cross(r1, v1);
        let ncp = mag(&cp);
        nvec = smul(1.0 / ncp, &cp);

        if ecc < 1.0 {
            let ptx = TWOPI * (p.powi(3) / (MU * (1.0 - ecc * ecc).powi(3))).sqrt();
            if dtsec % ptx > ptx * 0.5 {
                x1 = x1_abs;
            }
        }
    } else {
        // Common path
        let y2a = MU / p - x1 * dnu.sin() + a * dnu.cos();
        let y2b = MU / p + x1 * dnu.sin() + a * dnu.cos();
        if (MU / magr2 - y2b).abs() < (MU / magr2 - y2a).abs() {
            x1 = x1_abs;
        }

        let cp = cross(r1, r2);
        let ncp = mag(&cp);
        nvec = if dnu % TWOPI > PI {
            smul(-1.0 / ncp, &cp)
        } else {
            smul(1.0 / ncp, &cp)
        };
    }

    let sqrtmup = (MU * p).sqrt();
    let nr1 = cross(&nvec, r1);
    let v1t = smul(
        sqrtmup / magr1,
        &vadd(&smul(x1 / MU, r1), &smul(1.0 / magr1, &nr1)),
    );

    let x2 = x1 * dnu.cos() + a * dnu.sin();
    let nr2 = cross(&nvec, r2);
    let v2t = smul(
        sqrtmup / magr2,
        &vadd(&smul(x2 / MU, r2), &smul(1.0 / magr2, &nr2)),
    );

    (v1t, v2t)
}

/// Direction of motion for Lambert problem.
#[derive(Clone, Copy, PartialEq)]
pub enum DirectionOfMotion {
    Short,
    Long,
}

/// Direction of energy for Lambert problem.
#[derive(Clone, Copy, PartialEq)]
pub enum DirectionOfEnergy {
    Low,
    High,
}

/// Solve Lambert's problem using Battin's method.
///
/// Given two position vectors and time of flight, find the transfer velocities.
///
/// Algorithm 61, Vallado 2022, pp. 505-510.
pub fn battin(
    r1: &[f64; 3],
    r2: &[f64; 3],
    v1: &[f64; 3],
    dm: DirectionOfMotion,
    de: DirectionOfEnergy,
    nrev: i32,
    dtsec: f64,
) -> Result<([f64; 3], [f64; 3]), String> {
    let magr1 = mag(r1);
    let magr2 = mag(r2);

    // Angle between position vectors
    let cosdeltanu = dot(r1, r2) / (magr1 * magr2);
    let magrcrossr = mag(&cross(r1, r2));
    let sign = if dm == DirectionOfMotion::Short {
        1.0
    } else {
        -1.0
    };
    let sindeltanu = sign * magrcrossr / (magr1 * magr2);
    let mut dnu = sindeltanu.atan2(cosdeltanu);
    if dnu < 0.0 {
        dnu += TWOPI;
    }

    // Chord and semiperimeter
    let chord = (magr1 * magr1 + magr2 * magr2 - 2.0 * magr1 * magr2 * cosdeltanu).sqrt();
    let s = (magr1 + magr2 + chord) * 0.5;
    let eps = magr2 / magr1 - 1.0;

    // Lambda, L, m
    let lam = (magr1 * magr2).sqrt() / s * (dnu * 0.5).cos();
    let l_ = ((1.0 - lam) / (1.0 + lam)).powi(2);
    let m = 8.0 * MU * dtsec * dtsec / (s.powi(3) * (1.0 + lam).powi(6));

    // Initial guess
    let mut xn = if nrev > 0 { 1.0 + 4.0 * l_ } else { l_ };

    if de == DirectionOfEnergy::High && nrev > 0 {
        // High energy multi-rev case
        xn = 1e-20;
        let mut x = 10.0;
        for _ in 0..20 {
            if (xn - x).abs() < SMALL {
                break;
            }
            x = xn;
            let temp = 1.0 / (2.0 * (l_ - x * x));
            let temp1 = x.sqrt();
            let temp2 =
                (nrev as f64 * PI * 0.5 + temp1.atan()) / temp1;
            let h1 = temp * (l_ + x) * (1.0 + 2.0 * x + l_);
            let h2 = temp * m * temp1 * ((l_ - x * x) * temp2 - (l_ + x));

            let b = 0.25 * 27.0 * h2 / (temp1 * (1.0 + h1)).powi(3);
            let f = if b < 0.0 {
                2.0 * ((b + 1.0).sqrt().acos() / 3.0).cos()
            } else {
                let a_ = (b.sqrt() + (b + 1.0).sqrt()).powf(1.0 / 3.0);
                a_ + 1.0 / a_
            };

            let y = 2.0 / 3.0 * temp1 * (1.0 + h1) * ((b + 1.0).sqrt() / f + 1.0);
            xn = 0.5
                * ((m / (y * y) - (1.0 + l_))
                    - ((m / (y * y) - (1.0 + l_)).powi(2) - 4.0 * l_).sqrt());
        }

        let a_orbit = s * (1.0 + lam).powi(2) * (1.0 + xn) * (l_ + xn) / (8.0 * xn);
        let p =
            2.0 * magr1 * magr2 * (1.0 + xn) * (dnu * 0.5).sin().powi(2)
                / (s * (1.0 + lam).powi(2) * (l_ + xn));
        let ecc = (1.0 - p / a_orbit).abs().sqrt();
        Ok(hodograph(r1, r2, v1, p, ecc, dnu, dtsec))
    } else {
        // Standard / low energy case
        let mut x = 10.0;
        let max_loops = 30;
        let mut loops = 0;
        let mut y = 0.0;

        while (xn - x).abs() >= SMALL && loops < max_loops {
            x = xn;
            loops += 1;

            let (h1, h2) = if nrev > 0 {
                let temp = 1.0 / ((1.0 + 2.0 * x + l_) * (4.0 * x * x));
                let temp1 =
                    (nrev as f64 * PI * 0.5 + x.sqrt().atan()) / x.sqrt();
                let h1 = temp * (l_ + x).powi(2) * (3.0 * (1.0 + x).powi(2) * temp1 - (3.0 + 5.0 * x));
                let h2 = temp * m * ((x * x - x * (1.0 + l_) - 3.0 * l_) * temp1 + (3.0 * l_ + x));
                (h1, h2)
            } else {
                let tempx = seebatt(x);
                let denom = 1.0 / ((1.0 + 2.0 * x + l_) * (4.0 * x + tempx * (3.0 + x)));
                let h1 = (l_ + x).powi(2) * (1.0 + 3.0 * x + tempx) * denom;
                let h2 = m * (x - l_ + tempx) * denom;
                (h1, h2)
            };

            let b = 0.25 * 27.0 * h2 / (1.0 + h1).powi(3);
            let u = 0.5 * b / (1.0 + (1.0 + b).sqrt());
            let k2 = kbatt(u);
            y = (1.0 + h1) / 3.0 * (2.0 + (1.0 + b).sqrt() / (1.0 + 2.0 * u * k2 * k2));
            xn = (((1.0 - l_) * 0.5).powi(2) + m / (y * y)).sqrt() - (1.0 + l_) * 0.5;
        }

        if loops >= max_loops {
            return Err("Lambert did not converge".to_string());
        }

        let p = 2.0 * magr1 * magr2 * y * y * (1.0 + x).powi(2) * (dnu * 0.5).sin().powi(2)
            / (m * s * (1.0 + lam).powi(2));
        let ecc = (eps * eps
            + 4.0 * magr2 / magr1 * (dnu * 0.5).sin().powi(2) * ((l_ - x) / (l_ + x)).powi(2))
            / (eps * eps + 4.0 * magr2 / magr1 * (dnu * 0.5).sin().powi(2));
        let ecc = ecc.sqrt();

        Ok(hodograph(r1, r2, v1, p, ecc, dnu, dtsec))
    }
}

use rustler::NifResult;

type Vec3 = (f64, f64, f64);

pub(crate) fn lambert_battin_impl(
    r1: Vec3,
    r2: Vec3,
    v1: Vec3,
    dm: i32,
    de: i32,
    nrev: i32,
    dtsec: f64,
) -> NifResult<(Vec3, Vec3)> {
    let r1a = [r1.0, r1.1, r1.2];
    let r2a = [r2.0, r2.1, r2.2];
    let v1a = [v1.0, v1.1, v1.2];

    let dm = if dm == 0 {
        DirectionOfMotion::Short
    } else {
        DirectionOfMotion::Long
    };
    let de = if de == 0 {
        DirectionOfEnergy::Low
    } else {
        DirectionOfEnergy::High
    };

    match battin(&r1a, &r2a, &v1a, dm, de, nrev, dtsec) {
        Ok((v1t, v2t)) => Ok(((v1t[0], v1t[1], v1t[2]), (v2t[0], v2t[1], v2t[2]))),
        Err(_) => Ok(((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))),
    }
}
