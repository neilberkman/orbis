//! Conjunction assessment: find closest approach between two satellites.
//!
//! Uses coarse-fine search: scan at configurable step size, then refine
//! with golden section search within each candidate interval.

use rustler::{Encoder, Env, NifResult, Term};

/// Distance between two 3-vectors.
fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Propagate and return position, given precomputed Constants.
fn propagate(c: &sgp4::Constants, tsince: f64) -> Option<[f64; 3]> {
    c.propagate_afspc_compatibility_mode(sgp4::MinutesSinceEpoch(tsince))
        .ok()
        .map(|p| p.position)
}

/// Golden section search for minimum distance within [a, b] (minutes).
fn golden_search(
    c1: &sgp4::Constants,
    c2: &sgp4::Constants,
    epoch_offset2: f64,
    mut a: f64,
    mut b: f64,
) -> Option<(f64, f64)> {
    let gr = (5.0_f64.sqrt() + 1.0) / 2.0;
    let tol = 1.0 / 60.0; // 1 second

    let mut c = b - (b - a) / gr;
    let mut d = a + (b - a) / gr;

    for _ in 0..50 {
        if (b - a).abs() < tol {
            break;
        }

        let dc = dist_at(c1, c2, epoch_offset2, c)?;
        let dd = dist_at(c1, c2, epoch_offset2, d)?;

        if dc < dd {
            b = d;
        } else {
            a = c;
        }

        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    }

    let mid = (a + b) / 2.0;
    let d = dist_at(c1, c2, epoch_offset2, mid)?;
    Some((mid, d))
}

fn dist_at(
    c1: &sgp4::Constants,
    c2: &sgp4::Constants,
    epoch_offset2: f64,
    tsince: f64,
) -> Option<f64> {
    let p1 = propagate(c1, tsince)?;
    let p2 = propagate(c2, tsince - epoch_offset2)?;
    Some(dist(&p1, &p2))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn conjunction_impl<'a>(
    env: Env<'a>,
    line1_a: &str,
    line2_a: &str,
    line1_b: &str,
    line2_b: &str,
    start_min: f64,
    end_min: f64,
    step_min: f64,
    threshold_km: f64,
) -> NifResult<Term<'a>> {
    let ok = rustler::types::atom::Atom::from_str(env, "ok")?;
    let error = rustler::types::atom::Atom::from_str(env, "error")?;

    let tle1 = format!("{}\n{}\n", line1_a.trim(), line2_a.trim());
    let tle2 = format!("{}\n{}\n", line1_b.trim(), line2_b.trim());

    // Parse and initialize once
    let e1 = match sgp4::parse_2les(&tle1) {
        Ok(e) if !e.is_empty() => e,
        _ => return Ok((error, "failed to parse TLE 1").encode(env)),
    };
    let e2 = match sgp4::parse_2les(&tle2) {
        Ok(e) if !e.is_empty() => e,
        _ => return Ok((error, "failed to parse TLE 2").encode(env)),
    };

    let c1 = match sgp4::Constants::from_elements_afspc_compatibility_mode(&e1[0]) {
        Ok(c) => c,
        Err(e) => return Ok((error, format!("SGP4 init TLE1: {e}")).encode(env)),
    };
    let c2 = match sgp4::Constants::from_elements_afspc_compatibility_mode(&e2[0]) {
        Ok(c) => c,
        Err(e) => return Ok((error, format!("SGP4 init TLE2: {e}")).encode(env)),
    };

    // Epoch offset: TLE2 epoch - TLE1 epoch in minutes
    let epoch1 =
        sgp4::julian_years_since_j2000_afspc_compatibility_mode(&e1[0].datetime);
    let epoch2 =
        sgp4::julian_years_since_j2000_afspc_compatibility_mode(&e2[0].datetime);
    let epoch_offset2 = (epoch2 - epoch1) * 365.25 * 24.0 * 60.0;

    // Coarse scan + golden section refinement
    let n_steps = ((end_min - start_min) / step_min).ceil() as usize;
    let mut results: Vec<(f64, f64)> = Vec::new();
    let mut prev_dist = f64::MAX;
    let mut prev_t = start_min;
    let mut decreasing = false;

    for i in 0..=n_steps {
        let t = (start_min + i as f64 * step_min).min(end_min);

        let d = match dist_at(&c1, &c2, epoch_offset2, t) {
            Some(d) => d,
            None => {
                prev_dist = f64::MAX;
                decreasing = false;
                prev_t = t;
                continue;
            }
        };

        if d > prev_dist && decreasing {
            let search_start = (prev_t - step_min).max(start_min);
            if let Some((tca, tca_dist)) =
                golden_search(&c1, &c2, epoch_offset2, search_start, t)
            {
                if tca_dist < threshold_km {
                    results.push((tca, tca_dist));
                }
            }
        }

        decreasing = d < prev_dist;
        prev_dist = d;
        prev_t = t;
    }

    Ok((ok, results).encode(env))
}
