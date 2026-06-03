//! Rustler boundary for the `astrodynamics-gnss` single-point-positioning (SPP)
//! least-squares PVT solve.
//!
//! This module is **pure glue**: it decodes Erlang terms into the crate's
//! [`SolveInputs`], calls [`solve_spp`], and encodes the [`ReceiverSolution`]
//! back. No transmit-time iteration, no least-squares numerics, no atmospheric
//! model, and no frame conversion lives here — those are the crate's
//! responsibility. The SP3 product is reused from the [`Sp3Resource`] handle the
//! `sp3_parse/1` NIF already returns; this call never touches the filesystem.
//!
//! Boundary units: pseudoranges and the initial guess are meters, epoch scalars
//! are seconds (and a fractional day-of-year), pressure is hPa, temperature is
//! kelvin, relative humidity is a `[0, 1]` fraction. The returned position is
//! ITRF/IGS ECEF meters and the geodetic latitude/longitude are radians, exactly
//! as the crate produces them.

use astrodynamics_gnss::{
    solve_spp, Corrections, GnssSatelliteId, GnssSystem, KlobucharCoeffs, Observation,
    ReceiverSolution, RejectionReason, SolveInputs, SppError, SurfaceMet,
};
use rustler::types::atom;
use rustler::types::tuple::make_tuple;
use rustler::{Encoder, Env, Error, NifResult, ResourceArc, Term};

use crate::sp3::Sp3Resource;

/// Map a GNSS single-letter system identifier (e.g. `"G"`) onto the crate's
/// [`GnssSystem`]. Pure identifier translation; mirrors `sp3::system_from_letter`.
fn system_from_letter(letter: &str) -> NifResult<GnssSystem> {
    let c = letter
        .chars()
        .next()
        .ok_or_else(|| Error::Term(Box::new("empty GNSS system letter")))?;
    GnssSystem::from_letter(c)
        .ok_or_else(|| Error::Term(Box::new(format!("unknown GNSS system letter {letter:?}"))))
}

/// Translate an [`SppError`] into the atom the Elixir wrapper maps to a public
/// `{:error, reason}`. The offending satellite (when the variant carries one) is
/// rendered with the crate's canonical `Display` token (e.g. `"G01"`) so the
/// reason term stays informative without leaking crate internals.
fn spp_error_term<'a>(env: Env<'a>, e: &SppError) -> Term<'a> {
    match e {
        SppError::TooFewSatellites { used } => {
            (atom::error(), atom_from(env, "too_few_satellites"), *used as i64).encode(env)
        }
        SppError::Singular(_) => (atom::error(), atom_from(env, "singular_geometry")).encode(env),
        SppError::DuplicateObservation { satellite } => (
            atom::error(),
            atom_from(env, "duplicate_observation"),
            satellite.to_string(),
        )
            .encode(env),
        SppError::EphemerisLost { satellite } => (
            atom::error(),
            atom_from(env, "ephemeris_lost"),
            satellite.to_string(),
        )
            .encode(env),
    }
}

/// Intern a runtime atom. Glue helper so error reasons and rejection reasons are
/// encoded as atoms (idiomatic on the Elixir side) rather than strings.
fn atom_from<'a>(env: Env<'a>, name: &str) -> Term<'a> {
    atom::Atom::from_str(env, name)
        .map(|a| a.encode(env))
        .unwrap_or_else(|_| name.encode(env))
}

/// Encode the converged [`ReceiverSolution`] as the `{:ok, solution}` term the
/// Elixir wrapper destructures. The solution body is a fixed-arity tuple:
///
/// ```text
/// {{x_m, y_m, z_m},                      # ITRF/IGS ECEF position, meters
///  rx_clock_s,                           # receiver clock bias, seconds
///  {lat_rad, lon_rad, height_m} | nil,   # geodetic, when requested
///  {gdop, pdop, hdop, vdop, tdop} | nil, # DOP, when the geometry is full rank
///  [residual_m, ...],                    # post-fit residuals, used_sats order
///  ["G01", ...],                         # used satellites
///  [{"G07", :low_elevation}, ...],       # rejected satellites + reason atom
///  {iterations, converged, ionosphere_applied, troposphere_applied}}
/// ```
fn encode_solution<'a>(env: Env<'a>, sol: &ReceiverSolution) -> Term<'a> {
    let pos = sol.position.as_array();
    let position = (pos[0], pos[1], pos[2]);

    let geodetic: Term<'a> = match sol.geodetic {
        Some(g) => (g.lat_rad, g.lon_rad, g.height_m).encode(env),
        None => atom::nil().encode(env),
    };

    let dop: Term<'a> = match sol.dop {
        Some(d) => (d.gdop, d.pdop, d.hdop, d.vdop, d.tdop).encode(env),
        None => atom::nil().encode(env),
    };

    let used_sats: Vec<String> = sol.used_sats.iter().map(|s| s.to_string()).collect();

    let rejected_sats: Vec<(String, Term<'a>)> = sol
        .rejected_sats
        .iter()
        .map(|r| {
            let reason = match r.reason {
                RejectionReason::NoEphemeris => atom_from(env, "no_ephemeris"),
                RejectionReason::LowElevation => atom_from(env, "low_elevation"),
            };
            (r.satellite_id.to_string(), reason)
        })
        .collect();

    let metadata = (
        sol.metadata.iterations as i64,
        sol.metadata.converged,
        sol.metadata.ionosphere_applied,
        sol.metadata.troposphere_applied,
    );

    // The body has eight fields, past the arity of the blanket tuple `Encoder`,
    // so it is assembled with `make_tuple` over the already-encoded terms.
    let body = make_tuple(
        env,
        &[
            position.encode(env),
            sol.rx_clock_s.encode(env),
            geodetic,
            dop,
            sol.residuals_m.encode(env),
            used_sats.encode(env),
            rejected_sats.encode(env),
            metadata.encode(env),
        ],
    );

    (atom::ok(), body).encode(env)
}

/// Solve single-point positioning for one receive epoch against a loaded SP3
/// handle.
///
/// Dirty-CPU: the transmit-time iteration and trust-region least-squares solve
/// are unbounded relative to the 1 ms NIF budget. `observations` is a list of
/// `{sat_token, pseudorange_m}` pairs where `sat_token` is the canonical
/// SP3/RINEX id string (e.g. `"G01"`); the system letter and PRN are parsed via
/// [`GnssSystem::from_letter`]. The three epoch scalars, the four-element initial
/// guess `[x_m, y_m, z_m, b_m]`, the correction toggles, the Klobuchar
/// alpha/beta coefficient tuples, and the surface meteorology are forwarded
/// verbatim into [`SolveInputs`]; no domain math happens here.
///
/// Returns `{:ok, solution}` (see [`encode_solution`]) or `{:error, reason}`
/// where `reason` is the mapped [`SppError`] atom.
#[rustler::nif(schedule = "DirtyCpu")]
#[allow(clippy::too_many_arguments)]
fn spp_solve<'a>(
    env: Env<'a>,
    handle: ResourceArc<Sp3Resource>,
    observations: Vec<(String, f64)>,
    t_rx_j2000_s: f64,
    t_rx_second_of_day_s: f64,
    day_of_year: f64,
    initial_guess: (f64, f64, f64, f64),
    apply_iono: bool,
    apply_tropo: bool,
    alpha: (f64, f64, f64, f64),
    beta: (f64, f64, f64, f64),
    pressure_hpa: f64,
    temperature_k: f64,
    relative_humidity: f64,
    with_geodetic: bool,
) -> NifResult<Term<'a>> {
    let mut obs = Vec::with_capacity(observations.len());
    for (token, pseudorange_m) in &observations {
        let (letter, rest) = token.split_at(token.char_indices().nth(1).map_or(0, |(i, _)| i));
        let system = system_from_letter(letter)?;
        let prn: u8 = rest
            .parse()
            .map_err(|_| Error::Term(Box::new(format!("bad satellite token {token:?}"))))?;
        obs.push(Observation {
            satellite_id: GnssSatelliteId::new(system, prn),
            pseudorange_m: *pseudorange_m,
        });
    }

    let inputs = SolveInputs {
        observations: obs,
        t_rx_j2000_s,
        t_rx_second_of_day_s,
        day_of_year,
        initial_guess: [
            initial_guess.0,
            initial_guess.1,
            initial_guess.2,
            initial_guess.3,
        ],
        corrections: Corrections {
            ionosphere: apply_iono,
            troposphere: apply_tropo,
        },
        klobuchar: KlobucharCoeffs {
            alpha: [alpha.0, alpha.1, alpha.2, alpha.3],
            beta: [beta.0, beta.1, beta.2, beta.3],
        },
        met: SurfaceMet {
            pressure_hpa,
            temperature_k,
            relative_humidity,
        },
    };

    match solve_spp(&handle.sp3, &inputs, with_geodetic) {
        Ok(sol) => Ok(encode_solution(env, &sol)),
        Err(e) => Ok(spp_error_term(env, &e)),
    }
}
