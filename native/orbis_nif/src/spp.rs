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
    solve_spp, Corrections, EphemerisSource, GnssSatelliteId, GnssSystem, KlobucharCoeffs,
    Observation, ReceiverSolution, RejectionReason, SolveInputs, SppError, SurfaceMet,
};

use crate::broadcast::BroadcastResource;
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

/// The Elixir-facing reason for a failed solve, as a pure value with no `Env`
/// dependency, so the `SppError` → public-reason mapping is unit-testable
/// without the BEAM runtime. The satellite-carrying variants render the offender
/// with the crate's canonical `Display` token (e.g. `"G01"`) so the reason stays
/// informative without leaking crate internals. [`spp_error_term`] is the thin
/// encoder that turns this into the actual `{:error, ...}` term.
#[derive(Debug, Clone, PartialEq, Eq)]
enum SppErrorReason {
    TooFewSatellites { used: i64 },
    SingularGeometry,
    DuplicateObservation { satellite: String },
    EphemerisLost { satellite: String },
}

impl SppErrorReason {
    /// The atom name the Elixir wrapper destructures as the error reason. These
    /// strings are the public contract (`Orbis.PointPositioning.solve/4`), so a
    /// rename here is a breaking change.
    fn atom_name(&self) -> &'static str {
        match self {
            SppErrorReason::TooFewSatellites { .. } => "too_few_satellites",
            SppErrorReason::SingularGeometry => "singular_geometry",
            SppErrorReason::DuplicateObservation { .. } => "duplicate_observation",
            SppErrorReason::EphemerisLost { .. } => "ephemeris_lost",
        }
    }
}

/// Map an [`SppError`] onto its pure [`SppErrorReason`]. Total over the enum, so
/// every variant — including the defensive `Singular` / `EphemerisLost` paths
/// that real SP3 inputs do not naturally reach — has a tested mapping.
fn spp_error_reason(e: &SppError) -> SppErrorReason {
    match e {
        SppError::TooFewSatellites { used } => {
            SppErrorReason::TooFewSatellites { used: *used as i64 }
        }
        SppError::Singular(_) => SppErrorReason::SingularGeometry,
        SppError::DuplicateObservation { satellite } => SppErrorReason::DuplicateObservation {
            satellite: satellite.to_string(),
        },
        SppError::EphemerisLost { satellite } => SppErrorReason::EphemerisLost {
            satellite: satellite.to_string(),
        },
    }
}

/// Translate an [`SppError`] into the `{:error, reason}` term the Elixir wrapper
/// maps to a public reason. A thin `Env`-bound wrapper over [`spp_error_reason`].
fn spp_error_term<'a>(env: Env<'a>, e: &SppError) -> Term<'a> {
    let reason = spp_error_reason(e);
    let tag = atom_from(env, reason.atom_name());
    match reason {
        SppErrorReason::TooFewSatellites { used } => (atom::error(), tag, used).encode(env),
        SppErrorReason::SingularGeometry => (atom::error(), tag).encode(env),
        SppErrorReason::DuplicateObservation { satellite } => {
            (atom::error(), tag, satellite).encode(env)
        }
        SppErrorReason::EphemerisLost { satellite } => (atom::error(), tag, satellite).encode(env),
    }
}

/// The atom name for a solver termination [`Status`]. Pure (no `Env`) so the
/// status → atom mapping is unit-testable; the strings are the public contract
/// surfaced as `Solution.metadata.status`.
///
/// [`Status`]: astrodynamics::math::least_squares::Status
fn status_atom_name(status: astrodynamics::math::least_squares::Status) -> &'static str {
    use astrodynamics::math::least_squares::Status;
    match status {
        Status::GradientTolerance => "gradient_tolerance",
        Status::CostTolerance => "cost_tolerance",
        Status::StepTolerance => "step_tolerance",
        Status::MaxEvaluations => "max_evaluations",
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
///  {iterations, converged, status, ionosphere_applied, troposphere_applied}}
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

    let status = atom_from(env, status_atom_name(sol.metadata.status));
    let metadata = make_tuple(
        env,
        &[
            (sol.metadata.iterations as i64).encode(env),
            sol.metadata.converged.encode(env),
            status,
            sol.metadata.ionosphere_applied.encode(env),
            sol.metadata.troposphere_applied.encode(env),
        ],
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
            metadata,
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
/// Decode the common SPP term arguments into a [`SolveInputs`]. Shared by the
/// SP3-backed and broadcast-backed entry points, which differ only in the
/// ephemeris source they pass to the solver.
#[allow(clippy::too_many_arguments)]
fn build_solve_inputs(
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
) -> NifResult<SolveInputs> {
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

    Ok(SolveInputs {
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
    })
}

/// Run the solve against any ephemeris source and encode the result term.
fn solve_to_term<'a>(
    env: Env<'a>,
    eph: &dyn EphemerisSource,
    inputs: &SolveInputs,
    with_geodetic: bool,
) -> Term<'a> {
    match solve_spp(eph, inputs, with_geodetic) {
        Ok(sol) => encode_solution(env, &sol),
        Err(e) => spp_error_term(env, &e),
    }
}

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
    let inputs = build_solve_inputs(
        observations,
        t_rx_j2000_s,
        t_rx_second_of_day_s,
        day_of_year,
        initial_guess,
        apply_iono,
        apply_tropo,
        alpha,
        beta,
        pressure_hpa,
        temperature_k,
        relative_humidity,
    )?;
    Ok(solve_to_term(env, &handle.sp3, &inputs, with_geodetic))
}

/// As [`spp_solve`] but against a parsed broadcast-navigation product
/// ([`BroadcastResource`]) instead of an SP3 precise product.
#[rustler::nif(schedule = "DirtyCpu")]
#[allow(clippy::too_many_arguments)]
fn spp_solve_broadcast<'a>(
    env: Env<'a>,
    handle: ResourceArc<BroadcastResource>,
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
    let inputs = build_solve_inputs(
        observations,
        t_rx_j2000_s,
        t_rx_second_of_day_s,
        day_of_year,
        initial_guess,
        apply_iono,
        apply_tropo,
        alpha,
        beta,
        pressure_hpa,
        temperature_k,
        relative_humidity,
    )?;
    Ok(solve_to_term(env, &handle.store, &inputs, with_geodetic))
}

#[cfg(test)]
mod mapping_tests {
    //! Mechanical coverage of the boundary mappings that term encoding wraps.
    //! These exercise every `SppError` variant and every solver `Status`,
    //! including the defensive `Singular` / `EphemerisLost` paths that a real
    //! SP3 product does not naturally reach, so the advertised public reasons
    //! stay correct without depending on a physics fixture to trigger them.
    use super::*;
    use astrodynamics::math::least_squares::{SolveError, Status};

    fn gps(prn: u8) -> GnssSatelliteId {
        GnssSatelliteId::new(GnssSystem::Gps, prn)
    }

    #[test]
    fn spp_error_reason_is_total_over_every_variant() {
        assert_eq!(
            spp_error_reason(&SppError::TooFewSatellites { used: 3 }),
            SppErrorReason::TooFewSatellites { used: 3 }
        );
        assert_eq!(
            spp_error_reason(&SppError::Singular(SolveError::SingularJacobian)),
            SppErrorReason::SingularGeometry
        );
        assert_eq!(
            spp_error_reason(&SppError::DuplicateObservation { satellite: gps(7) }),
            SppErrorReason::DuplicateObservation {
                satellite: "G07".to_string()
            }
        );
        assert_eq!(
            spp_error_reason(&SppError::EphemerisLost { satellite: gps(12) }),
            SppErrorReason::EphemerisLost {
                satellite: "G12".to_string()
            }
        );
    }

    #[test]
    fn error_reason_atom_names_are_the_documented_public_reasons() {
        assert_eq!(
            SppErrorReason::TooFewSatellites { used: 0 }.atom_name(),
            "too_few_satellites"
        );
        assert_eq!(
            SppErrorReason::SingularGeometry.atom_name(),
            "singular_geometry"
        );
        assert_eq!(
            SppErrorReason::DuplicateObservation {
                satellite: String::new()
            }
            .atom_name(),
            "duplicate_observation"
        );
        assert_eq!(
            SppErrorReason::EphemerisLost {
                satellite: String::new()
            }
            .atom_name(),
            "ephemeris_lost"
        );
    }

    #[test]
    fn status_atom_names_cover_every_status() {
        assert_eq!(
            status_atom_name(Status::GradientTolerance),
            "gradient_tolerance"
        );
        assert_eq!(status_atom_name(Status::CostTolerance), "cost_tolerance");
        assert_eq!(status_atom_name(Status::StepTolerance), "step_tolerance");
        assert_eq!(status_atom_name(Status::MaxEvaluations), "max_evaluations");
    }
}
