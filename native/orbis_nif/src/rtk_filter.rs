//! Rustler boundary for the sequential RTK filter kernel.
//!
//! This is intentionally a traceable primitive, not the public RTK API: Elixir
//! still owns normalization/reporting while the kernel migration is gated. Terms
//! are plain tuples/lists so parity tests can feed the exact epoch/state stream
//! into Rust without introducing a second Elixir struct layer.

use astrodynamics_gnss::rtk_filter::{
    update_epoch, Epoch, FilterState, MeasModel, SatMeas, SearchOpts, StochasticModel, UpdateError,
    UpdateOpts,
};
use rustler::{Encoder, Env, NifResult, Term};
use std::collections::BTreeMap;

type Vec3 = (f64, f64, f64);
type SatIdsTerm = (String, String);
type SatObsTerm = (f64, f64, f64, f64);
type SatPosTerm = (Vec3, Vec3, Vec3);
type SatTerm = (SatIdsTerm, SatObsTerm, SatPosTerm);
type EpochTerm = (SatTerm, Vec<SatTerm>);
type StateHeaderTerm = (u16, String, Vec<String>, f64);
type StateTerm = (
    StateHeaderTerm,
    Vec3,
    Vec<f64>,
    Vec<f64>,
    Vec<(String, i64)>,
    Vec<(String, f64)>,
);
type ModelTerm = (f64, f64, String, bool, bool);
type UpdateOptsTerm = (f64, f64, f64, usize, f64);

mod atoms {
    rustler::atoms! {
        ok,
        error,
        invalid_stochastic_model,
        reference_changed,
        missing_ambiguity_column,
        missing_wavelength,
        missing_offset,
        singular_geometry,
        no_integer_candidates,
        too_many_integer_candidates,
        invalid_dimensions,
        non_finite_input,
        search_limit_exceeded
    }
}

#[rustler::nif(schedule = "DirtyCpu")]
#[allow(clippy::too_many_arguments)]
pub fn rtk_filter_update_epoch<'a>(
    env: Env<'a>,
    state_term: StateTerm,
    epoch_term: EpochTerm,
    base: Vec3,
    model_term: ModelTerm,
    wavelengths: Vec<(String, f64)>,
    offsets: Vec<(String, f64)>,
    opts_term: UpdateOptsTerm,
) -> NifResult<Term<'a>> {
    let Some(model) = decode_model(model_term) else {
        return Ok((atoms::error(), atoms::invalid_stochastic_model()).encode(env));
    };

    let update = match update_epoch(
        decode_state(state_term),
        &decode_epoch(epoch_term),
        vec3(base),
        &model,
        &wavelengths.into_iter().collect::<BTreeMap<_, _>>(),
        &offsets.into_iter().collect::<BTreeMap<_, _>>(),
        &decode_opts(opts_term),
    ) {
        Ok(update) => update,
        Err(err) => return Ok((atoms::error(), encode_update_error(env, err)).encode(env)),
    };

    Ok((
        atoms::ok(),
        (
            encode_state(update.state),
            update.integer_ratio,
            update.integer_fixed,
            update.newly_fixed,
            update.fixed_ids,
        ),
    )
        .encode(env))
}

fn encode_update_error<'a>(env: Env<'a>, err: UpdateError) -> Term<'a> {
    match err {
        UpdateError::ReferenceChanged { expected, actual } => {
            (atoms::reference_changed(), expected, actual).encode(env)
        }
        UpdateError::MissingAmbiguityColumn(id) => {
            (atoms::missing_ambiguity_column(), id).encode(env)
        }
        UpdateError::MissingWavelength(id) => (atoms::missing_wavelength(), id).encode(env),
        UpdateError::MissingOffset(id) => (atoms::missing_offset(), id).encode(env),
        UpdateError::SingularGeometry => atoms::singular_geometry().encode(env),
        UpdateError::Ils(err) => encode_ils_error(env, err),
    }
}

fn encode_ils_error<'a>(env: Env<'a>, err: astrodynamics_gnss::ils::IlsError) -> Term<'a> {
    match err {
        astrodynamics_gnss::ils::IlsError::Singular => atoms::singular_geometry().encode(env),
        astrodynamics_gnss::ils::IlsError::NoCandidates(n) => {
            (atoms::no_integer_candidates(), n).encode(env)
        }
        astrodynamics_gnss::ils::IlsError::TooManyCandidates { evaluated, limit } => {
            (atoms::too_many_integer_candidates(), evaluated, limit).encode(env)
        }
        astrodynamics_gnss::ils::IlsError::InvalidDimensions { n, rows } => {
            (atoms::invalid_dimensions(), n, rows).encode(env)
        }
        astrodynamics_gnss::ils::IlsError::NonFinite => atoms::non_finite_input().encode(env),
        astrodynamics_gnss::ils::IlsError::SearchLimitExceeded => {
            atoms::search_limit_exceeded().encode(env)
        }
    }
}

fn decode_state(term: StateTerm) -> FilterState {
    let (
        (version, reference_sat, sd_ambiguity_ids, ambiguity_prior_sigma_m),
        baseline_m,
        sd_ambiguities_m,
        information,
        fixed_cycles,
        fixed_m,
    ) = term;

    FilterState {
        version,
        reference_sat,
        sd_ambiguity_ids,
        baseline_m: vec3(baseline_m),
        sd_ambiguities_m,
        information,
        ambiguity_prior_sigma_m,
        fixed_cycles: fixed_cycles.into_iter().collect(),
        fixed_m: fixed_m.into_iter().collect(),
    }
}

fn encode_state(state: FilterState) -> StateTerm {
    (
        (
            state.version,
            state.reference_sat,
            state.sd_ambiguity_ids,
            state.ambiguity_prior_sigma_m,
        ),
        tuple3(state.baseline_m),
        state.sd_ambiguities_m,
        state.information,
        state.fixed_cycles.into_iter().collect(),
        state.fixed_m.into_iter().collect(),
    )
}

fn decode_epoch(term: EpochTerm) -> Epoch {
    let (reference, nonref) = term;
    Epoch {
        reference: decode_sat(reference),
        nonref: nonref.into_iter().map(decode_sat).collect(),
    }
}

fn decode_sat(term: SatTerm) -> SatMeas {
    let (
        (sat, sd_ambiguity_id),
        (base_code_m, base_phase_m, rover_code_m, rover_phase_m),
        (base_tx_pos, rover_tx_pos, pos),
    ) = term;

    SatMeas {
        sat,
        sd_ambiguity_id,
        base_code_m,
        base_phase_m,
        rover_code_m,
        rover_phase_m,
        base_tx_pos: vec3(base_tx_pos),
        rover_tx_pos: vec3(rover_tx_pos),
        pos: vec3(pos),
    }
}

fn decode_model(term: ModelTerm) -> Option<MeasModel> {
    let (code_sigma_m, phase_sigma_m, stochastic, elevation_weighting, sagnac) = term;
    let stochastic = match stochastic.as_str() {
        "simple" => StochasticModel::Simple {
            elevation_weighting,
        },
        "rtklib" => StochasticModel::Rtklib,
        _ => return None,
    };

    Some(MeasModel {
        code_sigma_m,
        phase_sigma_m,
        sagnac,
        stochastic,
    })
}

fn decode_opts(term: UpdateOptsTerm) -> UpdateOpts {
    let (hold_sigma_m, position_tol_m, ambiguity_tol_m, max_iterations, ratio_threshold) = term;
    UpdateOpts {
        hold_sigma_m,
        position_tol_m,
        ambiguity_tol_m,
        max_iterations,
        search: SearchOpts { ratio_threshold },
    }
}

fn vec3(v: Vec3) -> [f64; 3] {
    [v.0, v.1, v.2]
}

fn tuple3(v: [f64; 3]) -> Vec3 {
    (v[0], v[1], v[2])
}
