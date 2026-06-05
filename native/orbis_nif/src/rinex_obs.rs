//! Rustler boundary for the `astrodynamics-gnss` RINEX 3 observation product
//! and Hatanaka (CRINEX) decoder.
//!
//! Pure glue: it decodes Erlang terms, calls the crate's `rinex` public APIs,
//! holds the parsed product as a resource handle, and
//! encodes results back. No CRINEX grammar, RINEX parsing, or pseudorange
//! selection numerics live here — those are the crate's responsibility.
//!
//! - `crinex_decode/1` expands CRINEX text to plain RINEX text.
//! - `rinex_obs_parse/1` parses plain RINEX observation text into a handle.
//! - `crinex_obs_parse/1` decodes CRINEX then parses, in one dirty call, so a
//!   multi-megabyte expanded RINEX string is consumed inside Rust rather than
//!   marshalled across the BEAM boundary only to be passed straight back.
//! - the accessors expose the header, the epoch list, and per-epoch
//!   single-frequency pseudoranges as the `[{sat_token, range_m}]` shape the
//!   point-positioning solver consumes.

use astrodynamics_gnss::rinex::{
    decode_crinex,
    observations::{pseudoranges, ObsEpoch, RinexObs, SignalPolicy},
};
use astrodynamics_gnss::GnssSystem;
use rustler::{Encoder, Env, Error, NifResult, ResourceArc, Term};

/// Resource handle holding a parsed RINEX observation product across NIF calls.
pub struct RinexObsResource {
    pub obs: RinexObs,
}

#[rustler::resource_impl]
impl rustler::Resource for RinexObsResource {}

/// Decode CRINEX (Hatanaka) text into the plain RINEX observation text it
/// expands to.
///
/// Dirty-CPU: a daily file's expansion is unbounded relative to the 1 ms NIF
/// budget. Returns the decoded String, or the crate's parse-error reason.
#[rustler::nif(schedule = "DirtyCpu")]
fn crinex_decode(text: String) -> NifResult<String> {
    decode_crinex(&text).map_err(|e| Error::Term(Box::new(e.to_string())))
}

/// Parse plain RINEX 3 observation text into a resource handle.
///
/// Dirty-CPU: parsing a full daily file is unbounded relative to the NIF
/// budget. On a malformed file returns the parser's error as a term.
#[rustler::nif(schedule = "DirtyCpu")]
fn rinex_obs_parse(text: String) -> NifResult<ResourceArc<RinexObsResource>> {
    let obs = RinexObs::parse(&text).map_err(|e| Error::Term(Box::new(e.to_string())))?;
    Ok(ResourceArc::new(RinexObsResource { obs }))
}

/// Decode CRINEX text and parse the result in one dirty call.
///
/// The expanded RINEX text is consumed inside Rust, so only the compact typed
/// handle crosses back to the BEAM (the expansion is never marshalled).
#[rustler::nif(schedule = "DirtyCpu")]
fn crinex_obs_parse(text: String) -> NifResult<ResourceArc<RinexObsResource>> {
    let decoded = decode_crinex(&text).map_err(|e| Error::Term(Box::new(e.to_string())))?;
    let obs = RinexObs::parse(&decoded).map_err(|e| Error::Term(Box::new(e.to_string())))?;
    Ok(ResourceArc::new(RinexObsResource { obs }))
}

/// The surveyed a-priori receiver position `{x_m, y_m, z_m}` (ECEF meters), or
/// the atom `nil` when the file carries no `APPROX POSITION XYZ`.
#[rustler::nif]
fn rinex_obs_approx_position(env: Env<'_>, handle: ResourceArc<RinexObsResource>) -> Term<'_> {
    match handle.obs.header.approx_position_m {
        Some([x, y, z]) => (x, y, z).encode(env),
        None => rustler::types::atom::nil().encode(env),
    }
}

/// The per-constellation observation-code table as `[{"G", ["C1C", ...]}, ...]`
/// in declared order (system letter, then the code list).
#[rustler::nif]
fn rinex_obs_codes(handle: ResourceArc<RinexObsResource>) -> Vec<(String, Vec<String>)> {
    handle
        .obs
        .header
        .obs_codes
        .iter()
        .map(|(sys, codes)| (sys.letter().to_string(), codes.clone()))
        .collect()
}

/// The number of parsed epochs.
#[rustler::nif]
fn rinex_obs_epoch_count(handle: ResourceArc<RinexObsResource>) -> usize {
    handle.obs.epochs.len()
}

/// The epoch list as `[{ {{y,mo,d},{h,mi,second_float}}, flag, sat_count }]`, so
/// Elixir can index/select epochs without pulling every observation across the
/// boundary. The civil-time tuple is exactly the form `solve/4` accepts.
#[rustler::nif]
fn rinex_obs_epochs(env: Env<'_>, handle: ResourceArc<RinexObsResource>) -> Term<'_> {
    let list: Vec<Term> = handle
        .obs
        .epochs
        .iter()
        .map(|e| encode_epoch(env, e))
        .collect();
    list.encode(env)
}

/// Single-frequency pseudoranges for one epoch (by index), with an optional
/// per-system code override map `[{"G", ["C1C"]}, ...]` (an empty list uses the
/// crate's version-aware defaults).
///
/// Returns `{:ok, [{"G01", range_m}, ...]}` (exactly the solver's input shape)
/// or `{:error, :epoch_out_of_range}`.
#[rustler::nif(schedule = "DirtyCpu")]
fn rinex_obs_pseudoranges(
    env: Env<'_>,
    handle: ResourceArc<RinexObsResource>,
    epoch_index: usize,
    overrides: Vec<(String, Vec<String>)>,
) -> Term<'_> {
    let Some(epoch) = handle.obs.epochs.get(epoch_index) else {
        let reason = rustler::types::atom::Atom::from_str(env, "epoch_out_of_range")
            .map(|a| a.encode(env))
            .unwrap_or_else(|_| "epoch_out_of_range".encode(env));
        return (rustler::types::atom::error(), reason).encode(env);
    };

    // An empty override list uses the crate's version-aware defaults across all
    // systems; a non-empty override list defines the policy on its own (only the
    // listed systems are extracted), so a GPS-only request never pulls in, say,
    // GLONASS satellites that a later correction cannot model.
    let policy = if overrides.is_empty() {
        SignalPolicy::default_for(handle.obs.header.version)
    } else {
        let mut codes = std::collections::BTreeMap::new();
        for (letter, code_list) in overrides {
            if let Some(c) = letter.chars().next() {
                if let Some(system) = GnssSystem::from_letter(c) {
                    codes.insert(system, code_list);
                }
            }
        }
        SignalPolicy { codes }
    };

    let prs: Vec<(String, f64)> = pseudoranges(&handle.obs, epoch, &policy)
        .into_iter()
        .map(|(sat, range_m)| (sat.to_string(), range_m))
        .collect();

    (rustler::types::atom::ok(), prs).encode(env)
}

/// Encode one epoch as `{ {{y,mo,d},{h,mi,second_float}}, flag, sat_count }`.
fn encode_epoch<'a>(env: Env<'a>, epoch: &ObsEpoch) -> Term<'a> {
    let t = &epoch.epoch;
    let datetime = (
        (t.year, t.month as i32, t.day as i32),
        (t.hour as i32, t.minute as i32, t.second),
    );
    (datetime, epoch.flag, epoch.sats.len()).encode(env)
}
