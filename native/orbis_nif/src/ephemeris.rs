//! JPL SPK/BSP ephemeris file reader.
//!
//! Reads DAF (Double precision Array File) format SPK files containing
//! Type 2 segments (Chebyshev polynomials for position) as used by
//! DE421, DE440, etc.
//!
//! The DAF format is documented at:
//!   https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
//! The SPK Type 2 format is documented at:
//!   https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html

use rustler::NifResult;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

// ---------------------------------------------------------------------------
// DAF file record constants
// ---------------------------------------------------------------------------

/// DAF record size in bytes.
const RECORD_SIZE: usize = 1024;

/// Size of a double-precision value in bytes.
const F64_SIZE: usize = 8;

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

/// A segment summary from the DAF file. For SPK files (ND=2, NI=6):
///   - start_epoch, end_epoch: coverage window in seconds past J2000 TDB
///   - target: NAIF body code for the target
///   - center: NAIF body code for the center/observer
///   - frame: reference frame code (typically 1 = J2000)
///   - data_type: SPK segment type (we support type 2)
///   - start_addr, end_addr: 1-based double-precision element indices
#[derive(Debug, Clone)]
struct Segment {
    start_epoch: f64,
    end_epoch: f64,
    target: i32,
    center: i32,
    #[allow(dead_code)]
    frame: i32,
    data_type: i32,
    start_addr: i32,
    end_addr: i32,
}

/// Parsed SPK file with all segments indexed by (target, center) pair.
#[derive(Debug)]
pub(crate) struct SpkFile {
    path: String,
    /// Segments grouped by (target, center) for fast lookup.
    segments: HashMap<(i32, i32), Vec<Segment>>,
    /// Whether the file uses little-endian byte order.
    little_endian: bool,
}

// ---------------------------------------------------------------------------
// DAF reading
// ---------------------------------------------------------------------------

/// Read 8 bytes as an f64, respecting file endianness.
fn read_f64(buf: &[u8], offset: usize, little_endian: bool) -> f64 {
    let bytes: [u8; 8] = buf[offset..offset + 8].try_into().unwrap();
    if little_endian {
        f64::from_le_bytes(bytes)
    } else {
        f64::from_be_bytes(bytes)
    }
}

/// Read 4 bytes as an i32, respecting file endianness.
fn read_i32(buf: &[u8], offset: usize, little_endian: bool) -> i32 {
    let bytes: [u8; 4] = buf[offset..offset + 4].try_into().unwrap();
    if little_endian {
        i32::from_le_bytes(bytes)
    } else {
        i32::from_be_bytes(bytes)
    }
}

/// Read a full 1024-byte record (1-based record number).
fn read_record(file: &mut File, record_num: usize) -> Result<[u8; RECORD_SIZE], String> {
    let offset = (record_num - 1) as u64 * RECORD_SIZE as u64;
    file.seek(SeekFrom::Start(offset))
        .map_err(|e| format!("seek error: {e}"))?;
    let mut buf = [0u8; RECORD_SIZE];
    file.read_exact(&mut buf)
        .map_err(|e| format!("read error: {e}"))?;
    buf[..].try_into().map_err(|_| "buffer conversion".to_string())
}

/// Detect endianness from the DAF file record. The file record contains
/// an 8-byte string at offset 88 indicating byte order: "LTL-IEEE" or
/// "BIG-IEEE".
fn detect_endianness(file_record: &[u8; RECORD_SIZE]) -> Result<bool, String> {
    let id = std::str::from_utf8(&file_record[88..96])
        .map_err(|_| "invalid endian marker".to_string())?;
    if id.starts_with("LTL") {
        Ok(true)
    } else if id.starts_with("BIG") {
        Ok(false)
    } else {
        Err(format!("unknown endian marker: {id:?}"))
    }
}

/// Parse the DAF file record (record 1) to extract ND, NI, and the
/// forward pointer to the first summary record.
fn parse_file_record(
    file_record: &[u8; RECORD_SIZE],
    little_endian: bool,
) -> Result<(usize, usize, usize), String> {
    // Bytes 8..12: ND (number of double-precision components in each summary)
    let nd = read_i32(file_record, 8, little_endian) as usize;
    // Bytes 12..16: NI (number of integer components in each summary)
    let ni = read_i32(file_record, 12, little_endian) as usize;

    // For SPK files: ND=2, NI=6
    if nd != 2 || ni != 6 {
        return Err(format!("expected SPK (ND=2, NI=6), got ND={nd}, NI={ni}"));
    }

    // Bytes 76..80: forward pointer (record number of first summary record)
    let fward = read_i32(file_record, 76, little_endian) as usize;

    Ok((nd, ni, fward))
}

/// Parse all segment summaries from the DAF file by traversing the
/// linked list of summary records.
fn parse_summaries(file: &mut File, little_endian: bool, fward: usize) -> Result<Vec<Segment>, String> {
    let mut segments = Vec::new();
    let mut current_record = fward;

    // Each summary for SPK has ND=2 doubles + NI=6 ints.
    // Packed size: 2*8 + 6*4 = 40 bytes.
    // But DAF stores summaries as (ND + (NI+1)/2) doubles = 2 + 3 = 5 doubles = 40 bytes.
    let summary_size = 5; // in doubles (40 bytes)

    while current_record != 0 {
        let record = read_record(file, current_record)?;

        // First 8 bytes: next summary record pointer (as f64)
        let next = read_f64(&record, 0, little_endian) as usize;
        // Bytes 8..16: previous summary record pointer (as f64)
        // let _prev = read_f64(&record, 8, little_endian) as usize;
        // Bytes 16..24: number of summaries in this record (as f64)
        let n_summaries = read_f64(&record, 16, little_endian) as usize;

        // Summaries start at byte 24.
        for i in 0..n_summaries {
            let base = 24 + i * summary_size * F64_SIZE;

            // 2 doubles: start_epoch, end_epoch (seconds past J2000 TDB)
            let start_epoch = read_f64(&record, base, little_endian);
            let end_epoch = read_f64(&record, base + 8, little_endian);

            // 6 integers packed into remaining space (24 bytes = 6 * i32)
            let int_base = base + 16;
            let target = read_i32(&record, int_base, little_endian);
            let center = read_i32(&record, int_base + 4, little_endian);
            let frame = read_i32(&record, int_base + 8, little_endian);
            let data_type = read_i32(&record, int_base + 12, little_endian);
            let start_addr = read_i32(&record, int_base + 16, little_endian);
            let end_addr = read_i32(&record, int_base + 20, little_endian);

            segments.push(Segment {
                start_epoch,
                end_epoch,
                target,
                center,
                frame,
                data_type,
                start_addr,
                end_addr,
            });
        }

        current_record = next;
    }

    Ok(segments)
}

/// Open and parse an SPK file, returning the indexed structure.
pub(crate) fn open_spk(path: &str) -> Result<SpkFile, String> {
    let mut file = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;

    // Read and verify file record (record 1).
    let file_record = read_record(&mut file, 1)?;

    // Check magic: first 7 bytes should be "DAF/SPK"
    let magic = std::str::from_utf8(&file_record[0..7])
        .map_err(|_| "invalid file header".to_string())?;
    if magic != "DAF/SPK" {
        return Err(format!("not an SPK file: header is {magic:?}"));
    }

    let little_endian = detect_endianness(&file_record)?;
    let (_nd, _ni, fward) = parse_file_record(&file_record, little_endian)?;

    let segment_list = parse_summaries(&mut file, little_endian, fward)?;

    // Index segments by (target, center).
    let mut segments: HashMap<(i32, i32), Vec<Segment>> = HashMap::new();
    for seg in segment_list {
        segments
            .entry((seg.target, seg.center))
            .or_default()
            .push(seg);
    }

    Ok(SpkFile {
        path: path.to_string(),
        segments,
        little_endian,
    })
}

// ---------------------------------------------------------------------------
// SPK Type 2 evaluation
// ---------------------------------------------------------------------------

/// Read a range of f64 values from the file. Addresses are 1-based
/// element indices (each element is 8 bytes).
fn read_f64_array(
    file: &mut File,
    start_addr: i32,
    count: usize,
    little_endian: bool,
) -> Result<Vec<f64>, String> {
    let byte_offset = (start_addr as u64 - 1) * F64_SIZE as u64;
    file.seek(SeekFrom::Start(byte_offset))
        .map_err(|e| format!("seek error: {e}"))?;

    let mut buf = vec![0u8; count * F64_SIZE];
    file.read_exact(&mut buf)
        .map_err(|e| format!("read error: {e}"))?;

    let mut result = Vec::with_capacity(count);
    for i in 0..count {
        result.push(read_f64(&buf, i * F64_SIZE, little_endian));
    }
    Ok(result)
}

/// SPK Type 2 segment metadata, stored at the end of the segment's
/// data array:
///   - init: initial epoch of the first record (seconds past J2000 TDB)
///   - intlen: interval length per record (seconds)
///   - rsize: record size (number of f64 elements per record)
///   - n: number of records
struct Type2Meta {
    init: f64,
    intlen: f64,
    rsize: usize,
    n: usize,
}

/// Read the Type 2 metadata from the last 4 doubles of a segment.
fn read_type2_meta(
    file: &mut File,
    segment: &Segment,
    little_endian: bool,
) -> Result<Type2Meta, String> {
    // The last 4 doubles in the segment are: INIT, INTLEN, RSIZE, N
    let meta_addr = segment.end_addr - 3; // 4 values starting here
    let values = read_f64_array(file, meta_addr, 4, little_endian)?;

    Ok(Type2Meta {
        init: values[0],
        intlen: values[1],
        rsize: values[2] as usize,
        n: values[3] as usize,
    })
}

/// Evaluate a Chebyshev polynomial using the recurrence relation.
///
/// Given coefficients c[0..degree+1] and normalized time `x` in [-1, 1],
/// returns the value of the polynomial.
fn chebyshev_eval(coeffs: &[f64], x: f64) -> f64 {
    let n = coeffs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return coeffs[0];
    }

    // Clenshaw recurrence (evaluated from high to low order).
    let mut b_k1 = 0.0; // b_{k+1}
    let mut b_k2 = 0.0; // b_{k+2}
    let two_x = 2.0 * x;

    for i in (1..n).rev() {
        let b_k = coeffs[i] + (two_x * b_k1 - b_k2);
        b_k2 = b_k1;
        b_k1 = b_k;
    }

    coeffs[0] + (x * b_k1 - b_k2)
}

/// Evaluate the derivative of a Chebyshev polynomial.
///
/// Given coefficients c[0..degree+1] and normalized time `x` in [-1, 1],
/// returns dP/dx. To get dP/dt, divide by the half-interval length.
fn chebyshev_deriv(coeffs: &[f64], x: f64) -> f64 {
    let n = coeffs.len();
    if n <= 1 {
        return 0.0;
    }

    // Use recurrence for derivative of Chebyshev series:
    //   T'_0 = 0, T'_1 = 1, T'_n = 2*T_{n-1} + (n/(n-2))*T'_{n-2} for n >= 2
    // But it's simpler to use: d/dx sum(c_i T_i(x)) via Clenshaw on derivative coefficients.
    //
    // More directly: if P(x) = sum c_i T_i(x), then
    //   P'(x) = sum c_i T_i'(x)
    // where T_0' = 0, T_1' = 1, T_n'(x) = 2x T_{n-1}'(x) - T_{n-2}'(x) + 2 T_{n-1}(x)
    //
    // Simplest approach: compute T_i(x) and T_i'(x) simultaneously.
    let mut t_prev = 1.0; // T_0
    let mut t_curr = x; // T_1
    let mut dt_prev = 0.0; // T_0'
    let mut dt_curr = 1.0; // T_1'

    let mut result = coeffs[1] * dt_curr; // c_1 * T_1'

    for coeff in coeffs.iter().take(n).skip(2) {
        let t_next = 2.0 * x * t_curr - t_prev;
        let dt_next = 2.0 * t_curr + 2.0 * x * dt_curr - dt_prev;

        result += coeff * dt_next;

        t_prev = t_curr;
        t_curr = t_next;
        dt_prev = dt_curr;
        dt_curr = dt_next;
    }

    result
}

/// Compute position (and optionally velocity) for a single segment
/// using split JD TDB (whole, fraction) for maximum precision.
///
/// Mirrors Skyfield/jplephem's split-time Chebyshev argument computation:
///   index1, offset1 = divmod((tdb - T0) * S_PER_DAY - init, intlen)
///   index2, offset2 = divmod(tdb2 * S_PER_DAY, intlen)
///   index3, offset  = divmod(offset1 + offset2, intlen)
fn evaluate_type2_segment(
    file: &mut File,
    segment: &Segment,
    jd_whole: f64,
    jd_fraction: f64,
    little_endian: bool,
    compute_velocity: bool,
) -> Result<([f64; 3], Option<[f64; 3]>), String> {
    let meta = read_type2_meta(file, segment, little_endian)?;

    // Split-precision epoch computation matching jplephem exactly.
    let s_per_day = 86400.0;
    let t0 = 2451545.0;

    // Whole-day part: (tdb_whole - T0) * S_PER_DAY - init
    let whole_seconds = (jd_whole - t0) * s_per_day - meta.init;
    let (index1, offset1) = split_divmod(whole_seconds, meta.intlen);

    // Fractional-day part: tdb_fraction * S_PER_DAY
    let frac_seconds = jd_fraction * s_per_day;
    let (index2, offset2) = split_divmod(frac_seconds, meta.intlen);

    // Combine offsets
    let (index3, offset) = split_divmod(offset1 + offset2, meta.intlen);
    let mut record_index = (index1 + index2 + index3) as usize;

    // Clamp to valid range
    if record_index >= meta.n {
        record_index = meta.n - 1;
    }

    let n_coeffs = (meta.rsize - 2) / 3;
    let record_addr = segment.start_addr + (record_index * meta.rsize) as i32;
    let data = read_f64_array(file, record_addr, meta.rsize, little_endian)?;

    // Normalize offset to [-1, 1] for Chebyshev evaluation.
    // s = 2 * offset / intlen - 1
    let s = 2.0 * offset / meta.intlen - 1.0;


    let cx = &data[2..2 + n_coeffs];
    let cy = &data[2 + n_coeffs..2 + 2 * n_coeffs];
    let cz = &data[2 + 2 * n_coeffs..2 + 3 * n_coeffs];

    let position = [
        chebyshev_eval(cx, s),
        chebyshev_eval(cy, s),
        chebyshev_eval(cz, s),
    ];


    let velocity = if compute_velocity {
        let radius = meta.intlen / 2.0;
        Some([
            chebyshev_deriv(cx, s) / radius,
            chebyshev_deriv(cy, s) / radius,
            chebyshev_deriv(cz, s) / radius,
        ])
    } else {
        None
    };

    Ok((position, velocity))
}

/// Python-style divmod for floats: returns (quotient_floor, remainder).
fn split_divmod(value: f64, divisor: f64) -> (f64, f64) {
    let q = (value / divisor).floor();
    let r = value - q * divisor;
    (q, r)
}

// ---------------------------------------------------------------------------
// Body code lookups and chain resolution
// ---------------------------------------------------------------------------

/// NAIF body codes used in JPL ephemeris files.
pub(crate) fn body_name_to_code(name: &str) -> Result<i32, String> {
    match name {
        "ssb" | "solar_system_barycenter" => Ok(0),
        "mercury_barycenter" | "mercury" => Ok(1),
        "venus_barycenter" | "venus" => Ok(2),
        "earth_moon_barycenter" | "emb" => Ok(3),
        "mars_barycenter" | "mars" => Ok(4),
        "jupiter_barycenter" | "jupiter" => Ok(5),
        "saturn_barycenter" | "saturn" => Ok(6),
        "uranus_barycenter" | "uranus" => Ok(7),
        "neptune_barycenter" | "neptune" => Ok(8),
        "pluto_barycenter" | "pluto" => Ok(9),
        "sun" => Ok(10),
        "moon" => Ok(301),
        "earth" => Ok(399),
        _ => Err(format!("unknown body: {name}")),
    }
}

/// Convert a split Julian Date (TDB) to seconds past J2000.0 TDB.
/// J2000.0 = JD 2451545.0
/// Using the split form (whole, fraction) preserves precision: the
/// subtraction of J2000 from the integer part is exact, leaving the
/// full precision of the fractional day for the multiplication.
fn jd_to_spk_epoch_split(whole: f64, fraction: f64) -> f64 {
    (whole - 2451545.0 + fraction) * 86400.0
}

/// Single-JD convenience (used only when split form is unavailable).
#[allow(dead_code)]
fn jd_to_spk_epoch(jd: f64) -> f64 {
    (jd - 2451545.0) * 86400.0
}

/// Compute position by finding a chain of segments connecting the
/// target to the observer through the available segment pairs.
///
/// For example, to get Moon relative to Earth:
///   Moon(301) → Earth-Moon Barycenter(3) [from segment (301, 3)]
///   Earth(399) → Earth-Moon Barycenter(3) [from segment (399, 3)]
///   Result = pos(301, 3) - pos(399, 3)
///
/// For Sun relative to Earth:
///   Sun(10) → SSB(0) [from segment (10, 0)]
///   Earth-Moon Barycenter(3) → SSB(0) [from segment (3, 0)]
///   Earth(399) → EMB(3) [from segment (399, 3)]
///   Result = pos(10, 0) - pos(3, 0) - pos(399, 3)
fn compute_position_chain(
    spk: &SpkFile,
    target: i32,
    observer: i32,
    jd_whole: f64,
    jd_fraction: f64,
    epoch_s: f64,
    skyfield_compat: bool,
) -> Result<[f64; 3], String> {
    let observer_steps = build_chain_with_centers(spk, observer, jd_whole, jd_fraction, epoch_s)?;
    let target_steps = build_chain_with_centers(spk, target, jd_whole, jd_fraction, epoch_s)?;

    if skyfield_compat {
        // Replicate Skyfield's exact VectorSum._at() accumulation order:
        // reversed observer chain (negated) then forward target chain,
        // all in AU, going through SSB even when a common ancestor exists.
        // This introduces intermediate cancellation that changes FP rounding
        // but matches Skyfield's output at 0 ULP.
        const AU_KM: f64 = 149597870.700;

        let mut p = [0.0_f64; 3];

        for (_center, pos) in &observer_steps {
            p[0] += -pos[0] / AU_KM;
            p[1] += -pos[1] / AU_KM;
            p[2] += -pos[2] / AU_KM;
        }

        for (_center, pos) in target_steps.iter().rev() {
            p[0] += pos[0] / AU_KM;
            p[1] += pos[1] / AU_KM;
            p[2] += pos[2] / AU_KM;
        }

        Ok([p[0] * AU_KM, p[1] * AU_KM, p[2] * AU_KM])
    } else {
        // Precise mode: find the common ancestor and subtract directly,
        // avoiding the SSB round-trip that loses precision through
        // cancellation of large intermediate values.
        let target_centers: Vec<i32> = target_steps.iter().map(|(c, _)| *c).collect();

        let common = observer_steps.iter()
            .map(|(c, _)| *c)
            .find(|c| target_centers.contains(c))
            .unwrap_or(0);

        let mut target_sum = [0.0_f64; 3];
        for (center, pos) in &target_steps {
            target_sum[0] += pos[0];
            target_sum[1] += pos[1];
            target_sum[2] += pos[2];
            if *center == common { break; }
        }

        let mut observer_sum = [0.0_f64; 3];
        for (center, pos) in &observer_steps {
            observer_sum[0] += pos[0];
            observer_sum[1] += pos[1];
            observer_sum[2] += pos[2];
            if *center == common { break; }
        }

        Ok([
            target_sum[0] - observer_sum[0],
            target_sum[1] - observer_sum[1],
            target_sum[2] - observer_sum[2],
        ])
    }
}

fn build_chain_with_centers(
    spk: &SpkFile,
    body: i32,
    jd_whole: f64,
    jd_fraction: f64,
    epoch_s: f64,
) -> Result<Vec<(i32, [f64; 3])>, String> {
    let mut chain = Vec::new();
    let mut current = body;

    for _ in 0..10 {
        if current == 0 {
            break;
        }
        let (center, position) = find_body_as_target(spk, current, jd_whole, jd_fraction, epoch_s)?;
        chain.push((center, position));
        current = center;
    }

    if current != 0 {
        return Err(format!(
            "could not build chain from body {body} to SSB (stuck at {current})"
        ));
    }

    Ok(chain)
}

#[allow(dead_code)]
fn sum_positions(chain: &[(i32, [f64; 3])]) -> [f64; 3] {
    let mut sum = [0.0; 3];
    for (_, pos) in chain {
        sum[0] += pos[0];
        sum[1] += pos[1];
        sum[2] += pos[2];
    }
    sum
}


fn find_body_as_target(
    spk: &SpkFile,
    body: i32,
    jd_whole: f64,
    jd_fraction: f64,
    epoch_s: f64,
) -> Result<(i32, [f64; 3]), String> {
    for (&(target, center), segs) in &spk.segments {
        if target == body {
            // Use epoch_s for the range check (sufficient precision).
            if let Some(seg) = segs.iter().find(|s| epoch_s >= s.start_epoch && epoch_s <= s.end_epoch) {
                if seg.data_type != 2 {
                    return Err(format!(
                        "unsupported SPK type {} for body {body}",
                        seg.data_type
                    ));
                }

                let mut file = File::open(&spk.path)
                    .map_err(|e| format!("cannot reopen {}: {e}", spk.path))?;

                // Pass split JD for full-precision Chebyshev argument.
                let (position, _) =
                    evaluate_type2_segment(&mut file, seg, jd_whole, jd_fraction, spk.little_endian, false)?;
                return Ok((center, position));
            }
        }
    }
    Err(format!(
        "no segment found for body {body} at the requested epoch"
    ))
}


// ---------------------------------------------------------------------------
// NIF entry point
// ---------------------------------------------------------------------------

/// NIF: get_body_position(file_path, target_name, observer_name, jd_whole, jd_fraction)
///
/// Accepts a split Julian Date (TDB) for full precision. Returns {x, y, z} in km.
pub(crate) fn get_body_position_impl(
    file_path: String,
    target_name: String,
    observer_name: String,
    jd_whole: f64,
    jd_fraction: f64,
    skyfield_compat: bool,
) -> NifResult<(f64, f64, f64)> {
    let target_code = body_name_to_code(&target_name)
        .map_err(|e| rustler::Error::Term(Box::new(e)))?;
    let observer_code = body_name_to_code(&observer_name)
        .map_err(|e| rustler::Error::Term(Box::new(e)))?;

    let spk = open_spk(&file_path)
        .map_err(|e| rustler::Error::Term(Box::new(e)))?;

    let epoch_s = jd_to_spk_epoch_split(jd_whole, jd_fraction);

    let position = compute_position_chain(
        &spk, target_code, observer_code, jd_whole, jd_fraction, epoch_s, skyfield_compat,
    )
    .map_err(|e| rustler::Error::Term(Box::new(e)))?;

    Ok((position[0], position[1], position[2]))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chebyshev_eval_constant() {
        // T_0(x) = 1, so coeffs = [5.0] should give 5.0 for any x.
        assert_eq!(chebyshev_eval(&[5.0], 0.0), 5.0);
        assert_eq!(chebyshev_eval(&[5.0], 0.5), 5.0);
        assert_eq!(chebyshev_eval(&[5.0], -1.0), 5.0);
    }

    #[test]
    fn test_chebyshev_eval_linear() {
        // T_0(x) = 1, T_1(x) = x
        // coeffs = [3.0, 2.0] => 3.0 + 2.0*x
        let c = [3.0, 2.0];
        assert!((chebyshev_eval(&c, 0.0) - 3.0).abs() < 1e-15);
        assert!((chebyshev_eval(&c, 1.0) - 5.0).abs() < 1e-15);
        assert!((chebyshev_eval(&c, -1.0) - 1.0).abs() < 1e-15);
        assert!((chebyshev_eval(&c, 0.5) - 4.0).abs() < 1e-15);
    }

    #[test]
    fn test_chebyshev_eval_quadratic() {
        // T_2(x) = 2x^2 - 1
        // coeffs = [1.0, 0.0, 1.0] => 1*T_0 + 0*T_1 + 1*T_2 = 1 + (2x^2-1) = 2x^2
        let c = [1.0, 0.0, 1.0];
        assert!((chebyshev_eval(&c, 0.0) - 0.0).abs() < 1e-15);
        assert!((chebyshev_eval(&c, 1.0) - 2.0).abs() < 1e-15);
        assert!((chebyshev_eval(&c, -1.0) - 2.0).abs() < 1e-15);
        assert!((chebyshev_eval(&c, 0.5) - 0.5).abs() < 1e-15);
    }

    #[test]
    fn test_chebyshev_deriv_linear() {
        // P(x) = 3 + 2x => P'(x) = 2
        let c = [3.0, 2.0];
        assert!((chebyshev_deriv(&c, 0.0) - 2.0).abs() < 1e-15);
        assert!((chebyshev_deriv(&c, 0.5) - 2.0).abs() < 1e-15);
    }

    #[test]
    fn test_chebyshev_deriv_quadratic() {
        // P(x) = T_0 + T_2 = 1 + (2x^2-1) = 2x^2
        // P'(x) = 4x
        let c = [1.0, 0.0, 1.0];
        assert!((chebyshev_deriv(&c, 0.0) - 0.0).abs() < 1e-14);
        assert!((chebyshev_deriv(&c, 0.5) - 2.0).abs() < 1e-14);
        assert!((chebyshev_deriv(&c, 1.0) - 4.0).abs() < 1e-14);
    }

    #[test]
    fn test_body_name_to_code() {
        assert_eq!(body_name_to_code("sun").unwrap(), 10);
        assert_eq!(body_name_to_code("earth").unwrap(), 399);
        assert_eq!(body_name_to_code("moon").unwrap(), 301);
        assert_eq!(body_name_to_code("mars").unwrap(), 4);
        assert_eq!(body_name_to_code("ssb").unwrap(), 0);
        assert!(body_name_to_code("unknown_body").is_err());
    }

    #[test]
    fn test_jd_to_spk_epoch() {
        // J2000.0 should map to 0.0 seconds.
        assert!((jd_to_spk_epoch(2451545.0) - 0.0).abs() < 1e-10);
        // One day later should be 86400 seconds.
        assert!((jd_to_spk_epoch(2451546.0) - 86400.0).abs() < 1e-6);
    }
}
