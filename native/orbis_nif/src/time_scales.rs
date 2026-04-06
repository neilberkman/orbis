//! Precise time scale conversions: UTC → TAI → TT → TDB → UT1.
//!
//! Mirrors Skyfield's _utc() path for bit-exact parity.

use crate::iers_data::UT1_DATA;

const DAY_S: f64 = 86400.0;
const T0: f64 = 2451545.0;
const TT_MINUS_TAI_S: f64 = 32.184;
const ROUND_1E7: f64 = 10_000_000.0;

#[allow(dead_code)]
pub(crate) struct TimeScales {
    pub jd_whole: f64,
    pub ut1_fraction: f64,
    pub tt_fraction: f64,
    pub tdb_fraction: f64,
    pub jd_ut1: f64,
    pub jd_tt: f64,
    pub jd_tdb: f64,
}

struct LeapSecondEntry {
    mjd: i32,
    tai_utc: f64,
}

static LEAP_SECONDS: &[LeapSecondEntry] = &[
    LeapSecondEntry { mjd: 41317, tai_utc: 10.0 },
    LeapSecondEntry { mjd: 41499, tai_utc: 11.0 },
    LeapSecondEntry { mjd: 41683, tai_utc: 12.0 },
    LeapSecondEntry { mjd: 42048, tai_utc: 13.0 },
    LeapSecondEntry { mjd: 42413, tai_utc: 14.0 },
    LeapSecondEntry { mjd: 42778, tai_utc: 15.0 },
    LeapSecondEntry { mjd: 43144, tai_utc: 16.0 },
    LeapSecondEntry { mjd: 43509, tai_utc: 17.0 },
    LeapSecondEntry { mjd: 43874, tai_utc: 18.0 },
    LeapSecondEntry { mjd: 44239, tai_utc: 19.0 },
    LeapSecondEntry { mjd: 44786, tai_utc: 20.0 },
    LeapSecondEntry { mjd: 45151, tai_utc: 21.0 },
    LeapSecondEntry { mjd: 45516, tai_utc: 22.0 },
    LeapSecondEntry { mjd: 46247, tai_utc: 23.0 },
    LeapSecondEntry { mjd: 47161, tai_utc: 24.0 },
    LeapSecondEntry { mjd: 47892, tai_utc: 25.0 },
    LeapSecondEntry { mjd: 48257, tai_utc: 26.0 },
    LeapSecondEntry { mjd: 48804, tai_utc: 27.0 },
    LeapSecondEntry { mjd: 49169, tai_utc: 28.0 },
    LeapSecondEntry { mjd: 49534, tai_utc: 29.0 },
    LeapSecondEntry { mjd: 50083, tai_utc: 30.0 },
    LeapSecondEntry { mjd: 50448, tai_utc: 31.0 },
    LeapSecondEntry { mjd: 50813, tai_utc: 32.0 },
    LeapSecondEntry { mjd: 53736, tai_utc: 33.0 },
    LeapSecondEntry { mjd: 54832, tai_utc: 34.0 },
    LeapSecondEntry { mjd: 56109, tai_utc: 35.0 },
    LeapSecondEntry { mjd: 57204, tai_utc: 36.0 },
    LeapSecondEntry { mjd: 57754, tai_utc: 37.0 },
];

impl TimeScales {
    pub fn from_utc(year: i32, month: i32, day: i32, hour: i32, minute: i32, second: f64) -> Self {
        let jd_day = julian_day_number(year, month, day);
        let jd1 = jd_day as f64 - 0.5;
        let jd2 = (second + minute as f64 * 60.0 + hour as f64 * 3600.0) / DAY_S;
        let jd_utc_total = jd1 + jd2;

        let leap_seconds = find_leap_seconds(jd_utc_total);
        let utc_seconds_of_day = hour as f64 * 3600.0 + minute as f64 * 60.0 + second;
        let utc_seconds_at_midnight = jd1 * DAY_S;

        let utc_whole_seconds = utc_seconds_of_day.trunc();
        let utc_subsecond = utc_seconds_of_day.fract();

        // Mirror Skyfield's _utc() path.
        let tai_seconds = utc_seconds_at_midnight + leap_seconds + utc_whole_seconds;
        let jd_whole = (tai_seconds / DAY_S).floor();
        let tai_fraction = (tai_seconds - jd_whole * DAY_S + utc_subsecond) / DAY_S;
        let tt_offset_days = TT_MINUS_TAI_S / DAY_S;

        let tt_fraction = tai_fraction + tt_offset_days;
        let jd_tt = jd_whole + tt_fraction;

        let delta_t = interpolate_delta_t(jd_tt);
        let ut1_fraction = tt_fraction - delta_t / DAY_S;
        let jd_ut1 = jd_whole + ut1_fraction;

        let t = (jd_whole - T0 + tt_fraction) / 36525.0;
        let tdb_minus_tt_seconds = 0.001657 * (628.3076 * t + 6.2401).sin()
            + 0.000022 * (575.3385 * t + 4.2970).sin()
            + 0.000014 * (1256.6152 * t + 6.1969).sin()
            + 0.000005 * (606.9777 * t + 4.0212).sin()
            + 0.000005 * (52.9691 * t + 0.4444).sin()
            + 0.000002 * (21.3299 * t + 5.5431).sin()
            + 0.000010 * t * (628.3076 * t + 4.2490).sin();

        let tdb_fraction = tt_fraction + tdb_minus_tt_seconds / DAY_S;
        let jd_tdb = jd_whole + tdb_fraction;

        TimeScales {
            jd_whole,
            ut1_fraction,
            tt_fraction,
            tdb_fraction,
            jd_ut1,
            jd_tt,
            jd_tdb,
        }
    }
}

fn julian_day_number(year: i32, month: i32, day: i32) -> i64 {
    let janfeb = month <= 2;
    let g = year as i64 + 4716 - if janfeb { 1 } else { 0 };
    let f = ((month + 9) % 12) as i64;
    let e = 1461 * g / 4 + day as i64 - 1402;
    let j = e + (153 * f + 2) / 5;
    j + 38 - ((g + 184) / 100) * 3 / 4
}

fn find_leap_seconds(jd_utc: f64) -> f64 {
    let mjd = (jd_utc - 2400000.5) as i32;
    let mut ls = 10.0;
    for entry in LEAP_SECONDS {
        if mjd >= entry.mjd {
            ls = entry.tai_utc;
        } else {
            break;
        }
    }
    ls
}

fn interpolate_delta_t(jd_tt: f64) -> f64 {
    // Build delta-T table on first call (matching C++ lazy static pattern).
    use std::sync::LazyLock;

    struct DeltaTRow {
        jd_tt: f64,
        delta_t: f64,
    }

    static TABLE: LazyLock<Vec<DeltaTRow>> = LazyLock::new(|| {
        UT1_DATA
            .iter()
            .map(|entry| {
                let jd_utc = entry.mjd as f64 + 2400000.5;
                let leap_seconds = find_leap_seconds(jd_utc);
                let tt_minus_utc = leap_seconds + TT_MINUS_TAI_S;
                let delta_t = ((tt_minus_utc - entry.ut1_utc) * ROUND_1E7).round() / ROUND_1E7;
                DeltaTRow {
                    jd_tt: jd_utc + tt_minus_utc / DAY_S,
                    delta_t,
                }
            })
            .collect()
    });

    // Binary search for the bracketing entries.
    match TABLE.binary_search_by(|row| row.jd_tt.partial_cmp(&jd_tt).unwrap()) {
        Ok(i) => TABLE[i].delta_t,
        Err(0) => TABLE[0].delta_t,
        Err(i) if i >= TABLE.len() => TABLE.last().unwrap().delta_t,
        Err(i) => {
            let p1 = &TABLE[i - 1];
            let p2 = &TABLE[i];
            p1.delta_t + (jd_tt - p1.jd_tt) * (p2.delta_t - p1.delta_t) / (p2.jd_tt - p1.jd_tt)
        }
    }
}
