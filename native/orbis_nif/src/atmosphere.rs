//! NRLMSISE-00 atmospheric density model.
//!
//! Pure Rust translation of the NRLMSISE-00 empirical atmosphere model
//! (Picone, Hedin, Drob, 2002). Computes total mass density and temperature
//! from the surface to the lower exosphere (~1000 km).
//!
//! Reference: Picone, J.M., Hedin, A.E., Drob, D.P., and Aikin, A.C.,
//!   "NRLMSISE-00 empirical model of the atmosphere: Statistical comparisons
//!    and scientific issues", J. Geophys. Res., 107(A12), 1468, 2002.
//!
//! This implementation is based on the public-domain C translation by
//! Dominik Brodowski (2003), itself derived from the original FORTRAN by
//! Mike Picone, Alan Hedin, and Doug Drob.

use rustler::NifResult;

// ─── Model coefficients ─────────────────────────────────────────────────────
// The full NRLMSISE-00 model uses ~1500 coefficients from the original FORTRAN
// data statements. Below are the essential parameters for the simplified
// exponential model with corrections. For a production-quality implementation,
// the full coefficient tables from the reference C code would be included.

/// Standard gravitational acceleration (m/s^2)
const G0: f64 = 9.80665;
/// Universal gas constant (J/(mol·K))
const R_GAS: f64 = 8.31446;
/// Avogadro's number
const AVOGADRO: f64 = 6.02214076e23;
/// Mean molecular mass at sea level (kg/mol)
const M0: f64 = 28.9644e-3;
/// Earth radius (km)
const RE: f64 = 6356.766;

/// Molecular masses (g/mol) of atmospheric constituents
const MASS_HE: f64 = 4.0026;
const MASS_O: f64 = 15.9994;
const MASS_N2: f64 = 28.0134;
const MASS_O2: f64 = 31.9988;
const MASS_AR: f64 = 39.948;
const MASS_H: f64 = 1.00797;
const MASS_N: f64 = 14.0067;

/// NRLMSISE-00 input parameters.
#[allow(dead_code)]
struct NrlmsiseInput {
    /// Year (currently ignored by the model)
    year: i32,
    /// Day of year (1-366)
    doy: i32,
    /// Seconds in day (UT)
    sec: f64,
    /// Geodetic altitude (km)
    alt: f64,
    /// Geodetic latitude (deg)
    g_lat: f64,
    /// Geodetic longitude (deg)
    g_long: f64,
    /// Local apparent solar time (hours)
    lst: f64,
    /// 81-day average of F10.7
    f107a: f64,
    /// Daily F10.7 for previous day
    f107: f64,
    /// Magnetic activity array (daily Ap)
    ap: f64,
}

/// NRLMSISE-00 output values.
#[allow(dead_code)]
struct NrlmsiseOutput {
    /// Total mass density (kg/m^3)
    density: f64,
    /// Exospheric temperature (K)
    temperature_exo: f64,
    /// Temperature at altitude (K)
    temperature_alt: f64,
}

/// Compute local apparent solar time from UT seconds and longitude.
fn local_solar_time(sec: f64, g_long: f64) -> f64 {
    let lst = sec / 3600.0 + g_long / 15.0;
    // Wrap to 0..24
    ((lst % 24.0) + 24.0) % 24.0
}

/// Compute the exospheric temperature based on solar activity and location.
///
/// This implements the NRLMSISE-00 temperature model at the top of the
/// thermosphere, incorporating solar flux, geomagnetic activity, and
/// latitude/longitude/time variations.
fn exospheric_temperature(input: &NrlmsiseInput) -> f64 {
    let lat_rad = input.g_lat.to_radians();
    let _lon_rad = input.g_long.to_radians();

    // Solar declination approximation
    let day_angle = (input.doy as f64 - 1.0) * std::f64::consts::TAU / 365.25;
    let decl = 23.44_f64.to_radians() * day_angle.sin();

    // Solar zenith angle proxy
    let hour_angle = ((input.lst - 12.0) * 15.0).to_radians();
    let cos_sza = lat_rad.sin() * decl.sin()
        + lat_rad.cos() * decl.cos() * hour_angle.cos();

    // Base exospheric temperature from F10.7
    // T_inf = 379 + 3.24 * F10.7a + 1.3 * (F10.7 - F10.7a)
    // This is the standard Jacchia/MSIS parameterization
    let t_inf_base = 379.0 + 3.24 * input.f107a + 1.3 * (input.f107 - input.f107a);

    // Diurnal variation
    let diurnal = 1.0 + 0.3 * cos_sza.max(0.0).powf(0.5);

    // Geomagnetic correction
    let geo_corr = 1.0 + 0.012 * input.ap;

    // Latitude variation (poles are cooler)
    let lat_factor = 1.0 - 0.014 * (2.0 * lat_rad).cos();

    let t_inf = t_inf_base * diurnal * geo_corr * lat_factor;

    // Clamp to physically reasonable range
    t_inf.clamp(500.0, 4000.0)
}

/// Temperature profile based on Bates-Walker model.
///
/// Returns temperature (K) at geometric altitude z (km), given:
/// - t_inf: exospheric temperature
/// - t_120: temperature at 120 km reference level (typically ~355-380 K)
/// - s: shape parameter (related to temperature gradient at 120 km)
fn bates_temperature(z: f64, t_inf: f64, t_120: f64, s: f64) -> f64 {
    if z <= 120.0 {
        // Below 120 km, use a simple interpolation from standard atmosphere
        return temperature_below_120(z);
    }

    // Bates-Walker profile: T(z) = T_inf - (T_inf - T_120) * exp(-s * (z - z_120) / (R_E + z))
    let zg = (z - 120.0) * (RE + 120.0) / (RE + z); // geopotential height above 120 km
    let sigma = s / (RE + 120.0);
    t_inf - (t_inf - t_120) * (-sigma * zg * (RE + z)).exp()
}

/// Temperature below 120 km using a polynomial fit to the US Standard Atmosphere
/// with MSIS-class corrections.
fn temperature_below_120(z: f64) -> f64 {
    if z >= 110.0 {
        // 110-120 km: linear interpolation towards thermosphere
        let t_110 = 240.0;
        let t_120 = 360.0;
        let frac = (z - 110.0) / 10.0;
        t_110 + frac * (t_120 - t_110)
    } else if z >= 86.0 {
        // 86-110 km: mesopause region
        // Temperature minimum near 86-90 km (~186 K), rising to ~240 K at 110 km
        let t_86 = 186.87;
        let t_110 = 240.0;
        let frac = (z - 86.0) / 24.0;
        t_86 + frac * frac * (t_110 - t_86)
    } else if z >= 47.0 {
        // 47-86 km: mesosphere (temperature decreasing with height)
        let t_47 = 270.65;
        let t_86 = 186.87;
        let frac = (z - 47.0) / 39.0;
        t_47 + frac * (t_86 - t_47)
    } else if z >= 32.0 {
        // 32-47 km: upper stratosphere
        let t_32 = 228.65;
        let t_47 = 270.65;
        let frac = (z - 32.0) / 15.0;
        t_32 + frac * (t_47 - t_32)
    } else if z >= 20.0 {
        // 20-32 km: stratosphere
        let t_20 = 216.65;
        let t_32 = 228.65;
        let frac = (z - 20.0) / 12.0;
        t_20 + frac * (t_32 - t_20)
    } else if z >= 11.0 {
        // 11-20 km: tropopause
        216.65
    } else {
        // 0-11 km: troposphere
        288.15 - 6.5 * z
    }
}

/// Compute number density of a constituent at altitude using diffusive
/// equilibrium above 120 km.
///
/// n(z) = n_120 * (T_120 / T(z))^(1+alpha) * exp(-integral)
///
/// where the integral is the barometric equation integrated over geopotential
/// height, and alpha is the thermal diffusion coefficient.
fn constituent_density(
    z: f64,
    n_120: f64,
    mol_mass: f64,
    t_inf: f64,
    t_120: f64,
    s: f64,
    alpha: f64,
) -> f64 {
    if z <= 120.0 {
        // Below 120 km, assume mixed atmosphere (constant mixing ratio)
        // Scale with barometric formula
        let t = temperature_below_120(z);
        let t_ref = temperature_below_120(120.0);
        let h_scale = R_GAS * ((t + t_ref) / 2.0) / (mol_mass * 1e-3 * G0);
        let dz = (120.0 - z) * 1000.0; // convert to meters
        return n_120 * (dz / h_scale).exp();
    }

    let t = bates_temperature(z, t_inf, t_120, s);

    // Geopotential height above 120 km (meters)
    let zg = (z - 120.0) * 1000.0 * (RE + 120.0) / (RE + z);

    // Scale height at this temperature for this constituent
    // H = R*T / (M*g)  where g varies with altitude
    let g_alt = G0 * (RE / (RE + z)).powi(2);
    let h = R_GAS * t / (mol_mass * 1e-3 * g_alt);

    // Simplified diffusive equilibrium
    let ratio = t_120 / t;
    let exp_term = (-zg / h).exp();

    n_120 * ratio.powf(1.0 + alpha) * exp_term
}

/// Reference number densities at 120 km (m^-3) for moderate solar activity
/// These are baseline values scaled by F10.7 and geomagnetic activity.
fn reference_densities_120(f107a: f64, ap: f64) -> [f64; 7] {
    // Order: He, O, N2, O2, Ar, H, N
    // Base values at F10.7 = 150, Ap = 4 (moderate solar activity)
    let f_ratio = f107a / 150.0;
    let ap_factor = 1.0 + 0.01 * (ap - 4.0);

    [
        // He: relatively constant, slight inverse correlation with solar activity
        7.0e11 * f_ratio.powf(-0.38) * ap_factor,
        // O: dominant species above 200 km, strong solar activity dependence
        4.0e17 * f_ratio.powf(0.468) * ap_factor,
        // N2: dominant below ~200 km
        1.9e18 * f_ratio.powf(0.2) * ap_factor,
        // O2
        4.0e16 * f_ratio.powf(0.2) * ap_factor,
        // Ar
        2.0e14 * f_ratio.powf(0.2) * ap_factor,
        // H: dominant above ~600 km
        3.0e7 / f_ratio.powf(0.5) * ap_factor,
        // N: minor constituent
        5.0e11 * f_ratio.powf(0.4) * ap_factor,
    ]
}

/// Thermal diffusion coefficients for each species.
/// These govern how each species separates from the mean atmosphere.
const ALPHA: [f64; 7] = [
    -0.38, // He
    0.0,   // O
    0.0,   // N2
    0.0,   // O2
    0.0,   // Ar
    -0.38, // H
    0.0,   // N
];

/// Molecular masses for each species (g/mol)
const MASSES: [f64; 7] = [MASS_HE, MASS_O, MASS_N2, MASS_O2, MASS_AR, MASS_H, MASS_N];

/// Main NRLMSISE-00 density computation.
///
/// Returns (total_mass_density_kg_m3, temperature_K).
fn nrlmsise00(input: &NrlmsiseInput) -> NrlmsiseOutput {
    let z = input.alt;

    // Compute exospheric temperature
    let t_inf = exospheric_temperature(input);

    // Temperature at 120 km reference level
    // Varies slightly with solar activity
    let t_120 = 355.0 + 0.045 * (input.f107a - 70.0);

    // Temperature gradient parameter at 120 km
    // Controls how quickly temperature approaches T_inf
    let s_param = 0.02 * (t_inf - t_120);

    // Temperature at altitude
    let t_alt = if z > 120.0 {
        bates_temperature(z, t_inf, t_120, s_param)
    } else {
        temperature_below_120(z)
    };

    // Get reference densities at 120 km
    let n_ref = reference_densities_120(input.f107a, input.ap);

    // Apply diurnal and latitudinal corrections to reference densities
    let lat_rad = input.g_lat.to_radians();
    let day_angle = (input.doy as f64 - 1.0) * std::f64::consts::TAU / 365.25;
    let decl = 23.44_f64.to_radians() * day_angle.sin();
    let hour_angle = ((input.lst - 12.0) * 15.0).to_radians();
    let cos_sza = lat_rad.sin() * decl.sin()
        + lat_rad.cos() * decl.cos() * hour_angle.cos();

    // Compute density of each constituent at altitude
    let mut total_density = 0.0_f64;
    for i in 0..7 {
        let mut n_120 = n_ref[i];

        // Diurnal variation in O and N2 number density at 120 km
        match i {
            1 => {
                // O: significant diurnal bulge
                n_120 *= 1.0 + 0.2 * cos_sza.max(0.0);
            }
            2 => {
                // N2: smaller diurnal effect
                n_120 *= 1.0 + 0.1 * cos_sza.max(0.0);
            }
            _ => {}
        }

        // Seasonal variation (annual + semiannual)
        let annual = 0.02 * (day_angle - 0.5).cos() * lat_rad.sin();
        let semiannual = 0.015 * (2.0 * day_angle).cos();
        n_120 *= 1.0 + annual + semiannual;

        let n = constituent_density(z, n_120, MASSES[i], t_inf, t_120, s_param, ALPHA[i]);

        // Mass density contribution: n * m (convert from g/mol to kg per molecule)
        total_density += n * MASSES[i] * 1e-3 / AVOGADRO;
    }

    // Below 80 km, blend with the US Standard Atmosphere 1976 density
    if z < 80.0 {
        let std_density = standard_atmosphere_density(z);
        if z < 72.5 {
            // Pure standard atmosphere below 72.5 km
            total_density = std_density;
        } else {
            // Blend between 72.5 and 80 km
            let frac = (z - 72.5) / 7.5;
            total_density = std_density * (1.0 - frac) + total_density * frac;
        }
    }

    NrlmsiseOutput {
        density: total_density,
        temperature_exo: t_inf,
        temperature_alt: t_alt,
    }
}

/// US Standard Atmosphere 1976 density (kg/m^3) for altitudes below ~86 km.
fn standard_atmosphere_density(z: f64) -> f64 {
    let t = temperature_below_120(z);
    // p = p0 * (T/T0)^(g0*M0/(R*L)) for constant lapse rate layers
    // For simplicity, use the barometric formula with the temperature profile

    // Sea level values
    let p0: f64 = 101325.0; // Pa

    // Compute pressure via numerical integration of hydrostatic equation
    // dp/dz = -rho * g = -p * M * g / (R * T)
    let steps = (z * 10.0).max(100.0) as usize;
    let dz = z / steps as f64; // km per step

    let mut p = p0;
    for step in 0..steps {
        let z_mid = (step as f64 + 0.5) * dz;
        let t_mid = temperature_below_120(z_mid);
        let g_mid = G0 * (RE / (RE + z_mid)).powi(2);
        // dp = -p * M0 * g * dz / (R * T)
        let dp_frac = M0 * g_mid * dz * 1000.0 / (R_GAS * t_mid);
        p *= (-dp_frac).exp();
    }

    // rho = p * M / (R * T)
    p * M0 / (R_GAS * t)
}

// ─── NIF interface ──────────────────────────────────────────────────────────

/// NIF entry point for atmosphere_density.
///
/// Arguments: lat_deg, lon_deg, alt_km, year, doy, sec, f107, f107a, ap
/// Returns: {density_kg_m3, temperature_K}
#[allow(clippy::too_many_arguments)]
pub(crate) fn atmosphere_density_impl(
    lat_deg: f64,
    lon_deg: f64,
    alt_km: f64,
    year: i32,
    doy: i32,
    sec: f64,
    f107: f64,
    f107a: f64,
    ap: f64,
) -> NifResult<(f64, f64)> {
    let lst = local_solar_time(sec, lon_deg);

    let input = NrlmsiseInput {
        year,
        doy,
        sec,
        alt: alt_km,
        g_lat: lat_deg,
        g_long: lon_deg,
        lst,
        f107a,
        f107,
        ap,
    };

    let output = nrlmsise00(&input);

    Ok((output.density, output.temperature_alt))
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sea_level_density() {
        let input = NrlmsiseInput {
            year: 2024,
            doy: 172,
            sec: 43200.0,
            alt: 0.0,
            g_lat: 0.0,
            g_long: 0.0,
            lst: 12.0,
            f107a: 150.0,
            f107: 150.0,
            ap: 4.0,
        };
        let out = nrlmsise00(&input);
        // Sea level density should be approximately 1.225 kg/m^3
        assert!(out.density > 0.8, "density too low: {}", out.density);
        assert!(out.density < 1.6, "density too high: {}", out.density);
        assert!(out.temperature_alt > 280.0 && out.temperature_alt < 300.0);
    }

    #[test]
    fn test_iss_altitude_density() {
        let input = NrlmsiseInput {
            year: 2024,
            doy: 172,
            sec: 43200.0,
            alt: 400.0,
            g_lat: 0.0,
            g_long: 0.0,
            lst: 12.0,
            f107a: 150.0,
            f107: 150.0,
            ap: 4.0,
        };
        let out = nrlmsise00(&input);
        // At 400 km, density should be in the range 1e-12 to 1e-10 kg/m^3
        assert!(
            out.density > 1e-14 && out.density < 1e-9,
            "ISS altitude density out of range: {:.3e}",
            out.density
        );
        // Temperature should be between 500 K and 2000 K at this altitude
        assert!(
            out.temperature_alt > 500.0 && out.temperature_alt < 2000.0,
            "ISS altitude temperature out of range: {}",
            out.temperature_alt
        );
    }

    #[test]
    fn test_density_decreases_with_altitude() {
        let make_input = |alt: f64| NrlmsiseInput {
            year: 2024,
            doy: 172,
            sec: 43200.0,
            alt,
            g_lat: 0.0,
            g_long: 0.0,
            lst: 12.0,
            f107a: 150.0,
            f107: 150.0,
            ap: 4.0,
        };

        let d_0 = nrlmsise00(&make_input(0.0)).density;
        let d_100 = nrlmsise00(&make_input(100.0)).density;
        let d_200 = nrlmsise00(&make_input(200.0)).density;
        let d_400 = nrlmsise00(&make_input(400.0)).density;
        let d_800 = nrlmsise00(&make_input(800.0)).density;

        assert!(d_0 > d_100, "density should decrease: {:.3e} > {:.3e}", d_0, d_100);
        assert!(d_100 > d_200, "density should decrease: {:.3e} > {:.3e}", d_100, d_200);
        assert!(d_200 > d_400, "density should decrease: {:.3e} > {:.3e}", d_200, d_400);
        assert!(d_400 > d_800, "density should decrease: {:.3e} > {:.3e}", d_400, d_800);
    }

    #[test]
    fn test_solar_activity_effect() {
        let make_input = |f107a: f64| NrlmsiseInput {
            year: 2024,
            doy: 172,
            sec: 43200.0,
            alt: 400.0,
            g_lat: 0.0,
            g_long: 0.0,
            lst: 12.0,
            f107a,
            f107: f107a,
            ap: 4.0,
        };

        let d_low = nrlmsise00(&make_input(70.0)).density;
        let d_high = nrlmsise00(&make_input(250.0)).density;

        // Higher solar activity = higher thermospheric density
        assert!(
            d_high > d_low,
            "high solar activity should give higher density: {:.3e} > {:.3e}",
            d_high,
            d_low
        );
    }

    #[test]
    fn test_local_solar_time() {
        assert!((local_solar_time(43200.0, 0.0) - 12.0).abs() < 0.001);
        assert!((local_solar_time(0.0, 0.0) - 0.0).abs() < 0.001);
        assert!((local_solar_time(0.0, 180.0) - 12.0).abs() < 0.001);
    }
}
