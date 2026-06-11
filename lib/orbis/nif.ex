checksum_file = Path.expand("../../checksum-Elixir.Orbis.NIF.exs", __DIR__)
version = Mix.Project.config()[:version]
source_checkout? = File.exists?(Path.expand("../../.git", __DIR__))

checksum_current? =
  File.exists?(checksum_file) and
    checksum_file |> File.read!() |> String.contains?("-v#{version}-")

force_build =
  System.get_env("ORBIS_BUILD") in ["1", "true"] or source_checkout? or not checksum_current?

defmodule Orbis.NIF do
  @moduledoc false

  use RustlerPrecompiled,
    otp_app: :orbis,
    crate: "orbis_nif",
    base_url: "https://github.com/neilberkman/orbis/releases/download/v#{version}",
    force_build: force_build,
    nif_versions: ["2.15"],
    targets: [
      "aarch64-apple-darwin",
      "aarch64-unknown-linux-gnu",
      "x86_64-apple-darwin",
      "x86_64-pc-windows-msvc",
      "x86_64-unknown-linux-gnu"
    ],
    version: version

  def propagate_with_elements(_tle_map, _datetime_tuple), do: :erlang.nif_error(:nif_not_loaded)

  def teme_to_gcrs(_x, _y, _z, _vx, _vy, _vz, _datetime, _skyfield_compat),
    do: :erlang.nif_error(:nif_not_loaded)

  def gcrs_to_itrs(_x, _y, _z, _datetime, _skyfield_compat),
    do: :erlang.nif_error(:nif_not_loaded)

  def itrs_to_geodetic(_x, _y, _z), do: :erlang.nif_error(:nif_not_loaded)

  def gcrs_to_topocentric(_sat_x, _sat_y, _sat_z, _lat, _lon, _alt, _datetime, _skyfield_compat),
    do: :erlang.nif_error(:nif_not_loaded)

  def atmosphere_density(_lat, _lon, _alt, _year, _doy, _sec, _f107, _f107a, _ap),
    do: :erlang.nif_error(:nif_not_loaded)

  def get_body_position(
        _file_path,
        _target_name,
        _observer_name,
        _jd_whole,
        _jd_fraction,
        _skyfield_compat
      ), do: :erlang.nif_error(:nif_not_loaded)

  def utc_to_tdb_jd(_year, _month, _day, _hour, _minute, _second),
    do: :erlang.nif_error(:nif_not_loaded)

  def utc_to_tdb_jd_split(_year, _month, _day, _hour, _minute, _second),
    do: :erlang.nif_error(:nif_not_loaded)

  def doppler_compute(
        _sat_x,
        _sat_y,
        _sat_z,
        _sat_vx,
        _sat_vy,
        _sat_vz,
        _station_lat_deg,
        _station_lon_deg,
        _station_alt_km,
        _datetime_tuple
      ), do: :erlang.nif_error(:nif_not_loaded)

  def iod_gibbs(_r1, _r2, _r3), do: :erlang.nif_error(:nif_not_loaded)
  def iod_hgibbs(_r1, _r2, _r3, _jd1, _jd2, _jd3), do: :erlang.nif_error(:nif_not_loaded)

  def lambert_battin(_r1, _r2, _v1, _dm, _de, _nrev, _dtsec),
    do: :erlang.nif_error(:nif_not_loaded)

  def find_conjunctions(_l1a, _l2a, _l1b, _l2b, _start, _end, _step, _threshold),
    do: :erlang.nif_error(:nif_not_loaded)

  def sp3_parse(_bytes), do: :erlang.nif_error(:nif_not_loaded)

  def sp3_time_scale(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def sp3_satellite_ids(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def sp3_to_iodata(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def broadcast_parse(_text), do: :erlang.nif_error(:nif_not_loaded)

  def broadcast_position(_handle, _system_letter, _prn, _t_j2000_s),
    do: :erlang.nif_error(:nif_not_loaded)

  def sp3_position(_handle, _system_letter, _prn, _scale, _jd_whole, _jd_fraction),
    do: :erlang.nif_error(:nif_not_loaded)

  def sp3_clock_reference_offset(_reference, _other, _min_common),
    do: :erlang.nif_error(:nif_not_loaded)

  def sp3_align_clock_reference(_reference, _other, _min_common),
    do: :erlang.nif_error(:nif_not_loaded)

  def sp3_merge(
        _handles,
        _position_tolerance_m,
        _clock_tolerance_s,
        _min_agree,
        _clock_min_common,
        _combine,
        _target_epoch_interval_s,
        _system_letters
      ), do: :erlang.nif_error(:nif_not_loaded)

  def crinex_decode(_text), do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_parse(_text), do: :erlang.nif_error(:nif_not_loaded)

  def crinex_obs_parse(_text), do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_approx_position(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_antenna_delta_hen(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_phase_shifts(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_codes(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_glonass_slots(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_epoch_count(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_epochs(_handle), do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_pseudoranges(_handle, _epoch_index, _overrides),
    do: :erlang.nif_error(:nif_not_loaded)

  def rinex_obs_values(_handle, _epoch_index, _overrides), do: :erlang.nif_error(:nif_not_loaded)

  def reduced_orbit_fit(_samples, _scale, _model), do: :erlang.nif_error(:nif_not_loaded)

  def ils_search(_float_cycles, _covariance, _radius, _candidate_limit, _ratio_threshold),
    do: :erlang.nif_error(:nif_not_loaded)

  def ils_lambda_search(_float_cycles, _covariance, _ratio_threshold),
    do: :erlang.nif_error(:nif_not_loaded)

  def rtk_filter_update_epoch(_state, _epoch, _base, _model, _wavelengths, _offsets, _opts),
    do: :erlang.nif_error(:nif_not_loaded)

  def reduced_orbit_position(_epoch, _scale, _elements, _query, _frame),
    do: :erlang.nif_error(:nif_not_loaded)

  def reduced_orbit_position_velocity(_epoch, _scale, _elements, _query, _frame),
    do: :erlang.nif_error(:nif_not_loaded)

  def reduced_orbit_drift(_epoch, _scale, _elements, _truth, _threshold_m),
    do: :erlang.nif_error(:nif_not_loaded)

  def spp_solve(
        _handle,
        _observations,
        _t_rx_j2000_s,
        _t_rx_second_of_day_s,
        _day_of_year,
        _initial_guess,
        _apply_iono,
        _apply_tropo,
        _alpha,
        _beta,
        _pressure_hpa,
        _temperature_k,
        _relative_humidity,
        _with_geodetic
      ), do: :erlang.nif_error(:nif_not_loaded)

  def spp_solve_broadcast(
        _handle,
        _observations,
        _t_rx_j2000_s,
        _t_rx_second_of_day_s,
        _day_of_year,
        _initial_guess,
        _apply_iono,
        _apply_tropo,
        _alpha,
        _beta,
        _pressure_hpa,
        _temperature_k,
        _relative_humidity,
        _with_geodetic
      ), do: :erlang.nif_error(:nif_not_loaded)

  def klobuchar_delay(
        _lat_deg,
        _lon_deg,
        _azimuth_deg,
        _elevation_deg,
        _t_gps_s,
        _frequency_hz,
        _alpha,
        _beta
      ), do: :erlang.nif_error(:nif_not_loaded)

  def ionex_parse(_bytes), do: :erlang.nif_error(:nif_not_loaded)

  def ionex_slant(
        _handle,
        _lat_rad,
        _lon_rad,
        _elevation_rad,
        _azimuth_rad,
        _epoch_j2000_s,
        _frequency_hz
      ), do: :erlang.nif_error(:nif_not_loaded)

  def tropo_zenith_delay(_lat_rad, _height_m, _pressure_hpa, _temperature_k, _relative_humidity),
    do: :erlang.nif_error(:nif_not_loaded)

  def tropo_mapping_factors(_elevation_rad, _lat_rad, _height_m, _jd_whole, _jd_fraction),
    do: :erlang.nif_error(:nif_not_loaded)

  def tropo_slant_delay(
        _elevation_rad,
        _lat_rad,
        _lon_rad,
        _height_m,
        _pressure_hpa,
        _temperature_k,
        _relative_humidity,
        _jd_whole,
        _jd_fraction
      ), do: :erlang.nif_error(:nif_not_loaded)

  def iod_gauss(
        _decl1,
        _decl2,
        _decl3,
        _rtasc1,
        _rtasc2,
        _rtasc3,
        _jd1,
        _jdf1,
        _jd2,
        _jdf2,
        _jd3,
        _jdf3,
        _rseci1,
        _rseci2,
        _rseci3
      ), do: :erlang.nif_error(:nif_not_loaded)
end
