defmodule Orbis.NIF do
  @moduledoc false

  use Rustler,
    otp_app: :orbis,
    crate: "orbis_nif"

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
