defmodule Orbis.GNSS.Core.Constants do
  @moduledoc false

  @speed_of_light_m_s 299_792_458.0
  @earth_rotation_rate_rad_s 7.2921151467e-5
  @gps_l1_hz 1_575_420_000.0
  @gps_l2_hz 1_227_600_000.0
  @galileo_e1_hz 1_575_420_000.0
  @galileo_e5a_hz 1_176_450_000.0
  @beidou_b1i_hz 1_561_098_000.0
  @beidou_b3i_hz 1_268_520_000.0
  @ca_code_length 1023
  @ca_chip_rate_hz 1_023_000
  @time_scales ~w(UTC TAI TT TDB GPST GST BDT)

  def speed_of_light_m_s, do: @speed_of_light_m_s
  def earth_rotation_rate_rad_s, do: @earth_rotation_rate_rad_s
  def gps_l1_hz, do: @gps_l1_hz
  def gps_l2_hz, do: @gps_l2_hz
  def galileo_e1_hz, do: @galileo_e1_hz
  def galileo_e5a_hz, do: @galileo_e5a_hz
  def beidou_b1i_hz, do: @beidou_b1i_hz
  def beidou_b3i_hz, do: @beidou_b3i_hz
  def ca_code_length, do: @ca_code_length
  def ca_chip_rate_hz, do: @ca_chip_rate_hz
  def time_scales, do: @time_scales

  def carrier_frequencies_hz do
    %{
      "G" => %{l1: @gps_l1_hz, l2: @gps_l2_hz},
      "E" => %{e1: @galileo_e1_hz, e5a: @galileo_e5a_hz},
      "C" => %{b1i: @beidou_b1i_hz, b3i: @beidou_b3i_hz}
    }
  end
end
