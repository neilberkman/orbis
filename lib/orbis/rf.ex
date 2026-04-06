defmodule Orbis.RF do
  @moduledoc """
  RF link budget primitives.

  Pure physics calculations with no system-specific assumptions.
  Combine these with Orbis geometry outputs (slant range, elevation)
  to build a complete link budget for your specific system.

  ## Example: simple uplink budget

      # Get geometry from Orbis
      {:ok, elements} = Orbis.Format.TLE.parse(line1, line2)
      {:ok, look} = Orbis.look_angle(elements, datetime, station)

      # Compute path loss
      fspl = Orbis.RF.fspl(look.range_km, 1616.0)

      # Your system parameters
      eirp_dbw = 27.0 + 3.0 - 30.0  # tx power + antenna gain - 30
      gt_dbk = -12.0                  # satellite G/T
      other_losses = 3.0              # atmospheric, polarization, etc.

      # Link margin
      margin = Orbis.RF.link_margin(%{
        eirp_dbw: eirp_dbw,
        fspl_db: fspl,
        receiver_gt_dbk: gt_dbk,
        other_losses_db: other_losses,
        required_cn0_dbhz: 35.0
      })

  """

  # Boltzmann constant: -228.6 dBW/Hz/K
  @boltzmann_k_dbw_hz_k 228.6

  @doc """
  Free-space path loss in dB.

  This is the inverse square law — signal attenuation over distance
  in vacuum. The foundational calculation in any link budget.

      FSPL = 32.45 + 20·log₁₀(f_MHz) + 20·log₁₀(d_km)

  ## Parameters

    * `distance_km` — slant range to satellite in km
    * `frequency_mhz` — carrier frequency in MHz

  ## Examples

      iex> Orbis.RF.fspl(1200.0, 1616.0)
      158.20245204972383

  """
  @spec fspl(float(), float()) :: float()
  def fspl(distance_km, frequency_mhz) do
    32.45 + 20 * :math.log10(frequency_mhz) + 20 * :math.log10(distance_km)
  end

  @doc """
  Effective Isotropic Radiated Power in dBW.

      EIRP = P_tx(dBm) + G_tx(dBi) - 30

  ## Examples

      iex> Orbis.RF.eirp(27.0, 3.0)
      0.0

  """
  @spec eirp(float(), float()) :: float()
  def eirp(tx_power_dbm, tx_antenna_gain_dbi) do
    tx_power_dbm + tx_antenna_gain_dbi - 30.0
  end

  @doc """
  Carrier-to-noise-density ratio (C/N₀) in dB-Hz.

      C/N₀ = EIRP + G/T - FSPL + 228.6 - other_losses

  The 228.6 is the Boltzmann constant expressed as a positive number
  in the conventional link budget equation.

  ## Parameters

    * `eirp_dbw` — transmitter EIRP in dBW
    * `fspl_db` — free-space path loss in dB (from `fspl/2`)
    * `receiver_gt_dbk` — receiver figure of merit (G/T) in dB/K
    * `other_losses_db` — sum of all other losses (atmospheric, polarization, pointing, etc.)

  ## Examples

      iex> Orbis.RF.cn0(0.0, 165.0, -12.0, 3.0)
      48.599999999999994

  """
  @spec cn0(float(), float(), float(), float()) :: float()
  def cn0(eirp_dbw, fspl_db, receiver_gt_dbk, other_losses_db) do
    eirp_dbw + receiver_gt_dbk - fspl_db + @boltzmann_k_dbw_hz_k - other_losses_db
  end

  @doc """
  Link margin in dB.

  Positive margin means the link closes. Negative means it doesn't.

  Takes a map so parameters are self-documenting:

  ## Parameters

    * `:eirp_dbw` — transmitter EIRP in dBW
    * `:fspl_db` — free-space path loss in dB
    * `:receiver_gt_dbk` — receiver G/T in dB/K
    * `:other_losses_db` — sum of miscellaneous losses in dB
    * `:required_cn0_dbhz` — minimum C/N₀ for demodulation in dB-Hz

  ## Examples

      iex> Orbis.RF.link_margin(%{
      ...>   eirp_dbw: 0.0,
      ...>   fspl_db: 165.0,
      ...>   receiver_gt_dbk: -12.0,
      ...>   other_losses_db: 3.0,
      ...>   required_cn0_dbhz: 35.0
      ...> })
      13.599999999999994

  """
  @spec link_margin(map()) :: float()
  def link_margin(%{
        eirp_dbw: eirp_dbw,
        fspl_db: fspl_db,
        receiver_gt_dbk: gt,
        other_losses_db: losses,
        required_cn0_dbhz: required
      }) do
    cn0(eirp_dbw, fspl_db, gt, losses) - required
  end

  @doc """
  Wavelength in meters for a given frequency.

  ## Examples

      iex> Orbis.RF.wavelength(1616.0e6)
      0.1855151349009901

  """
  @spec wavelength(float()) :: float()
  def wavelength(frequency_hz) do
    299_792_458.0 / frequency_hz
  end

  @doc """
  Antenna gain in dBi for a parabolic dish.

      G = 10·log₁₀(η · (π·D/λ)²)

  ## Parameters

    * `diameter_m` — dish diameter in meters
    * `frequency_hz` — frequency in Hz
    * `efficiency` — aperture efficiency (default 0.55)

  ## Examples

      iex> Orbis.RF.dish_gain(1.0, 1616.0e6)
      21.97903741903791

  """
  @spec dish_gain(float(), float(), float()) :: float()
  def dish_gain(diameter_m, frequency_hz, efficiency \\ 0.55) do
    lambda = wavelength(frequency_hz)
    10 * :math.log10(efficiency * (:math.pi() * diameter_m / lambda) ** 2)
  end
end
