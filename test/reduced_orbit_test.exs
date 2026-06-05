defmodule Orbis.ReducedOrbitTest do
  use ExUnit.Case, async: true

  alias Orbis.ReducedOrbit

  @sp3_path Path.join(__DIR__, "fixtures/sp3/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3")

  # A near-circular Galileo MEO satellite: the circular_secular model is a good
  # approximation here (GPS, with e ~ 1e-2, is deliberately a poor fit — see the
  # module docs).
  @sat "E01"
  @t0 ~N[2020-06-24 00:00:00]
  @t1 ~N[2020-06-24 06:00:00]
  @day_end ~N[2020-06-25 00:00:00]

  setup do
    sp3 = Orbis.SP3.load!(@sp3_path)
    {:ok, model} = ReducedOrbit.fit(sp3, satellite_id: @sat, window: {@t0, @t1}, cadence_s: 900)
    {:ok, sp3: sp3, model: model}
  end

  describe "fit/2 from SP3" do
    test "fits a near-circular MEO satellite to a few-km residual", %{model: m} do
      assert m.model == "circular_secular"
      assert m.frame == "GCRS"
      assert m.version == 1
      # ~29 600 km Galileo semi-major axis, ~56 deg inclination.
      assert_in_delta m.a_m / 1000.0, 29_600.0, 100.0
      assert_in_delta m.i_rad * 180.0 / :math.pi(), 56.0, 1.0
      # In-window residual is a few km (short-period + unmodelled secular terms).
      assert m.fit.rms_m < 5_000.0
      assert m.fit.max_m < 10_000.0
      assert m.fit.n_samples >= 20
      assert m.fit.source == "sp3:#{@sat}"
      # Epochs are interpreted in the product's own scale, not silently UTC.
      assert m.time_scale == "GPST"
    end

    test "the nodal rate is fitted off the J2 seed, both finite", %{model: m} do
      assert is_float(m.raan_rate_rad_s)
      assert is_float(m.raan_rate_j2_rad_s)
      assert m.raan_rate_mode == "fitted_j2_seeded"
      # Fitted value tracks the real drift, distinct from the pure-J2 seed.
      assert m.raan_rate_rad_s != m.raan_rate_j2_rad_s
    end

    test "recovers the SP3 position across the window to a few km", %{sp3: sp3, model: m} do
      for offset <- [0, 3600, 10_800, 21_600] do
        epoch = NaiveDateTime.add(@t0, offset, :second)
        {:ok, truth} = Orbis.SP3.position(sp3, @sat, epoch)
        {:ok, pos} = ReducedOrbit.position(m, epoch)

        err =
          :math.sqrt(
            :math.pow(pos.x_m - truth.x_m, 2) +
              :math.pow(pos.y_m - truth.y_m, 2) +
              :math.pow(pos.z_m - truth.z_m, 2)
          )

        assert err < 10_000.0, "epoch +#{offset}s error #{Float.round(err / 1000, 2)} km"
      end
    end
  end

  describe "position/3 frames and velocity" do
    test "ecef and gcrs differ but share the orbital radius", %{model: m} do
      epoch = NaiveDateTime.add(@t0, 5400, :second)
      {:ok, ecef} = ReducedOrbit.position(m, epoch, frame: :ecef)
      {:ok, gcrs} = ReducedOrbit.position(m, epoch, frame: :gcrs)

      radius = fn p -> :math.sqrt(p.x_m ** 2 + p.y_m ** 2 + p.z_m ** 2) end
      assert_in_delta radius.(ecef), m.a_m, 1.0
      assert_in_delta radius.(gcrs), m.a_m, 1.0
      # Earth rotation separates the two frames by thousands of km.
      assert abs(ecef.x_m - gcrs.x_m) > 1_000_000.0
    end

    test "position_velocity returns a circular-orbit speed", %{model: m} do
      epoch = NaiveDateTime.add(@t0, 5400, :second)
      {:ok, %{position: p, velocity: v}} = ReducedOrbit.position_velocity(m, epoch, frame: :gcrs)
      assert_in_delta :math.sqrt(p.x_m ** 2 + p.y_m ** 2 + p.z_m ** 2), m.a_m, 1.0
      # A circular orbit's speed is exactly a*n (~3.67 km/s at the Galileo radius).
      speed = :math.sqrt(v.vx_m_s ** 2 + v.vy_m_s ** 2 + v.vz_m_s ** 2)
      assert_in_delta speed, m.a_m * m.mean_motion_rad_s, 1.0
    end
  end

  describe "drift/3 (source-backed)" do
    test "reports bounded error over a full day against SP3", %{sp3: sp3, model: m} do
      {:ok, d} =
        ReducedOrbit.drift(m, sp3,
          satellite_id: @sat,
          window: {@t0, @day_end},
          cadence_s: 1800,
          threshold_m: 100_000.0
        )

      assert length(d.per_epoch) >= 40
      assert is_float(d.max_m) and d.max_m > 0.0
      # Near-circular drift stays well under 20 km over the day.
      assert d.max_m < 20_000.0
      assert d.rms_m < 15_000.0
      # Never crosses 100 km, so no horizon.
      assert d.threshold_horizon == nil
      assert %{epoch: %NaiveDateTime{}, error_m: e} = hd(d.per_epoch)
      assert is_float(e)
    end

    test "reports the first epoch error crosses a tight threshold", %{sp3: sp3, model: m} do
      {:ok, d} =
        ReducedOrbit.drift(m, sp3,
          satellite_id: @sat,
          window: {@t0, @day_end},
          cadence_s: 1800,
          threshold_m: 1000.0
        )

      assert %NaiveDateTime{} = d.threshold_horizon
    end
  end

  describe "to_map/1 and from_map/1" do
    test "round-trips a fitted model", %{model: m} do
      map = ReducedOrbit.to_map(m)
      assert map["version"] == 1
      assert map["model"] == "circular_secular"
      assert map["frame"] == "GCRS"
      assert map["units"]["length"] == "m"
      assert map["elements"]["raan_rate_mode"] == "fitted_j2_seeded"
      assert is_float(map["elements"]["raan_rate_j2_rad_s"])

      assert {:ok, back} = ReducedOrbit.from_map(map)
      assert back.a_m == m.a_m
      assert back.i_rad == m.i_rad
      assert back.raan_rate_rad_s == m.raan_rate_rad_s
      assert back.epoch == m.epoch
      assert back.fit.rms_m == m.fit.rms_m
    end

    test "the map is JSON-encodable and survives a JSON round-trip", %{model: m} do
      # No tuples may leak out of to_map, and from_map must read string-keyed
      # JSON-decoded maps.
      decoded = m |> ReducedOrbit.to_map() |> Jason.encode!() |> Jason.decode!()
      assert {:ok, back} = ReducedOrbit.from_map(decoded)
      assert_in_delta back.a_m, m.a_m, 1.0e-6
      assert_in_delta back.raan_rate_rad_s, m.raan_rate_rad_s, 1.0e-18
      assert back.raan_rate_mode == "fitted_j2_seeded"
      assert back.epoch == m.epoch
    end

    test "rejects an unknown version and model" do
      assert {:error, {:unsupported_version, 99}} =
               ReducedOrbit.from_map(%{"version" => 99, "model" => "circular_secular"})

      assert {:error, {:unsupported_model, "wat"}} =
               ReducedOrbit.from_map(%{"version" => 1, "model" => "wat"})
    end
  end

  describe "tagged errors" do
    test "too few samples", %{sp3: sp3} do
      base = ~N[2020-06-24 00:00:00]

      samples =
        for k <- 0..2 do
          epoch = NaiveDateTime.add(base, k * 900, :second)
          {:ok, st} = Orbis.SP3.position(sp3, @sat, epoch)
          {epoch, {st.x_m, st.y_m, st.z_m}}
        end

      assert {:error, {:too_few_samples, 3, 4}} = ReducedOrbit.fit(samples)
    end

    test "inverted window", %{sp3: sp3} do
      assert {:error, :invalid_window} =
               ReducedOrbit.fit(sp3, satellite_id: @sat, window: {@t1, @t0})
    end

    test "missing satellite id", %{sp3: sp3} do
      assert {:error, :satellite_id_required} =
               ReducedOrbit.fit(sp3, window: {@t0, @t1})
    end

    test "unsupported source frame on a sample list" do
      samples = [{~N[2020-06-24 00:00:00], {1.0, 2.0, 3.0}}]

      assert {:error, {:unsupported_source_frame, :gcrs}} =
               ReducedOrbit.fit(samples, frame: :gcrs)
    end

    test "drift rejects a negative threshold", %{model: m} do
      truth = [{~N[2020-06-24 00:00:00], {1.0, 2.0, 3.0}}]
      assert {:error, :invalid_threshold} = ReducedOrbit.drift(m, truth, threshold_m: -1.0)
    end

    test "a window with no time span surfaces the crate's :invalid_window" do
      # Four samples at a single epoch span no time; the crate rejects the fit and
      # the atom-form error flows back through the NIF.
      ep = ~N[2020-06-24 00:00:00]
      samples = for _ <- 0..3, do: {ep, {26_560_000.0, 0.0, 0.0}}
      assert {:error, :invalid_window} = ReducedOrbit.fit(samples)
    end

    test "fit rejects a non-positive cadence", %{sp3: sp3} do
      assert {:error, :invalid_cadence} =
               ReducedOrbit.fit(sp3, satellite_id: @sat, window: {@t0, @t1}, cadence_s: 0)
    end

    test "drift rejects a model whose scale differs from the SP3 product", %{sp3: sp3} do
      # A model fit from a UTC sample list must not be drifted against a GPST SP3.
      samples = sample_list(sp3)
      {:ok, utc_model} = ReducedOrbit.fit(samples)
      assert utc_model.time_scale == "UTC"

      assert {:error, {:time_scale_mismatch, "UTC", "GPST"}} =
               ReducedOrbit.drift(utc_model, sp3, satellite_id: @sat, window: {@t0, @t1})
    end
  end

  describe "sample-list ordering" do
    test "the fitted model is independent of sample order", %{sp3: sp3} do
      samples = sample_list(sp3)
      {:ok, ordered} = ReducedOrbit.fit(samples)
      {:ok, reversed} = ReducedOrbit.fit(Enum.reverse(samples))

      # t0 is the earliest sample regardless of input order, so both agree.
      assert NaiveDateTime.compare(ordered.epoch, @t0) == :eq
      assert NaiveDateTime.compare(reversed.epoch, @t0) == :eq
      assert_in_delta ordered.a_m, reversed.a_m, 1.0e-6
      assert_in_delta ordered.raan_rate_rad_s, reversed.raan_rate_rad_s, 1.0e-18
    end
  end

  # A short list of real ECEF samples for the near-circular satellite.
  defp sample_list(sp3) do
    for k <- 0..10 do
      ep = NaiveDateTime.add(@t0, k * 900, :second)
      {:ok, st} = Orbis.SP3.position(sp3, @sat, ep)
      {ep, {st.x_m, st.y_m, st.z_m}}
    end
  end
end
