defmodule Orbis.CrinexTest do
  use ExUnit.Case, async: true

  alias Orbis.RinexObs

  @crx_path Path.join(__DIR__, "fixtures/obs/ESBC00DNK_R_20201770000_01D_30S_MO_trim.crx")
  @rnx_path Path.join(__DIR__, "fixtures/obs/ESBC00DNK_R_20201770000_01D_30S_MO_trim.rnx")

  # CRINEX 1.0 / RINEX 2 fixture (mixed GPS+GLONASS, 20 satellites, 8 obs types).
  @crx_v1_path Path.join(__DIR__, "fixtures/obs/algo0010_2015001_v1_trim.crx")
  @rnx_v1_path Path.join(__DIR__, "fixtures/obs/algo0010_2015001_v1_trim.rnx")

  describe "decode_crinex/1" do
    test "reproduces the crx2rnx reference decode byte-for-byte" do
      crx = File.read!(@crx_path)
      reference = File.read!(@rnx_path)

      assert {:ok, decoded} = RinexObs.decode_crinex(crx)

      # Compare line by line so a mismatch points at the offending record.
      decoded_lines = String.split(decoded, "\n", trim: false)
      reference_lines = String.split(reference, "\n", trim: false)

      # The reference file's trailing newline yields a final empty element; the
      # decoder collects newline-terminated lines, matching after trimming the
      # common trailing blank.
      decoded_lines = Enum.reverse(Enum.drop_while(Enum.reverse(decoded_lines), &(&1 == "")))
      reference_lines = Enum.reverse(Enum.drop_while(Enum.reverse(reference_lines), &(&1 == "")))

      assert length(decoded_lines) == length(reference_lines)

      for {d, r} <- Enum.zip(decoded_lines, reference_lines) do
        assert d == r
      end
    end

    test "decodes CRINEX 1.0 (RINEX 2) byte-for-byte" do
      crx = File.read!(@crx_v1_path)
      reference = File.read!(@rnx_v1_path)

      assert {:ok, decoded} = RinexObs.decode_crinex(crx)

      decoded_lines = String.split(decoded, "\n", trim: false)
      reference_lines = String.split(reference, "\n", trim: false)
      decoded_lines = Enum.reverse(Enum.drop_while(Enum.reverse(decoded_lines), &(&1 == "")))
      reference_lines = Enum.reverse(Enum.drop_while(Enum.reverse(reference_lines), &(&1 == "")))

      assert length(decoded_lines) == length(reference_lines)

      for {d, r} <- Enum.zip(decoded_lines, reference_lines) do
        assert d == r
      end
    end

    test "rejects non-CRINEX input" do
      assert {:error, _} = RinexObs.decode_crinex("not a crinex file\n")
    end
  end

  describe "parse_crinex/1 and parse_auto/1" do
    test "parses the CRINEX fixture into a handle" do
      assert {:ok, %RinexObs{}} = RinexObs.parse_crinex(File.read!(@crx_path))
    end

    test "auto-detects CRINEX vs plain RINEX" do
      assert {:ok, %RinexObs{} = from_crx} = RinexObs.parse_auto(File.read!(@crx_path))
      assert {:ok, %RinexObs{} = from_rnx} = RinexObs.parse_auto(File.read!(@rnx_path))

      # Both paths recover the same surveyed position.
      assert RinexObs.approx_position(from_crx) == RinexObs.approx_position(from_rnx)
    end

    test "load/1 decodes a .crx file from disk" do
      assert {:ok, %RinexObs{}} = RinexObs.load(@crx_path)
    end
  end

  describe "header accessors" do
    setup do
      {:ok, obs} = RinexObs.load(@rnx_path)
      {:ok, obs: obs}
    end

    test "approx_position/1 returns the surveyed ECEF position", %{obs: obs} do
      {x, y, z} = RinexObs.approx_position(obs)
      assert_in_delta x, 3_582_105.291, 1.0e-3
      assert_in_delta y, 532_589.7313, 1.0e-3
      assert_in_delta z, 5_232_754.8054, 1.0e-3
    end

    test "observation_codes/1 returns per-system code lists in order", %{obs: obs} do
      codes = RinexObs.observation_codes(obs)
      assert hd(codes["G"]) == "C1C"
      assert length(codes["G"]) == 18
      assert hd(codes["C"]) == "C2I"
    end

    test "epochs/1 lists the two epochs with indices and counts", %{obs: obs} do
      epochs = RinexObs.epochs(obs)
      assert length(epochs) == 2
      assert [%{index: 0, flag: 0, sat_count: 43} | _] = epochs
      assert {{2020, 6, 25}, {0, 0, _}} = hd(epochs).epoch
    end
  end
end
