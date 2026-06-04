defmodule Orbis.GnssData.CatalogTest do
  use ExUnit.Case, async: true

  alias Orbis.GnssData.Catalog

  doctest Orbis.GnssData.Catalog

  describe "GPS-week and day-of-year arithmetic" do
    test "worked example: 2020-06-24 is GPS week 2111, day 3, doy 176" do
      d = ~D[2020-06-24]
      assert Catalog.gps_week(d) == 2111
      assert Catalog.gps_day_of_week(d) == 3
      assert Catalog.day_of_year(d) == 176
    end

    test "GPS epoch day itself is week 0, day 0" do
      d = ~D[1980-01-06]
      assert Catalog.gps_week(d) == 0
      assert Catalog.gps_day_of_week(d) == 0
    end

    test "the day before the next GPS week rolls the day-of-week to 6" do
      assert Catalog.gps_day_of_week(~D[1980-01-12]) == 6
      assert Catalog.gps_week(~D[1980-01-12]) == 0
      assert Catalog.gps_week(~D[1980-01-13]) == 1
    end

    test "day-of-year handles the leap day and year boundaries" do
      assert Catalog.day_of_year(~D[2020-01-01]) == 1
      # 2020 is a leap year: Dec 31 is the 366th day.
      assert Catalog.day_of_year(~D[2020-12-31]) == 366
      # 2021 is not a leap year: Dec 31 is the 365th day.
      assert Catalog.day_of_year(~D[2021-12-31]) == 365
    end
  end

  describe "canonical_filename/4" do
    test "encodes the IGS long-name for several centers and content types" do
      assert {:ok, "GFZ0OPSRAP_20201760000_01D_15M_ORB.SP3"} =
               Catalog.canonical_filename(:gfz, :sp3, ~D[2020-06-24], "15M")

      assert {:ok, "COD0MGXFIN_20201760000_01D_05M_ORB.SP3"} =
               Catalog.canonical_filename(:cod, :sp3, ~D[2020-06-24], "05M")

      assert {:ok, "WUM0MGXFIN_20201760000_01D_30S_CLK.CLK"} =
               Catalog.canonical_filename(:wum, :clk, ~D[2020-06-24], "30S")
    end

    test "navigation uses the no-sample RINEX long-name (real BRDC00IGS product)" do
      # Matches the real merged-nav product served by IGS/CDDIS/ESA: no SMP
      # field, lowercase .rnx extension, single-char source code "R".
      assert {:ok, "BRDC00IGS_R_20201770000_01D_MN.rnx"} =
               Catalog.canonical_filename(:igs, :nav, ~D[2020-06-25], "01D")
    end

    test "IONEX uses the GIM content code (real combined product)" do
      assert {:ok, "IGS0OPSFIN_20241760000_01D_01H_GIM.INX"} =
               Catalog.canonical_filename(:igs, :ionex, ~D[2024-06-24], "01H")

      assert {:ok, "COD0OPSFIN_20241760000_01D_01H_GIM.INX"} =
               Catalog.canonical_filename(:cod, :ionex, ~D[2024-06-24], "01H")
    end

    test "a center that does not publish a content type is rejected" do
      # GFZ operational serves SP3/CLK but not broadcast nav or IONEX.
      assert {:error, {:unsupported_product, {:content_not_served, :nav}}} =
               Catalog.canonical_filename(:gfz, :nav, ~D[2020-06-25], "01D")

      # IGS serves nav/IONEX but not precise orbits here.
      assert {:error, {:unsupported_product, {:content_not_served, :sp3}}} =
               Catalog.canonical_filename(:igs, :sp3, ~D[2020-06-24], "05M")
    end

    test "zero-pads the day-of-year to three digits" do
      assert {:ok, name} = Catalog.canonical_filename(:cod, :sp3, ~D[2020-01-01], "05M")
      assert name =~ "_20200010000_"
    end

    test "rejects an unknown center" do
      assert {:error, {:unsupported_product, {:center, :nope}}} =
               Catalog.canonical_filename(:nope, :sp3, ~D[2020-06-24], "05M")
    end

    test "rejects an unknown content type" do
      assert {:error, {:unsupported_product, {:content, :bogus}}} =
               Catalog.canonical_filename(:cod, :bogus, ~D[2020-06-24], "05M")
    end

    test "rejects a malformed sampling code" do
      assert {:error, {:unsupported_product, {:sample, "5M"}}} =
               Catalog.canonical_filename(:cod, :sp3, ~D[2020-06-24], "5M")
    end
  end

  describe "archive_url/4" do
    test "builds the GFZ HTTPS URL with the class/week directory and .gz suffix" do
      assert {:ok,
              "https://isdc-data.gfz.de/gnss/products/rapid/w2111/GFZ0OPSRAP_20201760000_01D_15M_ORB.SP3.gz"} =
               Catalog.archive_url(:gfz, :sp3, ~D[2020-06-24], "15M")
    end

    test "builds an FTP URL for an MGEX precise product on the GSSC mirror" do
      assert {:ok, url} = Catalog.archive_url(:wum, :sp3, ~D[2020-06-24], "05M")

      assert url ==
               "ftp://gssc.esa.int/gnss/products/2111/WUM0MGXFIN_20201760000_01D_05M_ORB.SP3.gz"
    end

    test "builds the GSSC daily-tree URL for broadcast navigation" do
      assert {:ok,
              "ftp://gssc.esa.int/gnss/data/daily/2020/177/BRDC00IGS_R_20201770000_01D_MN.rnx.gz"} =
               Catalog.archive_url(:igs, :nav, ~D[2020-06-25], "01D")
    end

    test "builds the GSSC ionex-tree URL for a global ionosphere map" do
      assert {:ok,
              "ftp://gssc.esa.int/gnss/products/ionex/2024/176/IGS0OPSFIN_20241760000_01D_01H_GIM.INX.gz"} =
               Catalog.archive_url(:igs, :ionex, ~D[2024-06-24], "01H")
    end

    test "propagates an unsupported-product error" do
      assert {:error, {:unsupported_product, _}} =
               Catalog.archive_url(:nope, :sp3, ~D[2020-06-24], "05M")
    end
  end

  describe "protocol/1 and allowed_hosts/0" do
    test "maps centers to their transfer protocol" do
      assert {:ok, :https} = Catalog.protocol(:gfz)
      assert {:ok, :ftp} = Catalog.protocol(:wum)
      assert {:error, {:unsupported_product, _}} = Catalog.protocol(:nope)
    end

    test "allowed hosts contain every catalog host and nothing else slips in" do
      hosts = Catalog.allowed_hosts()
      assert MapSet.member?(hosts, "isdc-data.gfz.de")
      assert MapSet.member?(hosts, "gssc.esa.int")
      refute MapSet.member?(hosts, "evil.example.com")
    end
  end
end
