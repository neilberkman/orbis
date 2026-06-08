defmodule Orbis.GNSS.Data.ProductTest do
  use ExUnit.Case, async: true

  alias Orbis.GNSS.Data
  alias Orbis.GNSS.Data.Product

  doctest Orbis.GNSS.Data

  describe "new/4" do
    test "builds a valid product" do
      assert {:ok, %Product{center: :cod, content: :sp3, sample: "05M"} = p} =
               Product.new(:cod, :sp3, ~D[2020-06-24], "05M")

      assert p.date == ~D[2020-06-24]
    end

    test "rejects an invalid center/content/sample" do
      assert {:error, {:unsupported_product, _}} =
               Product.new(:nope, :sp3, ~D[2020-06-24], "05M")

      assert {:error, {:unsupported_product, _}} =
               Product.new(:cod, :bogus, ~D[2020-06-24], "05M")

      assert {:error, {:unsupported_product, _}} =
               Product.new(:cod, :sp3, ~D[2020-06-24], "bad")
    end
  end

  describe "resolution helpers" do
    setup do
      {:ok, p} = Product.new(:cod, :sp3, ~D[2020-06-24], "05M")
      {:ok, product: p}
    end

    test "canonical_filename/1, gps_week/1, day_of_year/1", %{product: p} do
      assert {:ok, "COD0MGXFIN_20201760000_01D_05M_ORB.SP3"} = Product.canonical_filename(p)
      assert Product.gps_week(p) == 2111
      assert Product.day_of_year(p) == 176
    end

    test "archive_url/1", %{product: p} do
      assert {:ok, url} = Product.archive_url(p)
      assert url =~ "/2111/COD0MGXFIN_20201760000_01D_05M_ORB.SP3.gz"
    end

    test "describe/1 returns the canonical name", %{product: p} do
      assert Product.describe(p) == "COD0MGXFIN_20201760000_01D_05M_ORB.SP3"
    end
  end

  describe "builders are pure and deterministic" do
    test "mgex_sp3/2 defaults to 5-minute sampling" do
      p = Data.mgex_sp3(:cod, ~D[2020-06-24])
      assert p == Data.mgex_sp3(:cod, ~D[2020-06-24])
      assert p.sample == "05M"
      assert p.content == :sp3
    end

    test "mgex_clk/2, mgex_nav/2, mgex_ionex/2 use sensible defaults" do
      assert Data.mgex_clk(:wum, ~D[2020-06-24]).sample == "30S"
      assert Data.mgex_nav(:igs, ~D[2020-06-25]).content == :nav
      assert Data.mgex_ionex(:igs, ~D[2024-06-24]).content == :ionex
      assert Data.mgex_ionex(:igs, ~D[2024-06-24]).sample == "01H"
    end

    test "ops_ultra_sp3/3 uses issue time and per-center sample defaults" do
      igs = Data.ops_ultra_sp3(:igs_ult, ~D[2024-09-03], issue: "0600")
      cod = Data.ops_ultra_sp3(:cod_ult, ~D[2024-09-03], issue: "0600")

      assert igs.sample == "15M"
      assert igs.issue == "0600"
      assert cod.sample == "05M"

      assert {:ok, "IGS0OPSULT_20242470600_02D_15M_ORB.SP3"} =
               Product.canonical_filename(igs)

      assert {:ok, "COD0OPSULT_20242470600_02D_05M_ORB.SP3"} =
               Product.canonical_filename(cod)
    end

    test "ops_ultra_sp3/3 resolves latest available issue for a target timestamp" do
      available = [{~D[2024-09-03], "0000"}, {~D[2024-09-03], "0600"}]
      p = Data.ops_ultra_sp3(:cod_ult, ~N[2024-09-03 13:00:00], available_issues: available)

      assert p.date == ~D[2024-09-03]
      assert p.issue == "0600"

      assert {:ok, "COD0OPSULT_20242470600_02D_05M_ORB.SP3"} =
               Product.canonical_filename(p)
    end

    test "ops_ultra_clk/3 builds the cataloged ultra clock product" do
      p = Data.ops_ultra_clk(:grg_ult, ~D[2024-09-03], issue: "0600")

      assert p.content == :clk
      assert p.sample == "05M"
      assert {:ok, "GRG0OPSULT_20242470600_02D_05M_CLK.CLK"} = Product.canonical_filename(p)
    end

    test "product/5 accepts product-specific issue options" do
      assert {:ok, p} = Data.product(:igs_ult, :sp3, ~D[2024-09-03], "15M", issue: "0600")
      assert p.issue == "0600"
      assert {:ok, "IGS0OPSULT_20242470600_02D_15M_ORB.SP3"} = Product.canonical_filename(p)
    end

    test "the sample override is honored" do
      assert Data.mgex_sp3(:cod, ~D[2020-06-24], sample: "15M").sample == "15M"

      assert Data.ops_ultra_sp3(:igs_ult, ~D[2024-09-03], issue: "0600", sample: "05M").sample ==
               "05M"
    end

    test "a bang builder raises on an invalid product" do
      assert_raise ArgumentError, fn -> Data.mgex_sp3(:nope, ~D[2020-06-24]) end
    end
  end
end
