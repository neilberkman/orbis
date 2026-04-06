defmodule Orbis.MixProject do
  use Mix.Project

  @version "0.5.0"
  @source_url "https://github.com/neilberkman/orbis"

  def project do
    [
      app: :orbis,
      version: @version,
      elixir: "~> 1.18",
      start_permanent: Mix.env() == :prod,
      deps: deps(),
      name: "Orbis",
      description: description(),
      package: package(),
      docs: docs(),
      source_url: @source_url
    ]
  end

  def application do
    [
      mod: {Orbis.Application, []},
      extra_applications: [:logger]
    ]
  end

  defp deps do
    [
      {:nx, "~> 0.7"},
      {:rustler, "~> 0.36"},
      {:req, "~> 0.5"},
      {:jason, "~> 1.4"},
      {:ex_doc, "~> 0.34", only: :dev, runtime: false},
      {:quokka, "~> 2.12", only: [:dev, :test], runtime: false}
    ]
  end

  defp description do
    """
    Satellite toolkit for Elixir — SGP4 propagation, coordinate transforms
    (0 ULP Skyfield parity), orbit determination, conjunction assessment,
    pass prediction, live TLE/OMM data, and real-time tracking. Rust NIF backend.
    """
  end

  defp package do
    [
      licenses: ["MIT"],
      links: %{"GitHub" => @source_url}
    ]
  end

  defp docs do
    [
      main: "Orbis",
      extras: [
        "README.md",
        "guides/track_iss.md",
        "guides/pass_prediction.md",
        "guides/conjunction_screening.md",
        "guides/accuracy.md",
        "guides/batch_analysis.md",
        "examples/iss_tracker.livemd"
      ],
      groups_for_extras: [
        Guides: Path.wildcard("guides/*.md")
      ],
      groups_for_modules: [
        Core: [Orbis, Orbis.Elements, Orbis.SGP4, Orbis.TemeState],
        Coordinates: [Orbis.Coordinates],
        "Ground Station": [Orbis.Passes, Orbis.Doppler, Orbis.RF, Orbis.Tracker],
        "Orbit Determination": [Orbis.IOD, Orbis.Lambert],
        "Space Environment": [Orbis.Eclipse, Orbis.Atmosphere, Orbis.Ephemeris, Orbis.Angles],
        Conjunction: [Orbis.Conjunction],
        "Data Sources": [Orbis.CelesTrak, Orbis.Constellation],
        "Batch Analysis": [
          Orbis.Nx,
          Orbis.Nx.Geometry,
          Orbis.Nx.Visibility,
          Orbis.Nx.RF,
          Orbis.Nx.Coverage
        ],
        Format: [Orbis.Format.TLE, Orbis.Format.OMM]
      ]
    ]
  end
end
