[![Hex.pm](https://img.shields.io/hexpm/v/orbis)](https://hex.pm/packages/orbis)
[![Hexdocs.pm](https://img.shields.io/badge/docs-hexdocs.pm-purple)](https://hexdocs.pm/orbis)
[![CI](https://github.com/neilberkman/orbis/actions/workflows/ci.yml/badge.svg)](https://github.com/neilberkman/orbis/actions)

# Orbis

Full-featured satellite toolkit for Elixir.

SGP4 propagation and high-accuracy coordinate transforms are handled by a Rust NIF. Everything else — pass prediction, orbit determination, conjunction assessment, constellation management, real-time tracking, batch analysis — is pure Elixir (and Nx for GPU workloads).

### Try it in Livebook

[![Run in Livebook](https://livebook.dev/badge/v1/blue.svg)](https://livebook.dev/run?url=https://github.com/neilberkman/orbis/blob/main/examples/iss_tracker.livemd)

## Features

| Category | What it does |
|----------|-------------|
| **Propagation** | SGP4/SDP4 via the [`sgp4`](https://crates.io/crates/sgp4) Rust crate (Rust NIF) |
| **Coordinate transforms** | TEME, GCRS, ITRS, geodetic, topocentric — 0 ULP Skyfield parity (Rust NIF) |
| **Ground station** | Pass prediction, look angles, Doppler shift, RF link budget |
| **Orbit determination** | Gibbs, Herrick-Gibbs, Gauss angles-only, Lambert/Battin (Rust NIF) |
| **Conjunction assessment** | Closest approach finder validated against the Iridium 33 / Cosmos 2251 collision |
| **Eclipse prediction** | Sunlit / penumbra / umbra with shadow fraction |
| **Atmospheric density** | NRLMSISE-00 model, surface to ~1000 km (Rust NIF) |
| **JPL ephemeris** | SPK/BSP reader for Sun, Moon, planets (Rust NIF) |
| **Live data** | CelesTrak TLE/OMM fetching, constellation loading, name search |
| **Real-time tracking** | GenServer with PubSub-compatible broadcasts |
| **RF primitives** | FSPL, EIRP, C/N₀, link margin, dish gain |
| **Batch analysis** | Nx-powered tensorized geometry, visibility, and RF (GPU-ready via EXLA/Torchx) |
| **Formats** | `Orbis.Elements` with TLE and OMM parsers/encoders |

## Installation

```elixir
def deps do
  [{:orbis, "~> 0.5.0"}]
end
```

Requires Rust for compiling the NIF.

## Quick Start

```elixir
# Fetch the ISS TLE from CelesTrak
{:ok, [iss]} = Orbis.CelesTrak.fetch_tle(25544)

# Propagate to now
{:ok, teme} = Orbis.propagate(iss, DateTime.utc_now())
teme.position   # {x, y, z} km
teme.velocity   # {vx, vy, vz} km/s

# Where is it over the Earth?
{:ok, geo} = Orbis.geodetic(iss, DateTime.utc_now())
# %{latitude: 23.4, longitude: -45.6, altitude_km: 420.1}
```

## Usage

### Parse from TLE or OMM

```elixir
# Two-Line Element format
{:ok, elements} = Orbis.Format.TLE.parse(line1, line2)

# OMM JSON (CelesTrak / Space-Track)
{:ok, elements} = Orbis.Format.OMM.parse(omm_json_map)

# Encode back to either format
{line1, line2} = Orbis.Format.TLE.encode(elements)
omm_map = Orbis.Format.OMM.encode(elements)
```

The `Orbis.Elements` struct is format-agnostic — parse from any source, serialize to any format.

### Coordinate Transforms

```elixir
gcrs = Orbis.teme_to_gcrs(teme, datetime)

station = %{latitude: 40.7128, longitude: -74.006, altitude_m: 10.0}
{:ok, geo} = Orbis.geodetic(elements, datetime)
{:ok, look} = Orbis.look_angle(elements, datetime, station)
# %{azimuth: 359.6, elevation: -41.9, range_km: 9130.5}
```

The GCRS transform includes IAU2000A nutation (1365 terms), IAU2006 precession, frame bias, and time scale conversions (UTC→TAI→TT→TDB→UT1).

### Constellation Management

```elixir
{:ok, constellation} = Orbis.Constellation.load("globalstar")
constellation.count  #=> 85

# Propagate all satellites in parallel
positions = Orbis.Constellation.propagate_all(constellation, DateTime.utc_now())

# Find visible satellites from a ground station
visible = Orbis.Constellation.visible_from(constellation, station, datetime)
```

### Real-Time Tracking

```elixir
{:ok, tracker} = Orbis.Tracker.start_link(elements, interval_ms: 1000)
Orbis.Tracker.subscribe(tracker)

receive do
  {:orbis_tracker, _pid, state} ->
    IO.puts("#{state.geodetic.latitude}, #{state.geodetic.longitude}")
end
```

### Pass Prediction

```elixir
passes = Orbis.Passes.predict(elements, station,
  ~U[2024-01-01 00:00:00Z], ~U[2024-01-02 00:00:00Z])

for pass <- passes do
  IO.puts("Rise: #{pass.rise} | Max el: #{pass.max_elevation}° | Set: #{pass.set}")
end
```

### Conjunction Assessment

```elixir
approaches = Orbis.Conjunction.find(elements1, elements2,
  end_min: 2880.0, step_min: 1.0, threshold_km: 50.0)

for {tca_min, distance_km} <- approaches do
  IO.puts("TCA: +#{Float.round(tca_min / 60, 1)}h, miss: #{Float.round(distance_km, 2)} km")
end
```

### RF Link Budget

```elixir
# Path loss from slant range
fspl = Orbis.RF.fspl(look.range_km, 1616.0)  # MHz

# Full link margin
margin = Orbis.RF.link_margin(%{
  eirp_dbw: Orbis.RF.eirp(27.0, 3.0),
  fspl_db: fspl,
  receiver_gt_dbk: -12.0,
  other_losses_db: 3.0,
  required_cn0_dbhz: 35.0
})
```

### Orbit Determination

```elixir
# Gibbs: 3 position vectors → velocity
{v2, theta12, theta23, copa} = Orbis.IOD.gibbs(r1, r2, r3)

# Gauss: 3 angular observations → full orbit
{r2, v2} = Orbis.IOD.gauss(
  decl1, decl2, decl3, rtasc1, rtasc2, rtasc3,
  jd1, jdf1, jd2, jdf2, jd3, jdf3,
  site1, site2, site3)

# Lambert: transfer orbit between two positions
{v1t, v2t} = Orbis.Lambert.solve(r1, r2, v1, dm, de, nrev, dtsec)
```

### Eclipse & Ephemeris

```elixir
eph = Orbis.Ephemeris.load("/path/to/de421.bsp")

{:ok, status} = Orbis.Eclipse.check(elements, datetime, eph)
# :sunlit | :penumbra | :umbra

mars = Orbis.Ephemeris.position(eph, :mars, :earth, datetime)
```

## Coordinate Frames

| Frame | Description |
|-------|-------------|
| TEME | True Equator Mean Equinox — SGP4 output frame |
| GCRS | Geocentric Celestial Reference System — inertial |
| ITRS | International Terrestrial Reference System — Earth-fixed |
| Geodetic | WGS84 latitude, longitude, altitude |
| Topocentric | Azimuth, elevation, range from a ground station |

## Accuracy

| Component | Reference | Accuracy |
|-----------|-----------|----------|
| TEME→GCRS→ITRS | Skyfield | 0 ULP (bit-identical) |
| SGP4 propagation | Skyfield | < 1 mm (via sgp4 crate) |
| Gibbs / Herrick-Gibbs | Vallado Python | 0 ULP |
| Gauss IOD | Vallado Python | 1e-12 relative |
| Lambert (Battin) | Vallado Python | 1e-12 relative |
| Conjunction | Iridium/Cosmos 2251 | < 2 km miss, < 1 min timing |
| RF (FSPL) | Analytical (inverse square law) | Exact |

## License

MIT
