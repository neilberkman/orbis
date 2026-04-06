# Collision Probability

Orbis provides a robust framework for assessing the risk of close approaches (conjunctions) between space objects. It supports the CCSDS Conjunction Data Message (CDM) standard and implements both analytical and numerical methods for computing the probability of collision (Pc).

## The Conjunction Data Message (CDM)

The CDM is the industry standard for communicating conjunction risk. It contains the states, covariances, and metadata for two objects at the predicted Time of Closest Approach (TCA).

```elixir
{:ok, cdm} = Orbis.CCSDS.CDM.parse(kvn_string)

cdm.tca                    # ~U[2026-04-10 12:34:56.789Z]
cdm.miss_distance_m        # 1250.0
cdm.collision_probability  # 1.23e-06
```

## Computing Pc

To compute Pc from a CDM, we first convert it to collision parameters:

```elixir
params = Orbis.CCSDS.CDM.to_collision_params(cdm)

# Compute Pc using the default 2D Foster equal-area method
{:ok, result} = Orbis.Collision.probability(params)
result.pc  #=> 1.23e-06
```

### Methods

Orbis supports multiple methods for Pc calculation:

1.  `:equal_area` (Default): The 2D Foster method using an equal-area square approximation. It is extremely fast and accurate for typical conjunction geometries.
2.  `:numerical`: Direct numerical integration of the 2D Gaussian over the hard-body radius circle. This is the gold standard for precision.

```elixir
{:ok, result} = Orbis.Collision.probability(params, method: :numerical)
```

## Catalog Screening

When managing a constellation, you need to screen thousands of pairs for potential risk. `Orbis.Screening` provides tools for catalog-scale analysis.

```elixir
# objects is a list of maps with :r, :v, :cov, and :hard_body_radius_km
results = Orbis.Screening.screen_catalog(objects, miss_threshold_km: 10.0)

# results is sorted by decreasing Pc
for res <- results do
  IO.puts("Risk: #{res.collision.pc} | Miss: #{res.candidate.miss_km} km")
end
```

## Frames and Covariances

Orbis handles the complex frame transforms required for collision analysis, including:
- **RTN (Radial, Transverse, Normal)** to **ECI (Earth-Centered Inertial)** covariance transforms.
- Construction of the **Encounter Frame** (centered at TCA, axes aligned with relative velocity and miss vector).
- Projection of the 3D covariance into the 2D **Encounter Plane**.

These primitives are available in `Orbis.Covariance` and `Orbis.Encounter`.
