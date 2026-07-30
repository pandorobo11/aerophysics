# aerophysics

`aerophysics` is a Python scientific-computing package for traceable,
vectorized atmospheric and aerodynamic physics models.

Version 0.3 includes perfect-gas properties, the U.S. Standard Atmosphere 1976,
isentropic compressible flow, normal and oblique shocks, Prandtl–Meyer
expansion, laminar and turbulent flat-plate boundary layers, flight conditions,
and explicit aviation unit conversions. Public calculation APIs use SI units;
angles use radians.

> The project is in its initial development phase. The public API may evolve
> under the documented pre-1.0 deprecation policy.

See [README.ja.md](README.ja.md) for a Japanese overview.

## Installation

The package requires Python 3.12 or newer.

```console
python -m pip install aerophysics
```

### Local GUI

Install the optional GUI dependencies to calculate and plot standard-atmosphere,
flight-condition, oblique-shock, and flat-plate boundary-layer results in a
local browser:

```console
python -m pip install "aerophysics[gui]"
aerophysics-gui
```

The GUI explicitly converts selected display units to the SI core API and can
export result CSV files and versioned calculation-configuration JSON files.

## Quick start

```python
from aerophysics import FlightCondition, standard_atmosphere

sea_level = standard_atmosphere(0.0)
print(sea_level.temperature)       # 288.15 K
print(sea_level.speed_of_sound)    # 340.294... m/s

condition = FlightCondition.from_mach(
    geometric_altitude=10_000.0,
    mach=0.8,
    characteristic_length=2.0,
)
print(condition.dynamic_pressure)  # Pa
print(condition.reynolds_number)
```

Calculation APIs use SI units. Use explicit functions from
`aerophysics.units` when converting aviation customary units.

## Models and references

- Atmospheric state and transport properties: *U.S. Standard Atmosphere,
  1976*, NOAA-S/T 76-1562 and NASA-TM-X-74335.
- Isentropic perfect-gas relations: NACA Report 1135, *Equations, Tables, and
  Charts for Compressible Flow*.
- Normal and oblique shocks, supersonic Pitot pressure, and Prandtl–Meyer
  expansion: NACA Report 1135.
- Smooth flat-plate boundary layers: Blasius and smooth-wall turbulent
  correlations, with Eckert and Van Driest II compressibility corrections.
- Direct drag of isolated boundary-layer protrusions using an
  effective-dynamic-pressure integral and optional Walz density profile.
- Unit conversion factors: NIST Special Publication 811, Appendix B.

The API documentation records assumptions, valid ranges, units, and ratio
conventions for each model.

## Development

```console
uv sync --all-groups --all-extras
uv run pytest
uv run ruff check .
uv run mypy
uv run sphinx-build -W -b html docs docs/_build/html
```

## Releases

Update the version in `pyproject.toml`, commit the change, and push a matching
`vX.Y.Z` tag:

```console
git tag v0.3.0
git push origin v0.3.0
```

The release workflow runs the full test suite and publishes the wheel, source
distribution, and `aerophysics-docs-X.Y.Z.zip` to the private repository's
GitHub Release. Extract the documentation archive and open
`aerophysics-docs-X.Y.Z/index.html` in a browser. Access to the release requires
read permission for this repository.

## License

MIT
