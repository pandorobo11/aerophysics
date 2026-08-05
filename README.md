# aerophysics

`aerophysics` is a Python scientific-computing package for traceable,
vectorized atmospheric and aerodynamic physics models.

Version 0.3 includes perfect-gas properties, the U.S. Standard Atmosphere 1976,
isentropic compressible flow, normal, oblique, and conical shocks, Prandtl–Meyer
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
flight-condition, isentropic-flow, normal-, oblique-, and conical-shock,
Prandtl-Meyer-expansion, and flat-plate boundary-layer results in a local
browser. Compressible boundary-layer profiles, boundary-layer-immersed
protrusion drag, Sutherland/Keyes/Blottner-Wilke viscosity comparisons, and
NASA7/NASA9 frozen-air thermochemistry are also available:

```console
python -m pip install "aerophysics[gui]"
aerophysics-gui
```

The GUI explicitly converts selected display units to the SI core API and can
export result CSV files and versioned calculation-configuration JSON files.
Turbulent flat-plate cases can be handed to the profile predictor, and generated
profiles can be handed to the protrusion calculator without leaving SI. Profile
and projected-shape CSV inputs are supported by the protrusion page.

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

Temperature-dependent isentropic flow uses the same API with a thermally
perfect gas and an explicit total temperature:

```python
from aerophysics import AIR_NASA9
from aerophysics.isentropic import isentropic_ratios

ratios = isentropic_ratios(
    2.0,
    AIR_NASA9,
    total_temperature=1000.0,
    allow_extrapolation=False,
)
print(ratios.total_pressure_ratio)  # 7.89467...
```

The full forward/inverse, area–Mach, critical-state, and choked-mass-flux
relations support frozen-composition NASA7/NASA9 gases. Their fitted range is
200–6000 K; isentropic calculations warn when their default extrapolation is
used outside it.

Harmonic-oscillator and Beattie–Bridgeman air presets use physics-based names:

```python
from aerophysics import AIR_BEATTIE_BRIDGEMAN
from aerophysics.isentropic import isentropic_state

state = isentropic_state(
    2.0,
    AIR_BEATTIE_BRIDGEMAN,
    total_temperature=1200.0,
    total_pressure=6.0e6,
    allow_extrapolation=False,
)
print(state.static_pressure, state.velocity)
```

`AIR_HARMONIC_OSCILLATOR` follows Kennard's frozen harmonic-vibration model.
`AIR_BEATTIE_BRIDGEMAN` uses the Beattie–Bridgeman equation of state with
Randall air constants. JAXA-RR-06-011 supplies the wind-tunnel implementation,
constants, and reference calculations; JAXA did not originate either model.
The documented reservoir ranges are 400–2000 K for both presets and 1–10 MPa
for Beattie–Bridgeman air. Neither includes dissociation, chemical equilibrium,
ionization, condensation, or phase change.

## Models and references

- Atmospheric state and transport properties: *U.S. Standard Atmosphere,
  1976*, Sutherland, Keyes, and frozen-composition Blottner/Wilke viscosity
  models.
- Calorically and thermally perfect isentropic relations, including
  NASA-polynomial enthalpy and entropy variation; constant-heat-capacity
  relations follow NACA Report 1135.
- Kennard harmonic-oscillator thermally perfect gas and the Beattie–Bridgeman
  dense-gas equation of state, using Randall/JAXA air constants and
  wind-tunnel calculation procedures.
- Normal and oblique shocks, supersonic Pitot pressure, and Prandtl–Meyer
  expansion: NACA Report 1135. Axisymmetric sharp-cone shocks use the
  Taylor–Maccoll model and NASA SP-3004 reference tables.
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
git tag v0.4.0
git push origin v0.4.0
```

The release workflow runs the full test suite and publishes the wheel, source
distribution, and `aerophysics-docs-X.Y.Z.zip` to the private repository's
GitHub Release. Extract the documentation archive and open
`aerophysics-docs-X.Y.Z/index.html` in a browser. Access to the release requires
read permission for this repository.

## License

MIT
