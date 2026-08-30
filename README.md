# aerophysics

`aerophysics` is a Python scientific-computing package for traceable,
vectorized atmospheric and aerodynamic physics models. It covers standard
atmosphere and gas properties, flight conditions, compressible flow, shocks
and expansions, and smooth flat-plate boundary layers. The documentation
records each model's assumptions, applicability, references, and verification.

Public calculation APIs use SI units and radians. The package provides local
engineering models rather than aircraft-level design or CFD analysis.

See [README.ja.md](README.ja.md) for the Japanese entry point.

## Installation

`aerophysics` requires Python 3.12 or newer. Until a PyPI distribution is
available, install the current wheel from the
[latest GitHub Release](https://github.com/pandorobo11/aerophysics/releases/latest):

```console
python -m pip install "https://github.com/pandorobo11/aerophysics/releases/download/v0.6.0/aerophysics-0.6.0-py3-none-any.whl"
```

To install and launch the optional local GUI:

```console
python -m pip install "aerophysics[gui] @ https://github.com/pandorobo11/aerophysics/releases/download/v0.6.0/aerophysics-0.6.0-py3-none-any.whl"
aerophysics-gui
```

Check the release page for a newer version before copying a versioned URL.

## Minimal example

```python
from aerophysics import FlightCondition, standard_atmosphere

sea_level = standard_atmosphere(0.0)
print(sea_level.temperature)       # K
print(sea_level.speed_of_sound)    # m/s

condition = FlightCondition.from_mach(
    geometric_altitude=10_000.0,
    mach=0.8,
    characteristic_length=2.0,
)
print(condition.dynamic_pressure)  # Pa
print(condition.reynolds_number)
```

## Documentation

- [Manual source and topic index](docs/index.rst)
- [Python quickstart](docs/getting_started/quickstart.rst)
- [Local GUI guide](docs/guides/gui.rst)
- [Verification](docs/verification/index.rst)
- [API reference](docs/api/index.rst)
- [Rendered HTML documentation and release downloads](https://github.com/pandorobo11/aerophysics/releases/latest)
- [Development and release workflow](DEVELOPMENT.md)

## License

[MIT](LICENSE)
