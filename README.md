# aerophysics

`aerophysics` is a Python scientific-computing package for traceable,
vectorized atmospheric and aerodynamic physics models.

Version 0.1 focuses on perfect-gas properties, the U.S. Standard Atmosphere
1976, isentropic compressible flow, flight conditions, and explicit aviation
unit conversions. Public calculation APIs use SI units.

> The project is in its initial development phase. The public API may evolve
> under the documented pre-1.0 deprecation policy.

See [README.ja.md](README.ja.md) for a Japanese overview.

## Development

```console
uv sync --all-groups
uv run pytest
uv run ruff check .
uv run mypy
```

## License

MIT

