# Changelog

All notable changes to this project are documented in this file. The project
uses Semantic Versioning, including during the pre-1.0 period.

## [Unreleased]

### Added

- NASA seven- and nine-coefficient thermochemistry, reusable ideal-gas species,
  and frozen-composition thermally perfect gas mixtures.
- Built-in NASA7 and NASA9 dry-air models with temperature-dependent heat
  capacities, heat-capacity ratio, enthalpy, entropy, and speed of sound from
  200 to 6000 K.
- Direct drag estimates for isolated boundary-layer-immersed protrusions using
  frontal-area effective dynamic pressure.
- Constant or height-dependent frontal width, provided velocity and density
  profiles, and a default turbulent one-seventh-power profile.
- Optional Walz temperature and ideal-gas density approximation for
  compressible turbulent boundary layers.

## [0.3.0] - 2026-07-23

### Added

- Smooth zero-pressure-gradient flat-plate boundary layers with explicit
  laminar, turbulent, and specified-transition regimes.
- Blasius laminar thickness and friction relations.
- Selectable one-fifth-power and Schlichting turbulent correlations, with
  Schlichting as the default and applicability warnings outside the nominal
  Reynolds-number range.
- Eckert reference-temperature and Reynolds-number-based Van Driest II
  compressibility corrections.
- Adiabatic-wall defaults, optional specified wall temperature, recovery
  temperature, wall shear, and one-sided drag per unit width.

## [0.2.0] - 2026-07-23

### Added

- Normal-shock state ratios and total-pressure loss.
- Oblique-shock theta–beta–Mach relations with explicit weak and strong
  branches, maximum attached deflection, and a dedicated no-attached-shock
  exception.
- Prandtl–Meyer angle, inverse Mach calculation, centered expansion state
  ratios, and the limiting expansion angle.
- Supersonic Pitot pressure relation.
- Vectorized scalar/array APIs with all angles expressed in radians.

## [0.1.0] - 2026-07-22

### Added

- Calorically perfect-gas properties and source-specific dry-air constants.
- Sutherland dynamic viscosity and U.S. Standard Atmosphere air thermal
  conductivity models.
- Vectorized U.S. Standard Atmosphere 1976 from -5 to 86 km geometric altitude.
- Explicit SI conversions for common aviation customary units.
- Isentropic total-to-static relations, inverse calculations, area-Mach
  branches, and choked mass flux.
- Integrated `FlightCondition` construction from Mach number or velocity.
- Typed public APIs, cross-platform CI, Sphinx documentation, and validated
  examples.

[0.3.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.3.0
[0.2.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.2.0
[0.1.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.1.0
