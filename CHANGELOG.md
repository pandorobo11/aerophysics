# Changelog

All notable changes to this project are documented in this file. The project
uses Semantic Versioning, including during the pre-1.0 period.

## [Unreleased]

### Changed

- Use one GUI selector for length and inverse-length display units so their
  dimensional bases stay consistent.
- Remove the pinned ``fluids`` comparison from required standard-atmosphere
  verification, separate numerical checks from the generated record, and
  regenerate that record before documentation builds.

## [0.6.0] - 2026-08-11

### Added

- Extend explicit unit conversions with inches, Celsius, Rankine, and
  pounds-force, and expand GUI display choices with millimetres, feet per
  second, engineering pressure units, and pounds mass per cubic foot.
- Add selectable area, force, and inverse-length display units, and apply unit
  preferences to Reynolds number per length, static density, direct drag, and
  frontal area in GUI tables and plots.

### Changed

- Add section landing pages and a complete top-level documentation map so the
  left sidebar consistently lists related pages instead of showing an empty
  section navigation panel.
- Force complete Sphinx rebuilds in the local validation gate so navigation
  hierarchy changes cannot leave stale sidebars in previously generated HTML.
- Replace the thermophysical CoolProp snapshot with U.S. Standard Atmosphere
  acceptance values and a reproducible, non-gating NIST dilute-air
  physical-accuracy assessment derived directly from Lemmon--Jacobsen.

## [0.5.0] - 2026-08-11

### Added

- Add reproducible verification records for compressible flow,
  thermophysical and real-gas properties, boundary layers, velocity profiles,
  and protrusion integration using primary references, pinned Cantera and
  CoolProp snapshots, generated figures, and offline invariant tests.
- Bundle the complete offline HTML manual in wheels and release archives, and
  add a local GUI documentation browser.

### Changed

- Reorganize the English manual into getting-started, task-guide,
  model-and-equation, verification, and domain API sections; update the local
  GUI and release artifacts to use the new document paths, simplify the
  English and Japanese README entry points, and move maintainer procedures to
  ``DEVELOPMENT.md``.
- Reuse Beattie--Bridgeman flow states in GUI calculations to reduce repeated
  nonlinear solves.

### Fixed

- Improve Beattie--Bridgeman isentropic branch solving and detached-shock
  numerical robustness near difficult geometries and solver limits.

## [0.4.0] - 2026-08-05

### Added

- Add a reproducible U.S. Standard Atmosphere 1976 verification record using
  official tables, a pinned ``fluids`` snapshot, physical invariants, and
  generated comparison figures.
- Add generic Kennard harmonic-oscillator and Beattie--Bridgeman gas models,
  physics-named dry-air presets, complete thermodynamic states, and real-gas
  isentropic forward/inverse, area--Mach, critical-state, and mass-flow APIs.
- Extend the isentropic GUI with harmonic-oscillator and Beattie--Bridgeman air,
  required real-gas reservoir pressure, absolute static properties, and
  applicability warnings.
- Add Taylor–Maccoll axisymmetric attached conical shocks, cone-surface state
  ratios, attached limits, vectorized APIs, and a dedicated GUI calculator.
- Extend all isentropic forward/inverse, area--Mach, critical-state, and mass
  flux relations to frozen-composition thermally perfect gases, and add
  NASA7/NASA9 selection to the local isentropic GUI.
- Separate transport properties into a dedicated module and add interchangeable
  Keyes and Blottner/Wilke dry-air viscosity models while preserving the
  Sutherland defaults and legacy imports.
- Add an optional Streamlit and Plotly local GUI prototype for flight
  conditions, oblique shocks, and flat-plate boundary layers.
- Add single-condition and parameter-sweep views, selectable display units,
  plots, flight-case handoff, CSV export, and versioned JSON settings.
- Extend the local GUI with isentropic forward and inverse relations, normal
  shocks and supersonic pitot pressure, and Prandtl-Meyer expansions.
- Add GUI pages for compressible turbulent boundary-layer profiles,
  boundary-layer-immersed protrusion drag, and NASA7/NASA9 frozen-air
  thermochemistry, including model comparison, validated CSV inputs, and
  SI-preserving handoff between viscous-flow pages.
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

[Unreleased]: https://github.com/pandorobo11/aerophysics/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/pandorobo11/aerophysics/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.5.0
[0.4.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.4.0
[0.3.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.3.0
[0.2.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.2.0
[0.1.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.1.0
