# Changelog

All notable changes to this project are documented in this file. The project
uses Semantic Versioning, including during the pre-1.0 period.

## [Unreleased]

### Added

- Add `transition_distance` to `marching_boundary_layer` with
  `regime=TRANSITIONAL` for a user-specified transition position.
  The resulting SA-noft2 transition length is not a physical prediction.

- Add `axisymmetric_body_flow` for the supersonic, zero-incidence flow around
  sharp, convex, piecewise linear bodies of revolution, such as cone-cylinders.
  It uses the rotational method of characteristics with a fitted leading shock.
  `boundary_layer_inputs` returns the body-streamline edge and matching wall
  radii for `marching_boundary_layer`; each convex corner becomes a short
  pressure ramp of length `corner_length`. Boundary-layer results within a few
  tens of `delta99` of a sharp corner depend on that length and are not reliable.
  Accuracy against experiment is not established.

- Add optional `wall_radius` to `marching_boundary_layer` for thin axisymmetric
  laminar and SA-noft2 flow from a sharp tip or leading edge. Radii use SI metres
  and linear interpolation in wall distance; Cf uses local edge dynamic pressure,
  heat flux is positive into the wall, and positive beta means acceleration.
  Perfect-gas, isentropic-edge assumptions remain; transverse/longitudinal
  curvature and normal pressure gradients are neglected. The new delta99/r_w
  warning at 0.1 is an uncalibrated indicator. Planar arrays remain bit-identical.
  Verify Mangler cone and variable-radius relations (the 1948 original was not
  obtained), balances, and refinement; add R&M 3633 measurements and shoulder
  stress tests, including locally resolved widths of 0.001 and 0.0002 m.
  Pressure-length warnings fire at 0.001 m for SA and 0.0002 m for both regimes.
  All 16 cases complete with zero Picard fallbacks; the recorded times are
  2.82--8.85 s. Shoulder-adjacent SA-variable oscillations reach 1.46404%, so
  applicability warnings do not establish local profile convergence. R&M maximum
  Cf difference is 0.0004443 and delta99/r_w reaches 0.4774; experimental accuracy,
  transition, and shoulder edge flow remain open.

- Extend `marching_boundary_layer` to planar SA-noft2 flow with varying
  isentropic edge conditions, prescribed or adiabatic wall temperature,
  conserved dimensional edge `nu_tilde`, and a Reynolds-scaled outer grid.
  Add `turbulent_prandtl_number`, `edge_sa_ratio`, read-only `nu_tilde_ratio`
  profiles, and a first-node `y+` warning. Laminar numerical results are
  unchanged. Verify ZPG regression,
  integral balances, refinement, and the NASA TMR bump; experimental accuracy
  and transition prediction are not established.
- Add `marching_boundary_layer` and `BoundaryLayerEdge` for laminar boundary
  layers with streamwise-varying edge conditions. The solver marches the
  compressible thin-layer equations from a sharp leading edge for a calorically
  perfect gas with an isentropic edge, an isothermal or adiabatic wall, and an
  optional wall-temperature distribution. Separation raises the new
  `BoundaryLayerSeparationError`, which holds the separation distance and the
  upstream result. A failed march raises the new
  `BoundaryLayerConvergenceError`.
- Add an independent compressible turbulent boundary-layer reference
  (`docs/scripts/build_turbulent_sa_reference.py`) that marches the
  zero-pressure-gradient thin-layer equations from the leading edge with the
  Spalart--Allmaras model. Its stored results for six Mach and wall-temperature
  conditions at Re_x = 1e6, 1e7 and 1e8 are compared with
  `VAN_DRIEST_II_WALL_TEMPERATURE` in a generated documentation table. The
  solver is a verification tool and is not part of the package.
- Add flat-plate wall heat transfer to `flat_plate_boundary_layer`: compressible
  results now include `stanton_number`, `heat_transfer_coefficient` and
  `wall_heat_flux` from the Reynolds analogy, with
  `reynolds_analogy_factor` selecting the turbulent factor (default 1.14).
- Add `flat_plate_boundary_layer_profile`, which returns laminar or turbulent
  velocity and temperature profiles at one flat-plate station. Laminar profiles
  come from the compressible similarity solution.
- Show the Stanton number, heat-transfer coefficient and wall heat flux on the
  GUI flat-plate page for compressible calculations, and add a GUI laminar
  boundary-layer profile page.

### Changed

- Use projected active-set Newton for the SA station's `nu_tilde >= 0` bound,
  returning wall-attached zero regions to Picard and preserving legacy Newton
  for explicitly zero SA seeds. Require the scaled local complementarity
  residual at convergence and set undamped active rows exactly to zero,
  preventing projection stagnation and BLAS-dependent cancellation remnants.
  On the current solver-verification tier,
  including legacy comparison calls, Picard fallbacks decrease 5,586 -> 637.
  Remeasured outer-front `nu_tilde_ratio` differences reach 6.03e-5 of the
  profile maximum; velocity/temperature differences stay below 8.2e-8 in the
  compared marches. Generated verification records are unchanged.
- Stop bundling documentation in editable installs; standard wheels still include it.
- Solve the Spalding--Coles thickness closure of
  `TurbulentModel.SCHLICHTING_LOG_WAKE` and
  `VAN_DRIEST_II_WALL_TEMPERATURE` for all points at once. A safeguarded
  Newton iteration with the analytic derivative replaces the per-point Brent
  solve on the same integration grid, making the closure about 15 to 30
  times faster (about 0.17 ms instead of 4.8 ms per incompressible point).
  delta99 and displacement thickness are now converged to about 1e-14
  instead of the former 1e-6 root tolerance, so they change by up to about
  3e-7 relative. Inadmissible points keep the previous NaN values, warnings,
  and errors.

### Fixed

- Keep the saved inputs when a settings JSON is loaded on the GUI compressible
  boundary-layer profile page. A file saved from a boundary-layer case was
  recalculated with whichever case the session held instead. The input source
  now starts as manual input after loading.

## [0.8.0] - 2026-11-02

### Added

- Select geometric or geopotential altitude in the atmosphere/flight GUI,
  including altitude sweeps, plot axes, and saved configuration replay.
  Support the same choice for standard-atmosphere flow-sequence initialization.

- Add scalar flow sequences with normal, oblique and conical shocks,
  isentropic transformations and centered expansions, for perfect and frozen
  thermally perfect gases. Support arbitrary static/total initial conditions
  and standard-atmosphere initialization. Add a GUI stage editor, per-stage
  properties and pressure recovery, plots, SI CSV export, configuration replay,
  and explicit partial results when a stage fails.

- Add `TurbulentModel.SCHLICHTING_LOG_WAKE` for `compressible=False` and its
  GUI option. It keeps the Schlichting friction and drag, sets
  theta = x * average_cf / 2, and finds delta99 and displacement thickness
  with the constant-property Spalding--Coles closure of
  `VAN_DRIEST_II_WALL_TEMPERATURE`. Points without an admissible wake
  parameter return NaN delta99 and displacement thickness with an
  `ApplicabilityWarning`. `SCHLICHTING`
  (one-fifth-power thickness) remains the default.

- Warn when `SCHLICHTING` one-fifth-power turbulent thicknesses are used
  above Re_x = 1e7, their nominal limit. Their delta99 is estimated at +40%
  against measurements at Re_x = 1e6 and -8% at 1e8. The documentation now
  gives the correlation's derivation and accuracy.

- Warn once for compressible laminar points with viscous-interaction
  parameter `chi_bar > 1`, reporting the maximum and affected point count.

- Add `TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE` and its GUI
  option: Schlichting friction with a VD2 wall-temperature correction,
  momentum-integral theta, and Van Driest/Walz profile thicknesses.

- Add GUI selection of Poll fully turbulent attachment-line heating and
  swept-cylinder laminar/turbulent comparison, with diagnostic columns and
  missing turbulent values outside the model domain. Existing saved
  configurations remain compatible.

- Add `swept_cylinder_turbulent_heat_flux` / `AttachmentLineHeatingModel.POLL`
  for explicitly selected fully turbulent smooth-wall attachment-line heating.
  Use the Poll reference-temperature equations reproduced in NASA CR-181893
  with modified-Newtonian edge properties; require M cos Λ > 1 and
  0 < Λ < π/2. Warn above T_0 = 2000 K and outside Bénard et al.'s turbulent
  database ranges for M_e and T_w/T_0. Existing laminar APIs remain unchanged;
  the separate Bénard correlation is not implemented.

- Add an aerodynamic heating GUI page for sphere stagnation-point heating
  using Fay--Riddell, Tauber equation 40, or a comparison of both, and swept
  cylinder stagnation-line heating using Tauber equation 41. Set freestream
  conditions from standard-atmosphere altitude or manual static temperature
  and pressure, and calculate single points, Mach sweeps, geometric-altitude
  sweeps (standard atmosphere only), or cylinder sweep-angle sweeps.
  Single-point cylinder results include a heat-flux versus sweep-angle plot. Add sidebar heat-flux display units of
  W/m², kW/m², and MW/m². Existing configuration JSON files with six or nine
  `display_units` fields still load, defaulting `heat_flux` to W/m²; newly
  exported configurations include `heat_flux`.

- Add `sphere_stagnation_heat_flux` with `StagnationHeatingModel.FAY_RIDDELL`
  and `TAUBER`, and `swept_cylinder_heat_flux` for laminar convective heating.
  The simplified Fay--Riddell model assumes Lewis number 1, no dissociation,
  and a Newtonian velocity gradient. Tauber's NASA TP-2914 equations 40 and
  41 cover a sphere and a swept infinite-cylinder stagnation line, respectively.
  These models require a supersonic perfect-gas freestream, use SI units, and
  return heat flux positive into the wall. They emit `ApplicabilityWarning`
  for Fay--Riddell when T_0 > 2000 K, and for Tauber when V > 14 km/s,
  T_w > 4000 K, or h_inf/(V^2/2) > 0.1 (approximately M < 7.07 for γ = 1.4).

### Changed

- Update the locked development and CI environment (certifi, filelock,
  msgpack, Streamlit, and Pillow) to clear the release-time dependency audit.
  The locked Streamlit no longer pulls in GitPython. Installed package
  requirements are unchanged.

- **Breaking:** After a prescribed transition, every turbulent model now
  continues the fully turbulent plate whose momentum thickness equals the
  laminar value at transition. Local friction, average friction, drag,
  thicknesses, and `effective_reynolds_number` use the Reynolds number
  `Re_eq = Re_0 + Re_x - Re_tr` of that plate, where `Re_0` is found from the
  model's own fully turbulent friction. The former Prandtl-type mixed average
  kept the leading-edge origin for friction, giving post-transition friction
  of a layer 2.8 times as thick as the actual one at Re_tr = 5e5. With
  Schlichting friction, the local coefficient just after transition rises by
  26% (Re_tr = 5e5) to 42% (5e6), and average friction by 6--15% near
  transition, falling below 1% beyond 100 Re_tr. Drag stays continuous for
  every model. Returned theta is continuous for `SCHLICHTING_LOG_WAKE` and
  `VAN_DRIEST_II_WALL_TEMPERATURE`; `SCHLICHTING` one-fifth-power thicknesses
  use the same origin, so theta is 1.3--2.2% above Blasius at transition.
  Local friction is the exact derivative d(x * average_cf)/dx except for
  `VAN_DRIEST_II`, whose Willems local and average equations are separate
  fits. `SCHLICHTING_LOG_WAKE` no longer returns NaN
  thicknesses just after a late transition. Transition Reynolds numbers too
  low for the turbulent correlation raise `ModelRangeError`. Fully turbulent
  and laminar results are unchanged.

- **Breaking:** Replace `compressibility_correction` / `CompressibilityCorrection`
  and `turbulent_correlation` / `TurbulentCorrelation` with `compressible` and
  `turbulent_model: TurbulentModel | None`, with no compatibility aliases.
  Incompatible or unused combinations raise `ValueError`. Numerical results
  for retained models are unchanged by this API reorganization. The GUI separates compressibility
  from the turbulent model and reads compatible legacy settings files; new exports
  use the new keys with schema version 1.

  | Former correction / correlation | compressible | turbulent_model |
  | --- | --- | --- |
  | NONE / SCHLICHTING | False | SCHLICHTING (default) |
  | NONE / POWER_LAW (turbulent or transitional) | — | Removed; select SCHLICHTING |
  | ECKERT / SCHLICHTING | True | ECKERT |
  | ECKERT / POWER_LAW (turbulent or transitional) | — | Removed; select ECKERT with Schlichting friction |
  | VAN_DRIEST_II / either | True | VAN_DRIEST_II |
  | VAN_DRIEST_II_WALL_TEMPERATURE / either | True | VAN_DRIEST_II_WALL_TEMPERATURE |
  | Laminar / NONE | False | None (Blasius) |
  | Laminar / any compressible correction | True | None (similarity) |

  Laminar calls must omit the turbulent model. Compressible turbulent and
  transitional calls require an explicit compressible model; incompressible
  calls accept SCHLICHTING or SCHLICHTING_LOG_WAKE.

- **Breaking:** Remove `TurbulentCorrelation.POWER_LAW` and its Eckert
  reference-state variant. The one-fifth-power local friction coefficient is
  about 19% below Schlichting at Re_x = 1e8 and 32% below at Re_x = 1e9,
  despite sharing the same nominal Reynolds warning range. Retain the
  incompressible one-fifth-power **thickness** correlation unchanged.
  Legacy GUI settings selecting `power_law` for turbulent or transitional
  NONE/ECKERT flow now fail on loading and ask the user to select Schlichting;
  they are never silently converted. Laminar and VD2 legacy settings still
  load because that friction choice was unused. `TurbulentModel` now contains
  only SCHLICHTING, SCHLICHTING_LOG_WAKE, ECKERT (Schlichting at the reference
  state), VAN_DRIEST_II, and VAN_DRIEST_II_WALL_TEMPERATURE.

- **Breaking:** Replace laminar cf, average friction, wall shear, accumulated
  drag, theta, displacement and delta99 for `compressible=True` with
  direct perfect-gas similarity solutions. Recovery temperature is the solved
  adiabatic wall temperature, omitted wall temperature imposes zero heat flux,
  and laminar `effective_reynolds_number` is the external Reynolds number.
  The laminar contribution to transitional average drag changes accordingly.
  Use the Blasius delta99 coefficient 4.9100 instead of 5.0 for `compressible=False`.

- **Breaking:** `TurbulentModel.ECKERT` and `TurbulentModel.VAN_DRIEST_II` return NaN for turbulent delta99,
  displacement thickness, and momentum thickness, including post-transition
  points. Their former compressible thickness estimates were inaccurate.
  Friction and drag are unchanged by the turbulent-thickness removal.
  The GUI displays unavailable metrics as a dash and prevents saving a
  boundary-layer case with unavailable thickness or wall shear.

- Revise the user manual for readers opening it from the GUI help page. Each
  model page now starts with a plain-language physical background, the
  conventions page adds common notation and a table of which gas models each
  calculation accepts, and the GUI guide uses the Japanese page names with a
  page-to-documentation map. Solver tolerances, settings-file fields, and
  other implementation details move to the new `devdocs/` directory.
- Reorganize the manual: the GUI manual is now its own section; shock waves
  are split into an overview and normal/oblique, conical, and detached pages;
  all gas models (perfect, NASA polynomial, harmonic oscillator, and
  Beattie--Bridgeman) are described on one gas-models page, with the standard
  atmosphere on its own page; verification and references become appendices,
  and regeneration commands move to `DEVELOPMENT.md`.
- Add a theory-page button to each GUI calculator page, and update the GUI
  documentation selector to the new page structure, including the
  aerodynamic heating, sequential flow, and conventions pages.

### Fixed

- Correct compressible laminar momentum thickness, previously too large by
  `T*/Te`, to `theta = x * cf_local` on the edge basis.

- Convert Eckert reference-temperature skin friction from the reference to the
  edge dynamic-pressure basis with `rho*/rho_e = T_e/T*`. Laminar and
  turbulent Eckert portions, including laminar portions of Van Driest II
  calculations, previously omitted this factor and overpredicted local and
  average friction, wall shear, and drag by `T*/T_e` (about 2.5 for a Mach 6
  cold wall). The Van Driest II turbulent path is unchanged by this friction
  correction.
- Solve Beattie--Bridgeman density for `rho / rho_ideal` with a normalized
  pressure residual and verify the EOS pressure closure before returning it.
  The fixed absolute `1e-12 kg/m^3` tolerance previously returned a density
  15.6% low at 1e-9 Pa and 300 K for a zero-correction model. States whose
  ideal-gas density is not a normal float64 now raise `ModelRangeError`.
- Raise `ShockConvergenceError` for unresolvable perfect-gas conical shocks
  instead of leaking `ValueError`, `RuntimeError`, or `RuntimeWarning` from
  the numerical solvers, and stop trial Taylor--Maccoll integrations at the
  singular line `a^2 = V_theta^2`. This applies to nonzero cone half-angles
  with `M - 1 < 1e-6`, which previously could hang or return a non-physical
  attached limit, and to very slender cones. A converged root must now match
  the requested cone half-angle to `1e-3` relative; slender cones below about
  `1e-3` rad previously could return the zero-angle Mach-wave state (unit
  surface ratios) or a root up to 24% off the requested angle. Resolvable
  results are unchanged.
- Recapture the NASA SP-3004 cone-table verification fixture from the report
  without solver values. The capture script previously read the wrong scan
  pages for 21 of 46 cases and silently replaced unmatched or mismatched cells
  with rounded `conical_shock` results, so 29 rows partly or wholly verified
  the solver against itself. Tables are now located from verified pages,
  cells are selected by table position, damaged OCR cells are transcribed
  from the scan with recorded locations, and every table passes internal
  consistency checks. Table numbers and printed pages are corrected. The
  solver still agrees with every cell within the `1e-4` tolerance.
- Keep Billig shock coordinates as owned, read-only snapshots so input-array
  changes cannot alter a previously computed shock shape.
- Reload flat-plate GUI settings with an unspecified Mach number when
  compressibility correction is disabled.
- Preserve validated sequence temperatures at thermal range boundaries when
  ratio restoration overshoots by exactly one ULP, while rejecting larger
  excursions and out-of-range user inputs. Bracket stagnation temperature
  locally so custom harmonic-oscillator gases without an applicability range
  do not evaluate properties at the maximum floating-point temperature.

## [0.7.0] - 2026-10-03

### Added

- Extend centered Prandtl--Meyer expansions to frozen thermally perfect gases
  with upstream static temperature, variable heat capacity, and strict static
  temperature bounds. Add thermal gas selection, temperature results,
  Mach/turn sweeps, and settings replay to the expansion GUI. Keep valid states
  when only the sonic reference for the absolute Prandtl--Meyer angle is outside
  the model range.

- Add single-point normal-shock and Prandtl--Meyer expansion GUI schematics
  with solved Mach labels, physical flow directions, a convex-wall expansion
  fan and angle arcs in the selected display unit.

- Add local gamma, full velocity ratios, sound speeds, dynamic-pressure ratios
  and entropy rise to shock GUI results, with optional cp/cv detail columns.
  Add optional static pressure, dry-air viscosity selection and local unit
  Reynolds numbers to normal, oblique and conical shocks; the isentropic page
  uses its existing total pressure. Optional length supplies Re_L. Preserve
  valid flow states when a viscosity correlation is outside its fitted range,
  and retain all output options through unit changes and settings replay.

- Add thermal gas selection, upstream static temperatures, Mach sweeps, and
  settings replay to the normal-shock GUI. Extend the supersonic pitot ratio
  to frozen thermally perfect gases; retain valid shock results when only the
  pitot stagnation state is outside the model range.

- Extend weak conical shocks and physical attached-cone limits to frozen
  thermally perfect gases using temperature-dependent Taylor--Maccoll flow.
  Enforce static temperature ranges through the cone surface and add thermal
  gas selection, temperatures, sweeps, and settings replay to the conical GUI.

- Extend normal and oblique shock relations to frozen thermally perfect gases
  using NASA7/NASA9 or harmonic-oscillator heat capacities and upstream static
  temperature, including weak/strong angle inversion and attached limits.
  Add gas selection, static temperatures, and settings replay to the oblique
  shock GUI while preserving the constant-gamma default.

- Add a fused ``isentropic_analysis`` API that reuses solved flow and critical
  states when ratios, area, mass-flow, and absolute-state results are needed
  together.

### Changed

- Build and validate release artifacts before publishing them from a separate
  job with repository write access. Verify wheel and source-distribution
  installations in clean environments.
- Clarify oblique, conical and detached shock GUI schematics with shaded bodies,
  horizontal angle references, labelled angle arcs and physical flow-direction
  arrows. Use bounded shock rays to keep near-normal shocks readable, and show
  detached-shock nose radius, standoff distance and upstream coordinate direction.
  Display flow from left to right consistently across all three schematics.
  Label oblique upstream/downstream and conical upstream/surface arrows with
  solved Mach numbers, matching the normal-shock and expansion schematics.
- Validate imported GUI settings against calculator-specific input, model, and
  sweep schemas before applying them, including finite and JSON-safe values.
- Use the fused isentropic analysis path in GUI sweeps so thermally perfect and
  Beattie--Bridgeman states are not solved repeatedly for each displayed
  quantity.
- Return vectorized calculation results as owned, read-only NumPy snapshots so
  later mutation of caller inputs cannot change previously computed states.
- Use one GUI selector for length and inverse-length display units so their
  dimensional bases stay consistent.
- Remove the pinned ``fluids`` comparison from required standard-atmosphere
  verification, separate numerical checks from the generated record, and
  regenerate that record before documentation builds.

### Fixed

- Update locked urllib3 to 2.8.0 to address streaming-response and HTTPS-proxy
  TLS vulnerabilities reported by the release dependency audit.
- Use the supplied boundary-layer edge viscosity consistently in Eckert and
  Van Driest II compressibility corrections.
- Bind the GUI launcher to the loopback interface by default.

- Keep the downstream Mach arrow between the convex wall and the final Mach
  line in expansion schematics, including narrow sectors at large turns.

- Render conical-shock sweep inputs and the calculation button without solving
  the attached limit first; use an editable 0--30 degree initial angle range.

- Preserve an in-range thermal conical attached limit when its temperature
  boundary is close to the physical maximum. Distinguish shock convergence
  failures from range and detachment errors, and retain failed GUI sweep points
  without discarding later valid results.

- Keep Beattie--Bridgeman density inversion on the lowest-density stable gas
  branch when custom coefficients permit multiple roots.
- Avoid spurious thermally perfect area-inversion applicability warnings caused
  only by intermediate root-solver probes outside the fitted temperature range.

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

[Unreleased]: https://github.com/pandorobo11/aerophysics/compare/v0.8.0...HEAD
[0.8.0]: https://github.com/pandorobo11/aerophysics/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/pandorobo11/aerophysics/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/pandorobo11/aerophysics/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.5.0
[0.4.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.4.0
[0.3.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.3.0
[0.2.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.2.0
[0.1.0]: https://github.com/pandorobo11/aerophysics/releases/tag/v0.1.0
