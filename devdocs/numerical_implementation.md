# Numerical implementation notes

This file collects solver details that were previously part of the user
manual. The governing equations, assumptions, and applicability ranges remain
in the corresponding `docs/models/*.rst` page; this file explains how the code
solves them and how it separates numerical failures from physical limits.

## Isentropic flow

### Fused analysis

`isentropic_analysis` solves each distinct requested state once and shares one
Mach-one critical state for each distinct reservoir condition. At Mach zero it
records the limiting area ratio as positive infinity, while the standalone
`area_ratio` still rejects Mach zero.

### Beattie--Bridgeman states

- For a specified `T0, p0, M`, the entropy solve is restricted to the
  low-density, mechanically stable gas branch: densities below the first
  `(dp/drho)_T = 0` spinodal.
- The density inverse solves for `rho / rho_ideal`, where
  `rho_ideal = p / (R T)`. Its accuracy is therefore relative to the state,
  not to a fixed density in kg/m³.
- The returned density must reproduce the requested pressure through the same
  equation of state to a relative tolerance of `1e-12`.
- A state whose `rho_ideal` is not a finite normal float64 value raises
  `ModelRangeError`.
- Non-positive compressibility, heat capacity, or sound speed, and absence of
  a stable gas root, raise `ModelRangeError`.

## Thermally perfect normal and oblique shocks

The temperature Hugoniot removes the zero-strength root using
`x = u_n2 / u_n1`, `tau = T2 / T1`, and `q = [h(T2) - h(T1)] / (R T1)`:

```text
x^2 + (2q + 1 - tau) x - tau = 0
u_n1^2 = R T1 (tau / x - 1) / (1 - x)
```

- The sonic endpoint is evaluated analytically.
- The positive-temperature compressive root is bracketed numerically.
- The solver maximizes the shock polar, then brackets the shock-angle root
  separately on the weak and strong branches.
- For custom caloric models, the branch construction assumes the usual
  single-maximum, convex-gas shock polar.
- A valid weak shock can be returned even when the normal shock or the polar
  maximum is outside the gas temperature range. An unavailable strong root or
  attached limit then raises `ModelRangeError`, not `NoAttachedShockError`.

### Thermal pitot ratio

`supersonic_pitot_pressure_ratio` with a thermal gas finds `T0` from
`h(T0) = h(T1) + u1^2 / 2`. It raises `ModelRangeError` when any required
static or stagnation temperature is outside the model range and never
extrapolates. This stricter requirement does not apply to `normal_shock` or
its `total_pressure_ratio`, which only need the static states.

## Conical shocks

### Calorically perfect gas

- The cone angle is integrated with the Taylor--Maccoll equations from each
  trial shock angle. A trial integration that reaches the singular line
  `a^2 = V_theta^2` is treated as unresolved.
- A converged shock-angle root is accepted only if its integrated cone
  half-angle matches the request to a relative tolerance of `1e-3`. Near the
  Mach wave the cone angle is extremely sensitive to the shock angle, so this
  rejects slender-cone roots that miss the requested cone or collapse onto the
  zero-angle Mach wave. In practice it affects only cone half-angles below
  about `1e-3` rad.
- For a nonzero cone half-angle with `M - 1 < 1e-6`, and for very slender
  cones whose trial shocks approach the Mach wave within floating-point
  resolution, the solver raises `ShockConvergenceError`.
- None of these cases substitute a normal, detached, or Mach-wave solution.

### Thermally perfect gas

- The thermal path scales velocities by the free-stream speed `U_inf` rather
  than by the vacuum limiting velocity.
- Enthalpy inversion uses the gas model's actual reference, including NASA
  polynomial enthalpy offsets. A local-gamma substitution in the calorically
  perfect energy relation is not used.
- A temperature-boundary event and a direct surface-enthalpy check keep every
  static state between the shock and the surface inside the inclusive gas
  temperature range, including integration steps that straddle the surface.
  Internal Runge--Kutta trial stages evaluate only in-range properties.
- The temperature-boundary classification uses the shock-angle separation
  between the optimized maximum and the boundary, accounting for the numerical
  maximizer's stopping precision. An angle-height difference alone is not
  reliable near the flat top of the cone-angle curve.
- The search assumes the usual single-maximum cone-angle curve of a convex
  caloric gas.
- Numerical integration, maximization, or root-resolution failure raises
  `ShockConvergenceError` (a `RuntimeError` subclass). It is distinct from
  physical detachment (`NoAttachedShockError`) and temperature-range failure
  (`ModelRangeError`). The solver does not widen its angular tolerance to
  return a state for very slender cones.

## Detached-shock correlations

- Billig's curvature product is evaluated in logarithmic form, and the
  hyperbola increment uses a cancellation-resistant form.
- Because Billig's fitted curvature diverges as `M -> 1+`, an admissible input
  can still exceed finite float64 representation. Unrepresentable curvature or
  non-finite shock coordinates raise `ValueError` instead of returning
  infinities.

## Thermally perfect Prandtl--Meyer expansion

- The turn-angle integral is evaluated with the substitution
  `x = sqrt(T / T1)` and is split at NASA polynomial region boundaries.
- The downstream temperature is solved directly inside the model range; no
  fixed `gamma(T1)` is substituted.
- Quadrature or root convergence failures raise `ExpansionConvergenceError`.

## Flat-plate boundary layer

- Van Driest II: the inverse-sine arguments are clipped to `[-1, 1]` only to
  suppress floating-point excursions at the endpoints.
- Both implicit Van Driest II friction equations are solved with Brent's
  method on a verified positive bracket. A missing bracket or failed solve
  raises `ModelRangeError`.
- At the incompressible adiabatic limit `M -> 0`, the implementation assigns
  `F_C = F_theta = F_x = 1` directly.

## Marching boundary layer

`marching_boundary_layer` is implemented in `_boundary_layer_march.py`.

- Each station solves `(f, u, g)` on a fixed `eta` grid with a damped Newton
  iteration. The block-tridiagonal Jacobian uses exact discrete derivatives for
  laminar flow and SA-noft2 with exact built-in `SutherlandModel` or `KeyesModel`
  types. SA solves `(f, u, g, n, Y)`. Other viscosity models (including subclasses)
  retain three-color finite differences. Both paths assemble directly into
  expanded `gbsv` bands and use the same banded LAPACK solve. Laminar convergence
  requires a maximum update below `1e-10` in `u` and in the relative `g`.
- The analytic path uses the model's viscosity value and validation at current
  temperatures, with a private derivative in Pa s/K. For Sutherland,
  `mu_T/mu = 1.5/T - 1/(T+S)`; for Keyes, with
  `a = A 10**(-B/T)`, `mu_T/mu = 1.5/T - (1+a ln(10) B/T**2)/(T+a)`.
  Thus `C_g = (T_e mu_T/mu_e - C)/g`. This derivative enters both face
  fluxes, dissipation, and the adiabatic half-cell balance. Convection and
  streamwise BDF terms are differentiated with edge and lagged states fixed.
  SA also differentiates the eddy transport, production, destruction, and
  stretched wall distance. Nonsmooth operations use the active branch, with
  zero derivatives inside floors/caps; exact ties use the right-sided derivative
  for a positive perturbation of each unknown (including the incoming sign).
  The residual, grids, step control, damping, and tolerances are unchanged.
  This changes no physical applicability range or coordinate convention:
  `g = T/T_e` and `eta` are dimensionless, temperatures are kelvin, positive
  `beta` accelerates the edge, and reported heat flux is positive into the wall.
- The adiabatic wall uses a half-cell energy balance, which keeps the Jacobian
  block tridiagonal: at the wall the convective terms vanish and only viscous
  dissipation remains.
- The grid has a first spacing of `1e-3`, a growth ratio of 1.015, a largest
  spacing of 0.05, and `eta_max = 16` (519 nodes). It is second order, and this
  resolution holds the Falkner--Skan wall shear within 0.05% down to
  `beta = -0.15`. Richardson estimates for an expansion from `M = 2` to `3` give
  wall-normal errors of about 0.01% in wall shear and 0.08% in wall heat flux.
  A station whose `u` or `g` differs from one by more than
  `1e-4` five nodes from the outer boundary raises
  `BoundaryLayerConvergenceError`.
- The march uses variable-step BDF2 in `xi`, with BDF1 on the first step. The
  first step is `1e-6` of the first output distance. Steps grow by at most 1.25
  and are limited to 3% of the distance from the leading edge, so the march
  always passes through the requested stations.
- A step is halved when Newton fails, the wall shear is not positive, the
  logarithm of the wall shear changes by more than 0.0125, the logarithm of `U_e`
  changes by more than 0.000625, or the wall temperature gradient changes too
  much. Without the temperature check, a wall-temperature change can converge
  in a single step while the heat flux stays unresolved.
- The temperature check limits the change in `g'_w` to 2% of the largest of its
  two magnitudes and `0.1 f''_w (1 + (gamma-1) M_e^2/2)`. The floor is needed
  because the heat flux can pass through zero. An adiabatic wall skips the check.
- When the step falls below `1e-10` of the final distance, the march reports
  separation if `beta < 0` and the wall shear is below a quarter of its peak.
  Otherwise it raises `BoundaryLayerConvergenceError`. The separation distance
  extrapolates the last six stations with a straight line in the squared wall
  shear, because the wall shear collapses as the square root of the remaining
  distance.
- Compared with fully refined marches at `M = 2` to `4`, the streamwise errors
  are at most about 0.4% in heat flux and 0.12% in skin friction through
  pressure changes as short as about 2 cm, and 0.3% in heat flux through wall-temperature
  changes of up to 100 K over 1 mm to 50 cm. The largest relative errors occur
  where the heat flux is small. One case takes about 1 to 2 s.
- `xi` is integrated with five-point Gauss--Legendre quadrature on 16 panels per
  input interval.

### SA-noft2 marching extension

- `Discretization` retains the laminar residual and numerical operations;
  `_boundary_layer_sa.TurbulentDiscretization` adds `(n, Y)` to `(f, u, g)`.
  Three representative laminar cases retain all 42 output arrays bit for bit.
- The SA grid uses geometric **spacings**, first `7e-4`, ratio `1.04`.
  Before the march, `Re_equivalent = max(R)^2/2` over the output stations sets
  `scale = max(1, (Re_equivalent/1e8)^0.42)`, extent `150 scale`, and outer
  spacing cap `scale`. The grid is unchanged at or below 1e8 (312 nodes;
  first positions 0.0007 and 0.001428), with a roundoff guard at that boundary.
  The estimate is not a guarantee for compression; the outer-boundary check
  still stops a layer that outgrows it.
  Refining to `(3.5e-4, 1.02, 0.5, 150)` changes expansion Cf by at most
  0.0547%, heat flux by 0.2024%, and delta99 by 0.0705%.
- The shared streamwise limits above are unchanged. In the M2 ZPG study there
  were 625 converged trial checks and zero shear-limit rejections. Typical V9
  cases took 3--10 s, with 0--234 Picard fallbacks; refined cases took 15--40 s.
  These are legacy positivity-damped Newton measurements. Strong expansions
  exceeded that cost, as recorded below.
- Newton uses 50 iterations and a `1e-9` applied-update tolerance. Interior SA
  rows solve `n >= 0`, `R <= 0`, `n R = 0`, with `R = residual[3]`.
  For negative finite `J_ii`, use `s_i = 1/abs(J_ii)` and select the active
  branch of `min(n_i, -s_i R_i)` when `n_i <= -s_i R_i`. Active rows become
  `delta n_i = -n_i` directly in the existing gbsv bands; inactive rows keep
  the finite-difference Newton equations. A nonnegative or nonfinite SA
  diagonal conservatively requests Picard. Velocity and temperature retain
  the base damping limits; SA is projected onto `n >= 0` after the update,
  with no SA positivity fraction. For an undamped step, active n values are
  assigned exactly zero rather than relying on LU cancellation; damped active
  rows retain the partial step toward zero. Exact boundary values are restored.
  Convergence also requires an unchanged active set in the last two
  linearizations and `max(abs(Phi))/max(1, max(n)) < 1e-9`. Phi and its scale
  use the already computed pre-update linearization, before active-row
  replacement. This is the same n scaling and tolerance as `update_size`,
  requiring no additional residual/Jacobian evaluation. If the applied update
  is small and the active set is stable but Phi fails this bound, Newton
  continues if the applied update is nonzero; a zero update requests Picard.
  The iteration limit still applies. Projection/damping stagnation is not
  convergence. At convergence, an active **or projected-zero** first interior
  SA node requests Picard: wall-attached zero regions select the trivial
  partially laminar branch, not the intended outer-front clipping. The actual
  projected n must also be checked because a residual-branch update can be
  projected to zero without selecting the n branch. Intermediate active sets
  may recover, so this guard is applied only at convergence. An explicitly
  zero `edge_sa_ratio` uses legacy Newton directly to preserve the intentional
  zero-SA/laminar limit; it does not enter active-set iteration.
  The private `_bound_constrained_newton = False` switch
  retains the old Newton, including its `1e-4` positivity-damping cutoff.
  Picard retains the frozen reference's
  constants, smooth S-tilde limiter, nonnegative clipping, 0.5 relaxation,
  4000-iteration cap, and `1e-9` tolerance. Beta and Lambda sinks are implicit
  only when they strengthen the diagonal. Failure of both methods halves the step.
- Picard solves a different discrete problem at clipped outer-front nodes.
  Its frozen tridiagonal solve is `A(n) z = b(n)`, with `R = A(n)n - b(n)`
  at a fixed point of the other fields. Projection gives
  `n = max(z, 0)`, equivalently `min(n, A(n)^(-1) R) = 0`.
  The clipping correction propagates through the off-diagonal entries of A,
  leaving nonzero R at positive neighbors of clipped nodes. Newton instead
  enforces zero R at every positive interior n and nonpositive R at zero n.
  Both use nondimensional `n = nu_tilde/nu_e`; the inverse Jacobian diagonal
  scales R to n units. No residual, physical model, coordinate convention, or
  applicability range has changed. Remeasured with the analytic Jacobian,
  the difference on one captured station was `1.26e-5` in n (`3.32e-8` of its
  maximum), but this is not a general bound. Remeasured with the wall guard at `Re=1e10`, M2, edge 220 K
  and wall 300 K, the largest n difference is `8.098729`, or `6.03083e-5` of
  maximum n, at node 301 (zero-based), eta `632.73694`, s=100 m. The next,
  last positive node before clipping (302, eta `639.65525`) changed from
  `1996.165999` to `2002.957444`, a difference of `6.791445`; node 303 is
  effectively zero in both solutions. Velocity and temperature profile
  differences are `6.03350e-8` and `8.18300e-8` of their maxima; the wall
  velocity slope changes by `1.53488e-7` relative, and delta99 by `7.22122e-7`.
  Complete M2--3 fixed-wall and adiabatic marches differed by at most
  `2.111e-6` of maximum n; other fields differed by at most `7.072e-8` of
  their respective maxima. The tested ZPG and variable-radius M2 cases were
  bit-identical. These comparisons do not establish physical validation or
  grid convergence of the outer SA front.
- Regression comparisons scale absolute errors by each legacy field maximum,
  including near-zero entries. The M2--3 bounds are `4e-6` for n and `2e-7`
  for other fields. The measured high-Re case uses `1e-4` for n, `2e-7` for
  velocity/temperature, `2e-6` for other fields, and `4e-7` relative for wall
  velocity slope. Complete marches and the captured-station comparison run
  in `solver_verification`. Six-node and synthetic station tests keep the
  Newton branches covered in the normal tier, including legacy dispatch,
  diagonal rejection, banded failures, exhaustion, projection and active-set
  convergence. The legacy-switch and laminar complete-march comparisons check
  bitwise equivalence.
- In the earlier finite-difference-Jacobian rounded-transition reproducer
  (transition at 0.1 m, wall 250 K), the unguarded Newton selected a zero region
  at nodes 1--82 on the 472nd SA station attempt, xi `3.4502924398659043e-4`. Its wall velocity slope
  was `0.469137`, versus Picard's `0.556645`. The guard rejects this state.
  The guarded reference march's Cf and heat flux match the legacy march
  exactly in the measured comparison, within the existing `1e-6` bound.
- After integrating PR #94, the **90-test `solver_verification` tier** was
  remeasured with the **analytic Jacobian**, using four workers on macOS arm64,
  Python 3.13.5, NumPy 2.5.1 / Accelerate and SciPy 1.18.0. Temporary in-memory
  counters changed no numerical decisions. Two runs set the class default to
  legacy or guarded active-set Newton; explicit test overrides were retained.
  Thus both totals include legacy comparison calls, and the legacy-default
  run includes explicitly requested active-set solves. Each run passed all
  90 tests. Counts include repeated worker-local fixtures and direct solver
  calls, so they are observations for this tier, not per-march guarantees.

  | Measurement | Previous finite-difference Jacobian | Merged analytic Jacobian |
  | --- | --- | --- |
  | Picard fallbacks, legacy -> guarded default | 5,586 -> 637 | 5,586 -> 634 |
  | Picard calls, legacy -> guarded default | 5,587 -> 638 | 5,587 -> 635 |
  | Picard sweeps, legacy -> guarded default | 806,235 -> 107,160 | 806,234 -> 105,432 |
  | Guarded-default active-set wall rejections | 32 | 30 |
  | Guarded-default active-set iteration exhaustion | 7 | 6 |
  | Guarded-default Phi-withheld iterations | 66 | 67 |

  Every Picard call converged. In the guarded-default run, invalid-diagonal,
  nonfinite-update, caught-exception and Phi-stagnation failures were all zero;
  legacy zero-damping failures remained 600. The Phi bound withheld convergence
  on 67 iterations; small nonzero updates still allow further Newton iterations.
  The high-Re n, u and g differences above are unchanged at their quoted precision; the
  relative wall velocity slope difference changed from `1.53615e-7` to
  `1.53488e-7`. The captured-station difference is `1.25960e-5` in n, or
  `3.32343e-8` of its maximum. Generated verification records and reference
  data were not modified; all 14 generated-asset checks passed.
- Both Jacobian paths passed verification: the default analytic path passed
  all 90 solver-tier tests, and forcing the private SA analytic switch off
  passed all 29 active-set and transition solver-tier tests. The normal SA
  and active-set tests passed 65 cases with 100% line and branch coverage of
  the SA module. A real six-node station compares the SA diagonals in gbsv
  storage at `rtol=2e-6` (forward-difference truncation/cancellation), and its
  converged active-set states at `rtol=atol=1e-12` (roundoff allowance below
  the `1e-9` stopping threshold); the maximum absolute state difference was
  `8.88e-16`. Synthetic-residual tests explicitly disable the analytic
  Jacobian so that derivatives match their substituted equations. The
  stagnation regression retains the real finite-difference assembly and
  banded solve for the indefinite middle block `[-1, 2; 2, -1]`: its negative
  step projects to zero but leaves positive R, so it must request Picard.
  Separate tests cover scaled Phi acceptance, continued damped progress,
  exact active zeros despite simulated one-ULP LU roundoff, and the captured
  station's branch criterion `n <= -R/abs(J_ii)`, including damped active rows.
- Edge dimensional `nu_tilde` is conserved. The viscosity exponent is a central
  difference in log temperature with half-width `1e-5`, valid for any viscosity
  protocol implementation. Isentropic density and temperature derivatives give
  `Lambda`. Fixed and adiabatic energy rows use the same residual in both
  solvers; the distinction above concerns the interior SA clipping.
- All 18 frozen ZPG rows satisfy 0.3%: maximum Cf 0.0553%, delta99 0.0645%,
  theta 0.0185%, and cold-wall heat coefficient 0.0283%. Recovery-temperature
  rows use an absolute heat difference scaled by Cf. Zero SA agrees with
  laminar Cf within 0.0142%. The reference's own grid uncertainty is about 0.1%.
- M2--3 momentum balance has median error 0.00437%. With 600 common output
  stations, energy balance errors are 0.2999% (default), 0.3386% (streamwise
  refinement), and 0.0873% (wall-normal refinement). Wall-normal discretization,
  including quadrature, dominates this defect. No experimental accuracy is inferred.
- The external NASA TMR bump comparison uses wall Cp as edge pressure, arc
  length, and assumed freestream coefficient normalization. Its measured
  maximum difference is 2.497%; the generated regression bound is 3.0%.
  The comparison and its exclusions are in the generated verification table.
- `nu_tilde_ratio` is copied and read-only in public results; laminar results
  use None. First-node y+ uses wall density and viscosity. The units, signs,
  applicability, V1 and V6--V11 measurements, and mutation tests
  are recorded in section 11 of `marching_boundary_layer_design.md`.

### SA outer-grid and transport limits

The Reynolds-scaled grid completes ZPG runs at Re_s = 3e8, 1e9, and 1e10
for M = 0.3, 2, and 6. At 1e9, comparison with an outer spacing cap of 1
and an extent of at least twice eta99 gives maximum Cf differences of
0.00120%, 0.00157%, and 0.00524%; heat-flux differences are 0.19503%,
0.00635%, and 0.00681%, respectively. The frozen reference ends at 1e8;
above it these are numerical comparisons only. At 1e10 the first-node y+
is 2.303, 1.911, and 0.993: the first two runs warn and do not establish
near-wall convergence. The grid estimate does not cover every compression.
M3 to 1.5 over a tanh width of 0.1 m, with local Re_s = 1e9 at s = 1 m,
outgrows its eta extent of 476.135 near s = 1.865 m. The outer-boundary
check stops the march there.

The added streamwise refinement checks use the same simultaneous tightening
of all limits as V9. M3 to 2, width 0.1 m and Tw = 300 K, changes Cf by
0.005283% and heat flux by 0.008101%. Adiabatic M2 to 3 changes Cf by
0.002497% and wall temperature by 0.001685%. Its momentum-balance median
error is 0.003724% (0.003723% refined). With zero wall heat flux, the energy
diagnostic is the absolute derivative of the enthalpy-deficit flux,
`|d/ds int rho u (H-H_e) dy|`, rather than a relative error divided by q_w.
Its median is 50.546 W/m² (53.895 W/m² refined), or 8.220e-7
(8.223e-7 refined) after division by `rho_e U_e cp T0`. These diagnostics
have no added acceptance threshold. The M3 to 1 separation position moves
from 1.04051134018 m to 1.04051326840 m, a change of 1.9282e-6 m.

The outer SA profile can oscillate even when wall quantities change little.
The following measurements use T0 = 450 K, p0 = 200 kPa, 201 input points
over 0--2 m, and the final station's n profile above node 30. The refined
outer grid retains the first spacing and ratio, reduces the cap to 0.25,
and retains at least the original extent. Relative amplitude is the
half-excursion between neighboring outer extrema divided by their mean,
`|a-b|/(a+b)`. At each extremum the smaller of the two adjacent excursions
is used, excluding the first and last extrema and the main profile peak,
so the monotone front drop is not counted as a wiggle. Sign-change counts
include changes at roundoff level; clipped interior nodes are counted
separately. Cf and delta99 differences are maxima over all output stations.

| Case, wall, initial SA ratio | Direction changes / clipped nodes, base → fine | Maximum relative amplitude, base → fine | Cf difference | delta99 difference |
| --- | --- | --- | ---: | ---: |
| M2→3, width 0.1 m, adiabatic, 3 | 14 / 0 → 4 / 0 | 11.469% → 0.00000275% | 0.005974% | 0.037897% |
| M1.2→4, width 0.05 m, adiabatic, 3 | 16 / 0 → 21 / 0 | 99.981% → 25.947% | 0.005003% | 0.050265% |
| M2→3, width 0.1 m, adiabatic, 0.001 | 46 / 1 → 48 / 1 | 73.420% → 61.113% | 0.005473% | 0.048579% |
| M2→3, width 0.1 m, 300 K, 0.001 | 48 / 2 → 48 / 1 | 100% → 59.305% | 0.006195% | 0.049224% |

The last two cases still clip an interior node after refinement. Local n
is not converged merely because Cf and delta99 are stable. The convection
scheme is unchanged; outer-front transport and clipping remain known limits.

The approximately 10 s performance target is not met for strong expansions.
On the final code, using the same thermodynamic conditions and input points:

| Expansion and tanh width | 300 K wall: seconds / Picard fallbacks | Adiabatic: seconds / Picard fallbacks |
| --- | ---: | ---: |
| M1.2→4, 0.05 m | 23.63 / 1004 | 24.97 / 1062 |
| M2→6, 0.1 m | 15.97 / 493 | 16.74 / 504 |
| M0.3→0.8, 0.1 m | 13.17 / 211 | 9.92 / 0 |

M1.2→4 and M2→6 also trigger the thin-layer pressure-gradient warning.
These timings and stable wall quantities do not establish physical accuracy
where pressure or wall-temperature changes occur on a boundary-layer scale.

## Aerodynamic heating: notes moved from the user manual

This literature check of the turbulent attachment-line correlations was part
of `docs/models/aerodynamic_heating.rst`.

Comparison with Poll (1981) itself remains incomplete;
identity with the original paper is not assumed. Bénard et al.'s separate
correlation, equations (11)--(12), is **not implemented**: the `e0` subscript
in its Stanton density `\rho_{e0}` is undefined in the nomenclature,
and the high-Reynolds-number Stanton plateau in figure 2, titled "Typical
transition case", reads approximately `St=2.0\times10^{-3}` at
`\overline R^*\approx330`--380. The figure conditions are
`\Lambda=60^\circ`, `M_\infty=7.14`,
`T_w/T_0\approx0.38`, and `M_e\approx3.28`.
Evaluating equation (11) at `\overline R^*=350`, assuming
`T_0=750` K within the paper's 650--800 K wind-tunnel range, and using
this library's air viscosity model to convert from the `T^*` to the
`\overline T` Reynolds-number basis, gives
`St\approx4.6\times10^{-3}` (about 2.3 times
the figure reading). This is inconsistent with the reported 11.7% standard
deviation. No inferred coefficient correction is used.

Tests verify analytical relations, algebraic equivalence with CR-181893
(5.12), direct edge-state calculations, and numerical regression values.
Independent comparison against experimental data has not yet been performed.
The existing laminar APIs and their angle conventions are unchanged.

## Flat-plate boundary layer: notes moved from the user manual

These paragraphs were part of `docs/models/flat_plate_boundary_layer.rst`
until that page was split and reduced to the models, their limits, and their
use. They are kept here unchanged in substance.

### Model-selection input handling

Set `compressible` independently of the turbulent model. It accepts Python
`bool` and NumPy `bool_` values. Thermal inputs (`mach`,
`edge_temperature`, and `wall_temperature`) require `compressible=True`;
Mach and edge temperature are mandatory in that case. Temperatures use K.

Other combinations raise `ValueError`; an irrelevant selection is never
silently ignored. Each member of
`aerophysics.boundary_layer.TurbulentModel` has a `compressible`
property indicating which setting it requires.

### One-fifth-power thickness after a specified transition

For `TRANSITIONAL`, `x` and `Re_x` in these equations are
`x_\mathrm{eq}` and `Re_\mathrm{eq}` of the equivalent turbulent
plate described under :ref:`flat-plate-specified-transition`, the same origin
used for friction. The power-law `\theta_m` is then 1.3--2.2% above the
Blasius value at transition for `Re_\mathrm{tr}=5\times10^5`--
`5\times10^6`, the power law's own mismatch with Schlichting friction
at `Re_\mathrm{eq}\approx1.4\times10^5`--`6\times10^5`.
`\delta_{99}` jumps there, because the assumed profile changes from
Blasius to the one-seventh-power shape.

### Specified transition: root finding

For Van Driest II with `w=1`, `k/s=F_\theta`, the momentum-thickness
Reynolds number transformation. The Schlichting root is found by Brent's
method in `\ln R`; the Willems average equation is inverted explicitly as
`R=S(\log_{10}S/0.242)^2` with `S=R\bar C_{f,i}`. Thickness
correlations, including the `SCHLICHTING` power law, use
`x_\mathrm{eq}` and `Re_\mathrm{eq}`, and
`effective_reynolds_number` returns `kRe_\mathrm{eq}` at turbulent
points (`k=1` for `compressible=False`). `Re_0` for
`Re_\mathrm{tr}=5\times10^5` is `1.417\times10^5` with
Schlichting friction.

### Specified transition: change from the earlier mixed average

Earlier versions kept the leading-edge origin for turbulent friction and
used the Prandtl-type mixed average
`\bar C_{f,t}(Re_x)+(Re_\mathrm{tr}/Re_x)[\bar C_{f,l}-\bar C_{f,t}]_{Re_\mathrm{tr}}`.
That average conserves drag, but its local friction belongs to a layer whose
momentum thickness is 2.8 times the laminar value at
`Re_\mathrm{tr}=5\times10^5`. With Schlichting friction, the local
coefficient just after transition is now 26% higher for
`Re_\mathrm{tr}=5\times10^5` and 42% higher for
`5\times10^6`. Average friction is 6--9% higher at
`Re_x/Re_\mathrm{tr}=1.2`--2 for `Re_\mathrm{tr}=5\times10^5`, up
to 15% for `5\times10^6`, and less than 1% higher beyond
`100Re_\mathrm{tr}`. Fully turbulent and laminar results are unchanged.

Drag continuity, the exact drag derivatives, the `SCHLICHTING` theta offset,
and the fully turbulent limits are verified analytically in the test suite. The equivalent plate has not been compared

### Compressible laminar similarity solver

The calculation groups exactly equal states (Mach, edge temperature, wall
condition, Pr, gamma) within each call, without rounding or a persistent
state cache. Conservative second-order finite differences and damped Newton
iteration operate in batches of at most 256 states. A stretched, fixed
161/321-node pair and Richardson extrapolation are followed by 321/641
nodes only for uncertain states. The solver checks residuals, boundaries,
positive temperatures and viscosities, grid differences, far-field decay,
and the independent momentum integral. Remaining states use
`scipy.integrate.solve_bvp(tol=1e-8)` with an enlarged-domain check.
Unresolved states raise `ModelRangeError` with their original input indices
and state. Blasius and Walz profiles supply initial guesses only; they are
not the returned physical approximation. There are no interpolation tables.


## Marching boundary layer: notes moved from the user manual

These paragraphs were part of `docs/models/marching_boundary_layer.rst` until
the user manual was reduced to the model, its limits, and its use. They are
kept here unchanged in substance. Some repeat material in the sections above.

### Station-solve Jacobians

The laminar and SA-noft2 station solves use analytic Jacobians for the exact
built-in Sutherland and Keyes viscosity model types. Other viscosity models,
including subclasses, use finite-difference Jacobians. SA differentiates the
active nonsmooth branch, with right-sided derivatives at exact ties. These
linearizations solve the same discrete equations with the same convergence
tolerances; they do not change the physical model or its applicability.

### SA unknowns, active-set Newton, Picard fallback and eta grid

The five unknowns per node are `(f, u, g, n, Y)`. Projected active-set Newton
solves the local interior SA conditions `n \geq 0`, `R \leq 0`,
and `nR=0`, where `R` is the discrete SA residual. Positive
interior values satisfy the SA equation; zero values admit a nonpositive
residual. Convergence requires a small applied update, the same active set in
two consecutive linearizations, and
`\max_i|\min(n_i,-R_i/|J_{ii}|)|/\max(1,\max n)<10^{-9}`
over interior SA rows. This uses the existing pre-update linearization and the
same `n` scale as the update tolerance. If Phi fails this bound, a nonzero
applied update permits further iteration within the iteration limit; a zero
update requests Picard. Projection or damping stagnation cannot count as
convergence. Undamped active rows are assigned exactly zero; damped rows
retain the partial step toward zero. Newton falls back to relaxed Picard
with nonnegative `n` clipping.
At convergence, an active or projected-zero first interior SA node also
requests Picard. Such wall-attached zero regions select the trivial partially
laminar branch; the intended clipped set is the outer turbulent front.
The explicitly zero edge-SA case retains legacy Newton for its intentional
zero-SA limit.
Picard projects the result of a coupled tridiagonal solve, satisfying
`\min(n,A(n)^{-1}R)=0` at its fixed point. This differs from Newton's
local complementarity condition: clipping can leave residuals at positive
neighbors of clipped nodes. Remeasured with the analytic Jacobian in M2--3
fixed-wall and adiabatic marches, the maximum difference in `n` was less than
`2.12\times10^{-6}` of maximum `n`; other output fields differed by less than
`7.08\times10^{-8}` of their respective maxima. These are not general
bounds: remeasured with the wall guard at `Re=10^{10}`, M2, edge 220 K
and wall 300 K, the largest `n` difference is 8.099
(`6.03\times10^{-5}` of the profile maximum), two nodes before clipping.
The last positive node changes from `n=1996.166` to `2002.957`.
Velocity and temperature profile differences are below
`8.2\times10^{-8}` of their maxima. These comparisons do not establish grid
convergence or physical validation of the outer SA front.
If both solvers fail the streamwise step halves.
The eta grid has first spacing 0.0007 and spacing
ratio 1.04. Before marching, define `Re_equivalent = max(R)^2/2` over the
output stations and `scale = max(1, (Re_equivalent/1e8)^0.42)`. The extent is
`150 scale` and the outer spacing cap is `scale`. Up to Re_equivalent=1e8
this is the original 312-node grid. These are geometric spacings, not
geometrically positioned nodes. The outer-boundary check remains active:
a compression can grow beyond this estimate and stop the march.

### Specified transition: implementation

SA-noft2 is solved over the entire domain on the same grid as `TURBULENT`.
Only `c_{b1}\widetilde S n` is multiplied by zero for
`s \le s_{tr}` and by one downstream. Destruction, diffusion, convection,
and the outer boundary condition remain active. The implementation applies
the multiplier in the residual, its numerical Newton Jacobian, and Picard's
lagged source and source-based relaxation. No profile or grid is replaced at
the switch. The march passes through the transition position internally; the
returned stations remain exactly `edge.distance[1:]`.

### Specified transition: upstream differences from the laminar solver

The upstream solution retains small eddy viscosity supplied from the edge.
For `edge_sa_ratio=3` the measured maximum upstream departures from the
laminar solver are 0.374% in Cf and 0.356% in wall heat flux over the four
uniform-edge verification cases (M0.3, M2 cooled/adiabatic, M6 cooled).
These are case-specific numerical measurements, not physical accuracy bounds.
Reducing the seed reduces eddy viscosity but does not remove differences
between the two existing discretizations: at M2, Te=220 K and Tw=300 K,
`edge_sa_ratio=1e-3` gives a 0.05304% upstream heat-flux difference from
`LAMINAR`. This difference is recorded as a discretization diagnostic, without
a pass/fail criterion. The laminar-limit check retains the V6 Cf target of 0.05%.
The transition implementation is checked separately against
`TURBULENT(edge_sa_ratio=0)` on the same SA grid and with the same operators,
using a maximum relative-difference limit of 1e-8% for both Cf and heat flux.

### Thin axisymmetric geometry: measured relations and synthetic shoulder

The measured cone relation errors are below 6e-12 percent. With 801 nodes,
the variable-radius laminar mapping differs by at most 0.000827 percent in Cf;
the interpolation is performed in different distance variables. Axisymmetric
momentum residual medians are below 0.05 percent. Streamwise and wall-normal
refinement are recorded separately in the development memo, section 12.
The R&M 3633 worked example reaches delta99/r_w=0.4774
and differs from measured Cf by up to 0.0004443 (measurement accuracy
+/-0.0002). This comparison does not establish a physical validity range.

The recorded synthetic cone-to-cylinder shoulder stress measurement (development
memo, section 12) used a tanh pressure drop; it is **not a physical shoulder edge
condition**. With input nodes locally spaced
at most width/5, the pressure-length warning fires for SA at width 0.001 m
(maximum delta99 times abs(d ln p/ds) 1.24051--1.26509), and for both regimes
at width 0.0002 m (1.33059--6.29451). Widths 0.02 and 0.005 m do not trigger
that warning. All 16 cases completed in 2.82--8.85 s in the recorded run,
with zero Picard fallbacks. This is a timing observation, not a performance bound.
The curvature maxima occur near the tip or downstream end, so that warning
alone does not identify a shoulder problem. At s=0.31 m, relative SA-variable
oscillations reach 1.46404%; the final station has no evaluable extrema triplets.
Some oscillations occur without a pressure-length warning. These applicability
warnings therefore do not establish local profile convergence or physical
accuracy; the need for an additional oscillation diagnostic remains open.

### Viscosity derivative, outer nu_tilde oscillation and solver cost

Turbulent flow requires the viscosity model to accept `T_e * exp(+/-1e-5)`
for an edge viscosity derivative, including when `edge_sa_ratio=0`. Near a
temperature-range limit, a model that raises outside its range can fail at
these temperatures, and a model that warns outside its range can warn.
Laminar flow does not evaluate this SA derivative.

The outer part of `nu_tilde_ratio` can oscillate and is not grid-converged
in strong expansions or with a small `edge_sa_ratio`. In the measured
cases, relative oscillation amplitudes reached about 100% on the default
grid and about 60% after outer-grid refinement, measured as half the excursion
between neighbouring extrema divided by their mean. Wall quantities and
delta99 were insensitive to that refinement in these cases: Cf changed by
less than 0.01% and delta99 by less than 0.06%. These measurements are not
general bounds.

With the legacy positivity-damped Newton, M1.2 to 4 took about 25 s on one
machine because of Picard fallbacks. Active-set Newton substantially reduced
fallbacks in the measured verification tier; timings remain machine-dependent.

## Flow sequence

- Temperature restoration from validated solver ratios repairs at most one
  ULP of boundary overshoot. Input temperature checks remain strict.
- Stagnation temperature is bracketed upwards from the static temperature,
  including for custom harmonic-oscillator gases without an applicability
  range.
