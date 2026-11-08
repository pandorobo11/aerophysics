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
  laminar flow with exact built-in `SutherlandModel` or `KeyesModel` types.
  Other viscosity models (including subclasses) and SA retain three-color
  finite differences. Both paths assemble directly into expanded `gbsv` bands
  and use the same banded LAPACK solve. Convergence requires a maximum
  update below `1e-10` in `u` and in the relative `g`.
- The analytic path uses the model's viscosity value and validation at current
  temperatures, with a private derivative in Pa s/K. For Sutherland,
  `mu_T/mu = 1.5/T - 1/(T+S)`; for Keyes, with
  `a = A 10**(-B/T)`, `mu_T/mu = 1.5/T - (1+a ln(10) B/T**2)/(T+a)`.
  Thus `C_g = (T_e mu_T/mu_e - C)/g`. This derivative enters both face
  fluxes, dissipation, and the adiabatic half-cell balance. Convection and
  streamwise BDF terms are differentiated with edge and lagged states fixed.
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
  with no SA positivity fraction. Exact boundary values are restored, and
  convergence also requires an unchanged active set in the last two
  linearizations. At convergence, an active **or projected-zero** first interior
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
  applicability range has changed. The observed difference on one captured
  station was `1.26e-5` in n (`3.32e-8` of its maximum), but this is not a
  general bound. Remeasured with the wall guard at `Re=1e10`, M2, edge 220 K
  and wall 300 K, the largest n difference is `8.098729`, or `6.03083e-5` of
  maximum n, at node 301 (zero-based), eta `632.73694`, s=100 m. The next,
  last positive node before clipping (302, eta `639.65525`) changed from
  `1996.165999` to `2002.957444`, a difference of `6.791445`; node 303 is
  effectively zero in both solutions. Velocity and temperature profile
  differences are `6.03350e-8` and `8.18300e-8` of their maxima; the wall
  velocity slope changes by `1.53615e-7` relative, and delta99 by `7.22122e-7`.
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
- On the rebased rounded-transition reproducer (transition at 0.1 m, wall
  250 K), the unguarded Newton selected a zero region at nodes 1--82 on the
  472nd SA station attempt, xi `3.4502924398659043e-4`. Its wall velocity slope
  was `0.469137`, versus Picard's `0.556645`. The guard rejects this state.
  The guarded reference march's Cf and heat flux match the legacy march
  exactly in the measured comparison, within the existing `1e-6` bound.
- Counts were remeasured on the **88-test `solver_verification` tier as of
  this change**, including transition and active-set comparison tests, using
  four workers. Two runs set the class default to legacy or guarded Newton;
  tests that explicitly toggle the private switch retained those overrides.
  Thus both totals include legacy comparison calls, and the legacy-default
  run also includes explicitly requested active-set solves. SA Newton attempts
  were `137,484` in both runs; successes were `131,896 -> 136,845`.
  Zero-damping failures were `5,586 -> 600` (first update `4,716 -> 526`, later
  `870 -> 74`), wall-guard rejections `2 -> 32`, and iteration exhaustion
  `0 -> 7`. Both runs had zero invalid-diagonal, nonfinite-update, and caught
  exception failures. Picard fallbacks were `5,586 -> 637`, total calls
  `5,587 -> 638`, and sweeps `806,235 -> 107,160`; every Picard call converged.
  Tests also call Newton/Picard directly, explaining why failure, fallback,
  and call counts differ. Counts include repeated worker-local fixtures.
- In the guarded-default run, the actual active-set path had `128,247` attempts,
  `128,208` successes, 32 wall-guard rejections and 7 exhausted solves. The
  legacy path (explicit comparisons plus zero-SA dispatch) had 9,237 attempts
  and 600 zero-damping failures. Of all returned active-set stations, 5,021
  had a nonempty active/bound set, each containing one node. Counting both
  selected active rows and projected zeros, **zero** stations had nodes
  disconnected from the outermost active region, and **zero** had active nodes
  at `u < 0.99`. No stricter region rule was added. Fallback stations retain
  Picard's different fixed point. Generated verification records remain
  unchanged at stored precision. The serial enthalpy-refinement test was
  remeasured on the rebased code: 40.39 s with legacy Newton and 32.50 s with
  guarded Newton (both passed). These elapsed times are machine-dependent.
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

## Flow sequence

- Temperature restoration from validated solver ratios repairs at most one
  ULP of boundary overshoot. Input temperature checks remain strict.
- Stagnation temperature is bracketed upwards from the static temperature,
  including for custom harmonic-oscillator gases without an applicability
  range.
