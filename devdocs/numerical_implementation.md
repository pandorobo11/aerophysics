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
  iteration. The Jacobian is block tridiagonal, assembled by three-color finite
  differences, and solved with banded LAPACK. Convergence requires a maximum
  update below `1e-10` in `u` and in the relative `g`.
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

## Flow sequence

- Temperature restoration from validated solver ratios repairs at most one
  ULP of boundary overshoot. Input temperature checks remain strict.
- Stagnation temperature is bracketed upwards from the static temperature,
  including for custom harmonic-oscillator gases without an applicability
  range.
