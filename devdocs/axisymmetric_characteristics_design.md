# Axisymmetric rotational characteristics for external body flow

Implementation and verification record, 2026-11-08. This work stays on
`feature/axisymmetric-characteristics`, based on `459b6bc`. No boundary-layer
module, shock public API, generated asset, reference data, CI configuration,
package version, commit or remote branch is changed.

## Contract and conventions

`axisymmetric_body_flow(x, r, *, mach, pressure, temperature, gas=AIR,
characteristic_count=88, fan_angle_step=pi/180)` returns
`AxisymmetricBodyFlowResult`. The requested names are retained. The result
also retains `gas: PerfectGas` so non-air flows produce consistent boundary-layer
edges. All result arrays are independent owned float64 arrays and read-only,
including when the result dataclass is constructed or copied separately.
The private `_solve(starting_distance=...)` option added during review starts
a genuine MOC net inside the cone; it is absent from the public API.

Lengths are metres, static pressure Pa, static and total temperature K. Angles
are radians from +x toward +r; a boattail has negative angle. Surface distance
is cumulative Euclidean meridian length from the tip. The sharp-tip state
is the cone surface limit, not a resolved singular tip cell. Upstream scalar
pressure and temperature must be positive and finite; Mach must exceed one.

The meridian starts at (0,0), has strictly increasing x, positive downstream
radii, and nonincreasing segment angle. Equal-angle vertices are retained as
output stations but do not create redundant wave fronts. The body is sharp,
axisymmetric, at zero incidence, with no swirl, and has an attached weak
nose shock. Gas is calorically perfect. No viscous displacement feedback,
concave-corner shocks, blunt noses, incidence, internal flow, chemistry or
shock-shock interaction is included.

At a convex corner all fan states share exactly the same physical coordinate.
The public array retains the **post-expansion limit only**; duplicating the
pre-expansion limit would violate `BoundaryLayerEdge`'s increasing-distance
contract. The result privately retains each vertex's pre-fan pressure, including
later corners where that pressure comes from the marched wall state.
`boundary_layer_inputs(*, corner_length=None, viscosity_model=AIR_VISCOSITY)`
replaces the unreleased `boundary_layer_edge` method, which is removed. It
returns the named tuple `BoundaryLayerInputs(edge, wall_radius)`, exported
from `aerophysics`. The added radius array is owned and read-only.

For each convex vertex at surface distance s_c the conversion inserts an
upstream-limit point at s_c-w, then a point at s_c-w/2 whose pressure is the
geometric mean of the upstream and downstream limits. The original vertex
retains its downstream state. The default w=0.001*s_c is set independently
at each corner; explicit `corner_length` is a positive finite scalar in metres
applied to every corner. If w reaches the preceding public surface station,
or floating-point arithmetic cannot distinguish the two inserted points,
`ValueError` is raised without deleting or moving any original stations.
The complete resulting radius array is interpolated linearly in surface
distance through the original geometry vertices, hence exactly represents
the piecewise linear meridian at every edge station. Constant body-streamline
total pressure and total temperature reconstruct all inserted states.

These inserted points control the pressure interpolation into a boundary
layer, not the inviscid MOC net. Public flow arrays and physical fan coordinates
are unchanged. The ramp is not a finite-radius shoulder or a resolved viscous
corner interaction. The latest measurements and warning limitations are in
**Approved corner-input API and revised default resolution** below.


## Equations and their derivation

Let q be speed, theta flow angle, mu=asin(1/M), and g=gamma. For a uniform
upstream total enthalpy, T0 is the same on every streamline. Store
P=ln(p/p_inf) and L=ln(p0/p_inf). Thermodynamic reconstruction is

```
T/T0 = exp((P-L)*(g-1)/g)
M² = 2/(g-1) * (exp((L-P)*(g-1)/g) - 1)
q² = 2 cp (T0-T)
s/R = constant - L                  (uniform T0)
```

With Ds along a streamline and Dn normal to it toward increasing angle,
normal momentum, streamline momentum, continuity and streamline entropy
conservation yield

```
Ds(theta) + Dn(p)/(rho*q²) = 0
Dn(theta) + (M²-1)*Ds(p)/(rho*q²) = -sin(theta)/r
Ds(p0) = 0
```

Combining the first two equations along directions theta +/- mu gives

```
dr/dx = tan(theta +/- mu)
d(theta) +/- K dP = H_+/- dx
K = sqrt(M²-1)/(g*M²)
H_+/- = -/+ sin(mu)*sin(theta)/(r*cos(theta +/- mu))
```

K is dimensionless and H has units 1/m. The private planar verification mode
sets H=0 and initializes with the uniform oblique-shock wedge state. It is
not a second public solver. These pressure-angle compatibility equations do
not assume constant entropy normal to streamlines. L is determined separately
by tracing a streamline; using one global p0 would be an incorrect
irrotational substitution for curved-shock flow.

## Initial line and analytic cone region

The entire initial conical portion, including collinear input vertices, uses
`conical_shock`'s analytic similarity state at the surface. This preserves
constant pressure/Mach without numerical drift and makes M1 an API/initialization
regression, not evidence of MOC accuracy. Likewise M5 at the exact vertex is
constructed from Prandtl--Meyer, and M7's constant-edge cone segment is an
API/initialization and coupling check. None of those construction-based
agreements is used to establish the MOC march's accuracy. The review-stage
forced-cone test below supplies that independent verification.

For the starting field, integrate the **existing** `_taylor_maccoll_rhs` from
shock polar angle beta to the cone angle, using the post-shock normalized
radial/polar velocities. No Taylor--Maccoll equation was duplicated and
`shocks.py` is unchanged. DOP853 uses rtol=1e-11, atol=1e-13, dense output.
The starting characteristic is the leading C+ ray through the first physical
corner, lying in the undisturbed conical field, rather than a vertical x line.
Its shape is obtained by integrating

```
d(log x)/d(alpha) = 1 / (lambda_plus*cos(alpha)² - sin(alpha)*cos(alpha))
r = x*tan(alpha)
```

from cone angle to shock angle, with the same tolerances. `characteristic_count`
is the number of equal polar-angle intervals on this line. This leading wave
can meet the shock beyond the body's end: for cone length 1 m the M3/9.46-degree
intersection is x=4.15866 m, and the M6/4.76-degree intersection
is x=4.48635 m. Conical surface samples use count+1 uniform axial positions
plus any input cone vertices. The first line ends on the fitted shock and is
retained privately for conservation tests.

## Unit operations and entropy transport

- **Interior:** C+ from the lower point and C- from the upper point intersect.
  A predictor uses upstream coefficients. An iterated trapezoidal corrector
  averages each endpoint's characteristic slope, K and H, then solves the two
  linear pressure-angle equations. Trace the streamline backward to the known
  edge of this cell. Its direction is the average of target and interpolated
  foot direction. Interpolate L there and carry it unchanged to the new point.
  Direction and thermodynamic reconstruction are corrected together. Linear
  interpolation of the foot entropy is second order for a smooth field.
- **Wall:** intersect the incoming C- with the straight meridian segment.
  Impose the segment angle and the constant body-streamline L. The one
  characteristic equation determines pressure. Endpoint averages use the
  same iterated trapezoidal corrector.
- **Specified vertex/end point:** solve inversely for the incoming foot on
  the first known C+ edge. Quadratic interpolation of geometry and state
  avoids a large endpoint interpolation error. A two-point linear stencil is
  used when only two outer-front points remain. Collinear vertices are sampled
  by this inverse operation without perturbing the marching front.
- **Shock:** solve for local beta using incoming C+ compatibility. For each
  trial beta, apply the local normal-shock density/pressure/total-pressure
  relations with upstream normal Mach M_inf*sin(beta), and turn the downstream
  velocity accordingly. Intersect C+ with the shock segment whose slope is
  the average of tan(beta_old) and tan(beta_new). A secant root solves
  compatibility to an angle tolerance of 1e-12 rad. Newly captured streamlines
  obtain this local post-shock L, not the cone-surface value.
- **Centered fan:** partition the turning angle into equal increments no
  larger than `fan_angle_step`. At the zero-length vertex the integrated
  axisymmetric source vanishes: use `prandtl_meyer_expansion` from the
  pre-corner surface state for each exact fan state. Each such state launches
  a C+ row; away from the vertex the full axisymmetric source applies.

Convergence of each geometric/thermodynamic corrector is required, with a
2e-10 tolerance on relative x/r and absolute theta/P/L updates; streamline
angle correction uses 2e-12 rad. At most 30 iterations are allowed. This is
an iterated second-order predictor-corrector, not a first-order explicit
march. There is no artificial constant-p0 reset in the interior. An exact
parallel shear-flow test (constant pressure/direction, transverse p0 gradient)
checks rotational transport independently of shock fitting.

A row progresses outward from the wall along C+, using incoming C- parents.
The next ordinary wall point uses the first interior C-; a fan adds outgoing
rays at the same wall coordinate. A tail fan ray need not ever reach the shock.
Once a row extends beyond the required axial region it is truncated, retaining
points for the inverse endpoint stencil. There is no imposed outflow pressure.
Returned shock stations can extend beyond the body end; they are actual fitted
characteristic intersections, not an interpolated or invented shock endpoint.
The implementation retains only the active row and the last complete
wall-to-shock row privately. It does not publish a field mesh API.

## Stops and numerical limitations

The boundary-layer pressure-length warning evaluates only its output stations.
PCHIP can have zero slope at both ends of a sharp pressure ramp, so endpoints
alone can miss its large interior gradient. This module adds a logarithmic
pressure midpoint so the gradient is sampled, and tests the public marcher's
warning for both laminar and SA default cone--cylinder inputs. This remains a
general limitation of the boundary-layer module, which is not modified here.
A warning marks loss of the thin-layer approximation; it does not repair it.


Invalid arrays/states/resolution raise ValueError. Concave corners, detached
or unresolved starting cone shocks, any sonic/subsonic internal state,
non-downstream characteristic direction, parallel/backward intersections,
streamline feet leaving their cells, ODE failures, corrector exhaustion,
more than 100000 surface stations or 10000 fan increments raise ModelRangeError.
Numerical errors carry x and r in metres. Failed calculations return no
partial success result. Failure-path tests cover these numerical safeguards.

The method does not automatically refine a poorly conditioned net. A backward
**converged** intersection can indicate either compression-wave coalescence or an
under-resolved characteristic cell; the diagnostic is not proof of a physical
secondary shock. For example, the coarse M3 body (0,0), (1,0.2), (10,1.1)
stops near (8.58015,3.16398) m; this case is retained as a failure regression.
Longer bodies, stronger boattails and near-sonic attached solutions require
individual resolution and applicability checks. The headline ten-diameter
cases are the acceptance domain measured here. No pressure clipping,
continued subsonic march, skipped tests or relaxed acceptance limits are used.

## Original verification procedure and acceptance rationale

This section records the pre-review count-24 setup. Current defaults and
corner-input verification are specified in the final section.

All air cases use gamma=1.4, p_inf=10000 Pa, T_inf=220 K. M3/M4 use cone length
1 m, r_shoulder=tan(cone angle), and cylinder length 20*r_shoulder (ten
diameters). Refinement doubles count and halves fan step with an exactly
nested fan division: count 24/48/96/192 and step angle/(count/2), with a
1e-7 relative roundoff allowance to avoid accidental extra fan subdivision.

M3 compares pressure on 501 fixed axial stations from the downstream shoulder
limit through the cylinder endpoint using cubic interpolation separately on
the cylinder. Do not interpolate through the shoulder jump when estimating
spatial order. Differences are maximum absolute pressure differences. The
fine-grid Richardson estimate is the last difference / (2²-1), not a
rigorous physical-error bound. The final observed order must be 1.8--2.2;
the coarsest difference / p_inf must be below 1e-4 and the finest estimate /
p_inf below 2e-6. These allowances cover interpolation/roundoff effects around
the measured order and have substantial margin over the measured differences,
without changing the requested approximately-second-order criterion.

M4 independently integrates the normal mass flux on the last complete
wall-to-shock characteristic polygon:

```
integral 2*pi*r*rho*(u dr - v dx)
```

and compares it with rho_inf*U_inf*pi*r_shock². Trapezoidal flux quadrature
is independent of the pressure compatibility residual. The 1e-4 relative
bound is about 2.4 times the largest measured default error. In the short
M6 case the leading shoulder wave reaches the shock beyond the body end,
so this complete line is still the analytic initial line. A separate longer
M6 calculation, x_end=5 m and r=0.1 m, count 48, fan step 0.25 degree, checks
nonuniform p0, curved-shock propagation and the same mass bound.

M1 uses the requested 1e-4 pressure/Mach bounds, with 1e-10 m for shock-ray
geometry in these metre-scale cases. M2 bounds post-corner planar Mach
relative error at 2e-6 (measured about 1.1e-6 at a 0.5-degree fan step); its
wedge state matches the existing oblique-shock API to 1e-12. M5 requires
2e-12 relative agreement with the exact local PM turn. Edge reconstruction
is checked at 2e-14 relative tolerance. M6 includes the requested inputs and
additional malformed arrays, unsupported gases and failed numerical processes.

M7 uses cone length 0.3 m, radius 0.05 m, cylinder length 1 m (ten diameters),
M_inf=3, wall temperature 300 K, and the default MOC resolution. Both laminar
and SA-noft2 must reach the last station with positive Cf. On cone stations
s<0.25 m, compare laminar Cf with an otherwise identical flat-plate march;
the relative error from the Mangler sqrt(3) factor must be below 5e-4, matching
the existing axisymmetric boundary-layer verification allowance. This checks
coupling, not corner interaction or physical validity of a thick downstream
SA layer.

## Original measurements (count 24, before the corner-input API)

The tests print M1--M7 records for reproducibility. Measurements on macOS
27.0.1 arm64 / CPython 3.13.5, with `.venv/bin/` executables:

| ID | Measured result | Decision |
| --- | --- | --- |
| M1, M2/10 degrees | Pressure relative error 0; Mach 1.11e-16; shock ray error 0 m | Pass |
| M1, M6/5 degrees | Pressure relative error 5.00e-15; Mach 8.88e-16; shock ray error 5.55e-16 m | Pass |
| M2, planar M3/10-degree wedge and 10-degree turn | Post-turn Mach maximum relative difference 1.08527e-6; wedge pressure/Mach agree to 1e-12 | Pass |
| M3, M3/9.46 degrees | Max pressure differences 0.695207, 0.136877, 0.0347727 Pa; orders 2.34456, 1.97686; finest error estimate 0.0115909 Pa | Pass |
| M3, M6/4.76 degrees | Max pressure differences 0.376064, 0.0967906, 0.0247265 Pa; orders 1.95804, 1.96881; finest error estimate 0.00824215 Pa | Pass |
| M4, M3 / M6 default nets | Signed relative mass differences +4.11483e-5 / -3.61816e-5, bound 1e-4 | Pass |
| M5, M3 / M6 | Shoulder downstream Mach 3.21571718768 / 6.29150033792; pressure 7253.94969350 / 7451.76242534 Pa; PM agreement within 2e-12 | Pass |
| M6 | Concave, detached, subsonic, repeated/decreasing x, negative radius and malformed input tests reject; numerical failures retain location | Pass |
| M7, laminar | Completes to s=1.30413812651 m; 1.135 s; Mangler Cf relative difference 9.33e-15; no warnings | Pass |
| M7, SA-noft2 | Completes to the same endpoint; 2.974 s; positive Cf; curvature applicability warning delta99/r=0.315703 at endpoint | Completion passes; thin-layer applicability is limited |

The longer rotational M6 test has a 0.0111856 rad range of fitted shock angles,
a 0.00269033 range in log(p0/p_inf) across the final complete characteristic,
and a +1.99957e-5 relative mass difference. These values test the genuinely
nonuniform-entropy implementation in addition to the initial conical field.

Default-resolution timings (count 24, step 1 degree), five calls in one fresh
Python process, excluding import time:

| Case | First call | Subsequent calls | Surface stations |
| --- | --- | --- | --- |
| M3/9.46 degrees | 0.06218 s | 0.02860--0.03061 s | 46 |
| M6/4.76 degrees | 0.04972 s | 0.02357--0.02444 s | 40 |

The cache in the existing conical-shock solver explains part of the first-call
difference. These timings are far below the several-seconds-per-condition
budget, not a portable performance guarantee. Refined count-192 solves took
0.481 s (M3) and 0.360 s (M6) without coverage instrumentation. The 51 ordinary
MOC tests before the final three error-propagation cases took 2.24 s serially;
the four marked tests contain refinement and the two boundary-layer marches.

### Commands and validation boundary

- `.venv/bin/ruff check .`: pass.
- `.venv/bin/ruff format --check .`: pass.
- `.venv/bin/mypy`: pass, 99 source files.
- `.venv/bin/pytest tests/test_characteristics.py -m "not generated_assets and not package_artifact" --no-cov -n auto`: pass; final count recorded below.
- `.venv/bin/pytest -m "not generated_assets and not package_artifact and not solver_verification" -n auto --cov-fail-under=0`: all but the existing local documentation server test pass; final counts below. `tests/gui/test_documentation.py::test_local_documentation_server` fails at `socket.bind` with `PermissionError: [Errno 1] Operation not permitted` in this sandbox. No escalation, exclusion, skip or xfail was used.
- `.venv/bin/python -m sphinx -W -b html docs docs/_build/html`: pass.
- `.venv/bin/python -m sphinx -W -b doctest docs docs/_build/doctest`: 227 doctests pass.
- `.venv/bin/coverage report --include='src/aerophysics/*.py' --fail-under=95`: repository numerical-core coverage gate passes; exact final metrics below.
- `git diff --check`: pass.

`uv`, `scripts/check.sh --full`, generated-asset regeneration/verification,
package build/clean install and the unrelated full slow solver tier were not
run, as requested. The reviewer owns the full gate. Thus a successful
unrestricted full repository gate is **not claimed**. No CI or PR was opened;
no commit, push or version change was made. The review-stage NACA TN 2250 comparison below uses computed literature
points, not experimental measurements.


### Validation totals before review

The original targeted suite passed 58 tests; the full ordinary suite passed
1320, with one sandbox-denied loopback socket bind. The core coverage.py
combined metric was 96.90%, satisfying the repository's >=95% gate. The
branch-only diagnostic was 92.56% and is **not** the gate criterion. The
review stage below independently measures the branch-only baseline at the
pre-change main commit, rather than estimating it by excluding this module.
Final review validation counts and coverage supersede these original counts.

Original implementation files (10; review additionally adds the supplied reference JSON):

- `src/aerophysics/characteristics.py`: public solver/result and private unit processes.
- `src/aerophysics/__init__.py`: public exports.
- `tests/test_characteristics.py`: analytical, rotational, conservation, refinement, input, failure and coupling tests.
- `devdocs/axisymmetric_characteristics_design.md`: this record.
- `docs/models/axisymmetric_characteristics.rst`: model, conventions and limits.
- `docs/models/index.rst`: model navigation.
- `docs/api/compressible_flow.rst`: API documentation.
- `docs/guides/boundary_layers.rst`: runnable cone--cylinder coupling example.
- `docs/references.rst`: checked bibliographic records and access limitations.
- `CHANGELOG.md`: Unreleased addition.

## Bibliographic provenance

The requested Zucrow and Hoffman title, volume 2, Wiley, copyright 1977,
ISBN 0471018066, were verified through WorldCat and KAIST library records:

- https://search.worldcat.org/title/Gas-dynamics.-Vol.-2-Multidimensional-flow/oclc/923378441
- https://library.kaist.ac.kr/search/ctlgSearch/posesn/view.do?bibctrlno=11136&se=b0&ty=B

Anderson, Modern Compressible Flow: With Historical Perspective, 4th edition,
copyright 2021, and Chapter 11, Numerical Techniques for Steady Supersonic
Flow, were verified against the publisher:

- https://www.mheducation.com/highered/product/Modern-Compressible-Flow-With-Historical-Perspective-Anderson.html
- https://www.mheducation.ca/higher-education/engineering-computer-science/mechanical-engineering/compressible-flow-gas-dynamics

The complete textbook chapters and their detailed rotational-MOC section or
equation numbers were **not** independently inspected; an accessible university
chapter PDF could not be fetched. The equations above were derived from the
stated Euler equations and checked analytically/numerically, not transcribed
with unverified page citations. The supplied NACA TN 2250 source metadata and
unaltered digitization support the computational comparison below. No
experimental comparison is claimed.


## First-review measurements (historical count 24 and old edge conversion)

The final section supersedes default/API and width-sensitivity statements in
this historical record; explicit grid studies remain applicable.

### Forced-cone march: M1-MOC

A private `starting_distance=0.05` inserts a zero-turn virtual vertex inside
a unit-length cone. The existing exact leading C+ initializer supplies data
there; every surface station with x>0.05 is then generated by the normal MOC
wall/interior/shock operations. The final x/start ratio is 20. Cone pressure
or Mach is never reset analytically during this march. A fast test explicitly
requires nonzero numerical truncation and multiple fitted shock updates, so
silently reverting to the analytic shortcut would fail.

This exposed a real numerical guard defect. In M2/10 degrees the first-order
predictor near the curved initial characteristic briefly lies upstream of a
parent. Previously the solver rejected it before the trapezoidal corrector
could converge. Causality and streamline-cell containment are now checked on
the **converged** cell; parallel/singular intersections and unphysical
thermodynamic trial states still stop. The final checks are not relaxed.
The same default initial-line spacing and unit processes are used by the
public solver; no cone-specific source or shock correction was introduced.

All errors below are maxima of absolute relative errors against
Taylor--Maccoll. Shock position error is abs(r_shock/(x_shock*tan(beta_TM))-1),
excluding the tip. The complete sampled shock trace is checked, including any
computed nodes just beyond the final wall x. Orders use log2(error_N/error_2N).

| Mach / half angle | N | Pressure | Surface Mach | Shock position |
| --- | ---: | ---: | ---: | ---: |
| 2 / 10 degrees | 24 | 1.07712e-2 | 3.80702e-3 | 1.78298e-3 |
| | 48 | 3.14631e-3 | 1.11597e-3 | 4.71870e-4 |
| | 96 | 8.40706e-4 | 2.98511e-4 | 1.18213e-4 |
| | 192 | 2.14588e-4 | 7.62164e-5 | 2.91458e-5 |
| 6 / 5 degrees | 24 | 7.14659e-4 | 1.18411e-4 | 7.02315e-4 |
| | 48 | 1.80492e-4 | 2.99143e-5 | 1.76151e-4 |
| | 96 | 4.52309e-5 | 7.49705e-6 | 4.40872e-5 |
| | 192 | 1.13155e-5 | 1.87559e-6 | 1.10257e-5 |

For M2 the three successive pressure orders are 1.7754, 1.9040, 1.9700;
Mach orders 1.7704, 1.9024, 1.9696; shock orders 1.9178, 1.9970, 2.0200.
For M6 they are 1.9853, 1.9966, 1.9990; 1.9849, 1.9964, 1.9990; and
1.9953, 1.9984, 1.9995. Final orders must lie in [1.8,2.2].

Measured default-error acceptance bounds (pressure, Mach, shock) are
(0.015,0.005,0.0025) for M2 and (0.001,0.00016,0.001) for M6: approximately
31--42% headroom above these measured errors. These are **new forced-march
bounds**, not changes to the original M1 API-initialization 1e-4 bound.
All pass. Default N=24 is not uniformly 1e-4 accurate for a long marched cone;
particularly the low-Mach cone needs refinement. The normal exact cone
shortcut remains the accurate and cheaper production choice.

### NACA TN 2250 Figure 1 comparison

The reviewer-supplied file was copied byte-for-byte to
`tests/reference_data/inviscid/naca_tn2250_cone_cylinder.json`.
Both source and copy have SHA-256
`1f7503299d821e7385ebb99c93aef73b15c81877af3288e19c21a8b955e0c78b`.
Source: Ehret, Rossow and Stevens, *An Analysis of the Applicability of the
Hypersonic Similarity Law to the Study of Flow About Bodies of Revolution at
Zero Angle of Attack*, NACA TN 2250, December 1950, Figure 1 (page 13), K=1.
https://ntrs.nasa.gov/citations/19930082890

Use gamma=1.4, cone length l=1 m, shoulder radius 1/(2*l/d) m, and x_end=3.4l.
Thus the cone angles are atan(1/6) and atan(1/12), not the rounded angles of
M3. The abscissa is 100*x/l **from the tip**. The comparison quantity is
Pi=(p-p_inf)/p_inf, not the conventional dynamic-pressure-normalized Cp.
The source calls the free-stream static pressure p0; that p0 is **not** the
stagnation pressure used elsewhere in this memo. Interpolation is restricted
to the downstream cylinder, and x/l=1 uses the post-expansion limit.

| Mach | 100 x/l | Supplied Pi | Rotational MOC Pi | Difference |
| --- | ---: | ---: | ---: | ---: |
| 3 | 100.0 | -0.273 | -0.274631 | -0.001631 |
| 3 | 119.7 | -0.207 | -0.206654 | +0.000346 |
| 3 | 151.6 | -0.130 | -0.135266 | -0.005266 |
| 3 | 199.0 | -0.077 | -0.077695 | -0.000695 |
| 3 | 260.7 | -0.040 | -0.042138 | -0.002138 |
| 3 | 339.0 | -0.012 | -0.022335 | -0.010335 |
| 6 | 100.0 | -0.246 | -0.254891 | -0.008891 |
| 6 | 112.9 | -0.199 | -0.211789 | -0.012789 |
| 6 | 136.1 | -0.141 | -0.154989 | -0.013989 |
| 6 | 177.5 | -0.088 | -0.094087 | -0.006087 |
| 6 | 236.4 | -0.044 | -0.051266 | -0.007266 |

Maximum absolute differences are **0.01033455 and 0.01398948**, both passing
the unchanged prescribed limit 0.02. The source uses irrotational graphical
MOC (mostly Sauer's method), with a curved bow shock but the rotation term
omitted. Its 3% criterion concerns scatter of surface Mach angle relative to
its total variation, not a quantitative 3% pressure bound. Reading uncertainty
is +/-0.005 in Pi and +/-1 percentage point in 100*x/l. The generous 0.02
comparison allowance accounts for digitization and the original graphical
method; agreement does not establish experimental accuracy or rotational
superiority. No reference values or uncertainty metadata were edited.

### M3 relative differences and default-resolution estimates

For the original nested count/fan-angle refinements, maximum relative
pressure differences use the finer solution's **local surface pressure** as
the denominator (not p_inf):

| Mach | 24 to 48 | 48 to 96 | 96 to 192 |
| --- | ---: | ---: | ---: |
| 3 | 7.06014e-5 | 1.38893e-5 | 3.52353e-6 |
| 6 | 3.95358e-5 | 1.01191e-5 | 2.57860e-6 |

The true default (count 24, fan step **1 degree**) is compared separately
with Richardson-extrapolated p_ref=p_192+(p_192-p_96)/3. The reference pair
uses the original exactly nested angular subdivisions. Near-shoulder means
1<=x/l<=1.1 (1001 fixed samples); downstream means 1.1<=x/l<=x_end/l (1501
samples). Cylinder-only cubic interpolation excludes the upstream jump.

| Mach | Region | Max default error estimate [Pa] | Max relative error estimate | Max fine-grid Richardson correction [Pa] |
| --- | --- | ---: | ---: | ---: |
| 3 | Near shoulder | 0.124460 | 1.64564e-5 | 0.001592 |
| 3 | Downstream | 0.848818 | 8.62003e-5 | 0.011591 |
| 6 | Near shoulder | 0.102977 | 1.34805e-5 | 0.000893 |
| 6 | Downstream | 0.245676 | 2.85162e-5 | 0.008242 |

These estimate inviscid discretization error, not the pressure-jump smoothing
error subsequently introduced by the boundary-layer interpolant. A 1.2e-4
relative regression bound provides about 39% margin over the largest measured
default estimate. The original Pa-based bounds remain unchanged.

### Historical shoulder-width sensitivity before API approval

Use the M7 M3 case: cone length 0.3 m, radius 0.05 m, cylinder length 1 m,
p_inf=10000 Pa, T_inf=220 K, wall T=300 K. Only insert one pre-shoulder
constant-cone pressure point for each narrowed case. Width is measured along
the surface, so the corresponding radius is s*sin(theta_cone). Default width
is 12.67242194 mm; alternatives are 1 mm and 0.2 mm. Edge inputs have 46/47/47
stations. The production `boundary_layer_edge` is **unchanged**.

The reference delta99 is computed on an independent constant-edge cone at
the shoulder: 0.5646614 mm for laminar and 2.5929231 mm for SA. Comparisons
are at shoulder + (5,20,100)*that upstream delta99, using the same stations
for all widths. Absolute surface distances are (0.3069614,0.3154314,0.3606043) m
for laminar and (0.3171027,0.3559966,0.5634304) m for SA. A common 93-station
output grid includes samples inside every pressure-ramp width to evaluate
warnings. The V20 fixed-output harness leaves each pressure interpolant
unchanged; pressure is not resampled onto the common output grid.

Entries below are signed percentages relative to the **0.2 mm** result:

| Regime | Width | Quantity | 5 delta99 | 20 delta99 | 100 delta99 |
| --- | --- | --- | ---: | ---: | ---: |
| Laminar | Default | Cf | -28.6400 | -10.8849 | -1.6467 |
| Laminar | Default | Heat flux | -12.3295 | -4.2419 | -0.0497 |
| Laminar | Default | delta99 | +3.1852 | +1.0886 | -0.1886 |
| Laminar | 1 mm | Cf | -4.2398 | -1.0123 | -0.1225 |
| Laminar | 1 mm | Heat flux | -2.5967 | -0.6360 | -0.0462 |
| Laminar | 1 mm | delta99 | +0.3363 | +0.0893 | -0.0134 |
| SA | Default | Cf | -1.9727 | -0.2782 | +0.1269 |
| SA | Default | Heat flux | +0.2181 | +0.3102 | +0.2683 |
| SA | Default | delta99 | -1.8178 | -1.5484 | -0.3718 |
| SA | 1 mm | Cf | -0.1243 | +0.0132 | +0.0111 |
| SA | 1 mm | Heat flux | +0.0508 | +0.0399 | +0.0200 |
| SA | 1 mm | delta99 | -0.1306 | -0.1145 | -0.0406 |

The narrowed reference is also checked with the common output intervals
halved and the relative marching-step cap halved (0.03 to 0.015). Maximum
relative changes at these comparison stations are 1.67e-5 / 7.39e-6 / 8.25e-7
for laminar Cf/heat flux/delta99, and 2.08e-4 / 7.47e-5 / 3.91e-4 for SA.
These are much smaller than the important width effects. The regression
limit for this step check is 0.002. This is not a wall-normal grid study or
a physical validation of the sharp-corner approximation.

| Regime | Width [mm] | Accepted marching steps | Max delta99 abs(d ln p/ds) | Pressure-length warning |
| --- | ---: | ---: | ---: | --- |
| Laminar | 12.6724 | 1079 | 0.05403 | No |
| Laminar | 1 | 1290 | 0.66258 | No |
| Laminar | 0.2 | 1483 | 3.27640 | Yes |
| SA | 12.6724 | 976 | 0.27047 | No |
| SA | 1 | 1047 | 3.48468 | Yes |
| SA | 0.2 | 1201 | 17.41804 | Yes |

The table's warning evaluation uses the common 93 output stations, including
points inside each ramp. A separate call to the unchanged public marcher with
only the one extra pre-shoulder edge station produces 46 downstream output
stations: at width 0.2 mm it emits no warning for laminar, and only the
transverse-curvature warning for SA. The pressure ramp endpoints have zero
PCHIP slope, so that sparse public warning sampling **misses** the large
interior pressure-length maximum. No warning suppression is applied. Any
proposed narrowing API must also arrange interior ramp diagnostics or output
sampling; adding just one edge point does not guarantee this warning fires.
The boundary-layer implementation was not changed to fix this separate
sampling limitation.

Accepted steps are counted by an observer of the existing step-resolution
predicate, without modifying discretization or acceptance decisions. All six
runs finish at s=1.30413812651 m. All SA cases additionally warn about
transverse curvature, with maximum delta99/r about 0.316 at the downstream
end. The dense-grid narrow-reference runs accept 2176 / 1878 steps. Serial
width-study times are approximately 1.2--2.3 s (laminar) and 3.2--4.6 s (SA).

**Decision: the effect is large**, particularly laminar Cf and heat flux near
the shoulder. Per the review instruction, stop at this proposal:

- Consider an explicit opt-in `pre_corner_width` (SI surface distance) on
  `boundary_layer_edge`, retaining the current default unless separately
  approved. For the cone--cylinder, insert the analytic cone state at
  s_shoulder-width while keeping the exact post-expansion shoulder point.
- An inserted point changes the edge length. Callers must interpolate the
  piecewise linear radius onto `edge.distance`, e.g. with
  `np.interp(edge.distance, result.surface_distance, result.radius)`;
  blindly passing `result.radius` would have the wrong shape. An API returning
  an aligned edge/radius pair is another design choice, not implemented here.
- Smaller width costs more marching steps and can explicitly violate the
  pressure-length applicability criterion. Include interior-ramp diagnostic
  sampling so that the existing endpoint-only warning blind spot is not
  mistaken for applicability. Even 1 mm versus 0.2 mm differs
  by 4.24% in laminar Cf at 5 delta99. There is no evidence here that 0.2 mm is
  a converged physical zero-width limit; a finite-radius shoulder or viscous
  interaction treatment may be required for local accuracy.

No production edge-conversion change was made, and no automatic choice of
width is endorsed by these completion tests.

### Total-pressure span at the terminal cross-section

Do not confuse `_solve`'s last complete C+ characteristic with a terminal
x=constant cut: its wall can still be at the shoulder and its shock well
beyond the body end. A test observer instead samples every active C+ row's
intersection with the final axial plane, then closes that section with the
body state and the local fitted shock state. Log total pressure is linearly
interpolated on row edges; the terminal shock angle is interpolated between
fitted shock stations and passed through the local shock relations. These
are sampled diagnostics, not a new public field-output API.

For the original ten-diameter M3 cases (cone length 1 m):

| Mach | Count / fan step | Terminal x [m] | Min p0 [Pa] | Max p0 [Pa] | Span [Pa] | Span/body p0 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 3 | 24 / 1 degree | 4.3325002 | 366987.4574 | 367013.7604 | 26.3030 | 7.1673e-5 |
| 3 | 192 / 0.125 degree | 4.3325002 | 366987.4574 | 367024.1128 | 36.6554 | 9.9882e-5 |
| 6 | 24 / 1 degree | 2.6653867 | 15774503.3712 | 15774503.3712 | 0 | 0 |
| 6 | 192 / 0.125 degree | 2.6653867 | 15774503.3712 | 15774503.3712 | 0 | 0 |

The coarse M3 shock trace has few points in this newly disturbed region, so
its small span is visibly resolution-sensitive. Report the refined value as
about 0.010% of body total pressure, not a highly accurate stagnation-pressure
loss estimate. M6's leading shoulder characteristic first reaches the shock
at x=4.48635 m, downstream of the terminal section, so no curved-shock entropy
variation has yet entered that section. In both NACA comparison geometries
x_end/l=3.4 is also upstream of the leading-wave/shock intersection, so the
terminal total-pressure span there is zero. Rotational effects are small or
absent in these short-body terminal sections; this does not justify dropping
entropy transport for longer bodies or more strongly curved shocks. The
existing longer M6 rotational test still measures a nonzero span.


### Review validation and independently measured main baseline

The baseline is the **pre-change main commit 459b6bcbaf9a8c9ff68050d17ca10961c9ccb4e4**,
which is this branch's HEAD/base and the available `origin/main` ref. The
local `main` ref still points to older 86fbaa2 and is not used as the feature
baseline. `git archive HEAD` produced a separate temporary source tree at
`/private/tmp/moc-review/main`; no checkout, branch, index or commit was
modified. The same existing .venv was used, with PYTHONPATH selecting that
archive's src and a separate COVERAGE_FILE. No dependency synchronization,
network access or permission elevation was attempted.

Both source trees ran exactly the ordinary-test selection
`-m "not generated_assets and not package_artifact and not solver_verification"`
with `-n auto --cov-fail-under=0`:

| Measurement | Main before change | Reviewed implementation |
| --- | ---: | ---: |
| Ordinary tests passed | 1266 | 1326 |
| Environment failures | 1 | 1 |
| Numerical-core covered statements | 5499/5622 | 5828/5951 |
| Numerical-core covered branch exits | 1220/1324 | 1312/1416 |
| Branch-only diagnostic | 92.1450% | 92.6554% |
| Repository `coverage report` metric | 96.7319% | 96.9187% |

The repository's **combined statement/branch coverage.report >=95% criterion
is satisfied**, as clarified in review. Branch-only coverage has increased
by 0.5104 percentage points. The new module is **100% statements, 100% branch
exits, 100% combined** (328 statements, 92 branch exits) in the ordinary suite.
This baseline is an actual executed run, not an exclusion estimate.

Both ordinary suites have the same single failure:
`tests/gui/test_documentation.py::test_local_documentation_server`, at
`socket.bind`, `PermissionError: [Errno 1] Operation not permitted`. The
baseline run took 32.32 s and the reviewed run 32.08 s, with 15 xdist workers.
No test was skipped, xfailed, excluded specially, or given a weaker threshold.

Final review checks:

- Ruff check and format --check: pass.
- Strict mypy: pass (99 source files).
- Sphinx HTML with -W: pass.
- Sphinx doctest with -W: 227 examples pass.
- `git diff --check`: pass.
- `.venv/bin/pytest tests/test_characteristics.py -m "not generated_assets and not package_artifact" --no-cov -n auto`: **70 passed in 22.92 s**, including ten marked slow tests.
- Supplied reference-data SHA-256 remains unchanged.

The previous full-gate boundary is retained: `scripts/check.sh --full`,
generated assets, package/clean-install checks and unrelated slow solver
verification remain for the reviewer. No commit, push, version change,
boundary-layer implementation change or edge-conversion API change occurred.


## Approved corner-input API and revised default resolution

This section supersedes the historical count-24/default-edge statements and
unimplemented proposal above. The user approved pre-corner insertion. No
boundary-layer module, shock routine, reference JSON, generated asset or CI
file was changed in this revision. No commit, push or privilege escalation
was performed.

### Implemented API and checks

`AxisymmetricBodyFlowResult.boundary_layer_inputs(*, corner_length=None,
viscosity_model=AIR_VISCOSITY)` returns the named tuple
`BoundaryLayerInputs(edge: BoundaryLayerEdge, wall_radius: FloatArray)`.
`boundary_layer_edge` has been removed. Original MOC surface arrays are
unchanged. The private immutable `_corner_states` records each convex
vertex's original surface index and its upstream limiting pressure. For
later corners this is the actual MOC wall state immediately before the fan,
not an extrapolated original-cone pressure. Total pressure and temperature
are constant on the body streamline and reconstruct all other edge states.

The conversion adds two stations per convex corner: the upstream pressure
at s_c-w and the geometric-mean pressure at s_c-w/2. At s_c it keeps the
original downstream state. Width w defaults to 0.001*s_c, or the positive
finite scalar `corner_length` in metres. Widths reaching the previous
output point and numerically indistinguishable inserted coordinates raise
`ValueError`. No points are silently removed. The returned radius array is
owned/read-only and aligns with every edge point by exact piecewise linear
meridian interpolation. `BoundaryLayerInputs` is exported at package level;
the guide and API/model documentation now use its paired fields.

Fast tests cover invalid widths (including no-corner inputs), overlap and
roundoff, both corners of a cone/cylinder/boattail, local pre-fan limits
compared with separately truncated bodies, logarithmic midpoint pressure,
unchanged original pressures, exact radius interpolation, ownership and
read-only arrays. M5 checks thermodynamic agreement at the original stations
in the augmented edge. M7 calls the public marcher with the paired radii and
requires a pressure-length warning for both laminar and SA default inputs.

### Selection of the new default count

The one-degree fan step is retained. For the M2/10-degree forced cone starting
at x=0.05 and ending at x=1, maximum relative surface-pressure errors near
the cutoff are:

| Count | Max relative pressure error | Below 0.001? |
|---:|---:|---|
| 80 | 0.00119736255 | No |
| 84 | 0.00108956092 | No |
| 87 | 0.00101767485 | No |
| 88 | 0.000995459513 | Yes |
| 89 | 0.000974078328 | Yes |

An exhaustive integer scan from 4 through 88 confirms this boundary: counts
4--6 stop with an invalid-shock-angle diagnostic, counts 7--87 complete but
exceed the error criterion, and only 88 passes. No unresolved run is treated
as an accurate result.

The minimum count reaching the requested criterion is **88**, now the public
and private solver default. A slow regression compares 87 with the actual
default (omitting the keyword), and verifies the public entry point uses the
same count. The 0.001 acceptance criterion is unchanged, with no rounding of
the computed value. This criterion leaves only 0.45% relative headroom to its
limit at N=88; it is a particular verification case, not a universal error
bound for every supported body.

At count 88 the maximum relative errors (pressure, Mach, shock-ray position)
are (9.9545951e-4, 3.5343445e-4, 1.3834825e-4) for M2/10 degrees and
(5.3818987e-5, 8.9204919e-6, 5.2867257e-5) for M6/5 degrees. The unchanged
24/48/96/192 study retains final orders (1.9700,1.9696,2.0200) and
(1.9990,1.9990,1.9995). Its old count-24 error limits describe that explicit
coarse grid, not the new default. M1/M5 and M7's conical edge remain
by-construction initialization/API checks, not MOC accuracy evidence.

For M3/9.46 degrees with a unit-length cone and a cylinder ten diameters
long, three serial N=88/one-degree calls took **0.10735, 0.08389, 0.08324 s**
(including initial-line integration). Thus the new default is well below
one second on this development machine and is not deferred. These are outer
flow timings, not boundary-layer timings.

With the same Richardson comparison as the first review, revised default
surface-pressure estimates are:

| Mach | Region in x/cone length | Max absolute estimate [Pa] | Max relative estimate |
|---:|---|---:|---:|
| 3 | Shoulder 1--1.1 | 0.124516 | 1.69807e-5 |
| 3 | Downstream 1.1--end | 0.068789 | 7.79875e-6 |
| 6 | Shoulder 1--1.1 | 0.211945 | 2.81856e-5 |
| 6 | Downstream 1.1--end | 0.571814 | 6.31338e-5 |

The fan remains at one degree; refining only the initial characteristic
count does not monotonically reduce every error (notably M6). The retained
coupled count/fan refinement test and its original limits still pass. At
the new default the NACA maximum absolute pressure-coefficient differences
are 0.0103484 and 0.0140098, both below the unchanged 0.02 limit. M4 relative
mass differences are +5.84908e-5 and -2.76578e-6, both below 1e-4.

### New shoulder-width measurement

Same physical conditions as the earlier study: M3, p_inf=10000 Pa,
T_inf=220 K, body (x,r)=(0,0),(0.3,0.05),(1.3,0.05) m, wall temperature
300 K, laminar and SA-noft2. The MOC count is now 88. The shoulder surface
distance is 0.3041381265 m, so the default width is **0.3041381265 mm**.
All three widths use the production API, including its pressure midpoint;
this differs from the old one-point experiment. All inputs have 166 stations.

A separately marched constant-cone edge gives upstream delta99=0.5646614 mm
(laminar) or 2.5917307 mm (SA). These fixed upstream thicknesses define the
sampling distances, rather than using the width-dependent local thickness.
Comparison uses the existing fixed-output harness without resampling edge
pressure: 230 common output stations, including 65 uniform points across
the widest (1 mm) ramp and the three downstream targets. At 5/20/100 delta99
the surface distances are (0.30696143,0.31543136,0.36060427) m for laminar
and (0.31709678,0.35597274,0.56331120) m for SA.

Entries below are **100*(default/reference - 1)**, in percent. The reference
width is explicit in each row.

| Regime | Reference width | Quantity | 5 delta99 [%] | 20 delta99 [%] | 100 delta99 [%] |
|---|---|---|---:|---:|---:|
| laminar | 0.2 mm | Cf | -0.61709 | -0.13466 | -0.01677 |
| laminar | 0.2 mm | heat_flux | -0.44232 | -0.10150 | -0.00915 |
| laminar | 0.2 mm | delta99 | +0.04616 | +0.01162 | -0.00177 |
| laminar | 1 mm | Cf | +3.82493 | +0.87692 | +0.10543 |
| laminar | 1 mm | heat_flux | +2.27129 | +0.55300 | +0.04255 |
| laminar | 1 mm | delta99 | -0.29455 | -0.07695 | +0.01146 |
| turbulent | 0.2 mm | Cf | -0.01555 | +0.00044 | +0.00159 |
| turbulent | 0.2 mm | heat_flux | +0.00583 | +0.00513 | +0.00274 |
| turbulent | 0.2 mm | delta99 | -0.01807 | -0.01549 | -0.00591 |
| turbulent | 1 mm | Cf | +0.09978 | +0.00263 | -0.00978 |
| turbulent | 1 mm | heat_flux | -0.03052 | -0.02879 | -0.01680 |
| turbulent | 1 mm | delta99 | +0.11426 | +0.09701 | +0.03214 |

The laminar width effect remains material versus 1 mm, including at 20
upstream thicknesses. The default is much closer to 0.2 mm, but neither
agreement nor successful completion makes the local thin-layer corner
solution physically reliable. At least the first several tens of upstream
thicknesses after the shoulder should not be treated as validated corner
predictions. No zero-width physical limit or wall-normal convergence is
claimed.

| Regime | Width [mm] | Accepted steps | Maximum pressure-length parameter on common grid |
|---|---:|---:|---:|
| Laminar | 0.304138 | 1536 | 2.03173 |
| Laminar | 1 | 1371 | 0.62432 |
| Laminar | 0.2 | 1589 | 3.08137 |
| SA | 0.304138 | 1281 | 10.89227 |
| SA | 1 | 1129 | 3.31313 |
| SA | 0.2 | 1326 | 16.55165 |

A direct public default-width run (its own 165 positive output stations)
completes for both regimes and emits pressure-length warnings with maxima
**1.42389** (laminar) and **7.47523** (SA). These differ from the densely
sampled common-grid maxima; the inserted midpoint guarantees a sampled steep
gradient here, not an exact maximum. Curvature warnings also occur: the finer
first cone station gives delta99/r=0.10594 for laminar, and SA reaches
0.316053 downstream. The M7 laminar cone/Mangler check remains within
4.89e-15 relative error; as stated above, its constant edge is by construction.

For the 0.2 mm reference, halve both the common output intervals and the
relative streamwise step from 0.03 to 0.015. Maximum changes over the three
targets in Cf/heat flux/delta99 are (2.55e-5,5.13e-5,6.19e-6) relative for
laminar and (9.25e-5,6.10e-6,8.57e-5) for SA, below the unchanged 0.002
streamwise-check limit. The refined runs use 2520 and 2267 accepted steps.
Some very small SA width differences are comparable to this discretization
uncertainty and their displayed last digits should not be interpreted as
resolved physical effects.

An initial attempt to union three separately spaced 17-point ramp output
grids failed in SA at s=0.304076 m. The final comparison uses one uniformly
spaced ramp grid common to all widths and the independent refinement check
above. No solver equation, tolerance, failure condition or acceptance limit
was changed, and no exception is suppressed. This observation is a limit
of solver robustness for prescribed output distributions, not evidence of
physical separation; its detailed nonlinear cause was not established.

### Validation after the approved changes

- Ruff check and format check: pass; mypy: pass (99 source files).
- Requested characteristics suite including solver verification: 78 passed
  in 28.39 s with `-n auto`; the fast subset has 67 tests (2.16 s serially).
- Requested ordinary suite: 1333 passed, one environment failure, 16 warnings,
  32.67 s. The sole failure is the same pre-existing
  `tests/gui/test_documentation.py::test_local_documentation_server` socket
  bind `PermissionError: [Errno 1] Operation not permitted` seen on baseline.
  No escalation, skip, xfail or reduced coverage limit was used.
- Numerical-core coverage (`src/aerophysics/*.py`): 5974 statements with 123
  missed, 1422 branches with 105 missed; combined **96.92%** passes 95%.
  Characteristics: 351/351 statements and 97/98 branch exits, **99.78%**
  combined. Branch-only core coverage is 1317/1422 = **92.6160%**, above the
  previously measured pre-change main 1220/1324 = 92.1450%.

- `.venv/bin/python -m sphinx -W -b html docs docs/_build/html`: pass.
- Additional Sphinx doctest build with `-W`: 228 passed.
- Final Ruff check: pass; format check: 117 files already formatted; mypy:
  99 source files pass. `git diff --check`: pass.

`scripts/check.sh --full`, generated-asset and package verification remain
assigned to the reviewer; no `uv` or elevated command was run.


## PR #89: roundoff-close output stations

The reviewed 30 m cone (`x=linspace(0,30,11)`, `r=0.4*x`, M3, default
count 88) retained automatic x=14.999999999999998 beside input x=15.
Before correction the minimum dx/ds was 1.77636e-15/3.55271e-15 m.
A downstream collinear vertex one ULP below a MOC wall point likewise
created duplicate-scale spacing; one ULP above it failed the inverse-wall
root bracket. Both paths now give exact input geometry priority.

For an input vertex v use the local overlap tolerance
`max(32*eps*min(h_left,h_right), 4*spacing(v))`, with the available adjacent
interval at an endpoint. The 4-ULP floor accounts for coordinate subtraction
when a neighboring interval is much shorter than the coordinate magnitude;
it is not a fixed absolute rounding scale. Only automatic output stations
are coalesced with input vertices. Two supplied vertices are never merged.
Cone sampling retains the exact input coordinates and radii. Downstream,
roundoff-close inverse-wall queries use the existing endpoint state and
output the exact input vertex, avoiding a nearly singular root bracket.
A full wall step reaching a segment end to roundoff is placed exactly there;
all intervening input vertices are retained, and the loop terminates at the
exact endpoint rather than skipping a small remainder by total body length.
Collinear output insertion otherwise leaves the characteristic net unchanged.
No global coordinate rounding, relaxed edge-distance validation, or changes
to boundary-layer modules, existing test limits, or reference data are used.

For the reported case the corrected minimum dx/ds is
**0.0681818182/0.0734340656 m**, with 97 surface stations and each of the 11
input vertices present exactly once. The paired edge/radius conversion
succeeds. At wall temperature 300 K, both laminar and SA reach surface
s=32.310988842807 m (initial isolated timings 0.515/6.721 s). Regression
checks also cover scale factors 1e-6/1/1e6, downstream +/-1 ULP vertices,
segment endpoints, explicitly adjacent input vertices, meaningful 1e-10 m
intervals and two distinct nearby corners. Collinear insertion preserves
cone states and the downstream pressure/Mach/shock arrays exactly in the
comparisons. Intentionally tiny *input* intervals are preserved, not promised
to be suitable boundary-layer marching intervals.

Commits `1306b42` and `69809ed` and their verification-gate changes are
retained without modification. In an isolated local clone at HEAD, using
`VERIFICATION_BASE=HEAD`, the pristine solver gate returns `no`. After copying
only the changed `tests/test_characteristics.py`, `git status` lists only
that file and `scripts/needs-verification.sh solver` returns `yes` with
`solver-relevant change: tests/test_characteristics.py`.


Final verification of this fix: Ruff check/format (117 files), mypy
(99 source files), `git diff --check`, and Sphinx HTML with `-W` all pass.
The requested characteristics command passes **86 tests in 31.35 s**,
including the unmodified cone, multiple-corner, refinement, mass, coupling
and NACA checks. The requested ordinary command gives **1339 passed,
1 failed, 16 warnings in 34.36 s**; the sole failure remains the sandbox's
`socket.bind` prohibition in `test_local_documentation_server`
(`PermissionError: [Errno 1] Operation not permitted`). Numerical-core
combined coverage is **96.95%** (6000 statements, 1434 branch exits);
`characteristics.py` is **100%** (377 statements, 110 branch exits).
The 95% core gate passes. No test limits were changed and no skip/xfail was
added. `scripts/check.sh --full` remains for the reviewer. No commit, push,
privilege escalation or generated-file edit was performed.


## Coarse-fan shock extrapolation failure (3-degree maximum fan step)

Reproduced M3, (x,r)=(0,0),(0.3,0.05),(1.3,0.05) m, p_inf=10000 Pa,
T_inf=220 K, count 88, fan_angle_step=3 degrees. The failing last interior
point was (1.2960647601,0.4729325710) m, theta=0.030357129 rad; the previous
shock point was (1.2470348168,0.4888831772) m, beta=21.4069926385 degrees.
The unrestricted secant trial sequence was 21.4069926385, 21.4072067084,
19.8709149865, then **19.3692770429 degrees**, below the free-stream Mach
angle 19.4712206345 degrees. The last trial, not a valid physical state,
caused the former `flow.shock` range error.

The review's geometric explanation is supported, but its approximate residual
root is not a root of the full discrete compatibility equation. With the
complete geometric corrector, at beta_previous the intersection is x=2.8760192 m
(outside x_limit=1.3 m), h*(x-a.x)=-0.03147338 and the full residual=+0.06788149.
At 20.7 degrees the intersection is x=3.5089212 m and the full residual is
+0.03797845. Sampling 2001 angles between mu_inf+1e-9 and beta_previous
found a minimum residual +0.01666599 near 19.89612 degrees, and every
intersection beyond x=2.8760192 m. Samples through 30 degrees also remain
positive. Thus merely bracketing between mu_inf and beta_previous does not
provide a sign change; the large, unnecessary downstream cell is not a
well-resolved shock continuation for this coarse fan. This is a numerical
extrapolation failure, not a detached shock or subsonic body flow.

Initial approach was **(b), applied as a fallback to preserve existing output**
(**the continuation-only authorization below is superseded by the PR #92
dependence and front-size checks**):
retain the existing arithmetic and successful secant solves, but guard every
trial angle before Rankine--Hugoniot evaluation. A private trial-angle
exception is distinct from physical `ModelRangeError`. If that numerical
trial leaves (mu_inf,pi/2), and the geometric continuation at the previous
shock tangent is beyond x_limit, return no fitted point. `_row` then stops
shock tracking with `has_shock=False`, retaining the computed characteristic
front, as in the existing downstream truncation path. At least three front
points are required. No artificial outer state is inserted. This initial condition restricted only
the prospective shock point, not the entire discarded outer region; it was
insufficient to authorize truncation on longer bodies. The corrected
wall-dependence criterion and mandatory front guards are specified below.
Otherwise report a location-bearing `ModelRangeError`. Subsonic states, characteristic crossing,
geometry-corrector failures and other `ModelRangeError`s are not caught or
silenced. The same failed unit process with x_limit=3 m is a regression check
that the failure remains visible when its continuation is in the requested
region. The angle guard is also checked by an observer that requires every
actual `flow.shock` call to receive an admissible angle.

Before editing, 59 representative solves were recorded: 55 successful and
four reproductions of the error. Recomputing all of them after the fix gives
**bit-identical arrays for all 55 previously successful cases**: maximum
relative wall-pressure difference **0**, unchanged surface coordinates and
unchanged shock coordinates/angles. They include the M3/9.46 and M6/4.76
default cone--cylinders, four coupled refinement levels, both forced-cone
marches at four counts, planar expansion, multiple corners, the rotational
entropy/mass case, NACA geometries and the boundary-layer coupling geometry.
For the unit-cone default cases the last shock x stays at 6.0480865521 m (M3)
and 4.4863467839 m (M6), respectively. The four previously failed solves now
complete. There is no new reference-data fixture or change to old thresholds.

For the reported 0.3 m cone, the successful 1-degree trace still ends at
x=1.8134098139 m. The newly successful 3-degree result stops its trace at the
last known fitted point, **x=1.2470348168 m**, before the body endpoint. It
does not supply the remaining shock segment, even though it supplies the
wall flow through x=1.3 m. Public signatures are unchanged. The model page
and solver docstring now explicitly allow either a successful trace beyond
the endpoint or an incomplete trace ending before it. Do not extrapolate
missing shock data to claim a terminal transverse total-pressure profile.
The previously reported default-resolution profiles remain unchanged.

The fast regression compares cylinder pressure on 1001 common axial stations,
using separate cylinder-only cubic interpolation, never across the shoulder
jump. Maximum 3-degree versus 1-degree difference is **0.6487273 Pa**, or
**8.84545098e-5 relative**. Its new acceptance bound is **1.2e-4**, about 35.7%
above the measurement. This measures finite-fan resolution difference, not
physical validation or a globally certified error bound. The slow-marked
sweep uses counts **8,24,88** and maximum steps **0.5,1,2,3,4,5 degrees** for
M3/9.46 and M6/4.76 unit cones followed by ten-diameter cylinders. All **36**
cases complete with positive pressure and supersonic surface states; isolated
run times were **0.712 s** (18 M3 cases) and **0.539 s** (18 M6 cases).


Final validation for this correction:

- Requested `uv run --no-sync pytest -q -p no:cacheprovider
  tests/test_characteristics.py`: **89 passed, 76.49 s**.
- Requested solver-verification-only command: **15 passed, 74 deselected,
  71.45 s**. Most time is the existing laminar/SA coupling verification.
  The new sweep measured with normal coverage instrumentation takes
  **1.810+1.202=3.012 s** for its 36 solves (3.99 s including pytest/coverage
  reporting). The reproduction also completes for every separately checked
  count 8,12,16,20,24,32,48,88.
- Requested Ruff check and format check: pass (117 files); mypy: pass
  (99 source files). Extra Sphinx HTML build with `-W` and `git diff --check`:
  pass.
- In the complete 89-test run, `characteristics.py` covers **391/391 statements
  and 118/118 branch exits: 100% combined**. No whole-repository coverage
  claim is made from these targeted commands; the full ordinary suite was
  not repeated in this correction.

The first literal uv invocation could not access
`/Users/satoshi/.cache/uv/sdists-v9/.git` (sandbox `Operation not permitted`).
All requested commands then succeeded with
`UV_CACHE_DIR=/private/tmp/moc-fan-fix/uv-cache`; this requires no escalation,
network synchronization, or environment change in the repository.
`scripts/check.sh --full` remains assigned to the reviewer. No commit, push,
privilege escalation, boundary-layer change, old threshold relaxation,
skip/xfail, reference-data or generated-file edit was made.


## PR #92: wall dependence and characteristic-front exhaustion

Work started on `fix/moc-coarse-fan-shock-truncation`, HEAD `5e30f70`.
The previous continuation-only condition was insufficient: discarding the
outer region can deprive later wall points of incoming C- data even though
the prospective shock intersection is farther downstream than the endpoint.
The reported case was reproduced: at count 88, 3-degree fan step, the retained
outer point has x=1.2960647601 m, r=0.4729325710 m and local C- slope
**-0.3311845151**. Its straight-ray cylinder foot is **x=2.5730947895 m**,
whereas the prospective shock intersection is x=2.8760191852 m. Terminations
at 2.5 and 2.8 m formerly exhausted the front and raised IndexError.

Truncation now needs BOTH conditions: the original shock-continuation test
must permit a candidate truncation, and the retained outermost point's local
C- ray must have an estimated wall foot at or beyond x_limit. Retaining the
original test is conservative: a failed shock fit inside the requested axial
extent is still reported, including the existing x_limit=3 m unit regression.
For outer point a, current wall point w, local C- slope lambda_minus and
actual body-segment slope b, the estimate is

```
gap = a.r - w.r - b * (a.x - w.x)
x_foot = a.x + gap / (b - lambda_minus)
```

Require positive gap and closing slope b-lambda_minus; if no forward foot
can be established, reject the truncation. The body-segment angle is supplied
explicitly during a centered fan: the intermediate flow angle in the fan is
not the wall tangent. On a convex piecewise linear body, future wall segments
lie below the current segment's extension, which is conservative with respect
to wall geometry. **The straight C- ray itself is only a local approximation,
not a guaranteed conservative bound for the curved characteristic.** In this
example it overestimates the available wall march: the count-88 finite front
runs out near x=2.49225 m, before its estimated foot at 2.57309 m.

Consequently front-size checks are mandatory even after the estimate passes:
`_solve` requires two points before indexing the incoming wall characteristic;
`_wall_at` requires two before building an inverse-wall stencil; and `_row`
requires a nonempty parent front before interior/shock operations. The row
entry check also covers `points[1:]` fan updates and `points[2:]` wall-step
consumption. A shortage raises a wall-position-bearing `ModelRangeError`
identifying front exhaustion and recommending **reduce fan_angle_step**.
There is no broad IndexError catch, invented outer state, missing-interval
extrapolation, or successful return before the requested wall endpoint.
Subsonic/geometry errors and the trial-angle guard retain their behavior.
These checks enforce the available discrete front, not a claim of a globally
exact continuum domain-of-dependence calculation.

Measured endpoint behavior at 3 degrees (M3, shoulder (0.3,0.05) m):

| End x [m] | Count 8 | Count 24 | Count 88 |
|---:|---|---|---|
| 1.3, 2.0, 2.4 | Complete | Complete | Complete |
| 2.5 | Front exhausted at wall x=2.49170635 | Front exhausted at wall x=2.49218854 | Front exhausted at wall x=2.49224752 |
| 2.8 | Reject at shoulder; foot=2.56863798 | Reject at shoulder; foot=2.57260786 | Reject at shoulder; foot=2.57309479 |
| 3.0 | Existing shock-angle diagnostic | Existing shock-angle diagnostic | Existing shock-angle diagnostic |

All failures are location-bearing `ModelRangeError`s; there are no IndexErrors.
The 2.5 m case demonstrates why the local foot test alone is not sufficient.
Long bodies with coarse fans may therefore still fail explicitly and need a
smaller `fan_angle_step`; they are not supplied with a fabricated continuation.
The 1.3 m case retains its pressure difference **8.84545098e-5** against
1 degree and the **unchanged 1.2e-4 bound**. Successful output is unchanged:
all arrays are bit-identical for **59 comparisons against HEAD 5e30f70**,
and also for the **55 saved successful cases from before the first coarse-fan
fix**. Maximum relative wall-pressure difference is **0**. The public API
is unchanged; the model page and docstrings state the revised guarantees.

New fast regressions cover ends 2.5/2.8/3.0 with counts 8 and 88, the estimated
foot one floating-point step before/equal to/after the endpoint, unavailable
forward intersections, and each of the row/inverse-wall/direct-wall point-count
guards. The foot test checks use of the real wall tangent during a fan. The
existing 36-case solver-verification sweep and its settings are unchanged;
a focused no-coverage run measured 0.785 s for M3 plus 0.596 s for M6.
No old thresholds, skip/xfail settings, boundary-layer modules, reference data
or generated files were changed.


Final verification on the corrected branch:

- Requested full `tests/test_characteristics.py` command: **103 passed,
  82.58 s**, including all unchanged numerical and boundary-layer checks.
- Requested `-m solver_verification` command: **15 passed, 88 deselected,
  86.97 s**. The 36-condition sweep remains unchanged and successful.
- Requested Ruff check/format: pass (121 files); mypy: pass (103 source files).
  Extra Sphinx HTML build with `-W` and `git diff --check`: pass.
- In the 103-test run, `characteristics.py` has **409/409 statements and
  128/128 branch exits covered: 100% combined**. These targeted commands
  do not measure whole-repository coverage or exercise the GUI socket server.

Commands used `UV_CACHE_DIR=/private/tmp/moc-domain-fix/uv-cache`, preserving
`uv run --no-sync` and avoiding the known sandbox restriction on the default
cache. No privilege escalation, commit, push, or environment synchronization
was performed. `scripts/check.sh --full --check-only` remains for the reviewer.


## Opt-in characteristic-net output

`include_characteristic_net=False` preserves the default numerical path: no
row snapshots or net arrays are collected and no additional Mach evaluations
are performed. When enabled,
`CharacteristicNet.rows` snapshots the exact initial line and every `_row`
result, including fans, inverse-wall steps and truncated rows. Owned read-only
SI/radian arrays carry static and streamline-dependent total pressure as well
as Mach. The optional result field follows `_corner_states` so positional
construction stays compatible; `last_complete` remains unchanged.

Each row is C+ from the wall outward. C- connectivity maps node k to previous
node k - 1 + s, with s=2 for regular/end wall steps and s=1 for fans and
inverse-wall steps. Shock nodes start new paths; coincident fan wall nodes
add no length. The inverse-wall foot is interpolated, not a stored node, so
its incoming segment is deliberately absent. The exact Taylor--Maccoll initial
C+ bounds an upstream similarity region, not a marched conical net.

Fast regression tests compare all existing arrays/scalars bit-for-bit and
`last_complete`, inspect every recorded `_row` state, check wall/shock
locations in order, and verify interior C+/C- trapezoidal chord slopes at
1e-8 relative tolerance. Cone-cylinder and second-corner/boattail tests also
compare every reconstructed C- edge, excluding the inverse-wall foot, and
check monotonicity, ownership, read-only flags, bool validation and public
exports. These tests verify output fidelity and connectivity, not convergence
or physical accuracy. The GUI opts in, displays bounded subsets and exports
the complete net, retaining existing applicability limits.

Runtime check (2026-11-10, macOS arm64, Python 3.13.5, existing locked
environment): solver-only `perf_counter` wall times, one warm-up per option
followed by seven measured calls each, alternating option order. No other
validation commands ran during these measurements. Both cases use M=3,
p=101325 Pa, T=288.15 K, 88 starting intervals, 1-degree fan increments,
L_n=0.3 m, R=0.05 m and L_cyl=1 m. The ogive has 32 equal-arc-angle chords.
These are local timings, not a portable performance guarantee; plotting and
CSV preparation are excluded.

| Geometry | Net off median | Net on median | Overhead | Rows / nodes |
| --- | ---: | ---: | ---: | ---: |
| Cone-cylinder | 0.084140 s | 0.087244 s | +3.69% | 86 / 4535 |
| N=32 tangent ogive-cylinder | 0.593883 s | 0.613761 s | +3.35% | 414 / 37535 |


Sandbox validation used `UV_CACHE_DIR=/tmp/aerophysics-uv-cache UV_NO_SYNC=1`
with `uv run --offline` because the default uv cache was inaccessible. Ruff
format/lint and mypy passed; solver plus verification-path tests passed
(116 tests), and the separate characteristic `solver_verification` tier
passed (13 tests). GUI tests passed 442 tests; the sole failure was the
unchanged `test_local_documentation_server`, whose loopback socket bind was
denied by the sandbox. The new body-flow tests all passed. Sphinx HTML built
with `-E -a -W --keep-going` successfully. `git diff --check` passed. No
permissions were elevated. The full `scripts/check.sh --full` gate remains
for the caller's outside-sandbox run, as requested; package and generated-asset
verification have not been claimed here.
