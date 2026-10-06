Shock waves
===========

These relations assume steady flow. State ratios are downstream over upstream,
and angles are in radians. Convert explicitly with
:func:`aerophysics.units.degrees_to_radians` and
:func:`aerophysics.units.radians_to_degrees`.

The normal, oblique, and conical solvers accept calorically or thermally perfect
gases. The Rayleigh--Pitot formula uses a calorically perfect gas.
Detached-shock engineering correlations provide standoff distance and shape;
they do not solve shock-layer thermodynamics. For task-oriented model
selection and examples, see :doc:`../guides/compressible_flow`.

Normal shocks
-------------

For upstream Mach number :math:`M_1>1`,

.. math::

   M_2^2
   =\frac{1+\frac{\gamma-1}{2}M_1^2}
          {\gamma M_1^2-\frac{\gamma-1}{2}},

.. math::

   \frac{p_2}{p_1}
   =1+\frac{2\gamma}{\gamma+1}(M_1^2-1),
   \qquad
   \frac{\rho_2}{\rho_1}
   =\frac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2},
   \qquad
   \frac{T_2}{T_1}=\frac{p_2/p_1}{\rho_2/\rho_1}.

The downstream-to-upstream total-pressure ratio is

.. math::

   \frac{p_{02}}{p_{01}}
   =\left[\frac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^
       {\gamma/(\gamma-1)}
    \left[\frac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^
       {1/(\gamma-1)}.

:func:`aerophysics.shocks.supersonic_pitot_pressure_ratio` returns the
Rayleigh--Pitot ratio

.. math::

   \frac{p_{02}}{p_1}
   =\left[\frac{(\gamma+1)M_1^2}{2}\right]^{\gamma/(\gamma-1)}
    \left[\frac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^
       {1/(\gamma-1)}.

Use :func:`aerophysics.shocks.normal_shock` for the complete result.
The one-dimensional relations follow
:ref:`NACA Report 1135 <ref-naca-report-1135>`.

Oblique shocks
--------------

The theta--beta--Mach relation is

.. math::

   \tan\theta
   =2\cot\beta\,
    \frac{M_1^2\sin^2\beta-1}
         {M_1^2(\gamma+\cos 2\beta)+2}.

The shock angle lies between the Mach angle
:math:`\mu_M=\sin^{-1}(1/M_1)` and :math:`\pi/2`. The normal component
reduces the state calculation to the normal-shock equations:

.. math::

   M_{n1}=M_1\sin\beta,
   \qquad
   M_2=\frac{M_{n2}}{\sin(\beta-\theta)}.

:func:`aerophysics.shocks.oblique_shock` defaults to the weak root. Select
:class:`~aerophysics.shocks.ShockBranch` explicitly when branch identity
matters. If :math:`\theta` exceeds the maximum attached-shock deflection,
:class:`~aerophysics.exceptions.NoAttachedShockError` is raised; the requested
solution is not silently replaced by a detached normal shock.
The theta--beta--Mach relation and branch convention follow
:ref:`NACA Report 1135 <ref-naca-report-1135>`.

.. list-table:: Shock symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - Unit
   * - :math:`M_1,M_2`
     - ``upstream_mach``, ``downstream_mach``
     - Upstream and downstream Mach numbers
     - dimensionless
   * - :math:`M_{n1},M_{n2}`
     - ``upstream_normal_mach``, ``downstream_normal_mach``
     - Mach components normal to the shock
     - dimensionless
   * - :math:`\theta`
     - ``deflection_angle``
     - Flow-deflection angle
     - rad
   * - :math:`\beta`
     - ``shock_angle``
     - Shock angle from the upstream velocity
     - rad
   * - :math:`p_{01},p_{02}`
     - total pressure
     - Upstream and downstream total pressures
     - Pa

>>> from aerophysics import ShockBranch, oblique_shock
>>> from aerophysics.units import degrees_to_radians, radians_to_degrees
>>> shock = oblique_shock(
...     2.0, degrees_to_radians(10.0), branch=ShockBranch.WEAK
... )
>>> round(radians_to_degrees(shock.shock_angle), 3)
39.314

Thermally perfect normal and oblique shocks
-----------------------------------------------

Pass a :class:`~aerophysics.thermochemistry.ThermallyPerfectGas` (including
``AIR_NASA7`` or ``AIR_NASA9``) or a
:class:`~aerophysics.real_gas.HarmonicOscillatorGas`, and provide
``upstream_temperature`` as the upstream **static** temperature in kelvin.
Mach number, angle, and temperature broadcast together. The default remains
the constant-:math:`\gamma` ``AIR`` model, with the original result fields and
downstream/upstream ratio convention.

The solver assumes steady, inviscid, adiabatic, two-dimensional ideal-gas
flow with fixed composition, thermal equilibrium, and :math:`p=\rho RT`.
The normal velocity components :math:`u_n` satisfy

.. math::

   \rho_1 u_{n1}=\rho_2 u_{n2},\qquad
   p_1+\rho_1u_{n1}^2=p_2+\rho_2u_{n2}^2,\qquad
   h(T_2)-h(T_1)=\frac{u_{n1}^2-u_{n2}^2}{2}.

Tangential velocity is unchanged. For a trial shock angle,
:math:`u_{n1}=M_1a(T_1)\sin\beta` and the conservation equations determine
:math:`T_2` and :math:`\rho_2/\rho_1`. The flow turn follows

.. math::

   \theta=\beta-\tan^{-1}\left[
      \frac{\rho_1}{\rho_2}\tan\beta\right],\qquad
   a(T)^2=\gamma(T)RT,\qquad \gamma(T)=\frac{c_p(T)}{c_p(T)-R}.

The temperature Hugoniot eliminates the zero-strength root using
:math:`x=u_{n2}/u_{n1}`, :math:`\tau=T_2/T_1`, and
:math:`q=[h(T_2)-h(T_1)]/(RT_1)`:

.. math::

   x^2+(2q+1-\tau)x-\tau=0,\qquad
   u_{n1}^2=RT_1\frac{\tau/x-1}{1-x}.

The sonic endpoint is evaluated analytically; the positive-temperature
compressive root is bracketed numerically. The solver then maximizes the
shock polar and brackets the angle root separately on the weak and strong
branches. ``theta_from_shock_angle``, ``shock_angle``,
``maximum_attached_deflection``, ``normal_shock``, and ``oblique_shock`` all
accept the same thermal gas and static-temperature arguments.

Since total enthalpy, composition, and total temperature are unchanged,
the total-pressure loss follows directly from entropy production:

.. math::

   \Delta s=s^\circ(T_2)-s^\circ(T_1)-R\ln(p_2/p_1),\qquad
   p_{02}/p_{01}=\exp(-\Delta s/R).

This does not require stagnation temperature to lie within the polynomial fit.
In contrast, the absolute pitot ratio
:func:`~aerophysics.shocks.supersonic_pitot_pressure_ratio` also needs the
stagnation state. For a thermal gas, supply ``upstream_temperature`` in kelvin;
Mach and temperature inputs broadcast together. The solver finds :math:`T_0`
from :math:`h(T_0)=h(T_1)+u_1^2/2`, then evaluates

.. math::

   \frac{p_{02}}{p_1}=\frac{p_{02}}{p_{01}}
   \exp\!\left(\frac{s^\circ(T_0)-s^\circ(T_1)}{R}\right).

The shock remains frozen, steady, and adiabatic. This calculation raises
:class:`~aerophysics.exceptions.ModelRangeError` when any required static or
stagnation temperature is outside the model range; it never extrapolates.
This stricter requirement does not change ``normal_shock`` or its
``total_pressure_ratio`` contract.

The governing thermal shock relations follow
:ref:`Tatum (1996), NASA CR-4749 <ref-tatum-1996>`.
Verification includes constant-heat-capacity limits, conservation residuals,
and an independent variable-heat-capacity solution; see
:doc:`../verification/compressible_flow`.

Temperature ranges are inclusive: both dry-air NASA presets currently share
200--6000 K, and ``AIR_HARMONIC_OSCILLATOR`` documents 400--2000 K. A custom
mixture uses the intersection of its species' fitted ranges. Required static
states outside these ranges raise :class:`~aerophysics.exceptions.ModelRangeError`;
the shock solver never extrapolates. A valid weak shock can still be returned
when the normal shock or polar maximum exceeds the range. An unavailable
strong root or attached limit raises ``ModelRangeError``, not a claim of
physical detachment. For custom caloric models, the branch construction
assumes the usual single-maximum, convex-gas shock polar.

These polynomial ranges describe thermodynamic data, not a certified physical
validity interval for chemically frozen air. Dissociation, ionization,
finite-rate chemistry, vibrational nonequilibrium, and dense-gas effects are
excluded; their importance depends on the actual pressure, composition, and
flow residence time. A frozen calculation at high temperature must be
interpreted with those assumptions.

>>> from aerophysics import AIR_NASA9
>>> thermal = oblique_shock(
...     3.0, degrees_to_radians(20.0), gas=AIR_NASA9,
...     upstream_temperature=300.0,
... )
>>> round(radians_to_degrees(thermal.shock_angle), 3)
37.68
>>> round(300.0 * thermal.static_temperature_ratio, 3)
466.118
>>> round(thermal.downstream_mach, 6)
2.006457
>>> round(thermal.total_pressure_ratio, 6)
0.796621

Conical shocks
--------------

For inviscid axisymmetric flow over a sharp circular cone at zero angle of
attack, the velocity between the shock and cone surface varies with polar
angle.  With radial and polar velocity components nondimensionalized by the
limiting velocity available from adiabatic expansion into a vacuum, the
Taylor--Maccoll equations are

.. math::

   \frac{dV_r}{d\theta}=V_\theta,

.. math::

   \frac{dV_\theta}{d\theta}
   =\frac{V_rV_\theta^2-a^2(2V_r+V_\theta\cot\theta)}
          {a^2-V_\theta^2},
   \qquad
   a^2=\frac{\gamma-1}{2}(1-V_r^2-V_\theta^2).

The Rankine--Hugoniot relations supply the velocity immediately behind a
trial shock angle :math:`\beta`.  Integration toward the axis locates the
cone surface where :math:`V_\theta=0`.  The weak attached solution is the
first shock angle above the Mach angle that produces the requested cone
half-angle :math:`\theta_c`.

:func:`aerophysics.shocks.conical_shock` returns the shock angle, Mach numbers
immediately behind the shock and at the cone surface, surface static-state
ratios over the free stream, and the post-shock/free-stream total-pressure
ratio.  :func:`aerophysics.shocks.maximum_attached_cone_angle` returns the
attached-shock limit.  A larger cone half-angle raises
:class:`~aerophysics.exceptions.NoAttachedShockError` rather than substituting
a detached-shock approximation.

The model assumes a sharp circular cone, zero angle
of attack, steady inviscid adiabatic flow, and an attached axisymmetric shock.
It does not model bluntness, viscosity, dense-gas effects, or asymmetric cone
flow. The default remains the calorically perfect ``AIR`` model.
Reference solutions for the Taylor--Maccoll model are tabulated by
:ref:`Sims (1964) <ref-sims-1964>`.

>>> from aerophysics import conical_shock
>>> cone = conical_shock(2.0, degrees_to_radians(10.0))
>>> round(radians_to_degrees(cone.shock_angle), 3)
31.206
>>> round(cone.surface_mach, 3)
1.834
>>> round(cone.surface_pressure_ratio, 3)
1.293

For a frozen thermally perfect gas, pass ``upstream_temperature`` in kelvin
to ``conical_shock`` or ``maximum_attached_cone_angle``. Mach number, cone
half-angle, and temperature broadcast together. The velocity equations have
the same form, but the thermal path scales velocities by the free-stream
speed :math:`U_\infty`, rather than a vacuum limiting velocity. It closes
the equations using

.. math::

   h(T)=h(T_2)+\frac{U_\infty^2}{2}
      \left(V_{r2}^2+V_{\theta2}^2-V_r^2-V_\theta^2\right),
   \qquad a^2=\frac{\gamma(T)RT}{U_\infty^2},
   \qquad \gamma(T)=\frac{c_p(T)}{c_p(T)-R}.

Here the subscript 2 denotes the state immediately behind the shock;
:math:`\theta` increases away from the cone axis, so :math:`V_\theta<0`
between the shock and the surface. The thermal normal-shock conservation
solver supplies the initial state, with tangential velocity unchanged.
Enthalpy inversion uses the gas model's actual reference, including NASA
polynomial enthalpy offsets. A local gamma substitution in the calorically
perfect energy relation is not used.

The shock fixes entropy. The smooth conical flow is isentropic, giving

.. math::

   \frac{p_s}{p_2}
   =\exp\left[\frac{s^\circ(T_s)-s^\circ(T_2)}{R}\right],
   \qquad \frac{\rho_s}{\rho_\infty}
   =\frac{p_s/p_\infty}{T_s/T_\infty}.

The total-pressure ratio is the shock loss, unchanged from shock to surface.
All **static** states along that trajectory must remain in the inclusive
gas temperature range. A temperature-boundary event and a direct surface
enthalpy check enforce this limit, including steps that straddle the surface.
Internal Runge--Kutta trial stages evaluate only in-range properties;
no thermodynamic property is extrapolated. A valid weak cone can be returned
even if its normal shock or physical cone-angle maximum exceeds the range.
An inaccessible requested solution or physical limit raises
:class:`~aerophysics.exceptions.ModelRangeError`; a cone beyond a resolved
physical maximum raises ``NoAttachedShockError``. Only the weak branch is
exposed, as in the original conical API. The search assumes the usual
single-maximum cone-angle curve for a convex caloric gas. Composition changes,
chemical reactions, and vibrational nonequilibrium are excluded.

The temperature-boundary classification uses the shock-angle separation
between the optimized maximum and the boundary, accounting for the numerical
maximizer's stopping precision. An angle-height difference alone is
insufficient near the flat top of the cone-angle curve. Numerical integration
or root-resolution failure raises
:class:`~aerophysics.exceptions.ShockConvergenceError`, distinct from physical
detachment or temperature-range failure. Very slender cones can reach this
limit because the shock angle approaches the Mach wave within floating-point
resolution; the solver does not inflate its angular tolerance to return a state.

>>> from aerophysics import AIR_HARMONIC_OSCILLATOR
>>> thermal_cone = conical_shock(
...     3.0, degrees_to_radians(10.0), gas=AIR_HARMONIC_OSCILLATOR,
...     upstream_temperature=1000.0,
... )
>>> round(radians_to_degrees(thermal_cone.shock_angle), 3)
21.626
>>> round(1000.0 * thermal_cone.surface_temperature_ratio, 3)
1110.946

.. _detached-shocks:

Detached shocks
---------------

The detached-shock correlations are implemented separately from the
Rankine--Hugoniot solvers in :mod:`aerophysics.detached_shock`. Two geometries
are explicit: an axisymmetric sphere or hemispherical nose, and a
two-dimensional cylindrical nose. In both cases :math:`R_n` is the nose
radius and :math:`\Delta` is the axial gap from the body vertex to the shock
vertex.

Ambrosio--Wortman standoff
^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`Ambrosio--Wortman correlations <ref-ambrosio-wortman-1962>` used by
:func:`aerophysics.detached_shock.shock_standoff_distance` are

.. math::

   \frac{\Delta}{R_n}=0.143\exp\left(\frac{3.24}{M^2}\right)
   \quad\text{(sphere or hemispherical nose)},

.. math::

   \frac{\Delta}{R_n}=0.386\exp\left(\frac{4.67}{M^2}\right)
   \quad\text{(two-dimensional cylindrical nose)}.

Billig shock shape
^^^^^^^^^^^^^^^^^^

:ref:`Billig (1967) <ref-billig-1967>` gives the shock-vertex curvature radius

.. math::

   \frac{R_c}{R_n}=1.143\exp\left[
       \frac{0.54}{(M-1)^{1.2}}\right]
   \quad\text{(sphere or hemispherical nose)},

.. math::

   \frac{R_c}{R_n}=1.386\exp\left[
       \frac{1.8}{(M-1)^{0.75}}\right]
   \quad\text{(two-dimensional cylindrical nose)}.

For a hemispherical or cylindrical nose followed by a parallel afterbody,
:func:`aerophysics.detached_shock.billig_shock_shape` uses
:math:`\beta=\sin^{-1}(1/M)` and the hyperbola

.. math::

   x=R_n+\Delta-R_c\cot^2\beta
   \left[
     \sqrt{1+\frac{y^2\tan^2\beta}{R_c^2}}-1
   \right].

The nose-curvature center is the origin, positive :math:`x` points upstream,
the body vertex is at :math:`x=R_n`, and the shock vertex is at
:math:`x=R_n+\Delta`. Billig shape calculations deliberately use the
Ambrosio--Wortman value of :math:`\Delta`; changing the displayed Seiff model
does not change that shape convention. A one-dimensional transverse
coordinate array is appended as the final output axis, so broadcast Mach and
radius cases retain their case axes.  The implementation evaluates the
curvature product in logarithmic form and uses a cancellation-resistant form
of the hyperbola increment.  Because Billig's fitted curvature diverges as
:math:`M\to1^+`, a mathematically admissible input can still exceed finite
``float64`` representation; unrepresentable curvature or non-finite shock
coordinates raise :class:`ValueError` rather than returning infinities.

Seiff density-ratio standoff
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a sphere, :ref:`Seiff's density-ratio relation <ref-seiff-1964>` is

.. math::

   \frac{\Delta}{R_n}=\frac{0.78}{\rho_2/\rho_1}.

:func:`aerophysics.detached_shock.seiff_standoff_distance` accepts
:math:`\rho_2/\rho_1` directly and therefore does not attach a Mach number to
its result. The convenience API
:func:`aerophysics.detached_shock.seiff_standoff_distance_from_mach` obtains
the density ratio from the existing calorically perfect-gas
:func:`aerophysics.shocks.normal_shock`. No cylindrical Seiff correlation is
claimed or accepted. :func:`aerophysics.detached_shock.compare_standoff_distances`
returns the Ambrosio--Wortman and Seiff sphere results together with their
signed and relative differences.

All functions require finite :math:`M>1` and :math:`R_n>0`; the low-level
Seiff function additionally requires finite :math:`\rho_2/\rho_1>1`. These
are physical and mathematical domains, not empirical fit limits, and invalid
values raise :class:`ValueError`. The cited original publications do not give
a sufficiently explicit numerical Mach fit range to justify inventing an
``ApplicabilityWarning`` or ``ModelRangeError`` boundary.
:ref:`NASA TN D-2780 <ref-inouye-1965>` independently compares
:math:`0.04<\rho_1/\rho_2<0.16`; this is recorded as a verification interval,
not as the full validity range of Seiff's correlation.
Its Table I air solution at :math:`M_\infty=8.949` independently tabulates
:math:`\rho_1/\rho_2=0.1253` and :math:`\Delta/R_b=0.0994`, where the table's
sphere nose radius :math:`R_b` is :math:`R_n` in this API.  The Seiff relation
gives ``0.097734`` for that printed density ratio, within the committed
``0.002`` absolute comparison tolerance.

These engineering correlations assume continuum, steady, low-temperature
flow. Their common use is for calorically perfect air; Billig's fitted curves
are primarily associated with :math:`\gamma=1.4`. They do not solve the
shock-layer thermodynamics. Real-gas Seiff models, Beattie--Bridgeman shock states,
rarefied-flow corrections, and shock fitting are outside this implementation.

>>> from aerophysics import DetachedShockGeometry, billig_shock_shape
>>> from aerophysics import compare_standoff_distances, shock_standoff_distance
>>> sphere = shock_standoff_distance(
...     4.0, 0.5, geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE
... )
>>> round(sphere.normalized_standoff_distance, 6)
0.175098
>>> shape = billig_shock_shape(
...     4.0, 0.5, [-1.0, 0.0, 1.0],
...     geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
... )
>>> bool(shape.shock_x[1] == sphere.nose_radius + sphere.standoff_distance)
True
>>> compare_standoff_distances(4.0, 0.5).seiff.density_ratio > 1.0
True

The formulas and geometry definitions follow
:ref:`Ambrosio and Wortman <ref-ambrosio-wortman-1962>`,
:ref:`Billig <ref-billig-1967>`, and :ref:`Seiff <ref-seiff-1964>`; the
independent Seiff cross-check follows :ref:`Inouye <ref-inouye-1965>`.
