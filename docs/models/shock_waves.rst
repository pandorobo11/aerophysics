Shock waves
===========

These relations assume steady flow. State ratios are downstream over upstream,
and angles are in radians. Convert explicitly with
:func:`aerophysics.units.degrees_to_radians` and
:func:`aerophysics.units.radians_to_degrees`.

The normal, oblique, and conical solvers use a calorically perfect gas.
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

As :math:`M_1\to1^+`, the Mach angle approaches :math:`\pi/2` and the
available attached-shock-angle interval collapses.  The zero-half-angle limit
remains the Mach wave.  For a positive cone half-angle, an interval or
Taylor--Maccoll root that is too degenerate to resolve raises
:class:`~aerophysics.exceptions.NoAttachedShockError`; low-level integration
and bracketing exceptions are not part of the public API.  Cone half-angle is
positive away from the axis, and a returned attached shock satisfies
:math:`\theta_c<\beta<\pi/2`.

The model assumes a calorically perfect gas, a sharp circular cone, zero angle
of attack, steady inviscid adiabatic flow, and an attached axisymmetric shock.
It does not model bluntness, viscosity, real-gas effects, or asymmetric cone
flow.
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
shock-layer thermodynamics. Real-gas Seiff models, NASA7/NASA9 or harmonic-
oscillator general normal-shock solvers, Beattie--Bridgeman shock states,
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
