.. _detached-shocks:

Detached shocks
===============

Lengths are in metres and Mach number is the free-stream value. The physical
meaning of a shock is explained in :doc:`shock_waves`.

Physical background
-------------------

A blunt nose cannot turn a supersonic flow through an attached shock.
Instead, a curved bow shock stands a short distance ahead of the body. Near
the axis it behaves like a normal shock; farther out it bends back and
weakens toward a Mach wave. The subsonic region between the shock and the
nose makes the full problem hard to solve analytically, so engineering
practice uses correlations fitted to experiments and numerical solutions.

The detached-shock correlations are implemented separately from the
Rankine--Hugoniot solvers in :mod:`aerophysics.detached_shock`. Two geometries
are available through :class:`~aerophysics.detached_shock.DetachedShockGeometry`:
an axisymmetric sphere or hemispherical nose, and a two-dimensional
cylindrical nose. In both cases :math:`R_n` is the nose radius and
:math:`\Delta` is the axial gap from the body vertex to the shock vertex
(the standoff distance).

Ambrosio--Wortman standoff
--------------------------

The :ref:`Ambrosio--Wortman correlations <ref-ambrosio-wortman-1962>` used by
:func:`aerophysics.detached_shock.shock_standoff_distance` are

.. math::

   \frac{\Delta}{R_n}=0.143\exp\left(\frac{3.24}{M^2}\right)
   \quad\text{(sphere or hemispherical nose)},

.. math::

   \frac{\Delta}{R_n}=0.386\exp\left(\frac{4.67}{M^2}\right)
   \quad\text{(two-dimensional cylindrical nose)}.

The standoff distance shrinks as Mach number increases and approaches a
constant fraction of the nose radius at high Mach number.

Billig shock shape
------------------

:ref:`Billig (1967) <ref-billig-1967>` gives the radius of curvature of the
shock at its vertex:

.. math::

   \frac{R_c}{R_n}=1.143\exp\left[
       \frac{0.54}{(M-1)^{1.2}}\right]
   \quad\text{(sphere or hemispherical nose)},

.. math::

   \frac{R_c}{R_n}=1.386\exp\left[
       \frac{1.8}{(M-1)^{0.75}}\right]
   \quad\text{(two-dimensional cylindrical nose)}.

For a hemispherical or cylindrical nose followed by a parallel afterbody,
:func:`aerophysics.detached_shock.billig_shock_shape` describes the shock as a
hyperbola. It starts at the shock vertex with curvature radius :math:`R_c`
and approaches the free-stream Mach angle :math:`\beta=\sin^{-1}(1/M)` far
from the axis:

.. math::

   x=R_n+\Delta-R_c\cot^2\beta
   \left[
     \sqrt{1+\frac{y^2\tan^2\beta}{R_c^2}}-1
   \right].

The origin is the centre of nose curvature, positive :math:`x` points
upstream, the body vertex is at :math:`x=R_n`, and the shock vertex is at
:math:`x=R_n+\Delta`. The shape always uses the Ambrosio--Wortman value of
:math:`\Delta`, even when a Seiff estimate is also displayed. The transverse
coordinate :math:`y` is a one-dimensional array and is appended as the last
output axis, so broadcast Mach and radius cases keep their own axes. Close to
:math:`M=1` the fitted curvature grows without bound. Inputs that would give
non-finite coordinates raise :class:`ValueError`.

Seiff density-ratio standoff
----------------------------

:ref:`Seiff's relation <ref-seiff-1964>` links the standoff distance of a
sphere to the density ratio across the normal part of the shock. A stronger
density rise compresses the shock layer, so the shock sits closer to the
body:

.. math::

   \frac{\Delta}{R_n}=\frac{0.78}{\rho_2/\rho_1}.

:func:`aerophysics.detached_shock.seiff_standoff_distance` accepts
:math:`\rho_2/\rho_1` directly, so it is not tied to a Mach number. The
convenience function
:func:`aerophysics.detached_shock.seiff_standoff_distance_from_mach` obtains
the density ratio from the calorically perfect
:func:`aerophysics.shocks.normal_shock`. The Seiff relation is for spheres
only; no cylindrical version is offered.
:func:`aerophysics.detached_shock.compare_standoff_distances` returns the
Ambrosio--Wortman and Seiff sphere results together with their absolute and
relative differences.

All functions require finite :math:`M>1` and :math:`R_n>0`; the Seiff
function also requires finite :math:`\rho_2/\rho_1>1`. Invalid values raise
:class:`ValueError`. The original publications do not state an explicit
Mach-number range for their fits, so no range warning is issued.
:ref:`NASA TN D-2780 <ref-inouye-1965>` compares Seiff's relation for
:math:`0.04<\rho_1/\rho_2<0.16`; this is a verification interval, not the full
range of validity. Its tabulated air solution at :math:`M_\infty=8.949`
(:math:`\rho_1/\rho_2=0.1253`, :math:`\Delta/R_n=0.0994`) agrees with the
Seiff value ``0.097734`` within the ``0.002`` verification tolerance.

These correlations assume continuum, steady, low-temperature flow. They are
normally used for calorically perfect air, and Billig's curves were fitted
mainly for :math:`\gamma=1.4`. They do not calculate the temperature or
pressure inside the shock layer. Real-gas Seiff models, Beattie--Bridgeman
shock states, rarefied-flow corrections, and shock fitting are outside this
implementation.

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
