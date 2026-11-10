Axisymmetric body characteristics
=================================

A cone solution stops describing the outer flow at a shoulder. Expansion
waves there interact with the axisymmetric flow between the surface and the
leading shock. :func:`~aerophysics.characteristics.axisymmetric_body_flow`
computes this inviscid outer flow, including entropy differences between
streamlines caused by the curved shock. It supplies pressure, temperature,
Mach number and flow direction on the surface, and a fitted shock trace.

Inputs and limits
-----------------

Supply a piecewise linear meridian in axial distance :math:`x` and radius
:math:`r`, both in metres. The first vertex must be :math:`(0,0)`, all later
radii must be positive, and axial distance must increase strictly. Segment
angles must decrease or remain constant: cones, cylinders and convex
boattails are supported. Concave corners require additional shocks and are
rejected. Free-stream pressure is in Pa, temperature in K and Mach number
must exceed one. Angles, including ``fan_angle_step``, are in radians.

The model assumes a calorically perfect gas, uniform free stream, attached
weak nose shock, zero incidence, no swirl, and inviscid adiabatic flow. It
rejects subsonic states, characteristics that cannot march downstream,
coalescing/crossed characteristics, and failed numerical unit processes with
:class:`~aerophysics.exceptions.ModelRangeError` and a location. It does not
model blunt noses, concave corners, internal flow, high-temperature chemistry,
viscous displacement or separation. A detected breakdown is not a prediction
of the flow beyond that location.

Characteristic equations
-------------------------

Let :math:`\theta` be flow angle from the positive x axis, :math:`q` speed,
:math:`\mu=\sin^{-1}(1/M)` the Mach angle, and
:math:`P=\ln(p/p_\infty)`. The two characteristic families satisfy

.. math::

   \frac{dr}{dx}=\tan(\theta\pm\mu),\qquad
   d\theta\pm K\,dP=H_\pm\,dx,

   K=\frac{\sqrt{M^2-1}}{\gamma M^2},\qquad
   H_\pm=\mp\frac{\sin\mu\sin\theta}{r\cos(\theta\pm\mu)}.

Total temperature is uniform. Along each streamline
:math:`dr/dx=\tan\theta`, total pressure is constant, but it can differ
between streamlines. Consequently this is a rotational formulation; it does
not impose a single isentropic relation across the entire shock layer.
At each new shock point, the local shock tangent and the incoming
characteristic compatibility equation are solved together with the
Rankine--Hugoniot relations. The local shock strength sets the total pressure
of newly entering fluid.

The initial conical region is exact Taylor--Maccoll flow. Its leading
outward characteristic through the first convex vertex supplies the starting
line. The implementation reuses the Taylor--Maccoll differential equation
from :mod:`aerophysics.shocks`. Internal points, wall points and shock points
use a predictor and iterated trapezoidal corrector. A centered fan uses the
local planar Prandtl--Meyer turn at its vertex; its finite-length waves then
include the axisymmetric source term.

Resolution and output
---------------------

The default ``characteristic_count=88`` divides the starting characteristic;
``fan_angle_step`` defaults to one degree expressed in radians. Double the
count and halve the fan step to refine. Surface stations are the computed
wall intersections, plus analytic cone samples and exact geometry vertices.
Result arrays are read-only. Surface distance follows the meridian from
the sharp tip, where the state is the cone-surface limit.

At a corner the surface state jumps from the value before the expansion to the
value after it. The result holds one station at the corner, with the state
**after** the expansion. See `Boundary-layer connection`_ for how the jump is
handled when the result is passed to the boundary-layer calculation.

The shock trace can extend beyond the end of the body, and it can also stop
before it: the trace ends at the last shock point that the calculation
actually reached. No values are extrapolated beyond that point, so do not use
the trace to construct a flow profile at the end of the body. The surface
result is always complete; if the surface cannot be completed, the function
raises ``ModelRangeError`` and names the location. On a long body with a
coarse fan the message may advise a smaller ``fan_angle_step``.

Accuracy
--------

The conical nose region uses the exact Taylor--Maccoll solution. For the rest
of the body, the numerical error of the default resolution was estimated by
grid refinement on M3 and M6 cone--cylinders ten diameters long: the surface
pressure error is below 0.01 percent. A separate test that marches a plain
cone with the same method agrees with the exact cone solution to within
0.1 percent in pressure at M2 and 10 degrees. These figures describe the
numerical method only, not the physical model.

Compared with the graphical characteristics solutions of
:ref:`NACA TN 2250, Figure 1 <ref-ehret-rossow-stevens-1950>`, the surface
pressure :math:`(p-p_\infty)/p_\infty` differs by at most 0.010 (M3) and 0.014
(M6). The figure can be read to about +/-0.005. That source is itself a
calculation that neglects entropy differences between streamlines, so
agreement with experiment has not been established.

Increasing only ``characteristic_count`` does not always reduce the error;
refine ``fan_angle_step`` as well. The default M3 cone--cylinder calculation
takes about 0.1 s. The full measurements are in
``devdocs/axisymmetric_characteristics_design.md``.

Boundary-layer connection
--------------------------

The body streamline crosses only the original conical shock. Its total
pressure is constant, so
:meth:`~aerophysics.characteristics.AxisymmetricBodyFlowResult.boundary_layer_inputs`
uses :meth:`~aerophysics.marching_boundary_layer.BoundaryLayerEdge.from_pressure`
with the post-conical-shock total pressure and the uniform total temperature.
It returns :class:`~aerophysics.characteristics.BoundaryLayerInputs`, a named
pair of ``edge`` and ``wall_radius``. The radius array matches
``edge.distance`` point for point. Always pass this array, rather than the
radius stored in the flow result, to the boundary-layer solver. See the runnable example in :doc:`../guides/marching_boundary_layers`.

A sharp corner would give the boundary layer two different outer states at
one position. The conversion therefore replaces the jump by a very short
pressure ramp that ends at the corner. Its length is ``0.001`` times the
surface distance of the corner by default; pass ``corner_length`` in metres to
``boundary_layer_inputs`` to choose another length. The ramp must fit between
the corner and the output point before it, otherwise ``ValueError`` is raised.

The ramp is artificial. A real boundary layer changes the outer flow near a
corner, which this one-way calculation does not include. The boundary-layer
march therefore warns about the steep pressure change, and its results near
the corner are not reliable. They also depend on the ramp length. For an M3
cone--cylinder (cone 0.3 m, radius 0.05 m, wall at 300 K), changing the ramp
from the default 0.3 mm to 1 mm changed the laminar Cf by 3.8 percent at
5 boundary-layer thicknesses downstream of the corner, 0.9 percent at 20, and
0.1 percent at 100. The turbulent Cf changed by 0.1 percent at 5 thicknesses.

The classical characteristic framework is described by
:ref:`Zucrow and Hoffman <ref-zucrow-hoffman-1977>` and
:ref:`Anderson <ref-anderson-moc>`. Bibliographic verification and the limits
of access to those chapters are recorded with the references.
