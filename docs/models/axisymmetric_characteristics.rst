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
Arrays are owned and read-only. Surface distance follows the meridian from
the sharp tip, where the state is the cone-surface limit.

A convex vertex has one public surface station, with the **downstream** fan
state. All fan states at that same coordinate are used internally. Keeping
both limits would violate the strictly increasing distance required by the
boundary-layer edge API. The boundary-layer conversion separately inserts
pre-corner and midpoint stations; the public inviscid arrays are unchanged.
This does not model viscous interaction at a sharp corner.

Shock stations may extend beyond the body's final axial coordinate because
a corner's leading wave can first reach the shock farther downstream.
Characteristics that do not reach the shock within the needed downstream
region are truncated without imposing an outflow state. Successful downstream
shock extrapolations are retained. If an angle iteration leaves the admissible
shock interval and the geometric continuation at the previous shock tangent
is beyond the body endpoint, that extrapolation is abandoned. The trace then
stops at the last fitted shock point, which can be before the body endpoint;
no terminal shock state or complete outer flow field is promised. Do not
extrapolate missing shock data to construct a terminal transverse profile.
Subsonic states, crossed characteristics and failures inside the requested
domain still raise location-bearing errors. The characteristic net is not
exposed as a public field solution.

For the M3/9.46-degree and M6/4.76-degree cone--cylinder verification cases,
the cylinder length is ten diameters. The design record
``devdocs/axisymmetric_characteristics_design.md`` records grid differences,
observed order, mass balance, runtime and boundary-layer coupling tests.
The original M1 cone-only result, M5 at the exact expansion vertex, and the
constant-edge cone part of M7 agree **by construction**: they check the API,
initialization and boundary-layer connection. They are not evidence for the
accuracy of the MOC march itself.

An additional private verification starts the characteristic net at 5 percent
of a cone's length with zero turn, and marches the remaining cone through the
same wall, interior and shock unit processes. For M2/10 degrees, maximum
relative errors at the default count 88 are 0.09955 percent in pressure,
0.03534 percent in Mach and 0.01383 percent in shock-ray position. Count 87
has pressure error 0.10177 percent; 88 is the smallest count reaching the
requested 0.1 percent criterion. For M6/5 degrees the errors are 0.005382,
0.000892 and 0.005287 percent. The retained 24/48/96/192 refinement study
gives final observed orders 1.970/1.970/2.020 and 1.999/1.999/1.999,
respectively. The normal public conical portion uses the exact similarity
solution. The default M3 cone--cylinder calculation takes about 0.08--0.11 s
on the development machine (cylinder length ten diameters).

For the cone--cylinder refinement cases, the estimated maximum relative
pressure errors of the default 88/1-degree net are 1.70e-5 (M3) and 2.82e-5
(M6) within the first 0.1 cone lengths after the shoulder, and 7.80e-6 / 6.31e-5
farther downstream. These are numerical Richardson estimates; they do not
include physical-model error. Increasing only the characteristic count does
not guarantee decreasing error: refine the fan angle step as well.

Comparison with the reviewer-supplied digitization of
:ref:`NACA TN 2250, Figure 1 <ref-ehret-rossow-stevens-1950>` gives maximum
absolute differences 0.01035 (M3, l/d = 3) and 0.01401 (M6, l/d = 6) in
(p - p_inf)/p_inf, within the prescribed bound 0.02. Cone length l is measured
from tip to shoulder and the cylinder ends at x/l = 3.4. The source solution
is **irrotational graphical MOC**, not a rotational solution or experiment.
Reading uncertainties are +/-0.005 vertically and +/-1 percentage point in
100*x/l. Experimental accuracy remains unestablished.

Boundary-layer connection
--------------------------

The body streamline crosses only the original conical shock. Its total
pressure is constant, so
:meth:`~aerophysics.characteristics.AxisymmetricBodyFlowResult.boundary_layer_inputs`
uses :meth:`~aerophysics.marching_boundary_layer.BoundaryLayerEdge.from_pressure`
with the post-conical-shock total pressure and the uniform total temperature.
It returns :class:`~aerophysics.characteristics.BoundaryLayerInputs`, a named
pair of ``edge`` and ``wall_radius``. The owned, read-only radius array matches
``edge.distance`` exactly by piecewise linear interpolation of the meridian.
Always pass this array, rather than the original flow-result radius, to the
boundary-layer solver. The original perfect-gas and viscosity models are
preserved. See the runnable example in :doc:`../guides/boundary_layers`.

At every convex corner at surface distance ``s_c``, conversion inserts the
pre-expansion wall pressure at ``s_c - w`` and a logarithmic-pressure midpoint
at ``s_c - w/2``. This uses the local marched pre-fan state at later corners,
not the original cone state. The corner itself retains the downstream state.
The default width is ``w = 0.001 * s_c``; an explicit ``corner_length`` sets
a common width in metres. A nonpositive or nonfinite width, a width reaching
the preceding output point, or unresolvable inserted coordinates raises
``ValueError``. No existing points are discarded.

The midpoint ensures that the public boundary-layer march samples the steep
pressure gradient and issues its pressure-length warning for the default
cone--cylinder case. The general boundary-layer warning samples only output
stations; without the midpoint the zero PCHIP endpoint slopes can hide a
large interior gradient. The boundary-layer module itself is unchanged.

Thin-layer predictions near the shoulder, including at least several tens
of upstream delta99 downstream, are not reliable local corner solutions.
They depend on the artificial width ``w``; the extra points do not resolve
viscous interaction or displacement feedback. Width-sensitivity measurements
for M3, cone length 0.3 m, radius 0.05 m, cylinder length 1 m and wall
temperature 300 K are recorded in the design memo. Here the default width
is 0.304138 mm. Relative to a 0.2 mm width, at 5 upstream thicknesses the
default changes laminar Cf by -0.617 percent, heat flux by -0.442 percent
and delta99 by +0.046 percent; the SA changes are -0.016, +0.006 and
-0.018 percent. Relative to 1 mm the laminar differences are +3.825,
+2.271 and -0.295 percent (SA +0.100, -0.031 and +0.114 percent).
At 20 thicknesses, default versus 1 mm still changes laminar Cf by
+0.877 percent and heat flux by +0.553 percent. At 100 thicknesses these
differences are +0.105 and +0.043 percent. Small width differences do not
establish physical validity; some SA differences are comparable to the
measured streamwise discretization uncertainty. The design memo gives all
three quantities at 5/20/100 thicknesses and the step-refinement check.

The classical characteristic framework is described by
:ref:`Zucrow and Hoffman <ref-zucrow-hoffman-1977>` and
:ref:`Anderson <ref-anderson-moc>`. Bibliographic verification and the limits
of access to those chapters are recorded with the references.
