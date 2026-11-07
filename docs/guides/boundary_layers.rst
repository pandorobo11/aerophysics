.. _boundary-layer-workflows:

Boundary-layer workflows
========================

Use these workflows to estimate a smooth flat-plate boundary layer, construct
a compressible turbulent mean profile, or account for the reduced dynamic
pressure acting on an immersed protrusion. All dimensional inputs and outputs
use SI units.

Choose a workflow
-----------------

.. list-table:: Boundary-layer tasks
   :header-rows: 1
   :widths: 23 30 27 20

   * - Task
     - Start with
     - Required state
     - Detailed model
   * - Estimate thickness, friction, and one-sided plate drag
     - ``flat_plate_boundary_layer``
     - Distance and edge velocity, density, and viscosity
     - :doc:`../models/flat_plate_boundary_layer`
   * - Follow a laminar or SA-noft2 layer through a changing edge pressure
     - ``marching_boundary_layer``
     - Edge pressure and total state; optional wall radii for thin axisymmetric flow
     - :doc:`../models/marching_boundary_layer`
   * - Transform a measured or computed compressible velocity profile
     - ``transform_compressible_velocity_profile``
     - Wall-normal velocity, density, and viscosity profiles plus wall shear
     - :doc:`../models/compressible_velocity_transformations`
   * - Construct a smooth turbulent ZPG mean profile
     - ``compressible_turbulent_boundary_layer_profile``
     - Edge state, :math:`\delta_{99}`, and wall shear
     - :doc:`../models/compressible_velocity_transformations`
   * - Estimate direct drag on one immersed protrusion
     - ``protrusion_drag``
     - Protrusion geometry, edge state, and either boundary-layer thickness or
       supplied profiles
     - :doc:`../models/protrusion_drag`

Estimate a flat-plate boundary layer
------------------------------------

Select the boundary-layer state explicitly. This laminar example uses SI
edge conditions and returns the one-sided drag per unit width:

>>> from aerophysics import BoundaryLayerRegime, flat_plate_boundary_layer
>>> layer = flat_plate_boundary_layer(
...     1.0,
...     edge_velocity=10.0,
...     edge_density=1.0,
...     edge_dynamic_viscosity=1e-5,
...     regime=BoundaryLayerRegime.LAMINAR,
... )
>>> round(layer.reynolds_number, 1)
1000000.0
>>> round(layer.boundary_layer_thickness, 6)
0.00491
>>> round(layer.drag_per_unit_width, 6)
0.0664

Use ``LAMINAR`` or ``TURBULENT`` when the state is known. ``TRANSITIONAL``
requires a caller-supplied transition Reynolds number; the API does not
predict natural transition. Downstream of it, the turbulent part continues
the fully turbulent plate with the same momentum thickness, so drag is
continuous while local friction jumps. Returned momentum thickness is
continuous for the log-wake and wall-temperature models; see the model page
for the exceptions. Compressibility
corrections are opt-in. See
:doc:`../models/flat_plate_boundary_layer` for the available correlations,
required thermodynamic inputs, and applicability ranges.

Construct or transform a mean profile
-------------------------------------

:func:`aerophysics.boundary_layer_profile.compressible_turbulent_boundary_layer_profile`
predicts a smooth-wall, zero-pressure-gradient mean profile from supplied edge
conditions, :math:`\delta_{99}`, and wall shear stress. It combines Spalding's
wall law and the Coles wake with either a Van Driest or Volpiani inverse
transformation. GRA is the default temperature--velocity closure; Walz is
selectable. The result includes velocity and thermodynamic arrays plus
profile-integrated displacement and momentum thicknesses.

The profile API accepts thickness and shear explicitly. For compressible
thicknesses, select ``VAN_DRIEST_II_WALL_TEMPERATURE`` in the flat-plate API;
``ECKERT`` and ``VAN_DRIEST_II`` do not supply turbulent thicknesses.
A typical sequence is:

#. Evaluate ``flat_plate_boundary_layer`` at the station of interest.
#. Create a wall-normal grid from zero through the returned
   ``boundary_layer_thickness``.
#. Pass that thickness and ``wall_shear_stress`` explicitly to
   ``compressible_turbulent_boundary_layer_profile``. To reproduce the new
   flat-plate closure, select ``VAN_DRIEST`` and ``WALZ`` and use the same
   wall temperature and thermophysical models.
#. Use the returned velocity, density, temperature, and viscosity arrays for
   analysis or as provided profiles for ``protrusion_drag``.

If profiles already come from measurements or simulation, use
``transform_compressible_velocity_profile`` instead. The model is limited to
fully turbulent, smooth, approximately zero-pressure-gradient flow. The
transformation equations, selection table, complete examples, and further
limitations are in
:doc:`../models/compressible_velocity_transformations`.

Estimate protrusion drag
------------------------

Scale a separately known free-stream drag coefficient by the dynamic pressure
available across a protrusion's frontal area:

>>> from aerophysics import protrusion_drag
>>> drag = protrusion_drag(
...     1.1,
...     height=0.02,
...     frontal_width=0.01,
...     edge_velocity=100.0,
...     edge_density=1.225,
...     boundary_layer_thickness=0.1,
... )
>>> round(drag.shielding_factor, 6)
0.491073
>>> round(drag.direct_drag, 6)
0.661721

The default uses a one-seventh-power turbulent velocity profile. For a
measured or computed profile, pass ``profile_height``, ``profile_velocity``,
and, when available, ``profile_density``. The profile arrays returned by
``compressible_turbulent_boundary_layer_profile`` can be supplied directly.
See :doc:`../models/protrusion_drag` for the integral definition, optional
compressible approximation, and excluded three-dimensional effects.

Follow a boundary layer through a pressure change
-------------------------------------------------

``flat_plate_boundary_layer`` evaluates each state independently.
``marching_boundary_layer`` instead marches from the leading edge, so each
station inherits the upstream history. Supply the edge pressure along the wall,
starting at the leading edge:

>>> import numpy as np
>>> from aerophysics import (
...     BoundaryLayerEdge,
...     BoundaryLayerRegime,
...     marching_boundary_layer,
... )
>>> distance = np.linspace(0.0, 1.0, 51)
>>> pressure = 2.0e4 * (1.0 - 0.2 * distance)
>>> edge = BoundaryLayerEdge.from_pressure(
...     distance,
...     pressure,
...     total_pressure=1.6e5,
...     total_temperature=500.0,
... )
>>> layer = marching_boundary_layer(
...     edge,
...     regime=BoundaryLayerRegime.LAMINAR,
...     wall_temperature=300.0,
... )
>>> layer.distance.size
50
>>> bool(np.all(layer.pressure_gradient_parameter > 0.0))
True

Use ``regime=BoundaryLayerRegime.TURBULENT`` for SA-noft2. The defaults are
``turbulent_prandtl_number=0.9`` and ``edge_sa_ratio=3.0`` at the first edge
station; dimensional ``nu_tilde`` is conserved along the outer streamline.
``edge_sa_ratio=0`` gives the laminar limit. Both options are validated in every
regime, and only the turbulent regime uses them. The result's read-only
``nu_tilde_ratio`` array has shape ``(N_s, N_eta)``; laminar results use None.
Omit ``wall_temperature`` for an adiabatic wall.

SA-noft2 develops its own laminar start and transition near the leading edge.
This is not a transition prediction. Accuracy against experiment is not
established; agreement with CFD using the same model does not establish it.

An adverse pressure gradient can separate the boundary layer. The march then
raises ``BoundaryLayerSeparationError``, which holds ``separation_distance`` and
the ``result`` at the stations upstream.

Thin axisymmetric boundary layers
---------------------------------

Supply ``wall_radius`` in metres at the edge nodes. Here a 20 degree cone
with constant edge conditions demonstrates the laminar Mangler relation:

>>> cone_edge = BoundaryLayerEdge.from_pressure(
...     [0.0, 1.0], [20000.0, 20000.0],
...     total_pressure=160000.0, total_temperature=500.0,
... )
>>> plate = marching_boundary_layer(cone_edge, regime=BoundaryLayerRegime.LAMINAR)
>>> cone = marching_boundary_layer(
...     cone_edge, regime=BoundaryLayerRegime.LAMINAR,
...     wall_radius=cone_edge.distance * np.sin(np.radians(20.0)),
... )
>>> round(float(cone.local_skin_friction_coefficient[0]
...             / plate.local_skin_friction_coefficient[0]), 6)
1.732051

The ratio is sqrt(3) at equal wall distance; the cone thickness is divided by
sqrt(3). This is an analytical verification example. Radius uses piecewise
linear interpolation, so curved bodies need sufficiently dense nodes.
``wall_radius=None`` keeps the planar calculation; a positive scalar describes
a constant-radius wall beginning at a sharp edge. Only the first array radius
may be zero. Results exclude the singular tip.

Transverse curvature is neglected. A warning reports delta99/r_w above 0.1,
an uncalibrated indicator; near-tip stations commonly trigger it. Longitudinal
curvature and normal pressure gradients are also neglected. Use SI units and
wall arc length, local-edge Cf, and positive-into-wall heat flux.

A cone-to-cylinder worked flow example is not supplied: downstream of its
shoulder, edge conditions require external Euler/CFD or measured data that are
not available here. ``flow_sequence`` cannot provide that distribution.
The synthetic shoulder verification is only a numerical stress test.

.. _rm3633-example:

Worked measured body: R&M 3633
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:ref:`Winter, Rotta and Smith <ref-winter-rotta-smith-1968>` measured a waisted
body of revolution, with a 20 degree conical nose, at M_inf=1.700 and
Re_l=10.08e6. Length is 60 in (1.524 m). The comparison below is generated,
not executed as a doctest. The generator contains the five Fig. 1 polynomials
for R/l versus X/l, with the final segment cubic, and checks all radius
junctions to 1e-5. Wall distance is the accurate quadrature of
``sqrt(1+(dR/dX)^2)`` using analytic polynomial derivatives.

The Table 1 phi=0 values are transcribed without alteration in
``tests/reference_data/viscous_flow/rm3633_m17.csv``. The total-pressure ratios
0.9926--0.9962 calculated from Cp and M_delta check transcription consistency.
Static edge
pressure uses ``p/p_inf = 1 + 0.7 M_inf**2 Cp`` and the edge total pressure is
the freestream total pressure times the Taylor--Maccoll loss for a 20 degree
cone at M1.7 (ratio 0.9957049344). At the added tip, pressure equals the X=3 in
value. Measured Cp is unsmoothed. PCHIP of ln(p) in wall distance supplies a
dense grid at approximately 0.25 in axial spacing, with measured stations and
geometry junctions included; radii on those nodes use the body polynomials.
The march ends at the last measured pressure station, X=59 in.

Total temperature and wall thermal condition are not stated. We assume
T0=310 K and an adiabatic wall, consistent with the paper's zero-heat-transfer
analysis, and scale freestream pressure to the measured Reynolds number.
The table includes 290 K and 330 K sensitivities and a separate substitution
of Taylor--Maccoll cone pressure at X<=6 in. ``edge_sa_ratio=3``.

The 0.005 in ballotini trip was at 1.5 in from the nose and was inadequate at
Mach numbers of 1.4 and above. X=3 in is transitional (Cf=0.00140 versus
0.00393 at 6 in); comparisons use X>=12 in only. Cf is based on local dynamic
pressure. SA-noft2 develops its own transition; prescribed transition is out
of scope. The source gives Cf accuracy +/-0.0002 and Cp accuracy +/-0.005.
A missing Cf at 15 in is not interpolated or treated as zero.

The calculated delta99/r_w is large near the waist, consistent with the
paper's warning that the layer is not small relative to radius there. Wall Cp
may also differ from edge Cp through normal pressure gradients. Thus the
comparison exposes limits of the thin-layer model; it is not an accuracy
certification or a fitted correlation.

.. include:: ../_generated/marching_axisymmetric_verification.rst

Check applicability before using the result
-------------------------------------------

These APIs are engineering models, and none of them is coupled to the outer
flow. The marching solver takes the edge state as an input and supports
laminar and SA-noft2 flow.
The flat-plate correlations exclude roughness, streamwise pressure gradients,
separation, suction, blowing, and natural-transition prediction. The profile
construction assumes a smooth, fully turbulent, approximately ZPG boundary
layer. The protrusion model integrates an undisturbed profile and does not
resolve wall interference, vortices, shocks, or interactions among multiple
elements.

Incompressible turbulent thicknesses
------------------------------------

With ``compressible=False``, turbulent thicknesses default to the
one-fifth-power correlation. Select ``TurbulentModel.SCHLICHTING_LOG_WAKE``
for the same friction with a momentum thickness consistent with the returned
drag and a log-wake profile delta99:

>>> from aerophysics import TurbulentModel
>>> layer = flat_plate_boundary_layer(
...     1.0,
...     edge_velocity=100.0,
...     edge_density=1.0,
...     edge_dynamic_viscosity=1e-5,
...     regime=BoundaryLayerRegime.TURBULENT,
...     turbulent_model=TurbulentModel.SCHLICHTING_LOG_WAKE,
... )
>>> round(layer.boundary_layer_thickness, 6)
0.014564
>>> round(layer.momentum_thickness, 6)
0.001502
>>> round(layer.distance * layer.average_skin_friction_coefficient / 2, 6)
0.001502

The one-fifth-power delta99 is about 40% thick at :math:`Re_x=10^6` and
becomes too thin above about :math:`5\times10^7`; it warns above
:math:`10^7`. See
:ref:`flat-plate-log-wake-thickness` for the closure, its accuracy, and the
points where it returns ``NaN``.

Compressible turbulent thicknesses
----------------------------------

Use ``compressible=True`` and
``turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE`` for the
wall-temperature-corrected Schlichting/VD2 friction pair and Van Driest/Walz
profile thicknesses. Momentum thickness follows the accumulated friction
integral; displacement thickness is integrated only through delta99.
Supply local edge conditions. See :doc:`../models/flat_plate_boundary_layer`
for the external validation envelope and limitations.

``ECKERT`` and ``VAN_DRIEST_II`` return ``NaN`` for all three turbulent
thicknesses, also after a prescribed transition. Their friction and drag remain available.


Compressible laminar thicknesses
--------------------------------

With ``compressible=True``, laminar points solve the perfect-gas similarity equations
for friction and thicknesses. Momentum thickness obeys ``theta = x * cf_local``
and average friction is twice local friction, on the edge dynamic-pressure
basis. Omitted wall temperature means zero wall heat flux; recovery temperature
is the solved adiabatic temperature even for a specified isothermal wall.
``effective_reynolds_number`` is the external Reynolds number at laminar points.
These rules also apply before transition, and to the laminar momentum
thickness that starts the turbulent part after transition. ``compressible=False`` retains the Blasius delta99 coefficient 4.9100.
See :doc:`../models/flat_plate_boundary_layer` for equations, numerical accuracy,
and the constant-cp, constant-Pr, zero-pressure-gradient assumptions.

An ``ApplicabilityWarning`` reports laminar viscous interaction when
``chi_bar = Mach**3 * sqrt(Cw / Re_x) > 1``. Leading-edge induced pressure
then invalidates the zero-pressure-gradient assumption. The warning excludes
``compressible=False`` and turbulent points; it does not correct the edge state.
