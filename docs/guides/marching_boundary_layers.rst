Boundary layers with a changing outer flow
==========================================

The flat-plate functions in :doc:`boundary_layers` assume that the flow outside
the boundary layer is the same all along the surface. On a real shape it is
not: the outer flow speeds up where the pressure falls and slows down where the
pressure rises. A boundary layer remembers this history. Falling pressure thins
it and raises the friction. Rising pressure thickens it and can make it
separate from the wall.

:func:`~aerophysics.marching_boundary_layer.marching_boundary_layer` follows
this history. It starts at a sharp leading edge and works downstream, one
position at a time. You supply the pressure along the wall. The function
returns the friction, heat flux, thicknesses, and profiles at every position.

The theory is in :doc:`../models/marching_boundary_layer`.

A first example
---------------

Two steps are needed. First describe the outer flow with
:meth:`BoundaryLayerEdge.from_pressure
<aerophysics.marching_boundary_layer.BoundaryLayerEdge.from_pressure>`. Then
run the calculation.

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
>>> round(float(layer.local_skin_friction_coefficient[-1]), 6)
0.000249
>>> round(float(layer.wall_heat_flux[-1]))
3231

Here the pressure falls by 20% over one metre, so the outer flow accelerates.

Describing the outer flow
-------------------------

``BoundaryLayerEdge.from_pressure`` needs:

* ``distance``: positions along the wall in metres, measured along the surface.
  The first value must be the leading edge, at zero.
* ``pressure``: the static pressure at each position, in Pa.
* ``total_pressure`` and ``total_temperature``: the total conditions of the
  outer flow, in Pa and K. The function works out the velocity and temperature
  at each position from these.

If the outer flow has passed through a shock before reaching the wall, use the
total pressure **behind** the shock. The total temperature does not change
across a shock.

Give a smooth pressure distribution with enough points to follow its shape. The
calculation depends on the slope of the pressure, so scatter in measured
pressures shows up directly in the result.

The wall
--------

* ``wall_temperature`` in K holds the wall at that temperature. It can be one
  number or one value per position.
* Leave it out for an insulated wall. The calculated temperature of the wall is
  then returned in ``layer.wall_temperature`` and the heat flux is zero.

Laminar and turbulent flow
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - ``regime``
     - Meaning
   * - ``LAMINAR``
     - Laminar everywhere.
   * - ``TURBULENT``
     - Turbulent, using the Spalart--Allmaras turbulence model (the "SA-noft2"
       form widely used in CFD).
   * - ``TRANSITIONAL``
     - Laminar up to a position you choose, then turbulent. Also pass
       ``transition_distance`` in metres from the leading edge.

>>> turbulent = marching_boundary_layer(
...     edge,
...     regime=BoundaryLayerRegime.TURBULENT,
...     wall_temperature=300.0,
... )
>>> round(float(turbulent.local_skin_friction_coefficient[-1]), 5)
0.00214

A turbulent calculation takes roughly a second for this case on a laptop; a
laminar one is several times faster.

For a transition at a known position::

    layer = marching_boundary_layer(
        edge,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_distance=0.3,
        wall_temperature=300.0,
    )

Keep three things in mind about transition:

* The package does not predict where transition happens. You choose the
  position.
* Just after the chosen position, the friction rises from the laminar to the
  turbulent level over some distance. That distance comes from the turbulence
  model and is **not** a prediction of the real transition length.
* With ``TURBULENT``, the model still begins with a short laminar-like stretch
  near the leading edge. Its length is likewise not a prediction.

The turbulent results agree with CFD codes that use the same turbulence model.
How well they match experiments has not been established in general.

Two further options, ``edge_sa_ratio`` and ``turbulent_prandtl_number``, adjust
the turbulence model. The defaults suit most cases; see the
:func:`API reference <aerophysics.marching_boundary_layer.marching_boundary_layer>`.

Reading the results
-------------------

The result holds one value per position, leaving out the leading edge itself,
where friction and heat flux are infinite:

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Field
     - Meaning
   * - ``distance``
     - Position along the wall [m]
   * - ``local_skin_friction_coefficient``, ``wall_shear_stress``
     - Friction, based on the local edge dynamic pressure
   * - ``wall_heat_flux``
     - Heat flux in W/m², positive into the wall
   * - ``wall_temperature``
     - Wall temperature [K]
   * - ``boundary_layer_thickness``, ``displacement_thickness``,
       ``momentum_thickness``
     - Thicknesses [m]
   * - ``pressure_gradient_parameter``
     - Positive where the outer flow accelerates, negative where it slows

Profiles through the layer are two-dimensional arrays with one row per
position: ``wall_distance`` in metres, ``velocity_ratio`` (velocity divided by
the edge velocity), and ``temperature_ratio``.

Separation
----------

If the pressure rises too steeply, the friction at the wall falls to zero and
the flow separates. The calculation cannot continue past that point and raises
:class:`~aerophysics.exceptions.BoundaryLayerSeparationError`. The exception
tells you where separation happened and keeps the results upstream of it:

>>> from aerophysics.exceptions import BoundaryLayerSeparationError
>>> rising = BoundaryLayerEdge.from_pressure(
...     distance,
...     2.0e4 * (1.0 + 2.0 * distance),
...     total_pressure=1.6e5,
...     total_temperature=500.0,
... )
>>> try:
...     marching_boundary_layer(
...         rising, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=300.0
...     )
... except BoundaryLayerSeparationError as error:
...     print(round(error.separation_distance, 2))
...     print(error.result.distance.size)
0.51
25

Bodies of revolution
--------------------

On a body of revolution such as a cone, the boundary layer spreads over a
growing circumference as it moves downstream, which makes it thinner. To
include this, pass ``wall_radius``: the distance from the axis to the wall at
each position, in metres.

This example compares a flat plate with a 20° cone that has the same outer
flow:

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

The friction on the cone is :math:`\sqrt{3}` times that on the plate, which is
the classical textbook result for laminar flow.

Notes on ``wall_radius``:

* Only the first radius may be zero, for a pointed tip.
* A single positive number describes a tube of constant radius.
* The radius is joined by straight lines between your positions, so give a
  curved body enough points.
* The method assumes that the boundary layer is thin compared with the radius.
  A warning appears when the thickness exceeds 10% of the radius. This always
  happens close to a pointed tip and can be ignored there.

Cone--cylinder: outer flow and boundary layer together
------------------------------------------------------

For a pointed body of revolution in supersonic flow, the package can also
calculate the outer flow, so you do not need a pressure distribution from
elsewhere. :func:`~aerophysics.characteristics.axisymmetric_body_flow` takes
the outline of the body and the free stream:

>>> from aerophysics import axisymmetric_body_flow
>>> outer = axisymmetric_body_flow(
...     [0.0, 0.3, 1.3],
...     [0.0, 0.05, 0.05],
...     mach=3.0,
...     pressure=10000.0,
...     temperature=220.0,
... )
>>> round(float(outer.mach[0]), 3)
2.731
>>> round(float(outer.pressure[0]))
15030

The two lists are the axial positions and radii of the corners of the outline,
in metres. This body is a cone 0.3 m long followed by a cylinder of radius
0.05 m. The result gives the pressure, temperature, Mach number, and flow
direction along the surface, and the position of the nose shock.

``boundary_layer_inputs()`` converts the result into what the boundary-layer
calculation needs:

>>> inputs = outer.boundary_layer_inputs()
>>> layer = marching_boundary_layer(
...     inputs.edge,
...     regime=BoundaryLayerRegime.LAMINAR,
...     wall_radius=inputs.wall_radius,
...     wall_temperature=300.0,
... )
>>> bool(layer.distance[-1] == outer.surface_distance[-1])
True

Always pass ``inputs.wall_radius`` together with ``inputs.edge``; the two
belong to the same set of positions.

Limits of this method:

* The body must be pointed, must face straight into a supersonic flow, and may
  only turn away from the flow at each corner (cone, then cylinder, then
  boattail). A blunt nose or a flare is rejected.
* At the sharp corner between cone and cylinder, the pressure drops almost
  instantly. A real boundary layer changes the outer flow there, which this
  one-way calculation ignores. Expect a warning about the steep pressure
  change, and do not rely on the results close to the corner. Farther
  downstream, at a hundred or so boundary-layer thicknesses, the effect of the
  corner treatment falls to about 0.1%.

The details are in :doc:`../models/axisymmetric_characteristics`.

Warnings and errors
-------------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Message
     - What to do
   * - ``BoundaryLayerSeparationError``
     - The flow separated. Use the upstream results in ``error.result``.
   * - ``BoundaryLayerConvergenceError``
     - The calculation could not be completed, most often in a very strong
       pressure rise at very high Reynolds number. Check the pressure input.
   * - Warning: the pressure changes too quickly
     - The pressure changes over a distance shorter than the boundary-layer
       thickness, as at a sharp corner. Results near that place are unreliable.
   * - Warning: thickness is large compared with the radius
     - The thin-layer assumption is doubtful for a body of revolution. Expected
       near a pointed tip.
   * - Warning: viscous interaction
     - At high Mach number and low Reynolds number near the leading edge, the
       boundary layer changes the pressure. The pressure you supplied may not
       include this.
   * - Warning: first grid point (``y+``)
     - At extremely high Reynolds number the turbulent calculation is too
       coarse near the wall. Treat friction and heat flux with caution.

What the calculation leaves out
-------------------------------

The method is for an attached boundary layer that starts at a sharp leading
edge or tip, in air with constant heat capacity. It does not include:

* any effect of the boundary layer on the outer flow;
* blunt noses, where the strong bow shock changes the flow near the wall;
* roughness, suction, or blowing;
* chemical reactions;
* flow after separation.

Comparisons with exact solutions, with CFD, and with one wind-tunnel
experiment are recorded in :doc:`../verification/viscous_flow`.
