.. _boundary-layer-workflows:

Flat-plate boundary layers
==========================

A boundary layer is the thin layer of slowed air next to a surface. It causes
skin-friction drag and, at high speed, heats the surface. This guide covers a
smooth flat plate with a steady outer flow, plus two related calculations:
velocity profiles through the layer, and the drag of a small object standing in
it.

If the outer flow changes along the surface, for example because the pressure
rises or falls, use :doc:`marching_boundary_layers` instead.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Task
     - Function
   * - Thickness, friction, drag, and heat transfer on a flat plate
     - :func:`~aerophysics.boundary_layer.flat_plate_boundary_layer`
   * - Velocity and temperature through the layer at one position
     - :func:`~aerophysics.boundary_layer.flat_plate_boundary_layer_profile`
   * - Drag of a small protrusion inside the layer
     - :func:`~aerophysics.protrusion.protrusion_drag`

Thickness, friction, and drag
-----------------------------

Give the distance from the leading edge and the flow just outside the boundary
layer, called the **edge** flow. On a flat plate at zero incidence the edge
flow is simply the free stream. You must also say whether the layer is laminar
or turbulent:

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
>>> round(layer.local_skin_friction_coefficient, 6)
0.000664
>>> round(layer.drag_per_unit_width, 6)
0.0664

``boundary_layer_thickness`` is the height at which the velocity reaches 99% of
the edge velocity, in metres. ``drag_per_unit_width`` is the friction drag of
one side of the plate from the leading edge to ``distance``, in newtons per
metre of width. The result also holds the displacement and momentum
thicknesses, the average friction coefficient, and the wall shear stress.

The edge flow can come straight from a flight condition:

>>> from aerophysics import FlightCondition
>>> flight = FlightCondition.from_mach(10_000.0, 0.8)
>>> air = flight.atmosphere
>>> layer = flat_plate_boundary_layer(
...     2.0,
...     edge_velocity=flight.velocity,
...     edge_density=air.density,
...     edge_dynamic_viscosity=air.dynamic_viscosity,
...     regime=BoundaryLayerRegime.TRANSITIONAL,
...     transition_reynolds=5e5,
... )
>>> round(layer.drag_per_unit_width, 2)
66.42

Laminar, turbulent, or both
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``LAMINAR`` and ``TURBULENT`` treat the whole plate as one or the other.
* ``TRANSITIONAL`` is laminar up to the Reynolds number you give in
  ``transition_reynolds`` and turbulent after it.

The package does not predict where transition happens; that depends on
roughness, noise, and other things outside these models. You choose it. A value
near :math:`5\times10^5` is a common assumption for a smooth plate in a quiet
flow.

Thickness at high Reynolds number
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default turbulent thickness formula is simple but loses accuracy above a
Reynolds number of about :math:`10^7`, where it emits a warning. For those
cases select the log-wake model, which gives the same friction and drag with a
better thickness:

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

High-speed flow and heat transfer
---------------------------------

At high Mach number, friction heats the air inside the boundary layer. This
changes its density and viscosity, and so the friction itself. To include the
effect, set ``compressible=True`` and add the Mach number and temperature at
the edge. For turbulent flow, also choose a compressible model:

>>> flight = FlightCondition.from_mach(20_000.0, 3.0)
>>> air = flight.atmosphere
>>> layer = flat_plate_boundary_layer(
...     1.0,
...     edge_velocity=flight.velocity,
...     edge_density=air.density,
...     edge_dynamic_viscosity=air.dynamic_viscosity,
...     regime=BoundaryLayerRegime.TURBULENT,
...     compressible=True,
...     turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
...     mach=3.0,
...     edge_temperature=air.temperature,
...     wall_temperature=400.0,
... )
>>> round(layer.local_skin_friction_coefficient, 6)
0.001904
>>> round(layer.recovery_temperature, 1)
566.2
>>> round(layer.wall_heat_flux)
14260

``recovery_temperature`` is the temperature the wall would reach if no heat
flowed through it. Here the wall is held at 400 K, which is cooler, so heat
flows into the wall: ``wall_heat_flux`` is positive, in W/m². The result also
holds the Stanton number and the heat-transfer coefficient.

If you leave out ``wall_temperature``, the wall is taken to be insulated. It
then sits at the recovery temperature and the heat flux is zero.

.. list-table:: Turbulent models for ``compressible=True``
   :header-rows: 1
   :widths: 40 60

   * - Model
     - Notes
   * - ``VAN_DRIEST_II_WALL_TEMPERATURE``
     - Recommended. Accounts for wall temperature and returns friction, heat
       transfer, and all three thicknesses.
   * - ``VAN_DRIEST_II``
     - Friction, drag, and heat transfer only. The thicknesses are ``NaN``.
   * - ``ECKERT``
     - A simpler reference-temperature method. Friction, drag, and heat
       transfer only. The thicknesses are ``NaN``.

Laminar flow needs no model choice. With ``compressible=True`` it is solved
directly from the laminar boundary-layer equations.

The equations, the ranges over which each model was checked, and the remaining
options are in :doc:`../models/flat_plate_boundary_layer`.

Profiles through the layer
--------------------------

:func:`~aerophysics.boundary_layer.flat_plate_boundary_layer_profile` returns
the velocity, temperature, density, viscosity, and local Mach number between
the wall and the edge of the layer at one position. It takes the edge velocity,
density, and temperature, and works out the Mach number and viscosity itself:

>>> from aerophysics import flat_plate_boundary_layer_profile
>>> profile = flat_plate_boundary_layer_profile(
...     1.0,
...     flight.velocity,
...     air.density,
...     air.temperature,
...     regime=BoundaryLayerRegime.TURBULENT,
...     wall_temperature=400.0,
... )
>>> profile.wall_distance.shape
(201,)
>>> round(float(profile.temperature[0]), 1)
400.0
>>> round(float(profile.temperature.max()), 1)
414.3

The temperature peaks inside the layer, above both the wall and the edge
values, because that is where friction heats the air. The matching thickness
and friction values are in ``profile.boundary_layer``.

Two more specialised functions are described in
:doc:`../models/compressible_velocity_transformations`:

* :func:`~aerophysics.boundary_layer_profile.compressible_turbulent_boundary_layer_profile`
  builds a turbulent profile from a thickness and wall shear stress that you
  supply, for example from a measurement.
* :func:`~aerophysics.boundary_layer_profile.transform_compressible_velocity_profile`
  converts a measured or computed high-speed profile into the form used to
  compare it with low-speed data.

Drag of a small protrusion
--------------------------

A rivet head, a sensor, or a step that sits inside the boundary layer meets
slower air than the free stream, so its drag is lower than its free-stream drag
coefficient suggests. :func:`~aerophysics.protrusion.protrusion_drag` makes
this correction. Give the drag coefficient the object would have in a uniform
stream, its height and width, the edge flow, and the boundary-layer thickness:

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

``shielding_factor`` is the fraction of the free-stream dynamic pressure that
the object actually meets, and ``direct_drag`` is the resulting force in
newtons.

By default the velocity in the layer follows a one-seventh power law, a common
approximation for a turbulent layer. To use a profile of your own, pass
``profile_height``, ``profile_velocity``, and optionally ``profile_density``.
The arrays from the profile functions above fit directly:

>>> drag = protrusion_drag(
...     1.1,
...     height=0.005,
...     frontal_width=0.01,
...     edge_velocity=flight.velocity,
...     edge_density=air.density,
...     boundary_layer_thickness=profile.boundary_layer.boundary_layer_thickness,
...     profile_height=profile.wall_distance,
...     profile_velocity=profile.velocity,
...     profile_density=profile.density,
... )
>>> 0.0 < drag.shielding_factor < 1.0
True

The model only accounts for the slower air. It does not describe how the
object disturbs the flow around it. See :doc:`../models/protrusion_drag`.

What these models leave out
---------------------------

These are quick engineering estimates for a smooth flat plate. They do not
include surface roughness, a pressure that changes along the plate, flow
separation, or suction and blowing through the wall. A warning tells you when
an input lies outside the range that a formula was fitted to; see
:doc:`vectorization_errors`.

At very high Mach number and low Reynolds number, the boundary layer near the
leading edge becomes thick enough to change the pressure on the plate. The
laminar compressible calculation warns when this is likely. The result is then
a rough guide only.
