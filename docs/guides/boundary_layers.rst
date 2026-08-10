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
0.005
>>> round(layer.drag_per_unit_width, 6)
0.0664

Use ``LAMINAR`` or ``TURBULENT`` when the state is known. ``TRANSITIONAL``
requires a caller-supplied transition Reynolds number; the API does not
predict natural transition. Compressibility corrections are opt-in. See
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

This profile API is independent of ``flat_plate_boundary_layer``. The latter
can estimate :math:`\delta_{99}` and wall shear, but the caller passes them
explicitly. A typical sequence is therefore:

#. Evaluate ``flat_plate_boundary_layer`` at the station of interest.
#. Create a wall-normal grid from zero through the returned
   ``boundary_layer_thickness``.
#. Pass that thickness and ``wall_shear_stress`` explicitly to
   ``compressible_turbulent_boundary_layer_profile``.
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

Check applicability before using the result
--------------------------------------------

These APIs are engineering models, not a coupled boundary-layer solver. The
flat-plate correlations exclude roughness, streamwise pressure gradients,
separation, suction, blowing, and natural-transition prediction. The profile
construction assumes a smooth, fully turbulent, approximately ZPG boundary
layer. The protrusion model integrates an undisturbed profile and does not
resolve wall interference, vortices, shocks, or interactions among multiple
elements.
