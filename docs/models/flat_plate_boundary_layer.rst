.. _flat-plate-boundary-layer-model:

Flat-plate boundary-layer model
===============================

These pages describe the models behind
:func:`~aerophysics.boundary_layer.flat_plate_boundary_layer` and
:func:`~aerophysics.boundary_layer.flat_plate_boundary_layer_profile`. This
page gives the background, the definitions shared by all the models, and a
guide to choosing one. The equations of each model are on the pages listed
under `Model pages`_. For worked examples, see
:doc:`../guides/boundary_layers`.

Physical background
-------------------

Air sticks to a solid surface: its velocity at the wall is zero (the no-slip
condition). The thin region in which the velocity rises from zero to the
outer-flow value :math:`U_e` is the **boundary layer**. Friction inside it
produces skin-friction drag, and it grows thicker with distance from the
leading edge.

A boundary layer starts **laminar**, with smooth, orderly layers. Farther
downstream, small disturbances grow and the flow becomes **turbulent**, with
strong mixing. Turbulent mixing brings fast fluid closer to the wall, so a
turbulent boundary layer has much higher skin friction and grows faster. The
change, called transition, typically happens between :math:`Re_x\approx
3\times10^5` and :math:`3\times10^6`, depending on surface finish, noise, and
free-stream turbulence. This package does not predict where it happens.

The flat plate with no pressure gradient is the simplest case and the
reference for many engineering estimates:

* The laminar case has an exact similarity solution by **Blasius**.
* Turbulent results are **empirical correlations** fitted to measurements.

At high Mach number, friction heats the gas near the wall. The hot gas has
lower density and higher viscosity, which changes skin friction. The
**compressibility corrections** account for this by evaluating the
incompressible formulas at suitably chosen reference conditions.

What is modelled
----------------

The model is one side of a smooth flat plate with a sharp leading edge. The
flow outside the boundary layer, called the edge flow, is the same all along
the plate, so there is no pressure gradient. Distance :math:`x` is measured
from the leading edge, and drag is reported per unit width of the plate.

The standard definitions and incompressible correlations are summarized by
:ref:`Schlichting and Gersten (2000) <ref-schlichting-gersten-2000>` and
:ref:`White (2006) <ref-white-2006>`.

Choosing a model
----------------

Three choices define a calculation:

* ``regime``: :class:`~aerophysics.boundary_layer.BoundaryLayerRegime`
  ``LAMINAR``, ``TURBULENT``, or ``TRANSITIONAL``. A transitional calculation
  needs a ``transition_reynolds`` value; the position of transition is never
  predicted.
* ``compressible``: ``False`` (the default) for low-speed flow, or ``True`` to
  include the heating of the gas inside the boundary layer. With ``True``,
  ``mach`` and ``edge_temperature`` are required, and ``wall_temperature`` may
  be given. Temperatures are in K.
* ``turbulent_model``: a member of
  :class:`~aerophysics.boundary_layer.TurbulentModel`, for turbulent or
  transitional flow only.

.. list-table:: Valid model selections
   :header-rows: 1
   :widths: 22 16 62

   * - Regime
     - ``compressible``
     - ``turbulent_model``
   * - ``LAMINAR``
     - False or True
     - Must be None: Blasius for False, compressible similarity for True
   * - ``TURBULENT`` / ``TRANSITIONAL``
     - False
     - None defaults to SCHLICHTING; SCHLICHTING or SCHLICHTING_LOG_WAKE
   * - ``TURBULENT`` / ``TRANSITIONAL``
     - True
     - Required: ECKERT, VAN_DRIEST_II,
       or VAN_DRIEST_II_WALL_TEMPERATURE

Any other combination raises ``ValueError``. A choice that does not apply is
never silently ignored.

.. list-table:: What each turbulent model provides
   :header-rows: 1
   :widths: 30 36 34

   * - Member
     - Skin friction
     - Thicknesses
   * - ``SCHLICHTING``
     - Schlichting's logarithmic formula
     - One-fifth-power formula; best below :math:`Re_x=10^7`
   * - ``SCHLICHTING_LOG_WAKE``
     - Same as ``SCHLICHTING``
     - From a log-wake velocity profile; consistent with the drag
   * - ``ECKERT``
     - Schlichting at Eckert's reference temperature
     - Not available (``NaN``)
   * - ``VAN_DRIEST_II``
     - Van Driest II transformation
     - Not available (``NaN``)
   * - ``VAN_DRIEST_II_WALL_TEMPERATURE``
     - Van Driest II with a wall-temperature correction
     - From a compressible velocity profile; consistent with the drag

For low-speed flow, ``SCHLICHTING_LOG_WAKE`` is recommended whenever a
thickness matters. For high-speed flow, ``VAN_DRIEST_II_WALL_TEMPERATURE`` is
recommended.

Model pages
-----------

.. toctree::
   :maxdepth: 1

   flat_plate_incompressible
   flat_plate_transition
   flat_plate_compressible
   flat_plate_heat_transfer

Comparisons of these models with CFD, DNS, and measurements are collected in
:doc:`../verification/flat_plate_boundary_layer`.

Common definitions
------------------

Reynolds number and dynamic pressure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   Re_x=\frac{\rho_eU_ex}{\mu_e},
   \qquad
   q_e=\frac{1}{2}\rho_eU_e^2.

.. list-table:: Common flat-plate symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`x`
     - ``distance``
     - Distance from the leading edge
     - m
   * - :math:`U_e`
     - ``edge_velocity``
     - Boundary-layer edge velocity
     - m/s
   * - :math:`\rho_e`
     - ``edge_density``
     - Edge density
     - kg/m³
   * - :math:`\mu_e`
     - ``edge_dynamic_viscosity``
     - Edge dynamic viscosity
     - Pa s
   * - :math:`Re_x`
     - ``reynolds_number``
     - Reynolds number based on :math:`x`
     - dimensionless
   * - :math:`q_e`
     - computed internally
     - Edge dynamic pressure
     - Pa

.. _flat-plate-skin-friction-coefficients:

Skin-friction coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~

The local coefficient relates wall shear stress to edge dynamic pressure:

.. math::

   C_{f,x}=\frac{\tau_w(x)}{q_e}.

The average coefficient from the leading edge to :math:`x` is

.. math::

   \bar C_f(x)
   =\frac{1}{x}\int_0^xC_{f,\xi}\,\mathrm{d}\xi
   =\frac{D'(x)}{q_ex},
   \qquad
   D'(x)=\frac{1}{2}\rho_eU_e^2x\bar C_f.

:math:`D'` is one-sided accumulated drag per unit width. These definitions
assume constant :math:`q_e` along the plate. They are flat-plate
*skin-friction coefficients*, not Darcy or Fanning pipe-friction factors.

.. list-table:: Skin-friction outputs
   :header-rows: 1
   :widths: 16 32 36 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`C_{f,x}`
     - ``local_skin_friction_coefficient``
     - Local skin-friction coefficient
     - dimensionless
   * - :math:`\bar C_f`
     - ``average_skin_friction_coefficient``
     - Leading-edge-averaged coefficient
     - dimensionless
   * - :math:`\tau_w`
     - ``wall_shear_stress``
     - Local wall shear stress
     - Pa
   * - :math:`D'`
     - ``drag_per_unit_width``
     - Accumulated one-sided drag per unit width
     - N/m

Thicknesses
~~~~~~~~~~~

.. list-table:: Thickness symbols
   :header-rows: 1
   :widths: 17 30 37 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`\delta_{99}`
     - ``boundary_layer_thickness``
     - Conventional 99-percent thickness
     - m
   * - :math:`\delta^*`
     - ``displacement_thickness``
     - Displacement thickness
     - m
   * - :math:`\theta_m`
     - ``momentum_thickness``
     - Momentum thickness
     - m

Worked example
--------------

For :math:`x=1\ \mathrm{m}`, :math:`U_e=100\ \mathrm{m/s}`,
:math:`\rho_e=1.225\ \mathrm{kg/m^3}`, and
:math:`\mu_e=1.7894\times10^{-5}\ \mathrm{Pa\,s}`, the Schlichting model
gives :math:`Re_x=6.84587\times10^6`, :math:`C_{f,x}=0.002670314`,
:math:`\bar C_f=0.003193859`, :math:`\tau_w=16.3557\ \mathrm{Pa}`, and
:math:`D'=19.5624\ \mathrm{N/m}`.

>>> from aerophysics import BoundaryLayerRegime, flat_plate_boundary_layer
>>> layer = flat_plate_boundary_layer(
...     1.0,
...     edge_velocity=100.0,
...     edge_density=1.225,
...     edge_dynamic_viscosity=1.7894e-5,
...     regime=BoundaryLayerRegime.TURBULENT,
... )
>>> f"{layer.reynolds_number:.5e}"
'6.84587e+06'
>>> round(layer.local_skin_friction_coefficient, 9)
0.002670314
>>> round(layer.average_skin_friction_coefficient, 9)
0.003193859
>>> round(layer.wall_shear_stress, 4)
16.3557
>>> round(layer.drag_per_unit_width, 4)
19.5624

The result object is
:class:`aerophysics.boundary_layer.FlatPlateBoundaryLayerResult`.

Applicability
-------------

The flat-plate model excludes surface roughness, streamwise pressure
gradients, separation, suction, blowing, and natural-transition prediction.
The gas has constant heat capacity, and the compressible models take their
viscosity from the supplied dynamic-viscosity model. The thickness models are
engineering estimates. Agreement with a CFD turbulence model, reported in the
verification appendix, is not an independent physical validation.
