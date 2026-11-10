Compressible flow
=================

This guide covers the basic building blocks of high-speed flow: isentropic
flow, shock waves, and expansion waves. To combine several of them, see
:doc:`flow_sequence`.

Two rules apply throughout:

* Angles are in **radians**. Use
  :func:`~aerophysics.units.degrees_to_radians` or ``math.radians`` if you
  think in degrees.
* The gas is ``AIR`` with constant heat capacity unless you pass another model.
  See `Hot flows`_ below.

Which function to use
---------------------

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Situation
     - Function
     - Model page
   * - Flow that speeds up or slows down smoothly without losses, such as in a
       nozzle
     - :func:`~aerophysics.isentropic.isentropic_ratios` and related functions
     - :doc:`../models/isentropic_flow`
   * - A shock standing straight across the flow
     - :func:`~aerophysics.shocks.normal_shock`
     - :doc:`../models/normal_oblique_shocks`
   * - A shock attached to a wedge or a compression corner
     - :func:`~aerophysics.shocks.oblique_shock`
     - :doc:`../models/normal_oblique_shocks`
   * - A shock attached to the tip of a cone
     - :func:`~aerophysics.shocks.conical_shock`
     - :doc:`../models/conical_shocks`
   * - A bow shock standing ahead of a blunt nose
     - :func:`~aerophysics.detached_shock.shock_standoff_distance`
     - :doc:`../models/detached_shocks`
   * - Supersonic flow turning away from itself around a corner
     - :func:`~aerophysics.expansion.prandtl_meyer_expansion`
     - :doc:`../models/expansion_waves`

Isentropic flow
---------------

``isentropic_ratios`` relates the total (stagnation) state to the static state
at a Mach number. Ratios are total divided by static:

>>> from aerophysics.isentropic import isentropic_ratios
>>> ratios = isentropic_ratios(2.0)
>>> round(ratios.total_temperature_ratio, 6)
1.8
>>> round(ratios.total_pressure_ratio, 6)
7.824449

When you know the total pressure and temperature,
:func:`~aerophysics.isentropic.isentropic_state` returns the actual static
values, speed, and dynamic pressure:

>>> from aerophysics.isentropic import isentropic_state
>>> state = isentropic_state(2.0, total_temperature=300.0, total_pressure=1.0e6)
>>> round(state.static_temperature, 2)
166.67
>>> round(state.static_pressure)
127805
>>> round(state.velocity, 1)
517.6

Nozzles
~~~~~~~

``area_ratio`` gives the flow area divided by the throat area for a Mach
number. Going the other way, one area ratio matches two Mach numbers, one
subsonic and one supersonic, so you must say which you want:

>>> from aerophysics.isentropic import MachBranch, area_ratio, mach_from_area_ratio
>>> round(area_ratio(2.0), 4)
1.6875
>>> round(mach_from_area_ratio(1.6875, MachBranch.SUPERSONIC), 4)
2.0
>>> round(mach_from_area_ratio(1.6875, MachBranch.SUBSONIC), 4)
0.3722

The mass flow per unit throat area of a choked nozzle follows from the total
conditions:

>>> from aerophysics.isentropic import choked_mass_flux
>>> round(choked_mass_flux(1.0e6, 300.0), 1)
2333.3

The value is in kg/(m² s). Multiply it by the throat area to get the mass flow.
If you need ratios, areas, and mass flow together,
:func:`~aerophysics.isentropic.isentropic_analysis` returns all of them in one
call.

Normal shocks
-------------

Shock results are downstream divided by upstream:

>>> from aerophysics import normal_shock
>>> shock = normal_shock(2.0)
>>> round(shock.downstream_mach, 4)
0.5774
>>> round(shock.static_pressure_ratio, 4)
4.5
>>> round(shock.total_pressure_ratio, 4)
0.7209

A pitot tube in supersonic flow measures the total pressure behind a normal
shock. :func:`~aerophysics.shocks.supersonic_pitot_pressure_ratio` gives that
pressure divided by the upstream static pressure:

>>> from aerophysics.shocks import supersonic_pitot_pressure_ratio
>>> round(supersonic_pitot_pressure_ratio(2.0), 4)
5.6404

Oblique shocks
--------------

Give the upstream Mach number and the angle through which the flow is turned:

>>> from math import degrees, radians
>>> from aerophysics import oblique_shock
>>> shock = oblique_shock(2.0, radians(10.0))
>>> round(degrees(shock.shock_angle), 3)
39.314
>>> round(shock.downstream_mach, 4)
1.6405

Two solutions exist for each turn angle. The weak one is the default and is the
one normally seen on a wedge. Pass ``branch=ShockBranch.STRONG`` for the other.

A shock can stay attached only up to a maximum turn angle. Beyond it, the
function raises :class:`~aerophysics.exceptions.NoAttachedShockError`; it never
quietly switches to another kind of shock. You can ask for the limit first:

>>> from aerophysics.shocks import maximum_attached_deflection
>>> round(degrees(maximum_attached_deflection(2.0).deflection_angle), 2)
22.97

Conical shocks
--------------

A cone disturbs the flow less than a wedge of the same angle. The result gives
the conditions on the cone surface, because the flow keeps changing between the
shock and the surface:

>>> from aerophysics import conical_shock
>>> cone = conical_shock(2.0, radians(10.0))
>>> round(degrees(cone.shock_angle), 2)
31.21
>>> round(cone.surface_mach, 4)
1.834
>>> round(cone.surface_pressure_ratio, 4)
1.2925

The cone must point straight into the flow. Like a wedge, it has a maximum
angle for an attached shock, given by
:func:`~aerophysics.shocks.maximum_attached_cone_angle`.

Bow shocks on blunt noses
-------------------------

Ahead of a blunt nose the shock stands off from the surface. The standoff
distance comes from a correlation and needs the Mach number, the nose radius in
metres, and the shape:

>>> from aerophysics import DetachedShockGeometry, shock_standoff_distance
>>> bow = shock_standoff_distance(
...     6.0, 0.1, geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE
... )
>>> round(bow.standoff_distance, 5)
0.01565

Use ``CYLINDRICAL_NOSE_2D`` for a two-dimensional cylinder.
:func:`~aerophysics.detached_shock.billig_shock_shape` returns the coordinates
of the whole shock curve, and
:func:`~aerophysics.detached_shock.compare_standoff_distances` compares two
standoff correlations for a sphere. These are quick estimates, not solutions of
the flow behind the shock.

Expansion waves
---------------

When supersonic flow turns away from itself, it speeds up through an expansion
fan:

>>> from aerophysics import prandtl_meyer_expansion
>>> expansion = prandtl_meyer_expansion(2.0, radians(10.0))
>>> round(expansion.downstream_mach, 6)
2.384887
>>> round(expansion.static_pressure_ratio, 4)
0.548

The turn angle must be zero or positive. For a turn into the flow, use
``oblique_shock`` instead.

Hot flows
---------

A strong shock heats the air enough that its heat capacity changes. To include
this, pass a temperature-dependent gas model and the upstream static
temperature:

>>> from aerophysics import AIR_NASA9
>>> hot = normal_shock(6.0, AIR_NASA9, upstream_temperature=300.0)
>>> round(hot.static_temperature_ratio, 2)
7.09
>>> round(normal_shock(6.0).static_temperature_ratio, 2)
7.94

The shock and expansion functions accept ``AIR``, ``AIR_NASA7``, ``AIR_NASA9``,
and ``AIR_HARMONIC_OSCILLATOR``. If any state in the calculation leaves the
temperature range of the gas model, the function raises
:class:`~aerophysics.exceptions.ModelRangeError`. See :doc:`gas_properties` for
the models and their ranges.

Next steps
----------

* :doc:`flow_sequence` chains shocks and expansions and keeps track of the
  state after each one.
* :doc:`../verification/compressible_flow` compares these functions with
  published tables.
