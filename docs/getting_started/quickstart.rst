Python quickstart
=================

This page walks through a few typical calculations so that you can see how the
package is used. Every input and output is in SI units (metres, kilograms,
seconds, kelvin, pascals), and angles are in radians.

Atmosphere and flight condition
-------------------------------

Ask for the standard atmosphere at sea level and read two of its properties:

>>> from aerophysics import FlightCondition, standard_atmosphere
>>> sea_level = standard_atmosphere(0.0)
>>> sea_level.temperature
288.15
>>> round(sea_level.speed_of_sound, 3)
340.294

Create a flight condition at 10 km and Mach 0.8. The characteristic length is
two metres, so the result also contains a dimensionless Reynolds number:

>>> condition = FlightCondition.from_mach(
...     geometric_altitude=10_000.0,
...     mach=0.8,
...     characteristic_length=2.0,
... )
>>> round(condition.velocity, 3)
239.625
>>> round(condition.dynamic_pressure, 1)
11872.0
>>> condition.reynolds_number is not None
True

Many values at once
-------------------

Pass an array instead of a single number to calculate many cases in one call.
The results are NumPy arrays of the same shape:

>>> import numpy as np
>>> altitudes = np.array([0.0, 5_000.0, 10_000.0])
>>> profile = standard_atmosphere(altitudes)
>>> profile.temperature.shape
(3,)
>>> np.round(profile.temperature, 3)
array([288.15 , 255.676, 223.252])
>>> bool(np.all(np.diff(profile.pressure) < 0.0))
True

A shock wave
------------

A wedge turns a Mach 2 flow through 10°. The functions take radians, so
convert the angle first:

>>> from math import degrees, radians
>>> from aerophysics import oblique_shock
>>> shock = oblique_shock(2.0, radians(10.0))
>>> round(degrees(shock.shock_angle), 2)
39.31
>>> round(shock.static_pressure_ratio, 3)
1.707

The shock stands at 39.3° to the flow, and the pressure behind it is 1.7 times
the pressure ahead of it.

Skin friction on a plate
------------------------

The flight condition from above supplies the flow over a flat plate. This
gives the thickness of a turbulent boundary layer 1 m from the leading edge,
and the friction drag of one side of the plate up to there, per metre of
width:

>>> from aerophysics import BoundaryLayerRegime, flat_plate_boundary_layer
>>> layer = flat_plate_boundary_layer(
...     1.0,
...     edge_velocity=condition.velocity,
...     edge_density=condition.atmosphere.density,
...     edge_dynamic_viscosity=condition.atmosphere.dynamic_viscosity,
...     regime=BoundaryLayerRegime.TURBULENT,
... )
>>> round(layer.boundary_layer_thickness, 4)
0.0159
>>> round(layer.drag_per_unit_width, 1)
38.0

How results are returned
------------------------

Every calculation returns a result object with named fields, such as
``shock.shock_angle`` above. In an interactive session, type the name of a
result to see all of its fields, or look the function up in
:doc:`../api/index`.

When a calculation cannot be done, the package raises an error that says why.
When it can be done but lies outside the range a formula was meant for, it
returns the value and emits a warning. See
:doc:`../guides/vectorization_errors`.

Where to go next
----------------

* :doc:`../guides/index` has a guide for each kind of calculation.
* :doc:`conventions` lists the symbols, the sign and ratio conventions, and
  tables for choosing a model.
* :doc:`../gui/index` describes the browser interface, if you prefer not to
  write Python.
