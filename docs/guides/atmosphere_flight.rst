Atmosphere and flight conditions
================================

Use :func:`~aerophysics.atmosphere.standard_atmosphere` when you need the
properties of the air at an altitude. Use
:class:`~aerophysics.flight.FlightCondition` when you also know how fast the
vehicle is flying and want quantities such as dynamic pressure or Reynolds
number.

Air properties at an altitude
-----------------------------

Pass the altitude in metres:

>>> from aerophysics import standard_atmosphere
>>> sea_level = standard_atmosphere(0.0)
>>> sea_level.temperature
288.15
>>> round(sea_level.pressure, 1)
101325.0
>>> round(sea_level.speed_of_sound, 3)
340.294

The result also holds density, viscosity, thermal conductivity, gravity, and
several other properties:

>>> state = standard_atmosphere(10_000.0)
>>> round(state.density, 4)
0.4135
>>> round(state.dynamic_viscosity, 9)
1.4577e-05

The altitude is the **geometric** altitude, which is the actual height above
sea level. The model covers -5 km to 86 km. It is a reference atmosphere and
does not describe the weather on a particular day. The full list of returned
properties and the equations are in :doc:`../models/standard_atmosphere`.

Some tables use **geopotential** altitude instead. Two helper functions convert
between the two:

>>> from aerophysics.atmosphere import geometric_to_geopotential
>>> round(geometric_to_geopotential(10_000.0), 1)
9984.3

Flight condition from a Mach number
-----------------------------------

Give the altitude and the Mach number. Add a characteristic length, such as a
chord or body length, if you want a Reynolds number:

>>> from aerophysics import FlightCondition
>>> condition = FlightCondition.from_mach(
...     0.0, 0.8, characteristic_length=1.5
... )
>>> round(condition.dynamic_pressure, 1)
45393.6
>>> round(condition.total_temperature, 4)
325.0332
>>> round(condition.reynolds_number)
27955292

Without a characteristic length, ``reynolds_number`` is ``None``. The Reynolds
number per metre is always available as ``reynolds_number_per_length``. The
atmospheric state used for the calculation is stored in
``condition.atmosphere``.

Flight condition from a speed
-----------------------------

When you know the speed in metres per second, use ``from_velocity``. The Mach
number is then calculated for you:

>>> condition = FlightCondition.from_velocity(10_000.0, 250.0)
>>> round(condition.mach, 4)
0.8346
>>> round(condition.dynamic_pressure, 1)
12922.2

The relations behind these values are in :doc:`../models/flight_conditions`.

Other units
-----------

Every calculation takes SI units. If your data are in feet, knots, or other
units, convert them first:

>>> from aerophysics.units import feet_to_meters, knots_to_meters_per_second
>>> feet_to_meters(10_000.0)
3048.0
>>> round(knots_to_meters_per_second(100.0), 6)
51.444444

Each conversion also has an inverse, such as ``meters_to_feet``. Angles are in
radians; use ``degrees_to_radians`` and ``radians_to_degrees``. All conversions
are listed in :doc:`../models/unit_conversions`.
