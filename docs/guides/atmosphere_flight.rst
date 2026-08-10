Atmosphere and flight-condition workflow
========================================

Use :func:`aerophysics.atmosphere.standard_atmosphere` when the atmospheric
state is the result you need. Use :class:`aerophysics.flight.FlightCondition`
when Mach or velocity must be combined with that state.

Atmospheric state
-----------------

>>> from aerophysics import standard_atmosphere
>>> sea_level = standard_atmosphere(0.0)
>>> sea_level.temperature
288.15
>>> round(sea_level.pressure, 1)
101325.0
>>> round(sea_level.speed_of_sound, 3)
340.294

The input is geometric altitude in metres. The implemented range is -5 to
86 km; this is a reference atmosphere, not a weather forecast. See
:doc:`../models/gas_and_atmosphere` for equations and all returned fields.

Mach-defined flight condition
-----------------------------

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

Use ``from_velocity`` instead when velocity, rather than Mach number, is the
independent input. Do not independently specify both. Details are in
:doc:`../models/flight_conditions`.

Explicit unit conversion
------------------------

Calculation APIs never infer customary units. Convert them at the package
boundary:

>>> from aerophysics.units import feet_to_meters, knots_to_meters_per_second
>>> feet_to_meters(10_000.0)
3048.0
>>> round(knots_to_meters_per_second(100.0), 6)
51.444444

All supported forward and inverse conversions are listed in
:doc:`../models/unit_conversions`.
