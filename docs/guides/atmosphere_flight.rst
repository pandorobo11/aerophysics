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

GUI altitude coordinate
-----------------------

The atmosphere/flight page accepts either geometric altitude ``h`` (the
default) or geopotential altitude ``H``. The selected coordinate applies to
the altitude input, altitude sweep bounds, and plot axes. Changing the
selection reinterprets the entered numbers in the selected coordinate.
Both altitudes remain available in the result table and CSV.

The existing conversion is ``H = R h / (R + h)``, with
``R = 6,356,766 m``; its inverse is ``h = R H / (R - H)``. Altitudes are
positive upwards from sea level. The valid range remains geometric
-5,000 to 86,000 m (approximately -5,003.936 to 84,852.046 geopotential m).
Values outside that physical range are rejected after conversion. The
standard-atmosphere model and Python API are unchanged.

Saved settings store the resolved geometric altitude in
``inputs_si.geometric_altitude`` and the selection in
``models.altitude_basis``. Sweep bounds use the selected coordinate in metres.
Older settings without ``altitude_basis`` use geometric altitude.
The flow-sequence page also offers this selection for standard-atmosphere
initialization; its ``inputs_si.initial.altitude`` is in the coordinate stored
in ``inputs_si.initial.altitude_basis`` (geometric when omitted).

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
