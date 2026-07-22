Quick start
===========

Atmosphere and flight condition
-------------------------------

Evaluate sea-level standard atmosphere:

>>> from aerophysics import FlightCondition, standard_atmosphere
>>> sea_level = standard_atmosphere(0.0)
>>> sea_level.temperature
288.15
>>> round(sea_level.speed_of_sound, 3)
340.294

Create a Mach-defined flight condition. Altitude and characteristic length are
in metres:

>>> condition = FlightCondition.from_mach(0.0, 0.8, characteristic_length=1.5)
>>> round(condition.dynamic_pressure, 1)
45393.6
>>> round(condition.total_temperature, 4)
325.0332

Isentropic flow
---------------

State ratios use the total-to-static convention:

>>> from aerophysics.isentropic import isentropic_ratios
>>> ratios = isentropic_ratios(2.0)
>>> round(ratios.total_temperature_ratio, 6)
1.8
>>> round(ratios.total_pressure_ratio, 6)
7.824449

Unit conversion
---------------

Customary units are converted explicitly rather than inferred:

>>> from aerophysics.units import feet_to_meters, knots_to_meters_per_second
>>> feet_to_meters(10_000.0)
3048.0
>>> round(knots_to_meters_per_second(100.0), 6)
51.444444
