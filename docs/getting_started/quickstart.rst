Python quickstart
=================

This example starts with the reference atmosphere, builds a consistent flight
condition, and then evaluates the atmosphere over an altitude array. All
inputs and outputs are SI unless an explicit conversion function is used.

Atmosphere and flight condition
-------------------------------

Evaluate sea level and inspect two returned properties:

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

Vectorized atmosphere profile
-----------------------------

Array-like inputs return float64 NumPy arrays with the broadcast input shape:

>>> import numpy as np
>>> altitudes = np.array([0.0, 5_000.0, 10_000.0])
>>> profile = standard_atmosphere(altitudes)
>>> profile.temperature.shape
(3,)
>>> np.round(profile.temperature, 3)
array([288.15 , 255.676, 223.252])
>>> bool(np.all(np.diff(profile.pressure) < 0.0))
True

Where to go next
----------------

* :doc:`../guides/atmosphere_flight` explains altitude, Mach, velocity,
  Reynolds number, and explicit unit conversion workflows.
* :doc:`conventions` explains model choice, ratio directions, radians,
  vectorization, warnings, and range errors.
* :doc:`../guides/compressible_flow` and :doc:`../guides/boundary_layers`
  provide task-oriented examples for the corresponding physics.
* :doc:`../api/index` lists the complete public API by subject.
