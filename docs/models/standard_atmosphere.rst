U.S. Standard Atmosphere 1976
=============================

:func:`aerophysics.atmosphere.standard_atmosphere` returns temperature,
pressure, density, speed of sound, gravity, viscosity, and thermal
conductivity at a given altitude. The layer definitions, constants, and
equations follow :ref:`U.S. Standard Atmosphere 1976
<ref-us-standard-atmosphere-1976>`. The air properties use the ``AIR`` gas
constants described in :doc:`gas_models` and the viscosity and conductivity
models in :doc:`transport_properties`.

Physical background
-------------------

A standard atmosphere is an agreed, idealized average of the real
atmosphere. It gives one temperature, pressure, and density for each
altitude, so that aircraft performance and test data can be compared on a
common basis. Actual weather can differ by tens of kelvin.

The model assumes dry air at rest that obeys the ideal-gas law and is in
hydrostatic balance: the pressure difference across a thin layer of air
supports that layer's weight,

.. math::

   \mathrm{d}p=-\rho g\,\mathrm{d}h=-\frac{pg}{RT}\,\mathrm{d}h.

Below 86 km the temperature profile is fixed as a series of layers in which
temperature either changes linearly with altitude or stays constant (for
example, falling by 6.5 K per km up to 11 km, then constant to 20 km).
Integrating the hydrostatic equation through each layer gives the pressure
formulas below.

Gravity weakens with height. To keep the integration simple, the standard
uses **geopotential altitude** :math:`H`: the height that would give the same
potential energy if gravity stayed at its sea-level value :math:`g_0`. With
:math:`g_0\,\mathrm{d}H=g\,\mathrm{d}h`, the hydrostatic equation becomes
:math:`\mathrm{d}p=-\rho g_0\,\mathrm{d}H`. Geopotential altitude is slightly
smaller than geometric altitude :math:`h`, the actual height above sea level;
at 10 km the difference is about 16 m.

Equations
---------

:func:`aerophysics.atmosphere.standard_atmosphere` accepts geometric altitude
:math:`h` and converts it to geopotential altitude :math:`H`. The conversions
are also available as :func:`~aerophysics.atmosphere.geometric_to_geopotential`
and :func:`~aerophysics.atmosphere.geopotential_to_geometric`:

.. math::

   H=\frac{r_0h}{r_0+h},
   \qquad
   h=\frac{r_0H}{r_0-H},

where :math:`r_0=6\,356\,766\ \mathrm{m}`. Within each atmosphere layer,

.. math::

   T=T_b+L_b(H-H_b).

For an isothermal layer (:math:`L_b=0`),

.. math::

   p=p_b\exp\left[-\frac{g_0(H-H_b)}{RT_b}\right].

For a nonzero lapse rate,

.. math::

   p=p_b\left(\frac{T_b}{T}\right)^{g_0/(RL_b)}.

The remaining properties are

.. math::

   \rho=\frac{p}{RT},
   \qquad
   g(h)=g_0\left(\frac{r_0}{r_0+h}\right)^2,
   \qquad
   \nu=\frac{\mu}{\rho},
   \qquad
   Pr=\frac{\mu c_p}{k},

where :math:`g_0=9.80665\ \mathrm{m/s^2}`.

.. list-table:: Atmosphere symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`h`
     - ``geometric_altitude``
     - Geometric altitude
     - m
   * - :math:`H`
     - ``geopotential_altitude``
     - Altitude used by layer hydrostatics
     - m
   * - :math:`H_b`
     - internal layer table
     - Active-layer base altitude
     - m
   * - :math:`T_b,p_b,L_b`
     - internal layer tables
     - Active-layer base state and lapse rate
     - K, Pa, K/m
   * - :math:`g`
     - ``gravity``
     - Local gravitational acceleration
     - m/s²
   * - :math:`\nu`
     - ``kinematic_viscosity``
     - Kinematic viscosity
     - m²/s

Applicability and example
-------------------------

The implemented geometric-altitude range is
:math:`-5\,000\le h\le86\,000\ \mathrm{m}`. Outside it,
:class:`~aerophysics.exceptions.ModelRangeError` is raised instead of
extrapolating. This is a reference atmosphere, not a weather forecast.

>>> from aerophysics import standard_atmosphere
>>> sea_level = standard_atmosphere(0.0)
>>> sea_level.temperature
288.15
>>> round(sea_level.pressure, 1)
101325.0
>>> round(sea_level.speed_of_sound, 3)
340.294

The returned :class:`~aerophysics.atmosphere.AtmosphereState` also contains
density, gravity, viscosity, conductivity, kinematic viscosity, and Prandtl
number.
