Gas properties and standard atmosphere
======================================

This page describes the constant-property gas and U.S. Standard Atmosphere
models. See :doc:`transport_properties` for viscosity and conductivity, or
:doc:`thermochemistry` when temperature-dependent heat capacity is required.

Calorically perfect gas
-----------------------

:class:`~aerophysics.gas.PerfectGas` holds the specific gas constant
:math:`R` and heat-capacity ratio :math:`\gamma` constant. The implemented
relations are

.. math::

   p=\rho RT,
   \qquad
   c_p=\frac{\gamma R}{\gamma-1},
   \qquad
   c_v=\frac{R}{\gamma-1},
   \qquad
   a=\sqrt{\gamma RT}.

``AIR`` uses :math:`R=8314.32/28.9644\ \mathrm{J/(kg\,K)}` and
:math:`\gamma=1.4`. All temperatures must be positive.
The constants follow :ref:`U.S. Standard Atmosphere 1976
<ref-us-standard-atmosphere-1976>`.

.. list-table:: Gas-property symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`p`
     - ``pressure``
     - Static pressure
     - Pa
   * - :math:`\rho`
     - ``density``
     - Density
     - kg/m³
   * - :math:`T`
     - ``temperature``
     - Absolute temperature
     - K
   * - :math:`R`
     - ``specific_gas_constant``
     - Specific gas constant
     - J/(kg K)
   * - :math:`\gamma`
     - ``heat_capacity_ratio``
     - Ratio :math:`c_p/c_v`
     - dimensionless
   * - :math:`c_p,c_v`
     - ``cp``, ``cv``
     - Specific heat capacities
     - J/(kg K)
   * - :math:`a`
     - ``speed_of_sound``
     - Ideal-gas speed of sound
     - m/s

Temperature-dependent properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`~aerophysics.thermochemistry.ThermallyPerfectGas` evaluates
temperature-dependent heat capacities from NASA-polynomial species data at
fixed composition. ``AIR_NASA7`` and ``AIR_NASA9`` cover 200--6000 K for
frozen dry air. See :doc:`thermochemistry` for their equations, reference
states, applicability, and examples.

These frozen-composition models include species heat-capacity variation but
not dissociation, chemical reactions, ionization, or thermodynamic
nonequilibrium.

U.S. Standard Atmosphere 1976
-----------------------------

The layer definitions, constants, and equations follow
:ref:`U.S. Standard Atmosphere 1976 <ref-us-standard-atmosphere-1976>`.

:func:`aerophysics.atmosphere.standard_atmosphere` accepts geometric altitude
:math:`h` and converts it to geopotential altitude :math:`H`:

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
^^^^^^^^^^^^^^^^^^^^^^^^^

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
