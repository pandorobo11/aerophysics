Gas properties
==============

This guide shows how to get the thermodynamic properties and the viscosity of
air, and how to choose a gas model for the flow calculations. The equations are
in :doc:`../models/gas_models` and :doc:`../models/transport_properties`.

Choosing a gas model
--------------------

A gas model tells a calculation how the heat capacity of the gas behaves. The
package provides five ready-made models for air:

.. list-table::
   :header-rows: 1
   :widths: 28 44 28

   * - Model
     - What it assumes
     - Temperature range
   * - ``AIR``
     - The heat capacity is constant (:math:`\gamma = 1.4`). This is the
       default everywhere and is accurate for air below about 800 K.
     - Any positive temperature
   * - ``AIR_NASA7``, ``AIR_NASA9``
     - The heat capacity changes with temperature, following NASA polynomial
       fits. The composition of the air stays fixed.
     - 200 K to 6000 K
   * - ``AIR_HARMONIC_OSCILLATOR``
     - The heat capacity changes with temperature, following a simple formula
       for molecular vibration.
     - 400 K to 2000 K
   * - ``AIR_BEATTIE_BRIDGEMAN``
     - As above, and the gas also departs from the ideal-gas law at high
       pressure.
     - About 39 K to 1222 K

Start with ``AIR``. Move to a temperature-dependent model when the flow becomes
hot enough for the heat capacity to change noticeably, for example behind a
strong shock. None of the models includes chemical reactions, so they do not
describe dissociated or ionized air.

Not every calculation accepts every model. The table in
:doc:`../getting_started/conventions` lists which combinations work.

Constant heat capacity
----------------------

``AIR`` is a :class:`~aerophysics.gas.PerfectGas`:

>>> from aerophysics import AIR
>>> AIR.heat_capacity_ratio
1.4
>>> round(AIR.cp, 2)
1004.69
>>> round(AIR.speed_of_sound(300.0), 2)
347.22

For another gas, create your own model from its gas constant in J/(kg K) and
its heat-capacity ratio. This example is nitrogen:

>>> from aerophysics import PerfectGas
>>> nitrogen = PerfectGas(specific_gas_constant=296.8, heat_capacity_ratio=1.4)
>>> round(nitrogen.speed_of_sound(300.0), 2)
353.07

Pass a gas model to a calculation through its ``gas`` argument, for example
``normal_shock(2.0, gas=nitrogen)``.

Heat capacity that changes with temperature
-------------------------------------------

For the temperature-dependent models, every property is a function of
temperature in kelvin:

>>> from aerophysics import AIR_NASA9
>>> round(AIR_NASA9.cp(300.0), 1)
1004.8
>>> round(AIR_NASA9.cp(2000.0), 1)
1250.3
>>> round(AIR_NASA9.heat_capacity_ratio(2000.0), 4)
1.298
>>> round(AIR_NASA9.speed_of_sound(2000.0), 1)
863.2

``sensible_enthalpy`` gives the enthalpy relative to a reference temperature,
which is 298.15 K unless you choose another:

>>> round(AIR_NASA9.sensible_enthalpy(2000.0) / 1e6, 3)
1.954

A temperature outside the model's range raises
:class:`~aerophysics.exceptions.ModelRangeError`:

>>> from aerophysics.exceptions import ModelRangeError
>>> try:
...     AIR_NASA9.cp(7000.0)
... except ModelRangeError as error:
...     print(error)
temperature must be within 200 K and 6000 K

``AIR_NASA7`` works the same way and gives almost the same values. The two use
different published fits of the same data.

High pressure
-------------

At high pressure, air is denser or thinner than the ideal-gas law predicts.
``AIR_BEATTIE_BRIDGEMAN`` takes both temperature and pressure and returns the
complete state:

>>> from aerophysics import AIR_BEATTIE_BRIDGEMAN
>>> state = AIR_BEATTIE_BRIDGEMAN.state(300.0, 2.0e7)
>>> round(state.density, 1)
226.9

The ideal-gas law would give 232.2 kg/m³ at this state. This model is accepted
only by the isentropic-flow functions.

Viscosity and thermal conductivity
----------------------------------

Three viscosity models for dry air are provided. Each one takes a temperature
in kelvin and returns the dynamic viscosity in Pa s:

>>> from aerophysics.transport import (
...     AIR_BLOTTNER_VISCOSITY,
...     AIR_KEYES_VISCOSITY,
...     AIR_VISCOSITY,
... )
>>> round(AIR_VISCOSITY.dynamic_viscosity(300.0), 8)
1.846e-05
>>> round(AIR_KEYES_VISCOSITY.dynamic_viscosity(300.0), 8)
1.852e-05
>>> round(AIR_BLOTTNER_VISCOSITY.dynamic_viscosity(3000.0), 8)
8.82e-05

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Model
     - Use
     - Fitted range
   * - ``AIR_VISCOSITY`` (Sutherland)
     - The default for atmosphere and boundary-layer work
     - No limit is enforced
   * - ``AIR_KEYES_VISCOSITY``
     - Low to moderate temperatures
     - 79 K to 1845 K
   * - ``AIR_BLOTTNER_VISCOSITY``
     - High temperatures
     - 1000 K to 30000 K

Outside its fitted range, the Keyes or Blottner model still returns a value but
emits an :class:`~aerophysics.exceptions.ApplicabilityWarning`. Functions that
need a viscosity, such as the boundary-layer and heating functions, accept any
of these through their ``viscosity_model`` argument.

Thermal conductivity of air is available from ``AIR_CONDUCTIVITY``:

>>> from aerophysics.transport import AIR_CONDUCTIVITY
>>> round(AIR_CONDUCTIVITY.thermal_conductivity(300.0), 5)
0.02625

At an altitude, :func:`~aerophysics.atmosphere.standard_atmosphere` already
returns the viscosity and conductivity, so you do not need to call these models
yourself.
