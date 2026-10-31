Conventions and model selection
===============================

Core conventions
----------------

* Calculation APIs use SI units. Angles are in radians.
* A scalar input returns a Python ``float``. Array-like inputs are broadcast
  together and return float64 NumPy arrays.
* Isentropic ratios are total over static, for example :math:`T_0/T`.
  Shock and expansion ratios are downstream over upstream, for example
  :math:`p_2/p_1`. The shock total-pressure ratio is :math:`p_{02}/p_{01}`.
* An invalid input (not finite, or physically impossible such as a negative
  temperature) makes the whole call fail. Outside a model's supported range,
  the model either raises :class:`~aerophysics.exceptions.ModelRangeError` or
  emits :class:`~aerophysics.exceptions.ApplicabilityWarning`. Each model page
  states which.

Notation
--------

The model pages use these common symbols and subscripts.

.. list-table:: Common notation
   :header-rows: 1
   :widths: 22 78

   * - Symbol
     - Meaning
   * - :math:`M`, :math:`V`, :math:`a`
     - Mach number, flow speed, and speed of sound, :math:`M=V/a`
   * - :math:`p`, :math:`T`, :math:`\rho`
     - Static pressure, temperature, and density: the values measured while
       moving with the flow
   * - :math:`h`, :math:`s`
     - Specific enthalpy and specific entropy (per unit mass)
   * - :math:`\gamma`, :math:`R`
     - Heat-capacity ratio :math:`c_p/c_v` and specific gas constant
   * - Subscript :math:`0`
     - Total (stagnation) value: the state reached by bringing the flow to
       rest without losses
   * - Superscript :math:`*`
     - Sonic (critical) value at :math:`M=1`, except :math:`T^*` in the
       Eckert method, which is a reference temperature
   * - Subscripts :math:`1`, :math:`2`
     - Upstream and downstream of a shock or expansion
   * - Subscript :math:`\infty`
     - Free stream, far ahead of the body
   * - Subscripts :math:`e`, :math:`w`
     - Boundary-layer edge and wall
   * - :math:`\theta`, :math:`\beta`
     - Flow-deflection angle and shock angle

Thermodynamic model selection
-----------------------------

.. list-table:: Thermodynamic models
   :header-rows: 1
   :widths: 17 22 24 20 17

   * - Model
     - Use when
     - Required state
     - Documented range
     - Not included
   * - ``AIR``
     - Constant heat capacity is good enough (air below about 800 K)
     - Mach number, or temperature
     - Any positive temperature
     - Heat-capacity variation and real-gas effects
   * - ``AIR_NASA7`` / ``AIR_NASA9``
     - Heat capacity must vary with temperature
     - Temperature (total temperature for isentropic flow)
     - 200--6000 K
     - Reactions, dissociation, ionization
   * - ``AIR_HARMONIC_OSCILLATOR``
     - A simple vibrational-excitation model is wanted
     - Temperature (total temperature for isentropic flow)
     - Reservoir 400--2000 K
     - Reactions and dense-gas corrections
   * - ``AIR_BEATTIE_BRIDGEMAN``
     - High reservoir pressure makes the ideal-gas law inaccurate
     - Temperature and pressure
     - Tabulated air range: 38.8889--1222.2222 K and 172.369 Pa--27.5790 MPa
     - Reactions, phase change, equilibrium chemistry

Not every calculation accepts every gas model:

.. list-table:: Gas models accepted by each calculation
   :header-rows: 1
   :widths: 34 12 18 18 18

   * - Calculation
     - ``AIR``
     - NASA7 / NASA9
     - Harmonic oscillator
     - Beattie--Bridgeman
   * - Isentropic flow
     - Yes
     - Yes
     - Yes
     - Yes
   * - Normal, oblique, and conical shocks
     - Yes
     - Yes
     - Yes
     - No
   * - Prandtl--Meyer expansion
     - Yes
     - Yes
     - Yes
     - No
   * - Sequential flow
     - Yes
     - Yes
     - Yes
     - No
   * - Detached-shock correlations, flat-plate boundary layer, aerodynamic
       heating
     - Yes (calorically perfect only)
     - No
     - No
     - No

See :doc:`../models/gas_models` and :doc:`../models/isentropic_flow` before
using a temperature-dependent or dense-gas model near the edge of its range.

Transport model selection
-------------------------

.. list-table:: Dynamic-viscosity models
   :header-rows: 1
   :widths: 19 25 21 18 17

   * - Model
     - Use
     - Required input
     - Nominal range
     - Outside the range
   * - Sutherland
     - Default for dry air in atmosphere and boundary-layer work
     - Temperature
     - Any positive temperature
     - Invalid temperatures raise ``ValueError``
   * - Keyes
     - Comparison at low to moderate temperature
     - Temperature
     - 79--1845 K
     - Warns and extrapolates
   * - Blottner/Wilke
     - Comparison at high temperature with fixed composition
     - Temperature
     - 1000--30000 K
     - Warns and extrapolates

These correlations do not model changes in composition. See
:doc:`../models/transport_properties`.

Choosing a calculation
----------------------

.. list-table:: Analysis choices
   :header-rows: 1
   :widths: 22 27 23 28

   * - Task
     - Entry point
     - Key choice
     - Main exclusion
   * - Total/static state, nozzle area, or mass flow
     - :mod:`aerophysics.isentropic`
     - Gas model and subsonic/supersonic branch
     - Shocks, friction, and heat transfer
   * - Normal, oblique, or conical shock
     - :mod:`aerophysics.shocks`
     - Geometry, and weak/strong branch for oblique shocks
     - Viscosity, chemistry, detached flow
   * - Bow-shock standoff and shape on a blunt nose
     - :mod:`aerophysics.detached_shock`
     - Sphere or two-dimensional cylinder
     - Shock-layer thermodynamics (correlation only)
   * - Expansion around a convex corner
     - :mod:`aerophysics.expansion`
     - Upstream Mach number and nonnegative turn angle
     - Losses and reacting flow
   * - Several waves in a row (for example, an inlet)
     - :mod:`aerophysics.flow_sequence`
     - Order and type of stages
     - Wave interactions, friction, and geometry
   * - Smooth flat plate
     - :mod:`aerophysics.boundary_layer`
     - Laminar/turbulent state, turbulent correlation, compressibility
       correction
     - Roughness, pressure gradient, separation, transition prediction
   * - Stagnation-point heating
     - :mod:`aerophysics.heating`
     - Sphere or swept cylinder; Fay--Riddell or Tauber
     - Turbulence, dissociation, radiation
   * - Compressible turbulent velocity profile
     - :mod:`aerophysics.boundary_layer_profile`
     - Velocity transformation and temperature relation
     - Three-dimensional or separated boundary layers
   * - Drag of a small protrusion inside a boundary layer
     - :mod:`aerophysics.protrusion`
     - Velocity profile and supplied drag coefficient
     - Flow around the protrusion itself

See :doc:`../guides/vectorization_errors` for array and error-handling
examples.
