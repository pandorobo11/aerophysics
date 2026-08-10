Conventions and model selection
===============================

Core conventions
----------------

* Calculation APIs use SI units. Angles are radians.
* Scalar inputs return Python ``float`` values; array-like inputs return
  float64 NumPy arrays after broadcasting.
* Isentropic ratios are total-to-static. Shock and expansion state ratios are
  downstream-to-upstream; shock total-pressure ratios are downstream total
  pressure divided by upstream total pressure.
* Non-finite and physically invalid inputs fail the whole call. Model range
  violations either raise :class:`~aerophysics.exceptions.ModelRangeError` or
  emit :class:`~aerophysics.exceptions.ApplicabilityWarning`, as documented by
  the selected model.

Thermodynamic model selection
-----------------------------

.. list-table:: Thermodynamic models
   :header-rows: 1
   :widths: 17 22 24 20 17

   * - Model
     - Start here when
     - Required state
     - Documented range
     - Excluded physics
   * - ``AIR``
     - Constant heat capacity is adequate
     - Temperature or Mach state
     - Positive temperature
     - Heat-capacity variation and real-gas effects
   * - ``AIR_NASA7`` / ``AIR_NASA9``
     - Frozen-air heat capacity must vary with temperature
     - Temperature; total temperature for isentropic flow
     - 200--6000 K
     - Reactions, dissociation, ionization
   * - ``AIR_HARMONIC_OSCILLATOR``
     - A Kennard frozen vibrational model is required
     - Temperature; total temperature for isentropic flow
     - Reservoir 400--2000 K
     - Reactions and dense-gas corrections
   * - ``AIR_BEATTIE_BRIDGEMAN``
     - The documented dense-gas engineering correction is required
     - Temperature and pressure
     - Tabulated air range: 38.8889--1222.2222 K and 172.369 Pa--27.5790 MPa
     - Reactions, phase change, equilibrium chemistry

See :doc:`../models/gas_and_atmosphere`,
:doc:`../models/thermochemistry`, and :doc:`../models/isentropic_flow` before
using a temperature-dependent or dense-gas model outside its documented
range.

Transport model selection
-------------------------

.. list-table:: Dynamic-viscosity models
   :header-rows: 1
   :widths: 19 25 21 18 17

   * - Model
     - Use
     - Required input
     - Nominal range
     - Range behavior
   * - Sutherland
     - Default dry-air atmosphere and boundary-layer work
     - Temperature
     - Positive temperature
     - Invalid temperatures raise ``ValueError``
   * - Keyes
     - Historical low-to-moderate-temperature comparison
     - Temperature
     - 79--1845 K
     - Extrapolation warns
   * - Blottner/Wilke
     - Frozen high-temperature species-mixture comparison
     - Temperature and fixed composition
     - 1000--30000 K
     - Extrapolation warns

These correlations do not model reacting composition. See
:doc:`../models/transport_properties`.

Flow and boundary-layer selection
---------------------------------

.. list-table:: Analysis choices
   :header-rows: 1
   :widths: 22 27 23 28

   * - Task
     - Recommended entry point
     - Key choice
     - Main exclusion
   * - Total/static state or area--Mach relation
     - :mod:`aerophysics.isentropic`
     - Gas model and Mach branch
     - Shocks and heat transfer
   * - Normal, oblique, conical, or detached shock
     - :doc:`../models/shock_waves`
     - Geometry and weak/strong branch where applicable
     - Most APIs assume a calorically perfect gas
   * - Centered expansion
     - :mod:`aerophysics.expansion`
     - Upstream Mach and nonnegative turn
     - Non-isentropic and reacting flow
   * - Smooth flat plate
     - :mod:`aerophysics.boundary_layer`
     - Regime, turbulent correlation, compressibility correction
     - Roughness, pressure gradient, separation, natural transition
   * - Compressible turbulent profile
     - :mod:`aerophysics.boundary_layer_profile`
     - Velocity transformation and temperature relation
     - General three-dimensional or separated boundary layers

See :doc:`../guides/vectorization_errors` for executable error-handling and
array examples.
