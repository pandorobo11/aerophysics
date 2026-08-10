aerophysics
===========

``aerophysics`` provides traceable, vectorized Python models for atmospheric
and aerodynamic physics. Calculation APIs use SI units and return Python
floats for scalar inputs or float64 NumPy arrays for array-like inputs.

Start with the :doc:`getting_started/quickstart` if you want to calculate an
atmospheric or flight condition. Use :doc:`getting_started/conventions` before
selecting a temperature-dependent gas, transport correlation, shock model, or
boundary-layer correction. The :doc:`guides/gui` guide covers the optional
local browser interface.

What the package covers
-----------------------

The package covers perfect and selected real-gas engineering models, U.S.
Standard Atmosphere 1976, isentropic compressible flow, shock and expansion
relations, smooth flat-plate boundary layers, integrated flight conditions,
and explicit unit conversions. Each model page records equations,
assumptions, valid or fitted ranges, and primary sources. Reproducible checks
are collected under :doc:`verification/index`.

What it does not cover
----------------------

``aerophysics`` is not an aircraft-design or CFD package. It does not predict
complete vehicle aerodynamics, wing lift, interference drag, separated flow,
chemical equilibrium, ionization, phase change, or weather. Individual model
pages document narrower exclusions.

.. toctree::
   :maxdepth: 1
   :caption: Getting started

   getting_started/installation
   getting_started/quickstart
   getting_started/conventions

.. toctree::
   :maxdepth: 1
   :caption: Task guides

   guides/atmosphere_flight
   guides/compressible_flow
   guides/boundary_layers
   guides/vectorization_errors
   guides/gui

.. toctree::
   :maxdepth: 1
   :caption: Models and equations

   models/gas_and_atmosphere
   models/transport_properties
   models/thermochemistry
   models/isentropic_flow
   models/shock_waves
   models/expansion_waves
   models/flat_plate_boundary_layer
   models/compressible_velocity_transformations
   models/protrusion_drag
   models/flight_conditions
   models/unit_conversions

.. toctree::
   :maxdepth: 1
   :caption: Verification

   verification/index

.. toctree::
   :maxdepth: 1
   :caption: API and references

   api/index
   references
