aerophysics
===========

``aerophysics`` provides traceable, vectorized engineering models for
atmospheric and aerodynamic physics. Calculation APIs use SI units and return
Python floats for scalar inputs or float64 NumPy arrays for array-like inputs.

The package covers calorically and thermally perfect gases, U.S. Standard
Atmosphere 1976, isentropic compressible flow, normal, oblique, and conical
shocks, Prandtl–Meyer expansion, smooth flat-plate boundary layers, integrated
flight conditions, and explicit unit conversions.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Thermodynamics and atmosphere

   gas_and_atmosphere
   transport_properties
   thermochemistry

.. toctree::
   :maxdepth: 2
   :caption: Compressible flow

   compressible_flow

.. toctree::
   :maxdepth: 2
   :caption: Boundary layers

   boundary_layers
   compressible_velocity_transformations

.. toctree::
   :maxdepth: 2
   :caption: Flight tools

   flight_conditions
   unit_conversions

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   references
