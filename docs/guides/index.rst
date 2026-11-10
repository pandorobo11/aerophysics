User guide
==========

Each guide shows how to do one kind of calculation in Python, with short
examples you can copy. The guides explain which function to call, what to pass
in, and what comes back. For the equations and the limits of each model, follow
the links to :doc:`../models/index`.

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Guide
     - What you can calculate
   * - :doc:`atmosphere_flight`
     - Air properties at an altitude, and the dynamic pressure, Reynolds
       number, and total conditions of a flight at a given Mach number or speed
   * - :doc:`gas_properties`
     - Heat capacity, speed of sound, enthalpy, viscosity, and thermal
       conductivity of air, and how to choose a gas model
   * - :doc:`compressible_flow`
     - Isentropic flow and nozzles, shock waves, bow-shock standoff, and
       expansion waves
   * - :doc:`flow_sequence`
     - Several shocks and expansions in a row, such as a supersonic inlet
   * - :doc:`boundary_layers`
     - Thickness, skin friction, drag, and heat transfer on a flat plate;
       velocity profiles; drag of a small protrusion
   * - :doc:`marching_boundary_layers`
     - A boundary layer on a wall where the outer flow changes, including
       bodies of revolution such as a cone--cylinder
   * - :doc:`aerodynamic_heating`
     - Heat flux at the nose of a sphere or the leading edge of a swept
       cylinder
   * - :doc:`vectorization_errors`
     - Calculating many cases at once with arrays, and what the warnings and
       errors mean

If you prefer a browser interface to Python, see :doc:`../gui/index`.

.. toctree::
   :maxdepth: 1
   :hidden:

   atmosphere_flight
   gas_properties
   compressible_flow
   flow_sequence
   boundary_layers
   marching_boundary_layers
   aerodynamic_heating
   vectorization_errors
