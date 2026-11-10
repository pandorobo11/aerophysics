aerophysics
===========

``aerophysics`` is a Python package for quick atmospheric and aerodynamic
calculations: the properties of the atmosphere, shock waves, boundary layers,
aerodynamic heating, and more. It uses well-known analytical and engineering
models, and this manual explains what each model assumes and where it can be
used.

New to the package? Start with the :doc:`getting_started/quickstart`. Then
pick the guide for your task from the :doc:`guides/index`. If you would rather
use a browser than write Python, see the :doc:`gui/index`.

What the package covers
-----------------------

* **Gas properties and atmosphere:** perfect, thermally perfect, and selected
  real-gas air models; viscosity and thermal conductivity; the U.S. Standard
  Atmosphere 1976; flight conditions such as dynamic pressure and Reynolds
  number.
* **Compressible flow:** isentropic flow and nozzles; normal, oblique, and
  conical shocks; bow-shock standoff and shape; Prandtl--Meyer expansions;
  and sequences of these waves.
* **Viscous flow and heating:** smooth flat-plate boundary layers with
  compressibility corrections; boundary layers that follow a changing outer
  flow, including bodies of revolution; compressible turbulent velocity
  profiles; stagnation-point aerodynamic heating; drag of small protrusions
  inside a boundary layer.
* **Unit conversion** between customary units and SI.

Reproducible checks against published tables and independent solutions are
collected under :doc:`verification/index`.

What it does not cover
----------------------

``aerophysics`` is not an aircraft-design or CFD package. It does not predict
complete vehicle aerodynamics, wing lift, interference drag, separated flow,
chemical equilibrium, ionization, phase change, or weather. Each model page
lists its own exclusions.

Documentation map
-----------------

The manual is arranged from "how do I use it" to "how does it work":

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Section
     - Read it when you want to
   * - :doc:`getting_started/index`
     - install the package, try it, and learn the shared conventions
   * - :doc:`guides/index`
     - carry out a particular calculation in Python, with examples
   * - :doc:`gui/index`
     - use the browser interface
   * - :doc:`models/index`
     - understand the physics, equations, assumptions, and sources of a model
   * - :doc:`api/index`
     - look up the exact inputs and outputs of a function
   * - :doc:`verification/index`
     - see how the results were checked against published data
   * - :doc:`references`
     - find the original papers and reports

The same sections are in the left sidebar. The headings of the current page
are in the right sidebar.

.. toctree::
   :maxdepth: 2
   :titlesonly:
   :hidden:

   getting_started/index
   guides/index
   gui/index
   models/index
   api/index

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :hidden:
   :caption: Appendices

   verification/index
   references
