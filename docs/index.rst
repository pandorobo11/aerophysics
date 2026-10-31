aerophysics
===========

``aerophysics`` is a Python package for atmospheric and aerodynamic
calculations based on analytical and engineering models. Every model page
explains the physics, the equations, the assumptions, and the range in which
the model can be used, with references to the original sources. Calculation
APIs use SI units and return Python floats for scalar inputs or float64 NumPy
arrays for array-like inputs.

Start with the :doc:`getting_started/quickstart` to calculate an atmospheric
or flight condition. Read :doc:`getting_started/conventions` for the notation
and for help choosing a gas, viscosity, shock, or boundary-layer model. The
:doc:`gui/index` covers the optional local browser interface.

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
  compressibility corrections; compressible turbulent velocity profiles;
  stagnation-point aerodynamic heating; drag of small protrusions inside a
  boundary layer.
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

The manual has five main sections and two appendices, also shown in the left
sidebar. Headings on the current page are listed in the right sidebar.

* :doc:`getting_started/index`: installation, a first example, notation, and
  model selection.
* :doc:`guides/index`: step-by-step Python workflows.
* :doc:`gui/index`: the local browser interface, page by page.
* :doc:`models/index`: physical background, equations, assumptions, and
  sources for each model.
* :doc:`api/index`: every public class and function.
* Appendices: :doc:`verification/index` (how the implementation was checked)
  and :doc:`references`.

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
