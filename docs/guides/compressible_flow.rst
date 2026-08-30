Compressible-flow workflows
===========================

Use this guide to choose a compressible-flow calculation and then follow the
linked model page for its equations, applicability limits, and complete symbol
definitions. All angles passed to the Python API are in radians.

Choose a workflow
-----------------

.. list-table::
   :header-rows: 1
   :widths: 24 27 27 22

   * - Task
     - Use
     - Required inputs
     - Important limit
   * - Relate a static state to its stagnation state or size a nozzle
     - :doc:`../models/isentropic_flow`
     - Mach number; reservoir temperature for a thermally perfect gas;
       reservoir pressure for Beattie--Bridgeman
     - Steady, adiabatic flow without entropy production
   * - Cross a normal, oblique, or attached conical shock
     - :doc:`../models/shock_waves`
     - Upstream Mach number; deflection or cone angle where applicable
     - Calorically perfect gas
   * - Estimate blunt-body shock standoff or shape
     - :ref:`detached-shocks` in :doc:`../models/shock_waves`
     - Mach number, nose radius, and explicit geometry
     - Engineering correlation, not a shock-layer solution
   * - Turn a supersonic stream through a centered expansion
     - :doc:`../models/expansion_waves`
     - Upstream Mach number and nonnegative turn angle
     - Calorically perfect gas; downstream angle below the limiting
       Prandtl--Meyer angle

Isentropic state ratios
-----------------------

State ratios use the total-to-static convention:

>>> from aerophysics.isentropic import isentropic_ratios
>>> ratios = isentropic_ratios(2.0)
>>> round(ratios.total_temperature_ratio, 6)
1.8
>>> round(ratios.total_pressure_ratio, 6)
7.824449

Use :class:`~aerophysics.isentropic.MachBranch` when inverting the
area--Mach relation because every area ratio above one has both a subsonic and
a supersonic solution. Thermally perfect and Beattie--Bridgeman workflows use
the same public API but require the reservoir inputs described in
:doc:`../models/isentropic_flow`.

When one workflow needs ratios, area, mass-flow, and absolute-state values
together, use :func:`~aerophysics.isentropic.isentropic_analysis`. It evaluates
the coupled results from one set of flow states and shares each Mach-one state
across matching reservoir conditions.

Shock and expansion turns
-------------------------

Angles are supplied in radians and shock branch selection is explicit:

>>> from aerophysics import ShockBranch, oblique_shock, prandtl_meyer_expansion
>>> from aerophysics.units import degrees_to_radians, radians_to_degrees
>>> shock = oblique_shock(
...     2.0, degrees_to_radians(10.0), branch=ShockBranch.WEAK
... )
>>> round(radians_to_degrees(shock.shock_angle), 3)
39.314
>>> expansion = prandtl_meyer_expansion(2.0, degrees_to_radians(10.0))
>>> round(expansion.downstream_mach, 6)
2.384887

An oblique or conical request above its attached-shock limit raises
:class:`~aerophysics.exceptions.NoAttachedShockError`; it is not silently
replaced by a normal-shock approximation. A Prandtl--Meyer turn must be an
expansion and must remain below the finite limiting angle.

Next steps
----------

- Read :doc:`../models/isentropic_flow` for perfect-gas, NASA-polynomial,
  harmonic-oscillator, and Beattie--Bridgeman state relations.
- Read :doc:`../models/shock_waves` for normal, oblique, conical, and
  detached-shock models.
- Read :doc:`../models/expansion_waves` for the Prandtl--Meyer equations and
  downstream state ratios.
- Check :doc:`../verification/compressible_flow` before using a model near its
  documented applicability limit.
