Several waves in a row
======================

Many practical flows pass through several waves, one after another. A
supersonic inlet, for example, may compress the air with one or two oblique
shocks, turn it through an expansion, and finish with a normal shock.

:mod:`aerophysics.flow_sequence` calculates such a path. You describe the
starting flow and a list of stages. Each stage uses the single-wave model from
:doc:`compressible_flow`, and the flow leaving one stage enters the next. The
GUI offers the same calculation on its **連続流れ計算** (sequential flow) page;
see :doc:`../gui/flight_and_sequence_pages`.

A first example
---------------

This inlet flies at Mach 3 at 10 km. It has an 8° compression ramp, a 3°
expansion corner, and a normal shock at the end:

>>> from math import radians
>>> from aerophysics.flow_sequence import (
...     ExpansionStep,
...     FlowState,
...     NormalShockStep,
...     ObliqueShockStep,
...     solve_flow_sequence,
... )
>>> inlet = FlowState.from_atmosphere(10_000.0, mach=3.0)
>>> result = solve_flow_sequence(
...     inlet,
...     [
...         ObliqueShockStep(radians(8.0), name="Compression ramp"),
...         ExpansionStep(radians(3.0), name="Expansion corner"),
...         NormalShockStep(name="Terminal shock"),
...     ],
... )
>>> result.complete
True
>>> for stage in result.steps:
...     print(
...         stage.step.name,
...         round(stage.properties.state.mach, 4),
...         round(stage.total_pressure_recovery, 4),
...     )
Compression ramp 2.6031 0.9799
Expansion corner 2.7402 0.9799
Terminal shock 0.4925 0.4014

The last column shows that about 40% of the free-stream total pressure is left
after the normal shock.

The starting flow
-----------------

A :class:`~aerophysics.flow_sequence.FlowState` holds the static pressure in
Pa, the static temperature in K, and the Mach number. Create it in whichever
way suits your data:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - You know
     - Use
   * - Static pressure, static temperature, and Mach number
     - ``FlowState(pressure, temperature, mach)``
   * - Static pressure, static temperature, and speed
     - ``FlowState.from_velocity(...)``
   * - Total pressure, total temperature, and Mach number
     - ``FlowState.from_total(...)``
   * - Altitude and Mach number
     - ``FlowState.from_atmosphere(geometric_altitude, mach)``

``from_atmosphere`` takes the pressure and temperature from the U.S. Standard
Atmosphere 1976 between -5 km and 86 km.

The stages
----------

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Stage
     - What it does
   * - ``NormalShockStep()``
     - A normal shock.
   * - ``ObliqueShockStep(angle)``
     - An oblique shock that turns the flow by ``angle``. The weak solution is
       used unless you pass ``branch=ShockBranch.STRONG``.
   * - ``ConicalShockStep(angle)``
     - A shock on a cone with half-angle ``angle``. The next stage starts from
       the flow on the cone surface.
   * - ``ExpansionStep(angle)``
     - A Prandtl--Meyer expansion that turns the flow by ``angle``.
   * - ``IsentropicStep(value, basis)``
     - A smooth change without losses, as in a nozzle or diffuser. ``value`` is
       the exit Mach number (``basis="mach"``, the default), the pressure ratio
       :math:`p_2/p_1` (``"pressure_ratio"``), or the area ratio
       :math:`A_2/A_1` (``"area_ratio"``).

Angles are in radians and are never negative: a shock stage always compresses
and an expansion stage always expands. The calculation follows only the state
of the gas. It does not keep track of the direction of the flow or the shape of
the duct.

With ``basis="area_ratio"``, also pass ``branch=MachBranch.SUBSONIC`` or
``MachBranch.SUPERSONIC`` to say which exit flow you want. Going from subsonic
to supersonic, or back, assumes that a throat lies between the two stations.
The stage fails if the exit area is smaller than the throat that the flow
needs.

Every stage accepts a ``name``, which is used only to label the results.

Reading the results
-------------------

``result.initial`` describes the starting flow, and ``result.steps`` holds one
entry per stage. For each stage:

* ``properties`` holds the flow leaving the stage: ``state`` (static pressure,
  temperature, and Mach number), density, velocity, speed of sound, dynamic
  pressure, mass flux, total enthalpy, total temperature, and total pressure.
  The mass flux is per unit area. No duct area is given, so a mass flow is not
  reported.
* ``pressure_ratio`` and ``temperature_ratio`` compare the stage's exit with
  its own inlet.
* ``total_pressure_recovery`` and ``entropy_change`` (in J/(kg K)) compare the
  stage's exit with the **starting** flow, so they accumulate along the
  sequence.
* Shock stages also report ``shock_angle``. A conical stage reports
  ``post_shock_mach``, the Mach number just behind the shock, in addition to
  the surface state.

When a stage cannot be calculated
---------------------------------

A stage can fail, for example because the turn angle is too large for an
attached shock. The calculation then stops at that stage, but it does not raise
an exception. Check ``result.complete`` and each stage's ``status``:

>>> failed = solve_flow_sequence(
...     inlet, [ObliqueShockStep(radians(40.0)), NormalShockStep()]
... )
>>> failed.complete
False
>>> [stage.status for stage in failed.steps]
['detached', 'not_computed']

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Status
     - Meaning
   * - ``ok``
     - The stage was calculated.
   * - ``invalid_input``
     - The stage cannot apply to its inlet flow, for example a shock in
       subsonic flow.
   * - ``out_of_range``
     - A temperature left the range of the gas model.
   * - ``detached``
     - The wedge or cone angle is too large for an attached shock.
   * - ``convergence_error``
     - The numerical solution did not converge.
   * - ``not_computed``
     - An earlier stage failed, so this one was skipped.

A failed stage keeps an explanation in ``message``, and the stages before it
keep their results. An invalid starting flow, such as a negative pressure,
does raise an exception.

Hot flows
---------

Pass a gas model as the third argument to ``solve_flow_sequence`` to let the
heat capacity change with temperature. ``AIR``, ``AIR_NASA7``, ``AIR_NASA9``,
and ``AIR_HARMONIC_OSCILLATOR`` are accepted. If you built the starting flow
with ``from_velocity`` or ``from_total``, give the same gas there.

With these models, the total temperature of a fast flow can lie above the
model's range even though the static temperature is inside it. The static
results are then still calculated, and ``total_temperature`` and
``total_pressure`` are ``None`` with a note in ``properties.warning``. An
``IsentropicStep`` needs the total conditions and fails in this situation.

What the calculation leaves out
-------------------------------

The sequence applies one ideal wave after another. It does not include:

* friction, heat addition, or flow separation;
* interactions between waves, such as a shock reflecting from a wall or meeting
  another shock;
* chemical reactions or changes in the composition of the gas;
* whether the duct can actually start or hold the flow you describe.

A cone followed by a wedge or a corner is therefore only a rough estimate. The
real flow around such a shape is axisymmetric throughout. For a cone followed
by a cylinder, :doc:`marching_boundary_layers` shows a more complete method.

The package's tests check each stage against the single-wave functions and
check that mass, momentum, and energy are conserved. They confirm that the
stages are combined correctly, not that the result matches a particular inlet
or nozzle.
