Sequential flow transformations
===============================

The **連続流れ計算** GUI page and :mod:`aerophysics.flow_sequence` apply an
ordered list of local transformations to a scalar flow state. Select the gas,
enter the initial condition, add stages, and press **計算**. Stages can be
named, duplicated, deleted, and moved up or down. Results include the initial
state and each exit, a stage plot, SI CSV export, and replayable settings JSON.
Input changes mark the previous result as stale until recalculated.

Python example
--------------

.. code-block:: python

   from math import radians
   from aerophysics.flow_sequence import (
       FlowState, ObliqueShockStep, ExpansionStep, NormalShockStep,
       solve_flow_sequence,
   )

   inlet = FlowState.from_atmosphere(10_000.0, mach=3.0)
   result = solve_flow_sequence(inlet, [
       ObliqueShockStep(radians(8), name="Compression ramp"),
       ExpansionStep(radians(3), name="Expansion corner"),
       NormalShockStep(name="Terminal shock"),
   ])
   assert result.complete
   for stage in result.steps:
       print(stage.step.name, stage.properties.state)
       print(stage.total_pressure_recovery, stage.entropy_change)

Initial conditions and units
----------------------------

``FlowState(pressure, temperature, mach)`` uses static pressure [Pa], static
absolute temperature [K], and non-negative Mach. ``from_velocity`` accepts a
speed magnitude [m/s] and gas; ``from_total`` accepts stagnation pressure,
stagnation temperature, Mach, and gas. Pass the same gas to initialization
and ``solve_flow_sequence``. ``from_atmosphere`` uses geometric altitude [m]
from -5000 to 86000 m and returns US Standard Atmosphere 1976 static p,T.
The GUI additionally accepts altitude and speed. For standard-atmosphere
initialization, altitude can be entered as geometric or geopotential altitude;
the selected coordinate is stored with the settings, while the Python API
continues to use geometric altitude. Standard-atmosphere p,T are initialization
data; density and sound speed use the selected gas model.

The API uses radians. GUI input and output honor the sidebar display units;
CSV and plots use SI. Flow speed is a magnitude, not a global direction.
Angles are non-negative local compression/expansion or cone half-angles.
The sequence does not track signed turns or global coordinates.

Stages and model boundaries
---------------------------

* ``NormalShockStep()`` applies the normal-shock equations.
* ``ObliqueShockStep(angle, branch)`` selects the weak (default) or strong
  attached solution.
* ``ConicalShockStep(angle)`` selects the weak attached Taylor--Maccoll
  solution and hands the **cone surface state** to the next stage. The
  immediate post-shock Mach and shock angle are also reported.
* ``ExpansionStep(angle)`` applies a centered Prandtl--Meyer expansion.
* ``IsentropicStep(value, basis, branch)`` specifies exit Mach, static
  pressure ratio ``p2/p1``, or area ratio ``A2/A1``. Area input requires an
  explicit ``MachBranch.SUBSONIC`` or ``MachBranch.SUPERSONIC`` exit branch.
  It computes ``A2/A* = (A2/A1)(A1/A*)`` at fixed total conditions and mass
  flow. Targets below the sonic minimum fail; an inlet at rest cannot
  define a finite area ratio. A branch change implies a suitable throat,
  whose geometry is not calculated.

The gas and composition remain fixed. The supported models are
``PerfectGas``, ``ThermallyPerfectGas`` (including NASA7/NASA9 air), and
``HarmonicOscillatorGas``. Thermal models use local heat capacity and
enthalpy, with their existing temperature limits and no extrapolation.
Temperature restoration from validated solver ratios repairs at most one ULP
of boundary overshoot; input temperature checks remain strict. Stagnation
temperature is bracketed upwards from the static temperature, including for
custom oscillator gases without an applicability range.

Chemical reactions, dissociation, heat addition, friction, flow separation,
wave interactions, duct startability and geometric realizability are outside
this calculation. In particular, chaining cone and planar stages is a local
state approximation, not a solution of a combined axisymmetric flow field.
See :doc:`../models/index` for each underlying model's sources and limits.

Results and failures
--------------------

Each successful stage supplies static p,T,M, density, velocity, sound speed,
dynamic pressure, mass flux, total enthalpy, total temperature and total
pressure. Pressure and temperature ratios refer to the immediately preceding
state; total-pressure recovery and entropy increase [J/(kg K)] refer to the
initial state. No cross-sectional area is specified, so mass flow is not
reported. Enthalpy uses the selected model's reference convention; compare
its conservation within the same model rather than across models.

If stagnation temperature alone exceeds a thermal model's range, static
properties remain available, while total pressure and temperature become
``None`` with a warning. Recovery and entropy still follow the stage losses.
An isentropic stage currently requires in-range stagnation conditions.

Invalid initial conditions raise an exception. A failed stage is classified
as ``invalid_input``, ``out_of_range``, ``detached``, or ``convergence_error``;
its message is retained. Earlier valid states remain, later stages have
``not_computed`` status, and ``result.complete`` is false. An empty sequence
returns just the initialized state. Missing values never imply a successful
calculation.

Verification
------------

The sequence tests compare normal shocks with independent constant-gamma
closed forms and mass/momentum conservation; compare cone and oblique
handoffs with their existing verified APIs; and verify total-enthalpy
conservation, shock entropy production, isentropic total-pressure
conservation, all three isentropic input modes, sonic/branch boundaries,
thermal limits, failure propagation and JSON replay. These checks establish
composition of the existing models, not CFD or experimental validation of
an arbitrary physical arrangement.
