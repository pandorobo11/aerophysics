**Overall status: Verified.**

NACA Report 1135 table comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Maximum difference divided by the adopted table tolerance
   :header-rows: 1

   * - Quantity
     - Maximum ratio
     - Mach
     - Result
   * - pressure
     - 0.4999
     - 0.78
     - Pass
   * - density
     - 0.4999
     - 25.4
     - Pass
   * - temperature
     - 0.4997
     - 4.36
     - Pass
   * - area
     - 0.4857
     - 1.09
     - Pass
   * - nu
     - 0.4853
     - 9.93
     - Pass
   * - m2
     - 0.4997
     - 4.93
     - Pass
   * - p2
     - 0.5
     - 80
     - Pass
   * - rho2
     - 0.4998
     - 4.36
     - Pass
   * - t2
     - 0.5
     - 2
     - Pass
   * - pt2
     - 0.4997
     - 14.4
     - Pass

NACA Charts 2--4 observation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The largest difference/tolerance ratio is ``0.4571`` at Mach ``2.5``, deflection ``15 deg`` (``Pass``).

NASA SP-3004 cone-table comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The maximum relative difference is ``8.191e-05`` at Mach ``2``, cone half-angle ``30 deg`` (``Pass``).

Physical and mathematical invariants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Dense-grid invariant checks
   :header-rows: 1

   * - Check
     - Criterion
     - Maximum or margin
     - Result
   * - Isentropic ratio inverse round trips
     - relative <= 1e-10
     - 2.275e-12
     - Pass
   * - Area-Mach branch round trips
     - relative <= 1e-10
     - 5.796e-12
     - Pass
   * - Prandtl-Meyer inverse round trip
     - relative <= 1e-10
     - 1.553e-13
     - Pass
   * - Mass-flow maximum at M=1
     - positive neighbour margin
     - 5.703e-07
     - Pass
   * - Normal-shock mass closure
     - relative <= 1e-12
     - 4.441e-16
     - Pass
   * - Normal-shock momentum closure
     - relative <= 1e-12
     - 4.441e-16
     - Pass
   * - Normal-shock energy closure
     - relative <= 1e-12
     - 5.551e-16
     - Pass
   * - Oblique normal-component closure
     - relative <= 1e-12
     - 4.441e-16
     - Pass
   * - Weak/strong attached-limit merger
     - absolute <= 1e-10 rad
     - 0
     - Pass
   * - Zero-angle cone Mach-wave limit
     - absolute <= 1e-12
     - 0
     - Pass
