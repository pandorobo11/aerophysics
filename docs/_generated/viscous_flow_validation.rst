**Overall status: Verified.**

Published-correlation comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Maximum relative differences from independently evaluated source equations
   :header-rows: 1

   * - Quantity
     - Maximum relative difference
     - Reynolds number
     - Result
   * - Blasius delta99
     - 4.441e-15
     - 2e+05
     - Pass
   * - Blasius displacement thickness
     - 8.882e-16
     - 5e+05
     - Pass
   * - Blasius momentum thickness
     - 1.998e-15
     - 1e+09
     - Pass
   * - Blasius local Cf
     - 1.998e-15
     - 1e+09
     - Pass
   * - Blasius average Cf
     - 6.661e-16
     - 1e+05
     - Pass
   * - power-law delta99
     - 1.998e-15
     - 1e+06
     - Pass
   * - power-law local Cf
     - 1.554e-15
     - 1e+07
     - Pass
   * - power-law average Cf
     - 3.109e-15
     - 1e+09
     - Pass
   * - Schlichting average Cf
     - 1.332e-15
     - 1e+07
     - Pass

NASA TN D-6945 chart comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For representative Figure 3(e) points, the largest absolute difference/tolerance ratio is ``0.3849`` at Mach ``5`` (``Pass``).

Physical and numerical invariants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Profile, limit, and integration checks
   :header-rows: 1

   * - Check
     - Criterion
     - Maximum difference
     - Result
   * - Schlichting local/average derivative
     - relative <= 1e-12
     - 0
     - Pass
   * - Van Driest II low-Mach limit
     - relative <= 1e-10
     - 0
     - Pass
   * - Van Driest constant-property transform
     - absolute <= 1e-10
     - 4.547e-13
     - Pass
   * - One-seventh-power protrusion integral
     - relative <= 1e-6
     - 6.059e-07
     - Pass
   * - Constant-profile protrusion integral
     - relative <= 1e-12
     - 2.22e-16
     - Pass
   * - Protrusion integration grid convergence
     - relative <= 1e-5
     - 8.173e-06
     - Pass
   * - Protrusion zero-height limit
     - shielding < 1e-3
     - 0.0002899
     - Pass
   * - Protrusion outside-layer limit
     - absolute distance from 1 <= 1e-4
     - 2.222e-05
     - Pass
