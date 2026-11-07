.. list-table:: R&M 3633, M_inf=1.7, Re_l=10.08e6: local Cf
   :header-rows: 1

   * - X [in]
     - Measured 10^3 cf
     - SA-noft2 10^3 cf
     - Difference (10^3 cf)
     - abs difference / 0.0002
     - delta99/r_w
   * - 12
     - 3.69
     - 3.29
     - -0.40
     - 2.0
     - 0.028
   * - 18
     - 2.58
     - 2.93
     - +0.35
     - 1.8
     - 0.039
   * - 21
     - 2.75
     - 2.67
     - -0.08
     - 0.4
     - 0.053
   * - 24
     - 2.34
     - 2.37
     - +0.03
     - 0.2
     - 0.073
   * - 27
     - 1.96
     - 2.04
     - +0.08
     - 0.4
     - 0.102
   * - 28.5
     - 1.90
     - 1.87
     - -0.03
     - 0.1
     - 0.120
   * - 30
     - 1.60
     - 1.74
     - +0.14
     - 0.7
     - 0.141
   * - 31.5
     - 1.67
     - 1.66
     - -0.01
     - 0.1
     - 0.165
   * - 33
     - 1.47
     - 1.58
     - +0.11
     - 0.5
     - 0.194
   * - 34.5
     - 1.38
     - 1.53
     - +0.15
     - 0.8
     - 0.229
   * - 36
     - 1.38
     - 1.49
     - +0.11
     - 0.5
     - 0.271
   * - 37.5
     - 1.20
     - 1.42
     - +0.22
     - 1.1
     - 0.320
   * - 40
     - 1.07
     - 1.21
     - +0.14
     - 0.7
     - 0.418
   * - 42
     - 0.83
     - 0.93
     - +0.10
     - 0.5
     - 0.474
   * - 44
     - 0.92
     - 0.99
     - +0.07
     - 0.4
     - 0.456
   * - 45.5
     - 1.18
     - 1.19
     - +0.01
     - 0.0
     - 0.406
   * - 47
     - 1.46
     - 1.40
     - -0.06
     - 0.3
     - 0.357
   * - 48.5
     - 1.75
     - 1.55
     - -0.20
     - 1.0
     - 0.317
   * - 50
     - 1.93
     - 1.65
     - -0.28
     - 1.4
     - 0.284
   * - 51.5
     - 1.94
     - 1.77
     - -0.17
     - 0.9
     - 0.255
   * - 53
     - 2.23
     - 1.82
     - -0.41
     - 2.0
     - 0.230
   * - 54.5
     - 2.38
     - 1.94
     - -0.44
     - 2.2
     - 0.208
   * - 56
     - 2.27
     - 1.94
     - -0.33
     - 1.6
     - 0.188
   * - 57.5
     - 2.38
     - 2.01
     - -0.37
     - 1.9
     - 0.171
   * - 59
     - 2.43
     - 2.06
     - -0.37
     - 1.8
     - 0.155

Measured Cf accuracy is +/-0.0002; Cp accuracy is +/-0.005 (section 3.6).
The X=15 in Cf is missing. X<12 in is excluded because of transition.
This is a wind-tunnel comparison, not a calibrated accuracy claim.

.. list-table:: R&M 3633 sensitivity to assumptions
   :header-rows: 1

   * - Assumption
     - Maximum absolute difference (10^3 cf)
     - Maximum delta99/r_w
     - Outcome
   * - T0=310 K
     - 0.44
     - 0.477
     - Completed
   * - T0=290 K
     - 0.44
     - 0.478
     - Completed
   * - T0=330 K
     - 0.45
     - 0.477
     - Completed
   * - T0=310 K, cone Cp replaced
     - 0.45
     - 0.478
     - Completed

Warnings from the baseline:

* ApplicabilityWarning: transverse curvature is neglected; maximum delta99/r_w = 0.477 at s = 1.105 m.

The regression bound is 0.0006: the measured maximum over these cases is
rounded up to 0.0005, plus 0.0001 margin. This is not a physical accuracy
bound. Input construction and limitations are in the boundary-layer guide.
