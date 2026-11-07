.. list-table:: SA-noft2 ZPG differences from the frozen reference [%]
   :header-rows: 1

   * - Case
     - Re_s
     - Cf
     - delta99
     - theta
     - Heat coefficient
   * - M0.3-AD
     - 1e+06
     - 0.0155
     - 0.0252
     - -0.0086
     - -0.0000
   * - M0.3-AD
     - 1e+07
     - 0.0264
     - -0.0028
     - -0.0088
     - 0.0000
   * - M0.3-AD
     - 1e+08
     - 0.0553
     - 0.0645
     - 0.0090
     - 0.0000
   * - M2-AD
     - 1e+06
     - 0.0139
     - 0.0163
     - -0.0082
     - -0.0006
   * - M2-AD
     - 1e+07
     - 0.0233
     - 0.0147
     - -0.0124
     - 0.0005
   * - M2-AD
     - 1e+08
     - 0.0425
     - 0.0039
     - -0.0011
     - 0.0009
   * - M5-TW300
     - 1e+06
     - 0.0083
     - 0.0037
     - -0.0028
     - 0.0059
   * - M5-TW300
     - 1e+07
     - 0.0173
     - 0.0100
     - -0.0149
     - 0.0158
   * - M5-TW300
     - 1e+08
     - 0.0269
     - -0.0110
     - -0.0129
     - 0.0240
   * - M6-0.76
     - 1e+06
     - 0.0089
     - -0.0058
     - -0.0054
     - -0.0104
   * - M6-0.76
     - 1e+07
     - 0.0156
     - -0.0047
     - -0.0135
     - 0.0138
   * - M6-0.76
     - 1e+08
     - 0.0240
     - -0.0079
     - -0.0157
     - 0.0282
   * - M6-0.25
     - 1e+06
     - 0.0061
     - -0.0007
     - -0.0038
     - 0.0032
   * - M6-0.25
     - 1e+07
     - 0.0151
     - 0.0034
     - -0.0142
     - 0.0139
   * - M6-0.25
     - 1e+08
     - 0.0211
     - 0.0138
     - -0.0166
     - 0.0202
   * - M10-0.20
     - 1e+06
     - -0.0013
     - -0.0052
     - 0.0103
     - -0.0089
   * - M10-0.20
     - 1e+07
     - 0.0106
     - -0.0091
     - -0.0092
     - 0.0077
   * - M10-0.20
     - 1e+08
     - 0.0160
     - -0.0126
     - -0.0184
     - 0.0138

The recovery-temperature rows (AD) use an absolute heat-coefficient
difference divided by the reference Cf. All other differences are relative
to the reference quantity. The limit is 0.3% per entry; reference grid
uncertainty is about 0.1%. The reference displacement integrates only to
delta99 and is excluded here. Unit Reynolds number is 1e7 per metre,
gamma is 1.4, Pr is 0.72, Prt is 0.9, and Sutherland S is 110.4 K.

.. list-table:: NASA TMR bump: local wall shear / freestream dynamic pressure
   :header-rows: 1

   * - Cp source / Cf reference
     - Region (x)
     - Median difference [%]
     - Maximum absolute difference [%]
     - Maximum code spread [%]
   * - CFL3D
     - Upstream [0.05, 0.3)
     - 0.125
     - 0.235
     - 0.020
   * - CFL3D
     - Windward [0.3, 0.675)
     - -0.030
     - 0.314
     - 0.296
   * - CFL3D
     - Crest [0.675, 0.825)
     - 0.917
     - 1.409
     - 0.086
   * - CFL3D
     - Leeward [0.825, 1.2)
     - 0.664
     - 2.420
     - 0.472
   * - CFL3D
     - Recovery [1.2, 1.45)
     - 0.067
     - 0.086
     - 0.022
   * - FUN3D
     - Upstream [0.05, 0.3)
     - 0.121
     - 0.221
     - 0.020
   * - FUN3D
     - Windward [0.3, 0.675)
     - -0.050
     - 0.289
     - 0.297
   * - FUN3D
     - Crest [0.675, 0.825)
     - 0.889
     - 1.327
     - 0.086
   * - FUN3D
     - Leeward [0.825, 1.2)
     - 0.654
     - 2.497
     - 0.470
   * - FUN3D
     - Recovery [1.2, 1.45)
     - 0.080
     - 0.087
     - 0.022

Each march uses that code's wall Cp as edge Cp. The comparison excludes
x < 0.05 and x >= 1.45; the march ends at the last source point below 1.48.
Cp is replaced by its linearly interpolated value at x = 0.02 for x < 0.02,
including the added leading edge x = 0. No other smoothing or fitting
is used.
The wall distance is the trapezoidal integral of sqrt(1 + (dy/dx)^2) on
the source x nodes. Regions are upstream [0.05, 0.3), windward [0.3, 0.675),
crest [0.675, 0.825), leeward [0.825, 1.2), recovery [1.2, 1.45).

M_inf = 0.2 and Re_inf = 3e6 per metre. Dimensional edge nu_tilde is
3 nu_inf, converted to the initial local edge ratio. T_inf = 300 K,
an adiabatic wall,
and freestream dynamic-pressure normalization for both Cp and Cf are
assumptions: the source does not state these three details. Repeating at
250 K and 350 K changes each regional maximum by less than 0.007 percentage
points; temperature has little effect at this Mach number. The small
upstream-flat discrepancy supports the normalization assumption. Wall Cp
need not equal boundary-layer edge Cp on the curved bump. The regression
limit is 3.0%, from the measured maximum of 2.50% plus 0.50 percentage
points of margin; it is not an accuracy claim. The same-model CFD spread
is smaller than our differences over much of the bump. Accuracy against
experiment is not established by this comparison.
