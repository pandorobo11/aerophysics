.. list-table:: R&M 3633, M_inf=1.7, Re_l=10.08e6: local Cf
   :header-rows: 1

   * - X [in]
     - Measured Cf
     - SA-noft2 Cf
     - Difference
     - abs difference / 0.0002
     - delta99/r_w
   * - 12
     - 0.003690
     - 0.003292
     - -0.000398
     - 1.99
     - 0.0275
   * - 18
     - 0.002580
     - 0.002933
     - +0.000353
     - 1.77
     - 0.0389
   * - 21
     - 0.002750
     - 0.002671
     - -0.000079
     - 0.40
     - 0.0528
   * - 24
     - 0.002340
     - 0.002373
     - +0.000033
     - 0.16
     - 0.0734
   * - 27
     - 0.001960
     - 0.002043
     - +0.000083
     - 0.42
     - 0.1024
   * - 28.5
     - 0.001900
     - 0.001874
     - -0.000026
     - 0.13
     - 0.1202
   * - 30
     - 0.001600
     - 0.001740
     - +0.000140
     - 0.70
     - 0.1406
   * - 31.5
     - 0.001670
     - 0.001657
     - -0.000013
     - 0.07
     - 0.1653
   * - 33
     - 0.001470
     - 0.001577
     - +0.000107
     - 0.54
     - 0.1944
   * - 34.5
     - 0.001380
     - 0.001534
     - +0.000154
     - 0.77
     - 0.2294
   * - 36
     - 0.001380
     - 0.001489
     - +0.000109
     - 0.54
     - 0.2712
   * - 37.5
     - 0.001200
     - 0.001423
     - +0.000223
     - 1.11
     - 0.3203
   * - 40
     - 0.001070
     - 0.001214
     - +0.000144
     - 0.72
     - 0.4182
   * - 42
     - 0.000830
     - 0.000926
     - +0.000096
     - 0.48
     - 0.4742
   * - 44
     - 0.000920
     - 0.000993
     - +0.000073
     - 0.37
     - 0.4561
   * - 45.5
     - 0.001180
     - 0.001185
     - +0.000005
     - 0.03
     - 0.4062
   * - 47
     - 0.001460
     - 0.001403
     - -0.000057
     - 0.29
     - 0.3571
   * - 48.5
     - 0.001750
     - 0.001545
     - -0.000205
     - 1.02
     - 0.3173
   * - 50
     - 0.001930
     - 0.001654
     - -0.000276
     - 1.38
     - 0.2836
   * - 51.5
     - 0.001940
     - 0.001767
     - -0.000173
     - 0.86
     - 0.2549
   * - 53
     - 0.002230
     - 0.001823
     - -0.000407
     - 2.04
     - 0.2298
   * - 54.5
     - 0.002380
     - 0.001936
     - -0.000444
     - 2.22
     - 0.2080
   * - 56
     - 0.002270
     - 0.001942
     - -0.000328
     - 1.64
     - 0.1881
   * - 57.5
     - 0.002380
     - 0.002006
     - -0.000374
     - 1.87
     - 0.1705
   * - 59
     - 0.002430
     - 0.002062
     - -0.000368
     - 1.84
     - 0.1548

Measured Cf accuracy is +/-0.0002; Cp accuracy is +/-0.005 (section 3.6).
The X=15 in Cf is missing. X<12 in is excluded because of transition.
This is a wind-tunnel comparison, not a calibrated accuracy claim.

.. list-table:: R&M 3633 sensitivity to assumptions
   :header-rows: 1

   * - Assumption
     - Maximum absolute Cf difference
     - Maximum delta99/r_w
     - Outcome
   * - T0=310 K
     - 0.000444
     - 0.4774
     - Completed
   * - T0=290 K
     - 0.000442
     - 0.4779
     - Completed
   * - T0=330 K
     - 0.000446
     - 0.4769
     - Completed
   * - T0=310 K, cone Cp replaced
     - 0.000445
     - 0.4779
     - Completed

Warnings from the baseline:

* transverse curvature is neglected but delta99/r_w = 0.477401 (maximum at s = 1.1049 m)

The regression bound is 0.0006: the measured maximum over these cases is
rounded up to 0.0005, plus 0.0001 margin. This is not a physical accuracy
bound. Input construction and limitations are in the boundary-layer guide.
