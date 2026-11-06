.. list-table:: Marching laminar boundary-layer verification
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - Case
     - Quantity
     - Computed
     - Reference
     - Difference [%]
   * - Falkner--Skan, beta = -0.15
     - ``f''(0)``
     - 0.2164
     - 0.2164
     - +0.03
   * - Falkner--Skan, beta = 0
     - ``f''(0)``
     - 0.4696
     - 0.4696
     - +0.01
   * - Falkner--Skan, beta = 0.5
     - ``f''(0)``
     - 0.9277
     - 0.9277
     - +0.00
   * - Similar, beta = -0.1, g_w = 0.2
     - ``f''(0)``
     - 0.4123
     - 0.4122
     - +0.01
   * - Similar, beta = -0.1, g_w = 0.2
     - ``g'(0)``
     - 0.365
     - 0.365
     - +0.01
   * - Similar, beta = 0.5, g_w = 2
     - ``f''(0)``
     - 1.235
     - 1.235
     - +0.00
   * - Similar, beta = 0.5, g_w = 2
     - ``g'(0)``
     - -0.5729
     - -0.5729
     - +0.01
   * - NACA 1293 Table II, beta = 1, S_w = 0
     - ``f''(0)``
     - 1.233
     - 1.233
     - +0.00
   * - NACA 1293 Table II, beta = 2, S_w = 0
     - ``f''(0)``
     - 1.687
     - 1.687
     - +0.01
   * - NACA 1293 Table II, beta = 0.5, S_w = 1
     - ``f''(0)``
     - 1.235
     - 1.235
     - -0.02
   * - NACA 1293 Table II, beta = 0.5, S_w = 1
     - ``g'(0)``
     - -0.5729
     - -0.5725
     - +0.07
   * - NACA 1293 Table II, beta = 2, S_w = 1
     - ``f''(0)``
     - 2.488
     - 2.488
     - -0.00
   * - NACA 1293 Table II, beta = 2, S_w = 1
     - ``g'(0)``
     - -0.6615
     - -0.6613
     - +0.03
   * - Howarth flow, M = 0.01
     - ``x_sep/L``
     - 0.1198
     - 0.1199
     - -0.07
   * - Howarth image, M0 = 3
     - ``x_sep/L``
     - 0.1198
     - 0.1198
     - -0.00
   * - Plate, M = 2, adiabatic
     - ``C_f``
     - 0.0001764
     - 0.0001763
     - +0.01
   * - Plate, M = 2, adiabatic
     - ``theta``
     - 0.0001763
     - 0.0001763
     - -0.02
   * - Plate, M = 10, T_w = 300 K
     - ``C_f``
     - 2.927e-05
     - 2.926e-05
     - +0.01
   * - Plate, M = 10, T_w = 300 K
     - ``theta``
     - 2.926e-05
     - 2.926e-05
     - -0.02
   * - Momentum integral, M = 2 to 3, T_w = 300 K
     - median residual
     - --
     - --
     - 0.02
   * - Total enthalpy integral, M = 2 to 3, T_w = 300 K
     - median residual
     - --
     - --
     - 0.07
