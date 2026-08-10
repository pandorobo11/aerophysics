**Overall status: Verified.**

Cantera 3.2.0 comparison
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Maximum relative differences from the pinned snapshot
   :header-rows: 1

   * - Quantity
     - Maximum relative difference
     - Temperature [K]
     - Result
   * - cp
     - 1.581e-09
     - 1000
     - Pass
   * - enthalpy
     - 4.83e-10
     - 1000
     - Pass
   * - entropy
     - 2.019e-10
     - 1000
     - Pass
   * - gamma
     - 5.315e-10
     - 1000
     - Pass
   * - speed of sound
     - 2.658e-10
     - 1000
     - Pass

Primary transport references
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Comparisons with U.S. Standard Atmosphere 1976
   :header-rows: 1

   * - Reference
     - Criterion
     - Maximum difference
     - Temperature [K]
     - Result
   * - Sutherland / USSA Table III
     - within 2 printed half-digits or 1e-4 relative
     - 4.202e-05
     - 255.676
     - Pass
   * - USSA conductivity / Equation (53) errata example
     - absolute <= 5e-7 W/(m K)
     - 1.157e-07
     - 288.15
     - Pass

Published transport equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Maximum relative differences from direct source-equation reproductions
   :header-rows: 1

   * - Model
     - Maximum relative difference
     - Temperature [K]
     - Result
   * - Sutherland viscosity
     - < 5e-13
     - —
     - Pass
   * - Keyes viscosity
     - < 5e-13
     - —
     - Pass
   * - Blottner/Wilke viscosity
     - < 5e-13
     - —
     - Pass
   * - USSA conductivity
     - < 5e-13
     - —
     - Pass

NIST physical-accuracy assessment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Lemmon--Jacobsen evaluated dilute-air correlation is an independent physical-accuracy reference, not an acceptance test for the intentionally simpler Sutherland and USSA correlations. Across 250--1500 K at zero density, the largest absolute relative differences are ``6.613%`` for viscosity and ``3.777%`` for conductivity. The NIST source estimates dilute-gas uncertainties of 1% and 2%, respectively, over this range.

Thermodynamic invariants
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Dense-grid checks
   :header-rows: 1

   * - Check
     - Criterion
     - Maximum or minimum
     - Result
   * - NASA7 heat-capacity identities
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - NASA7 dh/dT = cp
     - relative <= 1e-7
     - < 5e-08
     - Pass
   * - NASA7 1000 K region continuity
     - relative <= 1e-7
     - < 5e-08
     - Pass
   * - NASA9 heat-capacity identities
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - NASA9 dh/dT = cp
     - relative <= 1e-7
     - < 5e-08
     - Pass
   * - NASA9 1000 K region continuity
     - relative <= 1e-7
     - < 5e-08
     - Pass
   * - Harmonic-oscillator cp-cv=R
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - Beattie-Bridgeman stable density root
     - minimum dp/drho > 0
     - 1.153e+05
     - Pass
   * - Beattie-Bridgeman pressure closure
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - NASA7 isentropic total enthalpy
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - NASA7 isentropic entropy
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - NASA9 isentropic total enthalpy
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - NASA9 isentropic entropy
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - Harmonic oscillator isentropic total enthalpy
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - Harmonic oscillator isentropic entropy
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - Beattie-Bridgeman isentropic total enthalpy
     - relative <= 1e-12
     - < 5e-13
     - Pass
   * - Beattie-Bridgeman isentropic entropy
     - relative <= 1e-12
     - < 5e-13
     - Pass
