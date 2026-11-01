Isentropic flow
===============

These models assume steady, adiabatic flow without entropy production.
They may use constant heat capacities, frozen-composition
temperature-dependent heat capacities, or the Beattie--Bridgeman dense-gas
equation of state. State ratios are total over static, and angles are in
radians. Convert explicitly with
:func:`aerophysics.units.degrees_to_radians` and
:func:`aerophysics.units.radians_to_degrees`.

For task-oriented model selection and examples, see
:doc:`../guides/compressible_flow`.

Physical background
-------------------

Isentropic means "at constant entropy": no friction, no heat transfer, and no
shock waves. Real flows are close to this away from walls and shocks, for
example in a well-designed nozzle or outside a thin boundary layer.

**Total (stagnation) state.** Imagine bringing a moving gas to rest without
losses. The temperature, pressure, and density it would then have are the
total, or stagnation, values :math:`T_0`, :math:`p_0`, and :math:`\rho_0`. In
steady adiabatic flow the total enthalpy is constant along a streamline:

.. math::

   h_0=h+\frac{V^2}{2}.

For a calorically perfect gas, :math:`h=c_pT`. Dividing by :math:`c_pT` and
using :math:`a^2=\gamma RT` and :math:`c_p=\gamma R/(\gamma-1)` gives the
temperature ratio below. Because the process is isentropic, pressure and
density follow :math:`p\propto T^{\gamma/(\gamma-1)}` and
:math:`\rho\propto T^{1/(\gamma-1)}`.

**Area and Mach number.** In a duct, mass conservation requires
:math:`\rho VA` to stay constant. In subsonic flow, a narrowing duct speeds
the flow up. In supersonic flow the opposite happens: density falls faster
than velocity rises, so the duct must widen to keep accelerating the flow.
Mach one can therefore occur only at the narrowest section, the throat. The
throat area for sonic flow is written :math:`A^*`.

**Choking.** Once the throat is sonic, lowering the downstream pressure
further cannot increase the mass flow. The flow is then choked, and the mass
flow per unit throat area depends only on the reservoir state.

Calorically perfect gas
-----------------------

Define the temperature factor

.. math::

   F(M)=1+\frac{\gamma-1}{2}M^2.

The total-to-static ratios returned by
:func:`aerophysics.isentropic.isentropic_ratios` are

.. math::

   \frac{T_0}{T}=F(M),
   \qquad
   \frac{p_0}{p}=F(M)^{\gamma/(\gamma-1)},
   \qquad
   \frac{\rho_0}{\rho}=F(M)^{1/(\gamma-1)}.

The inverse functions ``mach_from_total_temperature_ratio``,
``mach_from_total_pressure_ratio``, and ``mach_from_total_density_ratio``
solve these expressions for :math:`M`. The values at :math:`M=1` are the
total-to-critical (sonic) ratios returned by
:func:`aerophysics.isentropic.critical_ratios`; for air with
:math:`\gamma=1.4`, :math:`T_0/T^*=1.2` and :math:`p_0/p^*\approx1.893`.

The quasi-one-dimensional area--Mach relation is

.. math::

   \frac{A}{A^*}
   =\frac{1}{M}
    \left[\frac{2}{\gamma+1}F(M)\right]^
    {(\gamma+1)/[2(\gamma-1)]}.

For :math:`A/A^*>1`, this relation has a subsonic and a supersonic root.
:func:`aerophysics.isentropic.mach_from_area_ratio` therefore requires
:class:`~aerophysics.isentropic.MachBranch` selection and solves only that
branch.

The dimensionless mass-flow parameter and mass flux are

.. math::

   \operatorname{MFP}(M)
   =\frac{\sqrt{\gamma}M}
          {F(M)^{(\gamma+1)/[2(\gamma-1)]}},
   \qquad
   \frac{\dot m}{A}
   =\frac{p_0}{\sqrt{RT_0}}\operatorname{MFP}(M).

:func:`~aerophysics.isentropic.mass_flow_parameter` returns
:math:`\operatorname{MFP}`, and :func:`~aerophysics.isentropic.choked_mass_flux`
returns the mass flux at :math:`M=1`, the largest value a duct can pass for
the given reservoir state.
These calorically perfect relations follow
:ref:`NACA Report 1135 <ref-naca-report-1135>`.

Fused multi-output analysis
---------------------------

Calling each relation separately is convenient when only one quantity is
needed. When a workflow needs ratios, area, mass flow, and absolute-state
values together, :func:`aerophysics.isentropic.isentropic_analysis` returns
them as one :class:`~aerophysics.isentropic.IsentropicAnalysis`. For gas
models that need a numerical solve, this avoids solving the same state more
than once.

>>> from aerophysics.isentropic import isentropic_analysis
>>> analysis = isentropic_analysis(
...     [0.0, 1.0, 2.0],
...     total_temperature=300.0,
...     total_pressure=101325.0,
... )
>>> [round(float(value), 6) for value in analysis.ratios.total_pressure_ratio]
[1.0, 1.892929, 7.824449]
>>> analysis.state is not None
True

The returned fields follow the broadcast shape of all supplied inputs.
Absolute state and mass-flux fields are ``None`` unless both reservoir
temperature and pressure are supplied. At Mach zero, the fused result records
the limiting area ratio as positive infinity; the standalone
:func:`~aerophysics.isentropic.area_ratio` continues to reject Mach zero.

Thermally perfect gas
---------------------

At high temperature the heat capacity of air rises, so the constant-:math:`\gamma`
formulas above lose accuracy. Passing a
:class:`~aerophysics.thermochemistry.ThermallyPerfectGas` to the same API uses
temperature-dependent enthalpy and entropy instead. The total temperature
:math:`T_0` must then be supplied, because the ratios now depend on the actual
temperature as well as on Mach number. The derivation is the same as above,
but :math:`h(T)` replaces :math:`c_pT`. Static temperature is the root of

.. math::

   h(T_0)-h(T)=\frac{M^2\gamma(T)RT}{2}.

At fixed frozen composition, the entropy and ideal-gas relations give

.. math::

   \frac{p_0}{p}
   =\exp\left[\frac{s^\circ(T_0)-s^\circ(T)}{R}\right],
   \qquad
   \frac{\rho_0}{\rho}=\frac{p_0}{p}\frac{T}{T_0}.

The temperature-dependent mass-flow parameter is

.. math::

   \operatorname{MFP}(M,T_0)
   =\frac{M\sqrt{\gamma(T)T_0/T}}{p_0/p}.

The sonic state at the same :math:`T_0` defines
:math:`\operatorname{MFP}^*`; consequently
:math:`A/A^*=\operatorname{MFP}^*/\operatorname{MFP}` and the choked mass
flux remains the value at :math:`M=1`. Ratio and area inverses use bounded
numerical roots and retain the explicit subsonic/supersonic branch choice.

>>> from aerophysics import AIR_NASA9
>>> from aerophysics.isentropic import isentropic_ratios
>>> thermal = isentropic_ratios(
...     2.0,
...     AIR_NASA9,
...     total_temperature=1000.0,
...     allow_extrapolation=False,
... )
>>> round(thermal.total_temperature_ratio, 6)
1.72214
>>> round(thermal.total_pressure_ratio, 6)
7.894673

``AIR_NASA7`` and ``AIR_NASA9`` have a fitted range of 200--6000 K. The
isentropic API extrapolates the nearest polynomial region by default and emits
one :class:`~aerophysics.exceptions.ApplicabilityWarning` per public call when
the total, static, or critical temperature is outside that range. Pass
``allow_extrapolation=False`` for strict range enforcement. Extrapolated
values, especially far below 200 K, should be treated cautiously.
The independent thermally perfect flow reference used for verification is
:ref:`Witte and Tatum (1994) <ref-witte-tatum-1994>`.

Harmonic-oscillator gas
-----------------------

:class:`~aerophysics.real_gas.HarmonicOscillatorGas` (see
:doc:`gas_models`) is also thermally perfect: its heat capacity depends on
temperature but not on pressure. The thermally perfect equations above
therefore apply unchanged, and no reservoir pressure is needed. The
``AIR_HARMONIC_OSCILLATOR`` preset is documented for reservoir temperatures of
400--2000 K.

.. _beattie-bridgeman-isentropic:

Beattie--Bridgeman real gas
---------------------------

For the Beattie--Bridgeman gas (see :ref:`beattie-bridgeman-gas`), enthalpy,
entropy, and sound speed depend on pressure as well as temperature. For a
specified :math:`T_0,p_0,M`, the solver uses the low-density gas solution of
the equation of state and solves

.. math::

   s(T,p)=s(T_0,p_0),\qquad
   h(T_0,p_0)-h(T,p)=\frac{M^2a_s(T,p)^2}{2}.

Because the properties now depend on pressure, total pressure is required in
addition to total temperature. With :math:`G=\rho M a_s`, the generalized
mass-flow parameter is

.. math::

   \Phi=\frac{G\sqrt{RT_0}}{p_0},\qquad
   \frac{A}{A^*}=\frac{G^*}{G}.

The sonic state for the same reservoir defines :math:`G^*`. If the equation
of state has no physically meaningful gas solution at the requested state
(for example, non-positive heat capacity or sound speed), the call raises
:class:`~aerophysics.exceptions.ModelRangeError`.

>>> from aerophysics import AIR_BEATTIE_BRIDGEMAN
>>> from aerophysics.isentropic import isentropic_state
>>> real = isentropic_state(
...     2.0,
...     AIR_BEATTIE_BRIDGEMAN,
...     total_temperature=1200.0,
...     total_pressure=6.0e6,
...     allow_extrapolation=False,
... )
>>> round(real.static_pressure)
765146
>>> round(real.velocity, 3)
1056.954

Outside the preset's tabulated range (38.8889--1222.2222 K and
172.369 Pa--27.5790 MPa), the isentropic API emits one
:class:`~aerophysics.exceptions.ApplicabilityWarning` per call, or raises
:class:`~aerophysics.exceptions.ModelRangeError` when
``allow_extrapolation=False``.

The GUI **等エントロピー流れ** (isentropic flow) page offers ``AIR``,
``NASA7``, ``NASA9``, ``HARMONIC_OSCILLATOR``, and ``BEATTIE_BRIDGEMAN``.
``BEATTIE_BRIDGEMAN`` always requires reservoir pressure. Exported rows
include total and static temperature, static pressure and density, velocity,
sound speed, and dynamic pressure.

.. list-table:: Isentropic-flow symbols
   :header-rows: 1
   :widths: 16 28 40 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`M`
     - ``mach``
     - Mach number
     - dimensionless
   * - :math:`T,p,\rho`
     - static state
     - Static temperature, pressure, and density
     - K, Pa, kg/m³
   * - :math:`T_0,p_0,\rho_0`
     - ``total_*``
     - Isentropic total state
     - K, Pa, kg/m³
   * - :math:`A`
     - area
     - Local flow area
     - m²
   * - :math:`A^*`
     - critical area
     - Sonic area for the same mass flow
     - m²
   * - :math:`\dot m/A`
     - ``mass_flux``
     - Mass flow per unit area
     - kg/(m² s)

>>> from aerophysics.isentropic import isentropic_ratios
>>> ratios = isentropic_ratios(2.0)
>>> round(ratios.total_temperature_ratio, 6)
1.8
>>> round(ratios.total_pressure_ratio, 6)
7.824449
