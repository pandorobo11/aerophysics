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

The inverse ratio functions solve these expressions algebraically for
:math:`M`. Their values at :math:`M=1` are also the total-to-critical ratios.

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

The choked mass flux is the value at :math:`M=1`.
These calorically perfect relations follow
:ref:`NACA Report 1135 <ref-naca-report-1135>`.

Fused multi-output analysis
---------------------------

Calling each relation independently is useful when only one quantity is
needed. When a workflow needs ratios, area, mass-flow, and absolute-state
values together, :func:`aerophysics.isentropic.isentropic_analysis` returns
them as one :class:`~aerophysics.isentropic.IsentropicAnalysis`. For numerical
gas models, the fused call solves each distinct requested state once and
shares one Mach-one critical state for each distinct reservoir condition.

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

Passing a :class:`~aerophysics.thermochemistry.ThermallyPerfectGas` to the
same API evaluates temperature-dependent enthalpy and entropy. The total
temperature :math:`T_0` must then be supplied because the ratios are no longer
functions of Mach number alone. Static temperature is the root of

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

Harmonic-oscillator thermally perfect gas
-----------------------------------------

:class:`~aerophysics.real_gas.HarmonicOscillatorGas` represents a
frozen-composition ideal gas with thermally excited harmonic-vibrational
modes.  Let :math:`w_j` be the weight of mode :math:`j`,
:math:`\theta_j` its characteristic temperature, and

.. math::

   x_j=\frac{\theta_j}{T}.

The three dimensionless vibrational functions used by the implementation are

.. math::

   \Psi_v(T)=\sum_j w_j
      \frac{x_j^2 e^{x_j}}{(e^{x_j}-1)^2},

.. math::

   \Omega_v(T)=\sum_j w_j\frac{x_j}{e^{x_j}-1},
   \qquad
   \Sigma_v(T)=\sum_j w_j\left[
      \frac{x_j}{e^{x_j}-1}-\ln\left(1-e^{-x_j}\right)\right].

Here :math:`\Psi_v`, :math:`\Omega_v`, and :math:`\Sigma_v` respectively
control vibrational heat capacity, energy, and entropy.  If
:math:`\gamma_b` is the base heat-capacity ratio before vibrational
excitation, the specific heats are

.. math::

   c_v(T)=\frac{R}{\gamma_b-1}+R\Psi_v(T),
   \qquad
   c_p(T)=c_v(T)+R,
   \qquad
   \gamma(T)=\frac{c_p(T)}{c_v(T)}.

With vibrational zero-point energy omitted, the internal energy and enthalpy
used by the model are

.. math::

   u(T)=\frac{RT}{\gamma_b-1}+RT\Omega_v(T),
   \qquad
   h(T)=u(T)+RT.

Entropy is defined up to an additive constant.  Relative to any reference
state :math:`T_r,p_r`, its change is

.. math::

   s(T,p)-s(T_r,p_r)
   =\frac{\gamma_bR}{\gamma_b-1}\ln\frac{T}{T_r}
    +R\left[\Sigma_v(T)-\Sigma_v(T_r)\right]
    -R\ln\frac{p}{p_r}.

The remaining state relations retain their ideal-gas form:

.. math::

   \rho=\frac{p}{RT},
   \qquad
   a_s(T)=\sqrt{\gamma(T)RT}.

Thus :math:`c_p`, :math:`c_v`, :math:`\gamma`, :math:`u`, :math:`h`, and
:math:`a_s` depend on temperature but not pressure.  Pressure affects density
and entropy only.  The thermally perfect isentropic equations above therefore
apply without a reservoir-pressure input.

>>> from aerophysics import AIR_HARMONIC_OSCILLATOR
>>> harmonic = AIR_HARMONIC_OSCILLATOR.state(1200.0, 6.0e6)
>>> round(harmonic.cp, 3), round(harmonic.cv, 3)
(1176.407, 889.354)
>>> round(harmonic.heat_capacity_ratio, 6)
1.322766
>>> round(harmonic.speed_of_sound, 3)
675.014

``AIR_HARMONIC_OSCILLATOR`` uses :math:`R=287.05287` J/(kg K), base
:math:`\gamma_b=1.4`, and one effective mode at 3055.56 K.  Its documented
reservoir range is 400--2000 K.  The model is named for its physics, not for
JAXA: :ref:`Kennard (1938) <ref-kennard-1938>` documents the
statistical-mechanics relation, while
:ref:`Watari's JAXA report <ref-watari-2007-report>` and
:ref:`proceedings paper <ref-watari-2007-proceedings>` record the effective
air constants and wind-tunnel use.

Beattie--Bridgeman real gas
---------------------------

:class:`~aerophysics.real_gas.BeattieBridgemanGas` adds density-dependent
corrections to the harmonic-oscillator caloric model.  Its equation of state
originates with :ref:`Beattie and Bridgeman (1928)
<ref-beattie-bridgeman-1928>` and is evaluated in the form

.. math::

   p=\rho RT\left(1+E_1\rho+E_2\rho^2+E_3\rho^3\right),

.. math::

   E_1=B_0-\frac{A_0}{RT}-\frac{c}{T^3},\quad
   E_2=\frac{A_0a}{RT}-B_0b-\frac{B_0c}{T^3},\quad
   E_3=\frac{B_0bc}{T^3}.

Heat capacities and sound speed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unlike a thermally perfect ideal gas, the Beattie--Bridgeman heat capacities
depend on both temperature and density.  Using :math:`\Psi_v(T)` and
:math:`\gamma_b` defined above, the constant-volume specific heat implemented
by :class:`~aerophysics.real_gas.BeattieBridgemanGas` is

.. math::

   c_v(T,\rho)
   =\frac{R}{\gamma_b-1}+R\Psi_v(T)
    +\frac{6Rc\rho}{T^3}
     \left(1+\frac{B_0\rho}{2}
              -\frac{B_0b\rho^2}{3}\right).

For compactness, define the isothermal pressure derivatives

.. math::

   p_\rho
   =\left(\frac{\partial p}{\partial\rho}\right)_T
   =RT\left(1+2E_1\rho+3E_2\rho^2+4E_3\rho^3\right),

.. math::

   p_T
   =\left(\frac{\partial p}{\partial T}\right)_\rho
   =\rho R\left(1+\frac{2c\rho}{T^3}\right)
      \left(1+B_0\rho-B_0b\rho^2\right).

The thermodynamic identity relating the two heat capacities then gives

.. math::

   c_p(T,\rho)
   =c_v(T,\rho)+\frac{T p_T^2}{\rho^2p_\rho},
   \qquad
   \gamma(T,\rho)=\frac{c_p(T,\rho)}{c_v(T,\rho)}.

The frozen-composition isentropic sound speed :math:`a_s` is

.. math::

   a_s^2
   =\left(\frac{\partial p}{\partial\rho}\right)_s
   =\frac{c_p}{c_v}p_\rho.

Thus pressure affects :math:`c_p`, :math:`c_v`, and :math:`a_s` through the
gas-branch density obtained from the equation of state.  Setting
:math:`A_0=B_0=a=b=c=0` gives :math:`p_\rho=RT`, :math:`p_T=\rho R`, and
:math:`c_p=c_v+R`; the model reduces exactly to the corresponding
harmonic-oscillator ideal gas.

>>> from aerophysics import AIR_BEATTIE_BRIDGEMAN
>>> reservoir = AIR_BEATTIE_BRIDGEMAN.state(1200.0, 6.0e6)
>>> round(reservoir.cp, 3), round(reservoir.cv, 3)
(1173.815, 882.346)
>>> round(reservoir.heat_capacity_ratio, 6)
1.330334
>>> round(reservoir.speed_of_sound, 3)
690.207

``AIR_BEATTIE_BRIDGEMAN`` uses the air constants reported by
:ref:`Randall (1957) <ref-randall-1957-air>` and applied by
:ref:`Watari (2007) <ref-watari-2007-report>`, with separate N2 and O2
vibrational modes.  The general caloric derivation is recorded in
:ref:`Randall's companion report <ref-randall-1957-gases>`. For a specified
:math:`T_0,p_0,M`, the solver
restricts the entropy solve to the low-density mechanically stable gas branch
(below the first :math:`(\partial p/\partial\rho)_T=0` spinodal) and solves

.. math::

   s(T,p)=s(T_0,p_0),\qquad
   h(T_0,p_0)-h(T,p)=\frac{M^2a_s(T,p)^2}{2}.

The total pressure is therefore required in addition to total temperature.
With :math:`G=\rho M a_s`, the generalized mass-flow parameter is

.. math::

   \Phi=\frac{G\sqrt{RT_0}}{p_0},\qquad
   \frac{A}{A^*}=\frac{G^*}{G}.

The sonic state for the same reservoir defines :math:`G^*`.  Non-positive
compressibility, heat capacity, or sound speed, and absence of a stable gas
root raise :class:`~aerophysics.exceptions.ModelRangeError`.

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

The harmonic-oscillator air preset's documented reservoir range is 400--2000 K.
The Beattie--Bridgeman air preset instead records the
:ref:`R. E. Randall, AEDC-TR-57-8 (1957) <ref-randall-1957-air>` tabulated
air-property range: 70--2200 degR
(38.8889--1222.2222 K) and 0.025--4000 psia (172.369 Pa--27.5790 MPa).
That is a tabulation range, not an inherent physical-validity limit of the
equation of state.  The isentropic API emits one
:class:`~aerophysics.exceptions.ApplicabilityWarning` per public call outside
that range, or raises in strict mode.  These are frozen-composition models:
dissociation, chemical equilibrium, ionisation, condensation, and phase
changes are outside their scope.

The GUI exposes ``AIR``, ``NASA7``, ``NASA9``, ``HARMONIC_OSCILLATOR``, and
``BEATTIE_BRIDGEMAN`` as mutually exclusive choices.  The latter always
requires reservoir pressure.  Exported rows include total/static temperature,
static pressure and density, velocity, sound speed, and dynamic pressure.

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
