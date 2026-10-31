Gas models
==========

A gas model supplies the equation of state, heat capacities, enthalpy,
entropy, and speed of sound used by the flow calculations. This page
describes every gas model in the package, from constant-property air to a
dense-gas equation of state. See :doc:`transport_properties` for viscosity
and conductivity, and :doc:`../getting_started/conventions` for which
calculations accept each model.

Physical background
-------------------

Gas models in compressible flow are usually described by three levels:

* A **calorically perfect** gas obeys :math:`p=\rho RT` and has constant
  :math:`c_p`, :math:`c_v`, and :math:`\gamma`. This is the ``AIR`` model.
  For air it is a good approximation below about 800 K, where molecular
  vibration is still largely inactive.
* A **thermally perfect** gas still obeys :math:`p=\rho RT`, but
  :math:`c_p` and :math:`c_v` depend on temperature. This captures the
  gradual activation of molecular vibration as air heats up.
* A **real** (imperfect) gas departs from :math:`p=\rho RT`, for example at
  high pressure. This is the Beattie--Bridgeman model below.

All models on this page keep a **frozen** composition: the mole fractions of N2, O2, Ar, and CO2 never change. In real air at about
one atmosphere, oxygen starts to dissociate above roughly 2000--2500 K, and
nitrogen above roughly 4000 K. Those chemical changes are not modelled.

Calorically perfect gas
-----------------------

:class:`~aerophysics.gas.PerfectGas` holds the specific gas constant
:math:`R` and heat-capacity ratio :math:`\gamma` constant. The implemented
relations are

.. math::

   p=\rho RT,
   \qquad
   c_p=\frac{\gamma R}{\gamma-1},
   \qquad
   c_v=\frac{R}{\gamma-1},
   \qquad
   a=\sqrt{\gamma RT}.

The speed of sound :math:`a` is the speed at which small pressure
disturbances travel through the gas. The Mach number :math:`M=V/a` compares
the flow speed with it. ``AIR`` uses
:math:`R=8314.32/28.9644\ \mathrm{J/(kg\,K)}` and
:math:`\gamma=1.4`. All temperatures must be positive.
The constants follow :ref:`U.S. Standard Atmosphere 1976
<ref-us-standard-atmosphere-1976>`.

.. list-table:: Gas-property symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`p`
     - ``pressure``
     - Static pressure
     - Pa
   * - :math:`\rho`
     - ``density``
     - Density
     - kg/m³
   * - :math:`T`
     - ``temperature``
     - Absolute temperature
     - K
   * - :math:`R`
     - ``specific_gas_constant``
     - Specific gas constant
     - J/(kg K)
   * - :math:`\gamma`
     - ``heat_capacity_ratio``
     - Ratio :math:`c_p/c_v`
     - dimensionless
   * - :math:`c_p,c_v`
     - ``cp``, ``cv``
     - Specific heat capacities
     - J/(kg K)
   * - :math:`a`
     - ``speed_of_sound``
     - Ideal-gas speed of sound
     - m/s

Thermally perfect gases (NASA polynomials)
------------------------------------------

``aerophysics`` provides NASA seven- and nine-coefficient polynomial
evaluators and a frozen-composition ideal-gas mixture model,
:class:`~aerophysics.thermochemistry.ThermallyPerfectGas`. These models make
heat capacity a function of temperature while keeping the ideal-gas equation
of state. NASA polynomials are curve fits of tabulated heat capacity,
enthalpy, and entropy for each species. They are widely used because they
are compact and accurate within their fitted temperature range.

The seven-coefficient form and data are documented by
:ref:`McBride, Gordon, and Reno (1993) <ref-mcbride-gordon-reno-1993>`; the
nine-coefficient form is documented by :ref:`McBride, Zehe, and Gordon (2002)
<ref-mcbride-zehe-gordon-2002>`.

NASA polynomial forms
^^^^^^^^^^^^^^^^^^^^^

For the NASA seven-coefficient form, each temperature region stores
:math:`a_1,\ldots,a_7` and evaluates

.. math::

   \frac{\bar c_p^\circ}{R_u}
   &= a_1+a_2T+a_3T^2+a_4T^3+a_5T^4,\\
   \frac{\bar h^\circ}{R_uT}
   &= a_1+\frac{a_2T}{2}+\frac{a_3T^2}{3}
      +\frac{a_4T^3}{4}+\frac{a_5T^4}{5}+\frac{a_6}{T},\\
   \frac{\bar s^\circ}{R_u}
   &= a_1\ln T+a_2T+\frac{a_3T^2}{2}
      +\frac{a_4T^3}{3}+\frac{a_5T^4}{4}+a_7.

The NASA nine-coefficient form adds inverse-temperature terms:

.. math::

   \frac{\bar c_p^\circ}{R_u}
   &= a_1T^{-2}+a_2T^{-1}+a_3+a_4T+a_5T^2+a_6T^3+a_7T^4,\\
   \frac{\bar h^\circ}{R_uT}
   &= -a_1T^{-2}+a_2\frac{\ln T}{T}+a_3+\frac{a_4T}{2}
      +\frac{a_5T^2}{3}+\frac{a_6T^3}{4}
      +\frac{a_7T^4}{5}+\frac{a_8}{T},\\
   \frac{\bar s^\circ}{R_u}
   &= -\frac{a_1T^{-2}}{2}-a_2T^{-1}+a_3\ln T+a_4T
      +\frac{a_5T^2}{2}+\frac{a_6T^3}{3}
      +\frac{a_7T^4}{4}+a_9.

``NASA7Polynomial`` and ``NASA9Polynomial`` return these three dimensionless
quantities through ``cp_over_r``, ``h_over_rt``, and ``s_over_r``. A shared
temperature boundary belongs to the lower-temperature region. Inputs outside
the fitted range raise :class:`~aerophysics.exceptions.ModelRangeError`.
Passing ``allow_extrapolation=True`` uses the nearest end region and emits
:class:`~aerophysics.exceptions.ApplicabilityWarning`.

Species and mixture properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An :class:`~aerophysics.thermochemistry.IdealGasSpecies` combines a polynomial
with its molar mass and standard-state pressure. Species methods return molar
SI properties. For fixed mole fractions :math:`x_i`, the mixture relations are

.. math::

   M &= \sum_i x_i M_i,
   &R &= \frac{R_u}{M},\\
   c_p &= \frac{\sum_i x_i\bar c_{p,i}^\circ}{M},
   &c_v &= c_p-R,\\
   \gamma &= \frac{c_p}{c_v},
   &a &= \sqrt{\gamma RT}.

The ideal-mixture entropy at pressure :math:`p` is

.. math::

   s(T,p)=\frac{1}{M}\sum_i x_i
   \left[
      \bar s_i^\circ(T)
      -R_u\ln\left(\frac{x_i p}{p_i^\circ}\right)
   \right].

NASA standard enthalpy includes formation enthalpy. It is therefore distinct
from the sensible enthalpy used for a temperature change. The
``standard_enthalpy`` and ``standard_internal_energy`` methods preserve the
NASA reference, while ``sensible_enthalpy`` and
``sensible_internal_energy`` subtract the value at an explicit reference
temperature.

Built-in frozen dry air
^^^^^^^^^^^^^^^^^^^^^^^

``AIR_NASA7`` and ``AIR_NASA9`` use normalized U.S. Standard Atmosphere dry-air
mole fractions for N2, O2, Ar, and CO2. Both coefficient sets have a common
fitted range of 200--6000 K.

``AIR_NASA7`` uses the :ref:`Cantera NASA gas data
<ref-cantera-nasa-gas-data>` and ``AIR_NASA9`` uses the
:ref:`NASA CEA data <ref-nasa-cea-data>` pinned by the package.

>>> from aerophysics import AIR_NASA9
>>> round(AIR_NASA9.cp(300.0), 3)
1004.829
>>> round(AIR_NASA9.cv(1000.0), 3)
853.978
>>> round(AIR_NASA9.heat_capacity_ratio(2000.0), 6)
1.297997
>>> round(AIR_NASA9.sensible_enthalpy(1000.0), 1)
747891.1

The two models are useful for comparing polynomial parameterizations. The
isentropic API accepts them together with total temperature. Normal, oblique,
and conical shocks and Prandtl--Meyer expansions accept them together with
upstream static temperature. All of these solve the flow with the gas's own
enthalpy and entropy. Substituting a local :math:`\gamma(T)` into a
constant-:math:`\gamma` formula is not equivalent and is not used.

At high temperature, real air dissociates and eventually ionizes. The built-in
models deliberately keep N2, O2, Ar, and CO2 mole fractions fixed through
6000 K, so their high-temperature results describe a frozen-composition
ideal-gas model rather than equilibrium air.

Harmonic-oscillator gas
-----------------------

:class:`~aerophysics.real_gas.HarmonicOscillatorGas` represents a
frozen-composition ideal gas whose molecules can vibrate. At room temperature
air molecules store energy only in translation and rotation, which gives
:math:`\gamma=1.4`. As temperature rises, molecular vibration starts to absorb
energy too, so the heat capacity increases and :math:`\gamma` falls. The model
treats each vibration as a quantum harmonic oscillator with a characteristic
temperature :math:`\theta_j`; a mode is mostly inactive for
:math:`T\ll\theta_j` and fully active for :math:`T\gg\theta_j`.

Let :math:`w_j` be the weight of mode :math:`j`,
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
reservoir range is 400--2000 K. :ref:`Kennard (1938) <ref-kennard-1938>`
documents the statistical-mechanics relation, while
:ref:`Watari's JAXA report <ref-watari-2007-report>` and
:ref:`proceedings paper <ref-watari-2007-proceedings>` record the effective
air constants and wind-tunnel use.

.. _beattie-bridgeman-gas:

Beattie--Bridgeman real gas
---------------------------

At high pressure, such as in a wind-tunnel reservoir at several megapascals,
gas molecules are close enough that their size and mutual attraction matter.
The ideal-gas law :math:`p=\rho RT` then loses accuracy.
:class:`~aerophysics.real_gas.BeattieBridgemanGas` adds density-dependent
corrections for these effects to the harmonic-oscillator model.  Its equation
of state
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
:ref:`Randall's companion report <ref-randall-1957-gases>`. The preset records
the :ref:`R. E. Randall, AEDC-TR-57-8 (1957) <ref-randall-1957-air>`
tabulated air-property range: 70--2200 degR (38.8889--1222.2222 K) and
0.025--4000 psia (172.369 Pa--27.5790 MPa). That is the range of the source
tables, not a strict physical limit of the equation of state. Like the other
models on this page, it keeps the composition frozen: dissociation, chemical
equilibrium, ionisation, condensation, and phase changes are outside its
scope. Its use in nozzle and reservoir calculations is described in
:ref:`beattie-bridgeman-isentropic`.
