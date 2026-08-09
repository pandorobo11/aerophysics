Compressible-flow relations
===========================

This page groups the ideal-gas relations used for isentropic flow, shock
waves, and Prandtl--Meyer expansions. Isentropic calculations support both
calorically and thermally perfect gases; the shock and expansion models remain
calorically perfect. State-ratio directions and angle conventions are stated
with each model so results can be compared without consulting a separate
conventions page.

Common assumptions and conventions
----------------------------------

The models assume steady flow. Isentropic relations additionally require
adiabatic flow without entropy production and may use constant heat capacities,
frozen-composition temperature-dependent heat capacities, or the
Beattie--Bridgeman dense-gas equation of state.
Shock and expansion state ratios are downstream over upstream; isentropic
state ratios are total over static. Angles are in radians. Convert explicitly
with :func:`aerophysics.units.degrees_to_radians` and
:func:`aerophysics.units.radians_to_degrees`.

Isentropic flow
---------------

Calorically perfect gas
^^^^^^^^^^^^^^^^^^^^^^^

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

Thermally perfect gas
^^^^^^^^^^^^^^^^^^^^^

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

Harmonic-oscillator thermally perfect gas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
JAXA: Kennard documents the statistical-mechanics relation, while
JAXA-RR-06-011 records the effective air constants and wind-tunnel use.

Beattie--Bridgeman real gas
^^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`~aerophysics.real_gas.BeattieBridgemanGas` adds density-dependent
corrections to the harmonic-oscillator caloric model.  Its equation of state
is evaluated in the form

.. math::

   p=\rho RT\left(1+E_1\rho+E_2\rho^2+E_3\rho^3\right),

.. math::

   E_1=B_0-\frac{A_0}{RT}-\frac{c}{T^3},\quad
   E_2=\frac{A_0a}{RT}-B_0b-\frac{B_0c}{T^3},\quad
   E_3=\frac{B_0bc}{T^3}.

Heat capacities and sound speed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

``AIR_BEATTIE_BRIDGEMAN`` uses the Randall/JAXA air constants and separate
N2 and O2 vibrational modes.  For a specified :math:`T_0,p_0,M`, the solver
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
The Beattie--Bridgeman air preset instead records the R. E. Randall,
*AEDC-TR-57-8* (1957) tabulated air-property range: 70--2200 degR
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

Normal shocks
-------------

For upstream Mach number :math:`M_1>1`,

.. math::

   M_2^2
   =\frac{1+\frac{\gamma-1}{2}M_1^2}
          {\gamma M_1^2-\frac{\gamma-1}{2}},

.. math::

   \frac{p_2}{p_1}
   =1+\frac{2\gamma}{\gamma+1}(M_1^2-1),
   \qquad
   \frac{\rho_2}{\rho_1}
   =\frac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2},
   \qquad
   \frac{T_2}{T_1}=\frac{p_2/p_1}{\rho_2/\rho_1}.

The downstream-to-upstream total-pressure ratio is

.. math::

   \frac{p_{02}}{p_{01}}
   =\left[\frac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^
       {\gamma/(\gamma-1)}
    \left[\frac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^
       {1/(\gamma-1)}.

:func:`aerophysics.shocks.supersonic_pitot_pressure_ratio` returns the
Rayleigh--Pitot ratio

.. math::

   \frac{p_{02}}{p_1}
   =\left[\frac{(\gamma+1)M_1^2}{2}\right]^{\gamma/(\gamma-1)}
    \left[\frac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^
       {1/(\gamma-1)}.

Use :func:`aerophysics.shocks.normal_shock` for the complete result.

Oblique shocks
--------------

The theta--beta--Mach relation is

.. math::

   \tan\theta
   =2\cot\beta\,
    \frac{M_1^2\sin^2\beta-1}
         {M_1^2(\gamma+\cos 2\beta)+2}.

The shock angle lies between the Mach angle
:math:`\mu_M=\sin^{-1}(1/M_1)` and :math:`\pi/2`. The normal component
reduces the state calculation to the normal-shock equations:

.. math::

   M_{n1}=M_1\sin\beta,
   \qquad
   M_2=\frac{M_{n2}}{\sin(\beta-\theta)}.

:func:`aerophysics.shocks.oblique_shock` defaults to the weak root. Select
:class:`~aerophysics.shocks.ShockBranch` explicitly when branch identity
matters. If :math:`\theta` exceeds the maximum attached-shock deflection,
:class:`~aerophysics.exceptions.NoAttachedShockError` is raised; the requested
solution is not silently replaced by a detached normal shock.

.. list-table:: Shock symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - Unit
   * - :math:`M_1,M_2`
     - ``upstream_mach``, ``downstream_mach``
     - Upstream and downstream Mach numbers
     - dimensionless
   * - :math:`M_{n1},M_{n2}`
     - ``upstream_normal_mach``, ``downstream_normal_mach``
     - Mach components normal to the shock
     - dimensionless
   * - :math:`\theta`
     - ``deflection_angle``
     - Flow-deflection angle
     - rad
   * - :math:`\beta`
     - ``shock_angle``
     - Shock angle from the upstream velocity
     - rad
   * - :math:`p_{01},p_{02}`
     - total pressure
     - Upstream and downstream total pressures
     - Pa

>>> from aerophysics import ShockBranch, oblique_shock
>>> from aerophysics.units import degrees_to_radians, radians_to_degrees
>>> shock = oblique_shock(
...     2.0, degrees_to_radians(10.0), branch=ShockBranch.WEAK
... )
>>> round(radians_to_degrees(shock.shock_angle), 3)
39.314

Conical shocks
--------------

For inviscid axisymmetric flow over a sharp circular cone at zero angle of
attack, the velocity between the shock and cone surface varies with polar
angle.  With radial and polar velocity components nondimensionalized by the
limiting velocity available from adiabatic expansion into a vacuum, the
Taylor--Maccoll equations are

.. math::

   \frac{dV_r}{d\theta}=V_\theta,

.. math::

   \frac{dV_\theta}{d\theta}
   =\frac{V_rV_\theta^2-a^2(2V_r+V_\theta\cot\theta)}
          {a^2-V_\theta^2},
   \qquad
   a^2=\frac{\gamma-1}{2}(1-V_r^2-V_\theta^2).

The Rankine--Hugoniot relations supply the velocity immediately behind a
trial shock angle :math:`\beta`.  Integration toward the axis locates the
cone surface where :math:`V_\theta=0`.  The weak attached solution is the
first shock angle above the Mach angle that produces the requested cone
half-angle :math:`\theta_c`.

:func:`aerophysics.shocks.conical_shock` returns the shock angle, Mach numbers
immediately behind the shock and at the cone surface, surface static-state
ratios over the free stream, and the post-shock/free-stream total-pressure
ratio.  :func:`aerophysics.shocks.maximum_attached_cone_angle` returns the
attached-shock limit.  A larger cone half-angle raises
:class:`~aerophysics.exceptions.NoAttachedShockError` rather than substituting
a detached-shock approximation.

The model assumes a calorically perfect gas, a sharp circular cone, zero angle
of attack, steady inviscid adiabatic flow, and an attached axisymmetric shock.
It does not model bluntness, viscosity, real-gas effects, or asymmetric cone
flow.

>>> from aerophysics import conical_shock
>>> cone = conical_shock(2.0, degrees_to_radians(10.0))
>>> round(radians_to_degrees(cone.shock_angle), 3)
31.206
>>> round(cone.surface_mach, 3)
1.834
>>> round(cone.surface_pressure_ratio, 3)
1.293

Detached shocks
---------------

The detached-shock correlations are implemented separately from the
Rankine--Hugoniot solvers in :mod:`aerophysics.detached_shock`. Two geometries
are explicit: an axisymmetric sphere or hemispherical nose, and a
two-dimensional cylindrical nose. In both cases :math:`R_n` is the nose
radius and :math:`\Delta` is the axial gap from the body vertex to the shock
vertex.

Ambrosio--Wortman standoff
^^^^^^^^^^^^^^^^^^^^^^^^^^

The Ambrosio--Wortman correlations used by
:func:`aerophysics.detached_shock.shock_standoff_distance` are

.. math::

   \frac{\Delta}{R_n}=0.143\exp\left(\frac{3.24}{M^2}\right)
   \quad\text{(sphere or hemispherical nose)},

.. math::

   \frac{\Delta}{R_n}=0.386\exp\left(\frac{4.67}{M^2}\right)
   \quad\text{(two-dimensional cylindrical nose)}.

Billig shock shape
^^^^^^^^^^^^^^^^^^

Billig gives the shock-vertex curvature radius

.. math::

   \frac{R_c}{R_n}=1.143\exp\left[
       \frac{0.54}{(M-1)^{1.2}}\right]
   \quad\text{(sphere or hemispherical nose)},

.. math::

   \frac{R_c}{R_n}=1.386\exp\left[
       \frac{1.8}{(M-1)^{0.75}}\right]
   \quad\text{(two-dimensional cylindrical nose)}.

For a hemispherical or cylindrical nose followed by a parallel afterbody,
:func:`aerophysics.detached_shock.billig_shock_shape` uses
:math:`\beta=\sin^{-1}(1/M)` and the hyperbola

.. math::

   x=R_n+\Delta-R_c\cot^2\beta
   \left[
     \sqrt{1+\frac{y^2\tan^2\beta}{R_c^2}}-1
   \right].

The nose-curvature center is the origin, positive :math:`x` points upstream,
the body vertex is at :math:`x=R_n`, and the shock vertex is at
:math:`x=R_n+\Delta`. Billig shape calculations deliberately use the
Ambrosio--Wortman value of :math:`\Delta`; changing the displayed Seiff model
does not change that shape convention. A one-dimensional transverse
coordinate array is appended as the final output axis, so broadcast Mach and
radius cases retain their case axes.

Seiff density-ratio standoff
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a sphere, Seiff's density-ratio relation is

.. math::

   \frac{\Delta}{R_n}=\frac{0.78}{\rho_2/\rho_1}.

:func:`aerophysics.detached_shock.seiff_standoff_distance` accepts
:math:`\rho_2/\rho_1` directly and therefore does not attach a Mach number to
its result. The convenience API
:func:`aerophysics.detached_shock.seiff_standoff_distance_from_mach` obtains
the density ratio from the existing calorically perfect-gas
:func:`aerophysics.shocks.normal_shock`. No cylindrical Seiff correlation is
claimed or accepted. :func:`aerophysics.detached_shock.compare_standoff_distances`
returns the Ambrosio--Wortman and Seiff sphere results together with their
signed and relative differences.

All functions require finite :math:`M>1` and :math:`R_n>0`; the low-level
Seiff function additionally requires finite :math:`\rho_2/\rho_1>1`. These
are physical and mathematical domains, not empirical fit limits, and invalid
values raise :class:`ValueError`. The cited original publications do not give
a sufficiently explicit numerical Mach fit range to justify inventing an
``ApplicabilityWarning`` or ``ModelRangeError`` boundary. NASA TN D-2780's
independent comparison over :math:`0.04<\rho_1/\rho_2<0.16` is recorded as a
verification interval, not as the full validity range of Seiff's correlation.

These engineering correlations assume continuum, steady, low-temperature
flow. Their common use is for calorically perfect air; Billig's fitted curves
are primarily associated with :math:`\gamma=1.4`. They do not solve the
shock-layer thermodynamics. Real-gas Seiff models, NASA7/NASA9 or harmonic-
oscillator general normal-shock solvers, Beattie--Bridgeman shock states,
rarefied-flow corrections, and shock fitting are outside this implementation.

>>> from aerophysics import DetachedShockGeometry, billig_shock_shape
>>> from aerophysics import compare_standoff_distances, shock_standoff_distance
>>> sphere = shock_standoff_distance(
...     4.0, 0.5, geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE
... )
>>> round(sphere.normalized_standoff_distance, 6)
0.175098
>>> shape = billig_shock_shape(
...     4.0, 0.5, [-1.0, 0.0, 1.0],
...     geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
... )
>>> shape.shock_x_coordinates[1] == sphere.nose_radius + sphere.standoff_distance
True
>>> compare_standoff_distances(4.0, 0.5).seiff.density_ratio > 1.0
True

The formulas and geometry definitions follow Ambrosio and Wortman
(DOI ``10.2514/8.5988``), Billig (DOI ``10.2514/3.28969``), and Seiff
(NASA SP-24); the Seiff cross-check follows Inouye, NASA TN D-2780.

Prandtl--Meyer expansions
-------------------------

For :math:`M\ge1`, the Prandtl--Meyer function is

.. math::

   \nu(M)
   =\sqrt{\frac{\gamma+1}{\gamma-1}}
    \tan^{-1}\left[
      \sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)}
    \right]
    -\tan^{-1}\left(\sqrt{M^2-1}\right).

Its finite limiting value is

.. math::

   \nu_\max=\frac{\pi}{2}
   \left(\sqrt{\frac{\gamma+1}{\gamma-1}}-1\right).

For a centered expansion through turn angle :math:`\delta`,

.. math::

   \nu(M_2)=\nu(M_1)+\delta.

:func:`aerophysics.expansion.mach_from_prandtl_meyer` solves this equation
numerically. Total temperature and total pressure remain constant, while

.. math::

   \frac{T_2}{T_1}=\frac{F(M_1)}{F(M_2)},
   \qquad
   \frac{p_2}{p_1}
   =\left(\frac{T_2}{T_1}\right)^{\gamma/(\gamma-1)},
   \qquad
   \frac{\rho_2}{\rho_1}
   =\left(\frac{T_2}{T_1}\right)^{1/(\gamma-1)}.

.. list-table:: Expansion symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - Unit
   * - :math:`\nu`
     - ``prandtl_meyer_angle``
     - Prandtl--Meyer angle
     - rad
   * - :math:`\delta`
     - ``turn_angle``
     - Flow turning angle
     - rad
   * - :math:`M_1,M_2`
     - ``upstream_mach``, ``downstream_mach``
     - Upstream and downstream Mach numbers
     - dimensionless

The input turn must be nonnegative and keep the downstream angle below
:math:`\nu_\max`. The complete state change is returned by
:func:`aerophysics.expansion.prandtl_meyer_expansion`. These governing
relations and reference values follow NACA Report 1135; see :doc:`references`.

>>> from aerophysics import prandtl_meyer_expansion
>>> expansion = prandtl_meyer_expansion(2.0, degrees_to_radians(10.0))
>>> round(expansion.downstream_mach, 6)
2.384887
