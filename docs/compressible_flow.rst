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

The models assume steady ideal-gas flow. Isentropic relations additionally
require adiabatic flow without entropy production and may use either constant
heat capacities or frozen-composition temperature-dependent heat capacities.
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
The local GUI exposes ``AIR``, ``NASA7``, and ``NASA9`` as mutually exclusive
gas-model choices and includes total and static temperature in exported rows.

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
