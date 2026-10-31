Normal and oblique shocks
=========================

State ratios are downstream over upstream, and angles are in radians. The
conservation laws and the physical meaning of a shock are explained in
:doc:`shock_waves`.

Normal shocks
-------------

For a calorically perfect gas (constant :math:`\gamma`), the conservation
equations can be solved in closed form. For upstream Mach number
:math:`M_1>1`,

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

The density ratio approaches :math:`(\gamma+1)/(\gamma-1)=6` for air as
:math:`M_1\to\infty`, while the pressure and temperature ratios keep growing.
The downstream-to-upstream total-pressure ratio is

.. math::

   \frac{p_{02}}{p_{01}}
   =\left[\frac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^
       {\gamma/(\gamma-1)}
    \left[\frac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^
       {1/(\gamma-1)}.

It equals one at :math:`M_1=1` and falls as the shock strengthens. The
entropy rise is :math:`\Delta s/R=-\ln(p_{02}/p_{01})`.

A pitot tube in supersonic flow has a normal shock in front of its opening.
The tube then measures the total pressure behind that shock, not the
free-stream total pressure. :func:`aerophysics.shocks.supersonic_pitot_pressure_ratio`
returns this Rayleigh--Pitot ratio, measured against upstream static
pressure:

.. math::

   \frac{p_{02}}{p_1}
   =\left[\frac{(\gamma+1)M_1^2}{2}\right]^{\gamma/(\gamma-1)}
    \left[\frac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^
       {1/(\gamma-1)}.

Use :func:`aerophysics.shocks.normal_shock` for the complete result.
The one-dimensional relations follow
:ref:`NACA Report 1135 <ref-naca-report-1135>`.

Oblique shocks
--------------

When a wall turns into the flow by the deflection angle :math:`\theta`, an
attached shock forms at the shock angle :math:`\beta`, measured from the
upstream velocity. The theta--beta--Mach relation is

.. math::

   \tan\theta
   =2\cot\beta\,
    \frac{M_1^2\sin^2\beta-1}
         {M_1^2(\gamma+\cos 2\beta)+2}.

The shock angle lies between the Mach angle
:math:`\mu_M=\sin^{-1}(1/M_1)` (an infinitely weak wave) and :math:`\pi/2`
(a normal shock). The velocity component along the shock does not change, so
the normal-shock equations apply to the normal Mach number:

.. math::

   M_{n1}=M_1\sin\beta,
   \qquad
   M_2=\frac{M_{n2}}{\sin(\beta-\theta)}.

For each deflection below a maximum value :math:`\theta_\max`, two shock
angles satisfy the relation:

* The **weak** solution has the smaller shock angle. It is the one normally
  seen on wedges and ramps, and the downstream flow is usually still
  supersonic.
* The **strong** solution has the larger shock angle and a subsonic
  downstream flow. It appears only when high downstream pressure forces it.

:func:`aerophysics.shocks.oblique_shock` defaults to the weak solution.
Select :class:`~aerophysics.shocks.ShockBranch` explicitly when the branch
matters. Above :math:`\theta_\max` no attached shock exists and the real
shock detaches from the body. In that case
:class:`~aerophysics.exceptions.NoAttachedShockError` is raised; the request is
not silently replaced by a normal shock.
:func:`aerophysics.shocks.maximum_attached_deflection` returns
:math:`\theta_\max`.
The theta--beta--Mach relation and branch convention follow
:ref:`NACA Report 1135 <ref-naca-report-1135>`.

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

Thermally perfect normal and oblique shocks
-------------------------------------------

At high temperature, the heat capacity of air increases as molecular
vibration becomes active. A constant :math:`\gamma` then overestimates the
temperature rise behind a strong shock. The solver can instead use a gas
whose heat capacity depends on temperature (see :doc:`gas_models`).

Pass a :class:`~aerophysics.thermochemistry.ThermallyPerfectGas` (including
``AIR_NASA7`` or ``AIR_NASA9``) or a
:class:`~aerophysics.real_gas.HarmonicOscillatorGas`, and provide
``upstream_temperature`` as the upstream **static** temperature in kelvin.
Temperature is needed because the ratios now depend on the actual
temperature, not on Mach number alone. Mach number, angle, and temperature
broadcast together. Without a gas argument, the constant-:math:`\gamma`
``AIR`` model is used.

The model assumes steady, inviscid, adiabatic, two-dimensional ideal-gas
flow with fixed composition, thermal equilibrium, and :math:`p=\rho RT`.
The Rankine--Hugoniot conditions above are solved with the gas's own
enthalpy :math:`h(T)`; the tangential velocity is unchanged. For a trial
shock angle, :math:`u_{n1}=M_1a(T_1)\sin\beta`, the conservation equations
give :math:`T_2` and :math:`\rho_2/\rho_1`, and the flow turn is

.. math::

   \theta=\beta-\tan^{-1}\left[
      \frac{\rho_1}{\rho_2}\tan\beta\right],\qquad
   a(T)^2=\gamma(T)RT,\qquad \gamma(T)=\frac{c_p(T)}{c_p(T)-R}.

``theta_from_shock_angle``, ``shock_angle``,
``maximum_attached_deflection``, ``normal_shock``, and ``oblique_shock`` all
accept the same thermal gas and static-temperature arguments.

Total enthalpy, composition, and total temperature do not change across the
shock, so the total-pressure loss follows directly from the entropy rise:

.. math::

   \Delta s=s^\circ(T_2)-s^\circ(T_1)-R\ln(p_2/p_1),\qquad
   p_{02}/p_{01}=\exp(-\Delta s/R).

This needs only the static temperatures on both sides, so the stagnation
temperature does not have to lie inside the gas model's range. The pitot
ratio :func:`~aerophysics.shocks.supersonic_pitot_pressure_ratio` is
different: it also needs the stagnation state. For a thermal gas, supply
``upstream_temperature``. The total temperature follows from
:math:`h(T_0)=h(T_1)+u_1^2/2`, and

.. math::

   \frac{p_{02}}{p_1}=\frac{p_{02}}{p_{01}}
   \exp\!\left(\frac{s^\circ(T_0)-s^\circ(T_1)}{R}\right).

The governing thermal shock relations follow
:ref:`Tatum (1996), NASA CR-4749 <ref-tatum-1996>`.
Verification includes constant-heat-capacity limits, conservation residuals,
and an independent variable-heat-capacity solution; see
:doc:`../verification/compressible_flow`.

Temperature ranges are inclusive. Both dry-air NASA presets cover
200--6000 K, and ``AIR_HARMONIC_OSCILLATOR`` covers 400--2000 K. A custom
mixture uses the overlap of its species' fitted ranges. The shock solver
never extrapolates: a required static state outside the range raises
:class:`~aerophysics.exceptions.ModelRangeError`. A weak shock can still be
returned when only the strong solution or the attached-shock limit is out of
range. In that case the unavailable value raises ``ModelRangeError``, which
means "outside the gas data", not "the shock has detached".

These ranges describe the thermodynamic data, not a guarantee that frozen
(non-reacting) air is physically realistic. Dissociation, ionization,
finite-rate chemistry, vibrational nonequilibrium, and dense-gas effects are
excluded. Whether they matter depends on pressure, composition, and how long
the gas stays hot, so interpret high-temperature results with these
assumptions in mind.

>>> from aerophysics import AIR_NASA9
>>> thermal = oblique_shock(
...     3.0, degrees_to_radians(20.0), gas=AIR_NASA9,
...     upstream_temperature=300.0,
... )
>>> round(radians_to_degrees(thermal.shock_angle), 3)
37.68
>>> round(300.0 * thermal.static_temperature_ratio, 3)
466.118
>>> round(thermal.downstream_mach, 6)
2.006457
>>> round(thermal.total_pressure_ratio, 6)
0.796621
