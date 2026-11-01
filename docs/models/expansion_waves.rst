Prandtl--Meyer expansion waves
==============================

This model assumes steady, inviscid, adiabatic, two-dimensional centered
expansion of an ideal gas with frozen composition. State ratios are downstream
over upstream, and angles are in
radians. Convert explicitly with
:func:`aerophysics.units.degrees_to_radians` and
:func:`aerophysics.units.radians_to_degrees`.

Physical background
-------------------

When a supersonic stream flows around a convex corner (a wall that turns
away from the flow), it does not turn suddenly. It turns gradually through a
fan of very weak waves called Mach waves, centred on the corner. Each Mach
wave is inclined to the local flow at the Mach angle

.. math::

   \mu=\sin^{-1}\frac{1}{M}.

Across each wave the flow turns slightly and speeds up, while pressure,
temperature, and density fall. Unlike a shock, this process is smooth and
loss-free, so total pressure and total temperature stay constant.

The **Prandtl--Meyer function** :math:`\nu(M)` is the angle through which a
sonic stream (:math:`M=1`) must turn to reach Mach number :math:`M`. For any
centred expansion, the turn angle is simply the difference between the
downstream and upstream values of :math:`\nu`. There is a largest possible
turn, :math:`\nu_\max`, at which the gas would expand to zero pressure.

Compression around a concave corner instead produces a shock; see
:doc:`shock_waves`.

Calorically perfect gas
-----------------------

For :math:`M\ge1`, the Prandtl--Meyer function is

.. math::

   \nu(M)
   =\sqrt{\frac{\gamma+1}{\gamma-1}}
    \tan^{-1}\left[
      \sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)}
    \right]
    -\tan^{-1}\left(\sqrt{M^2-1}\right).

Its finite limiting value, reached as :math:`M\to\infty`, is

.. math::

   \nu_\max=\frac{\pi}{2}
   \left(\sqrt{\frac{\gamma+1}{\gamma-1}}-1\right),

about 130.45 degrees for :math:`\gamma=1.4`.

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
relations and reference values follow
:ref:`NACA Report 1135 <ref-naca-report-1135>`.

>>> from aerophysics import prandtl_meyer_expansion
>>> from aerophysics.units import degrees_to_radians
>>> expansion = prandtl_meyer_expansion(2.0, degrees_to_radians(10.0))
>>> round(expansion.downstream_mach, 6)
2.384887

Thermally perfect gas
---------------------

Pass ``AIR_NASA7``, ``AIR_NASA9``, ``AIR_HARMONIC_OSCILLATOR``, or a custom
:class:`~aerophysics.thermochemistry.ThermallyPerfectGas` or
:class:`~aerophysics.real_gas.HarmonicOscillatorGas`, together with
``upstream_temperature`` in K. Mach, positive expansion turn in radians, and
temperature broadcast together. Without a gas argument, the constant-gamma
``AIR`` model is used. These models include variable heat capacity; they do
not include chemistry, relaxation, dissociation, or real-gas pressure effects.

When :math:`\gamma` varies with temperature, :math:`\nu` is no longer a
function of Mach number alone, so the turn must be integrated through the
fan. Along a Mach wave, the flow-turning relation and energy conservation give

.. math::

   d\theta=\sqrt{M^2-1}\,\frac{du}{u},\qquad
   h_0=h(T_1)+\frac{u_1^2}{2},\qquad
   u_1=M_1 a(T_1),

   u^2(T)=u_1^2+2[h(T_1)-h(T)],\qquad
   a^2(T)=\gamma(T)RT,\qquad M(T)=\frac{u(T)}{a(T)},

   \delta=\int_{T_2}^{T_1}
     \frac{c_p(T)}{u^2(T)}\sqrt{M^2(T)-1}\,dT.

The solver evaluates this integral numerically and finds the downstream
temperature :math:`T_2` that gives the requested turn :math:`\delta`. It does
not replace the variable :math:`\gamma(T)` by a fixed upstream value.
Pressure follows from constant entropy, with :math:`s^\circ(T)` the entropy at
a fixed reference pressure:

.. math::

   \frac{p_2}{p_1}=\exp\left[
       \frac{s^\circ(T_2)-s^\circ(T_1)}{R}\right],\qquad
   \frac{\rho_2}{\rho_1}=\frac{p_2/p_1}{T_2/T_1}.

Total enthalpy and total pressure remain constant. Variable-specific-heat
Prandtl--Meyer flow is treated in
:ref:`NACA TN 2125 <ref-noyes-1950>`.

The NASA presets cover **200--6000 K**, and the harmonic-oscillator air preset
covers **400--2000 K**. Every static state in the fan must stay in the gas's
``temperature_range``; no extrapolation is performed. A requested turn that
requires a downstream temperature below the lower bound raises
:class:`~aerophysics.exceptions.ModelRangeError`. Numerical quadrature or root
convergence failures raise
:class:`~aerophysics.exceptions.ExpansionConvergenceError` instead.

``available_turn_angle`` is the largest turn from this upstream state before
the gas reaches the model's minimum temperature. It is a limit of the gas
data, not the physical expansion-to-vacuum limit. For a calorically perfect
gas, or a custom harmonic-oscillator gas without a temperature range, it is
the vacuum limit (which itself cannot be reached).

The absolute angles :math:`\nu_1,\nu_2` are measured from the sonic state with
the same total enthalpy :math:`h_0`, so they depend on the upstream
temperature as well as on Mach number. That sonic state is hotter than the
upstream flow. If it lies above the gas's temperature range, the expansion
itself is still returned, but both absolute angles are ``None`` for a scalar,
or ``NaN`` at the affected array entries. The Mach-only helpers ``prandtl_meyer_angle``,
``mach_from_prandtl_meyer``, and ``maximum_prandtl_meyer_angle`` remain limited
to :class:`~aerophysics.gas.PerfectGas`.

>>> from aerophysics import AIR_HARMONIC_OSCILLATOR
>>> thermal = prandtl_meyer_expansion(
...     3.0, degrees_to_radians(10.0), AIR_HARMONIC_OSCILLATOR,
...     upstream_temperature=1000.0,
... )
>>> round(thermal.downstream_mach, 6)
3.495015
>>> round(1000.0 * thermal.static_temperature_ratio, 6)
815.57947
>>> thermal.upstream_prandtl_meyer_angle is None
True

The tests check the constant-heat-capacity limit for three values of
:math:`\gamma`, compare with an independent numerical solution for a test gas
whose heat capacity varies linearly with temperature, and confirm that energy
and entropy are conserved across NASA polynomial boundaries.
