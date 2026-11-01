Conical shocks
==============

State ratios are downstream over upstream, and angles are in radians. The
conservation laws and the physical meaning of a shock are explained in
:doc:`shock_waves`.

Sharp cone at zero incidence
----------------------------

A sharp circular cone at zero angle of attack carries an attached conical
shock. Unlike the wedge case, the flow behind a conical shock is not
uniform. Every streamline must keep turning to follow the cone, so the flow
continues to compress smoothly (and without further loss) between the shock
and the surface. The surface pressure is therefore higher than the pressure
just behind the shock.

The flow depends only on the polar angle :math:`\theta` measured from the
cone axis. With radial and polar velocity components divided by the limiting
velocity (the speed reached by expanding the gas into a vacuum), the
Taylor--Maccoll equations are

.. math::

   \frac{dV_r}{d\theta}=V_\theta,

.. math::

   \frac{dV_\theta}{d\theta}
   =\frac{V_rV_\theta^2-a^2(2V_r+V_\theta\cot\theta)}
          {a^2-V_\theta^2},
   \qquad
   a^2=\frac{\gamma-1}{2}(1-V_r^2-V_\theta^2).

The calculation picks a trial shock angle :math:`\beta`, uses the
Rankine--Hugoniot relations to find the velocity just behind the shock, and
integrates toward the axis. The cone surface is where the flow is parallel to
the surface, :math:`V_\theta=0`. The weak attached solution is the first
shock angle above the Mach angle that produces the requested cone half-angle
:math:`\theta_c`.

:func:`aerophysics.shocks.conical_shock` returns the shock angle, the Mach
numbers just behind the shock and at the cone surface, surface static-state
ratios over the free stream, and the total-pressure ratio across the shock.
:func:`aerophysics.shocks.maximum_attached_cone_angle` returns the largest
cone half-angle with an attached shock. A larger cone half-angle raises
:class:`~aerophysics.exceptions.NoAttachedShockError` rather than substituting
a detached-shock approximation.

Inputs are the upstream Mach number :math:`M_\infty>1` and the cone
half-angle in radians, :math:`0\le\theta_c<\pi/2`. For :math:`\theta_c>0` the
shock angle satisfies :math:`\theta_c<\beta<\pi/2`. For :math:`\theta_c=0`
it equals the Mach angle. As :math:`M_\infty\to1^+`, the attached cone-angle
limit approaches :math:`0.7792\sqrt{M_\infty-1}`. For extremely slender
cones (below about :math:`10^{-3}` rad) or Mach numbers within about
:math:`10^{-6}` of one, the shock lies so close to the Mach wave that double
precision may not resolve it. Such cases raise
:class:`~aerophysics.exceptions.ShockConvergenceError` instead of returning an
inaccurate state.

The model assumes a sharp circular cone, zero angle of attack, steady
inviscid adiabatic flow, and an attached axisymmetric shock. It does not
model nose bluntness, viscosity, dense-gas effects, or asymmetric cone flow.
Without a gas argument, the calorically perfect ``AIR`` model is used.
Reference solutions are tabulated by :ref:`Sims (1964) <ref-sims-1964>`.

>>> from aerophysics import conical_shock
>>> from aerophysics.units import degrees_to_radians, radians_to_degrees
>>> cone = conical_shock(2.0, degrees_to_radians(10.0))
>>> round(radians_to_degrees(cone.shock_angle), 3)
31.206
>>> round(cone.surface_mach, 3)
1.834
>>> round(cone.surface_pressure_ratio, 3)
1.293

Thermally perfect conical shocks
--------------------------------

For a frozen thermally perfect gas, pass ``upstream_temperature`` in kelvin
to ``conical_shock`` or ``maximum_attached_cone_angle``. Mach number, cone
half-angle, and temperature broadcast together. The velocity equations keep
the same form, but velocities are divided by the free-stream speed
:math:`U_\infty`, and temperature comes from the energy equation with the
gas's own enthalpy:

.. math::

   h(T)=h(T_2)+\frac{U_\infty^2}{2}
      \left(V_{r2}^2+V_{\theta2}^2-V_r^2-V_\theta^2\right),
   \qquad a^2=\frac{\gamma(T)RT}{U_\infty^2},
   \qquad \gamma(T)=\frac{c_p(T)}{c_p(T)-R}.

Here the subscript 2 denotes the state just behind the shock. The angle
:math:`\theta` increases away from the cone axis, so :math:`V_\theta<0`
between the shock and the surface.

The shock sets the entropy, and the flow between the shock and the surface is
isentropic:

.. math::

   \frac{p_s}{p_2}
   =\exp\left[\frac{s^\circ(T_s)-s^\circ(T_2)}{R}\right],
   \qquad \frac{\rho_s}{\rho_\infty}
   =\frac{p_s/p_\infty}{T_s/T_\infty}.

The total-pressure ratio is the shock loss and does not change between the
shock and the surface. Every **static** state from the shock to the surface
must stay inside the gas temperature range; nothing is extrapolated. A valid
weak cone solution can still be returned when the normal shock or the
attached-limit angle is out of range. A requested solution or limit that
cannot be reached inside the range raises
:class:`~aerophysics.exceptions.ModelRangeError`. A cone beyond a resolved
physical maximum raises ``NoAttachedShockError``, and a numerical failure
raises :class:`~aerophysics.exceptions.ShockConvergenceError`. Only the weak
solution is available. Composition changes, chemical reactions, and
vibrational nonequilibrium are excluded.

>>> from aerophysics import AIR_HARMONIC_OSCILLATOR
>>> thermal_cone = conical_shock(
...     3.0, degrees_to_radians(10.0), gas=AIR_HARMONIC_OSCILLATOR,
...     upstream_temperature=1000.0,
... )
>>> round(radians_to_degrees(thermal_cone.shock_angle), 3)
21.626
>>> round(1000.0 * thermal_cone.surface_temperature_ratio, 3)
1110.946
