.. _marching-boundary-layer-model:

Marching boundary-layer model
=============================

This page describes the laminar and SA-noft2 boundary-layer solver for edge
conditions that change along a wall, including its equations, assumptions,
and API mapping. For a task-oriented example, see
:doc:`../guides/boundary_layers`.

Physical background
-------------------

The flat-plate models in :doc:`flat_plate_boundary_layer` assume that the flow
outside the boundary layer does not change along the wall. On most real shapes
it does: the outer flow speeds up through an expansion and slows down through a
compression.

A boundary layer remembers what happened upstream. Accelerating flow thins it
and raises the wall shear, while decelerating flow thickens it and lowers the
wall shear. If the deceleration is strong enough, the wall shear falls to zero
and the flow **separates**. Two stations with the same local outer flow can
therefore have different skin friction, because their upstream histories differ.

This solver follows that history. It starts at a sharp leading edge and marches
downstream, using each solved station as the upstream state of the next one.

The outer flow is an input. The solver does not calculate it or feed the
boundary layer back into it. Optional wall radii describe thin axisymmetric
geometry. The equations below first show the planar form.

Governing equations
-------------------

Distance :math:`s` is measured along the wall from the leading edge. The thin
boundary-layer equations are written in Levy--Lees variables
(:ref:`Levy, 1954 <ref-levy-1954>`; :ref:`Lees, 1956 <ref-lees-1956>`):

.. math::

   \xi(s) = \int_0^s \rho_e \mu_e U_e \, ds, \qquad
   \eta = \frac{U_e}{\sqrt{2\xi}} \int_0^y \rho \, dy, \qquad
   f' = \frac{u}{U_e}, \qquad g = \frac{T}{T_e}

With :math:`C = \rho\mu/(\rho_e\mu_e)`, momentum and energy become

.. math::

   \left(C f''\right)' + f f'' + \beta\left(g - f'^2\right)
   = 2\xi\left(f' f'_\xi - f'' f_\xi\right)

.. math::

   \left(\frac{C}{Pr} g'\right)' + f g' + (\gamma - 1) M_e^2 C f''^2
   = 2\xi\left(f' g_\xi - g' f_\xi\right)

where primes denote :math:`\partial/\partial\eta` and

.. math::

   \beta = 2\xi \frac{d \ln U_e}{d\xi}

is positive for accelerating flow. The energy equation has no explicit
:math:`\beta` term. For an isentropic edge the pressure work cancels the change
in edge enthalpy, in the same way that a fluid particle expanding isentropically
keeps a constant :math:`T/T_e`.

At the wall, :math:`f = f' = 0`, and :math:`g` is either the specified wall
temperature ratio or follows from a zero wall heat flux. At the outer boundary,
:math:`f' = g = 1`.

The laminar station solve uses an analytic Jacobian for the built-in Sutherland
and Keyes viscosity models. Other viscosity models and the SA equations use
finite-difference Jacobians. These linearizations solve the same discrete
equations with the same convergence tolerances; they do not change the physical
model or its applicability.

SA-noft2 transport
~~~~~~~~~~~~~~~~~~

For turbulent flow define :math:`n = \widetilde\nu/\nu_e`,
:math:`Y = \int_0^\eta g\,d\eta`, :math:`\chi = n/(C g^2)`, and
:math:`\epsilon = \chi f_{v1}`. Replace the momentum diffusivity by
:math:`C(1+\epsilon)` and the energy diffusivity by
:math:`C/Pr\,(1+\epsilon Pr/Pr_t)`. Dissipation uses :math:`C(1+\epsilon)` too.
The additional equation is

.. math::

   \frac{((Cg+n/g)n')'}{\sigma g}
   + \frac{c_{b2} n'^2}{\sigma g^2} + f n'
   + c_{b1}\widetilde S n - c_{w1} f_w\frac{n^2}{Y^2}
   - \Lambda f'n
   = 2\xi(f'n_\xi - n'f_\xi),

.. math::

   R = \frac{\sqrt{2\xi}}{\mu_e},\qquad
   \Lambda = 2\xi\frac{d\ln\nu_e}{d\xi},\qquad
   \widetilde S = \frac{R|f''|}{g} + \frac{n f_{v2}}{\kappa^2Y^2}.

The original model is due to
:ref:`Spalart and Allmaras, 1994 <ref-spalart-allmaras-1994>`.
The constants and smooth modified-vorticity limiter follow
:ref:`Allmaras, Johnson and Spalart (2012) <ref-allmaras-johnson-spalart-2012>`.
The trip term ft2 is omitted. The wall has :math:`n=0`; the outer boundary has
:math:`n_e(s)=n_e(0)\nu_e(0)/\nu_e(s)`, so dimensional modified kinematic
viscosity is conserved.
The adiabatic wall uses a half-cell energy balance with zero wall flux.

The five unknowns per node are ``(f, u, g, n, Y)``. A damped Newton iteration
falls back to relaxed Picard with nonnegative ``n`` clipping. If both fail the
streamwise step halves. The eta grid has first spacing 0.0007 and spacing
ratio 1.04. Before marching, define ``Re_equivalent = max(R)^2/2`` over the
output stations and ``scale = max(1, (Re_equivalent/1e8)^0.42)``. The extent is
``150 scale`` and the outer spacing cap is ``scale``. Up to Re_equivalent=1e8
this is the original 312-node grid. These are geometric spacings, not
geometrically positioned nodes. The outer-boundary check remains active:
a compression can grow beyond this estimate and stop the march.

Specified transition
~~~~~~~~~~~~~~~~~~~~

Use ``regime=BoundaryLayerRegime.TRANSITIONAL`` with an explicit
``transition_distance`` in metres along the wall from the leading edge.
The distance must be a finite scalar with
``0 < transition_distance < edge.distance[-1]``. There is no default
transition position; specifying it in another regime raises ``ValueError``.
``edge_sa_ratio`` must be positive: production is proportional to
:math:`\widetilde\nu`, so an initially zero field cannot grow.

SA-noft2 is solved over the entire domain on the same grid as ``TURBULENT``.
Only :math:`c_{b1}\widetilde S n` is multiplied by zero for
:math:`s \le s_{tr}` and by one downstream. Destruction, diffusion, convection,
and the outer boundary condition remain active. The implementation applies
the multiplier in the residual, its numerical Newton Jacobian, and Picard's
lagged source and source-based relaxation. No profile or grid is replaced at
the switch. The march passes through the transition position internally; the
returned stations remain exactly ``edge.distance[1:]``.

This follows the production suppression described for FUN3D by
:ref:`Hildebrand et al. (2022) <ref-hildebrand-2022>`. It supplies a user-chosen
onset, not natural-transition prediction, a roughness trip, or an intermittency
model. **The transition length produced by SA-noft2 is a numerical model
response, not a prediction of physical transition length.** The limitations of
SA activation discussed by :ref:`Rumsey and Spalart (2009)
<ref-rumsey-spalart-2009>` still apply.

The upstream solution retains small eddy viscosity supplied from the edge.
For ``edge_sa_ratio=3`` the measured maximum upstream departures from the
laminar solver are 0.374% in Cf and 0.356% in wall heat flux over the four
uniform-edge verification cases (M0.3, M2 cooled/adiabatic, M6 cooled).
These are case-specific numerical measurements, not physical accuracy bounds.
Reducing the seed reduces eddy viscosity but does not remove differences
between the two existing discretizations: at M2, Te=220 K and Tw=300 K,
``edge_sa_ratio=1e-3`` gives a 0.05304% upstream heat-flux difference from
``LAMINAR``. This difference is recorded as a discretization diagnostic, without
a pass/fail criterion. The laminar-limit check retains the V6 Cf target of 0.05%.
The transition implementation is checked separately against
``TURBULENT(edge_sa_ratio=0)`` on the same SA grid and with the same operators,
using a maximum relative-difference limit of 1e-8% for both Cf and heat flux.

The result includes read-only ``nu_tilde_ratio`` profiles. The ``y+ > 1``
warning considers only output stations strictly downstream of the specified
transition; all other applicability checks remain active throughout the flow.
Planar and thin axisymmetric geometry, perfect-gas and isentropic-edge
assumptions, and the existing separation and outer-boundary limits still apply.

Results
~~~~~~~

.. math::

   C_f = \frac{\tau_w}{\tfrac12 \rho_e U_e^2}
       = \frac{2 C_w \mu_e f''_w}{\sqrt{2\xi}}, \qquad
   q_w = \frac{C_w \rho_e \mu_e U_e c_p T_e}{Pr \sqrt{2\xi}} \, g'_w

Heat flux into the wall is positive. :math:`q_w` comes from the computed wall
temperature gradient, not from the Reynolds analogy that
``flat_plate_boundary_layer`` uses. With :math:`Y = \int_0^\eta g \, d\eta`,

.. math::

   y = \frac{\sqrt{2\xi}}{\rho_e U_e} Y, \qquad
   \delta^* = \frac{\sqrt{2\xi}}{\rho_e U_e} \int_0^{\eta_{max}} (g - f') \, d\eta,
   \qquad
   \theta = \frac{\sqrt{2\xi}}{\rho_e U_e} \int_0^{\eta_{max}} f' (1 - f') \, d\eta

:math:`\delta_{99}` is the first wall distance where :math:`u/U_e = 0.99`.

Thin axisymmetric geometry
--------------------------

``wall_radius=None`` retains the planar computation bit for bit. A positive
scalar radius represents a cylinder beginning at a sharp leading edge, such as
a thin-walled tube. An array matching ``edge.distance`` gives the local wall
distance from the symmetry axis, in metres. Its first value may be zero for a
sharp tip; all subsequent values must be strictly positive and finite.
The input is copied. The result dataclass is unchanged.

Radius is interpolated **piecewise linearly** in wall distance, preserving a
cone-to-cylinder shoulder when its station is included. Curved bodies require
nodes dense enough to resolve their radius with linear segments. Pressure still
uses PCHIP in log pressure. The quadrature panels split at every input station.
No explicit derivative of radius is used.

For :math:`r=r_w(s)` and :math:`j=1`, the Mangler / Levy--Lees variables are

.. math::

   \xi(s) = \int_0^s \rho_e\mu_e U_e r^2\,ds, \qquad
   \eta = \frac{U_e r}{\sqrt{2\xi}}\int_0^y\rho\,dy,

.. math::

   \beta = 2\xi\frac{d\ln U_e}{d\xi}
     = \frac{2\xi\,d\ln U_e/ds}{\rho_e\mu_e U_e r^2},\qquad
   \Lambda = \frac{2\xi\,d\ln\nu_e/ds}{\rho_e\mu_e U_e r^2},

.. math::

   R = \frac{\sqrt{2\xi}}{r\mu_e},\qquad
   y = \frac{\sqrt{2\xi}Y}{U_e r\rho_e},\qquad
   C_f = \frac{2 C_w f''_w}{R},\qquad
   q_w = \frac{C_w\rho_e U_e c_p T_e g'_w}{Pr R}.

The transformed PDEs above are unchanged. Only at :math:`\xi=0` we use the
sharp-tip limits :math:`\beta=\Lambda=R=0`, without dividing by zero radius.
The SA outer grid uses this axisymmetric :math:`R` in
``Re_equivalent = max(R)^2/2``. For a constant-edge cone it equals
:math:`Re_s/3`, not :math:`Re_s`.

Thicknesses retain their planar thin-layer definitions, without a
:math:`1+y\cos(\alpha)/r` factor:

.. math::

   \delta^* = \int\left(1-\frac{\rho u}{\rho_e U_e}\right)dy,\qquad
   \theta = \int\frac{\rho u}{\rho_e U_e}\left(1-\frac{u}{U_e}\right)dy.

The corresponding integral balances are

.. math::

   \frac{d(r\rho_e U_e^2\theta)}{ds}
     + r\rho_e U_e\delta^*\frac{dU_e}{ds}=r\tau_w,\qquad
   \frac{d(rF_H)}{ds}=-r q_w,\qquad
   F_H=\int\rho u(H-H_e)dy.

All quantities use SI. Distance follows the wall from the tip or sharp leading
edge; :math:`y` points outward normal to the wall. Cf uses local edge dynamic
pressure, heat flux into the wall is positive, and positive beta means
acceleration. This is a thin-layer, perfect-gas, isentropic-edge model. It
neglects transverse curvature, longitudinal curvature and normal pressure
gradients; it does not cover blunt-nose entropy layers or separated flow.

An ``ApplicabilityWarning`` reports the largest :math:`\delta_{99}/r_w` and
its output station when the ratio exceeds 0.1. This threshold is an
**uncalibrated order-of-magnitude indicator**, not an accuracy limit.
Near a sharp tip the ratio grows without bound, so warnings are expected for
stations near the tip. The :math:`\delta_{99}|d\ln p/ds|` thin-layer warning
and all other existing warnings remain active.

The laminar constant-edge cone has local Cf and heat flux equal to
:math:`\sqrt{3}` times their plate values at the same wall distance; all
thicknesses are divided by :math:`\sqrt{3}`. Adiabatic wall temperature is
unchanged. The transformation itself yields these relations analytically.
The original :ref:`Mangler (1948) <ref-mangler-1948>` paper was not obtained.
The relation is confirmed through :ref:`Tauber, NASA TP-2914, equation (47),
page 9 <ref-tauber-1989>`, :ref:`Seiff, NACA TN 3284, equation (C13), page 24
<ref-seiff-1954>` (mean-friction ratio :math:`2/\sqrt{3}`), and
:ref:`Moore, NACA TN 2521 <ref-moore-1951>`.

For variable radius, the laminar planar equivalent uses
:math:`X=\int_0^s(r/L)^2ds`. Thus :math:`C_{f,a}=C_{f,p}r/L`,
:math:`q_a=q_pr/L` and :math:`\theta_a=\theta_pL/r`.
SA is not invariant under this mapping: its source contains :math:`R`, and
:math:`R_a=(L/r)R_p`. No exact turbulent cone-to-plate ratio is imposed.

The measured cone relation errors are below 6e-12 percent. With 801 nodes,
the variable-radius laminar mapping differs by at most 0.000827 percent in Cf;
the interpolation is performed in different distance variables. Axisymmetric
momentum residual medians are below 0.05 percent. Streamwise and wall-normal
refinement are recorded separately in the development memo, section 12.
The :ref:`R&M 3633 worked example <rm3633-example>` reaches delta99/r_w=0.4774
and differs from measured Cf by up to 0.0004443 (measurement accuracy
+/-0.0002). This comparison does not establish a physical validity range.

The synthetic cone-to-cylinder shoulder stress test uses a tanh pressure drop;
it is **not a physical shoulder edge condition**. With input nodes locally spaced
at most width/5, the pressure-length warning fires for SA at width 0.001 m
(maximum delta99 times abs(d ln p/ds) 1.24051--1.26509), and for both regimes
at width 0.0002 m (1.33059--6.29451). Widths 0.02 and 0.005 m do not trigger
that warning. All 16 cases completed in 2.82--8.85 s in the recorded run,
with zero Picard fallbacks. This is a timing observation, not a performance bound.
The curvature maxima occur near the tip or downstream end, so that warning
alone does not identify a shoulder problem. At s=0.31 m, relative SA-variable
oscillations reach 1.46404%; the final station has no evaluable extrema triplets.
Some oscillations occur without a pressure-length warning. These applicability
warnings therefore do not establish local profile convergence or physical
accuracy; the need for an additional oscillation diagnostic remains open.

Edge state
----------

The equations require :math:`dp_e = -\rho_e U_e \, dU_e` and a constant total
enthalpy along the edge, so pressure, temperature, and velocity are not
independent. ``BoundaryLayerEdge.from_pressure`` takes the edge pressure and the
edge total pressure and temperature, then calculates the rest from isentropic
relations. Downstream of a shock, the total pressure is the post-shock value.

Between the supplied stations, :math:`\ln p_e` follows a monotone piecewise cubic
(PCHIP) interpolation. Supply a smooth, well-resolved pressure distribution:
:math:`\beta` depends on the pressure gradient, so noise in the pressure appears
directly in the solution.

Assumptions and range
---------------------

* Steady planar or thin axisymmetric laminar or SA-noft2 attached flow from a
  sharp tip or sharp leading edge.
* Calorically perfect gas with constant :math:`Pr`, and pressure constant across
  the layer.
* Isentropic edge flow with constant total enthalpy. Entropy layers from blunt
  noses are not represented.
* No specified transition, roughness, suction, blowing, transverse or longitudinal curvature terms, or chemical
  reactions, and no interaction with the outer flow.
* Wall temperature changes over distances much longer than the boundary-layer
  thickness. Unlike a pressure change, a shorter change does not raise a warning.

``BoundaryLayerRegime.TRANSITIONAL`` raises ``ModelRangeError``.
SA-noft2 develops its own laminar start and transition near the leading edge.
This is not a transition prediction. Accuracy against experiment is not
established. Agreement with CFD that uses the same model does not establish it.
Wall-temperature changes as short as delta99 can be numerically resolved, but
neglect of streamwise diffusion then has no established accuracy.

Turbulent flow requires the viscosity model to accept ``T_e * exp(+/-1e-5)``
for an edge viscosity derivative, including when ``edge_sa_ratio=0``. Near a
temperature-range limit, a model that raises outside its range can fail at
these temperatures, and a model that warns outside its range can warn.
Laminar flow does not evaluate this SA derivative.

The outer part of ``nu_tilde_ratio`` can oscillate and is not grid-converged
in strong expansions or with a small ``edge_sa_ratio``. In the measured
cases, relative oscillation amplitudes reached about 100% on the default
grid and about 60% after outer-grid refinement, measured as half the excursion
between neighbouring extrema divided by their mean. Wall quantities and
delta99 were insensitive to that refinement in these cases: Cf changed by
less than 0.01% and delta99 by less than 0.06%. These measurements are not
general bounds.

A strong expansion can take tens of seconds because of Picard fallbacks;
M1.2 to 4 took about 25 s in a measurement on one machine.

.. list-table:: Reported conditions
   :header-rows: 1
   :widths: 40 60

   * - Condition
     - Report
   * - Invalid input or an edge state that is not isentropic
     - ``ValueError``
   * - The wall shear collapses in an adverse pressure gradient
     - ``BoundaryLayerSeparationError``, holding ``separation_distance`` and the
       ``result`` at the stations upstream. The march cannot continue because
       the thin-layer equations are not valid past separation.
   * - The solution cannot be resolved or the layer reaches the outer grid boundary
     - ``BoundaryLayerConvergenceError``. The outer grid is sized from the
       Reynolds number before the march. Strong compression at high Reynolds
       number can therefore stop with this error, as measured for M3 to 1.5
       at :math:`Re_s \simeq 10^9`.
   * - :math:`\delta_{99}\,|d \ln p_e/ds| > 1`
     - ``ApplicabilityWarning``. The pressure changes too quickly for the
       thin-layer approximation, as near a sharp corner, where viscous
       interaction with the outer flow matters.
   * - Turbulent first-node :math:`y^+ = y_1\sqrt{\rho_w\tau_w}/\mu_w > 1`
     - ``ApplicabilityWarning`` using local wall density and viscosity.
   * - :math:`\bar\chi = M_e^3 \sqrt{C_w/Re_s} > 1`
     - ``ApplicabilityWarning``. Leading-edge viscous interaction changes the
       pressure, which the supplied edge state may not include.

Verification
------------

:doc:`../verification/viscous_flow` records the comparisons. They use the
Falkner--Skan wedge flows (:ref:`Falkner and Skan, 1931
<ref-falkner-skan-1931>`), similar solutions with wall heat transfer in the form
of :ref:`Cohen and Reshotko, 1956 <ref-cohen-reshotko-1956>`, including their
Table II, Howarth's linearly retarded flow (:ref:`Howarth, 1938
<ref-howarth-1938>`, with the separation position of :ref:`Hartree, 1939
<ref-hartree-1939>`), and its compressible image under the transformation of :ref:`Illingworth, 1949
<ref-illingworth-1949>` and :ref:`Stewartson, 1949 <ref-stewartson-1949>`.

The frozen SA reference covers constant-edge Mach numbers 0.3--10 and
Reynolds numbers 1e6--1e8. ZPG marches at M0.3, M2, and M6 also complete at
Re_s=3e8, 1e9, and 1e10. At Re_s=1e9, comparison with an extended grid capped
at unit spacing changes Cf by less than 0.1% and heat flux by less than 0.5%.
Above 1e8 these comparisons establish numerical convergence only, without a
frozen reference or experimental validation. Completion at 1e10 does not
establish wall-normal convergence: the first-node y+ warning fires in the
M0.3 and M2 cases. M2--3 and M3--4 expansions, compression and separation,
and a low-Mach bump are also checked. These are numerical verification
conditions, not an experimentally established validity range.

API mapping
-----------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Symbol
     - Field of ``MarchingBoundaryLayerResult``
   * - :math:`s`
     - ``distance``
   * - :math:`\rho_e U_e s/\mu_e`
     - ``reynolds_number``
   * - :math:`\delta_{99},\ \delta^*,\ \theta`
     - ``boundary_layer_thickness``, ``displacement_thickness``,
       ``momentum_thickness``
   * - :math:`C_f,\ \tau_w`
     - ``local_skin_friction_coefficient``, ``wall_shear_stress``
   * - :math:`T_w,\ q_w`
     - ``wall_temperature``, ``wall_heat_flux``
   * - :math:`\beta`
     - ``pressure_gradient_parameter``
   * - :math:`n = \widetilde\nu/\nu_e`
     - ``nu_tilde_ratio`` (None for laminar flow)
   * - :math:`\eta,\ y,\ u/U_e,\ T/T_e`
     - ``eta``, ``wall_distance``, ``velocity_ratio``, ``temperature_ratio``

The leading edge is excluded from the result, because :math:`\tau_w` and
:math:`q_w` are singular there.
