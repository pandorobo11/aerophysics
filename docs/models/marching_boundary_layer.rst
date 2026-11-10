.. _marching-boundary-layer-model:

Marching boundary-layer model
=============================

This page describes the laminar and SA-noft2 boundary-layer solver for edge
conditions that change along a wall, including its equations, assumptions,
and API mapping. For worked examples, see
:doc:`../guides/marching_boundary_layers`.

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

The wall-normal grid is sized from the Reynolds number before the march
starts. A strong compression can thicken the layer beyond that estimate; the
march then stops with ``BoundaryLayerConvergenceError`` instead of returning a
truncated layer. The numerical method is recorded in the developer notes
(``devdocs/numerical_implementation.md``).

Specified transition
~~~~~~~~~~~~~~~~~~~~

Use ``regime=BoundaryLayerRegime.TRANSITIONAL`` with an explicit
``transition_distance`` in metres along the wall from the leading edge.
The distance must be a finite scalar with
``0 < transition_distance < edge.distance[-1]``. There is no default
transition position; specifying it in another regime raises ``ValueError``.
``edge_sa_ratio`` must be positive: production is proportional to
:math:`\widetilde\nu`, so an initially zero field cannot grow.

SA-noft2 is solved over the entire wall, exactly as for ``TURBULENT``. The
only change is that the production term :math:`c_{b1}\widetilde S n` is switched
off for :math:`s \le s_{tr}`. Destruction, diffusion, convection, and the outer
boundary condition remain active. The returned stations are the same as in the
other regimes.

This follows the production suppression described for FUN3D by
:ref:`Hildebrand et al. (2022) <ref-hildebrand-2022>`. It supplies a user-chosen
onset, not natural-transition prediction, a roughness trip, or an intermittency
model. **The transition length produced by SA-noft2 is a numerical model
response, not a prediction of physical transition length.** The limitations of
SA activation discussed by :ref:`Rumsey and Spalart (2009)
<ref-rumsey-spalart-2009>` still apply.

Upstream of the transition position, a little eddy viscosity still enters from
the edge, so the solution there is close to, but not identical to, the
``LAMINAR`` result. With the default ``edge_sa_ratio=3``, the largest
differences seen in the verification cases were about 0.4% in both Cf and wall
heat flux. These are measurements for those cases, not accuracy bounds.

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

``wall_radius=None`` gives the planar calculation. A positive
scalar radius represents a cylinder beginning at a sharp leading edge, such as
a thin-walled tube. An array matching ``edge.distance`` gives the local wall
distance from the symmetry axis, in metres. Its first value may be zero for a
sharp tip; all subsequent values must be strictly positive and finite.
The result has the same fields as in the planar case.

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

The cone relation is reproduced to within 6e-12 percent, and the
variable-radius laminar mapping to within 0.001 percent in Cf. The
:ref:`R&M 3633 comparison <rm3633-example>` reaches delta99/r_w=0.4774 and
differs from measured Cf by up to 0.0004443 (measurement accuracy +/-0.0002).
This comparison does not establish a physical validity range.

A sharp cone-to-cylinder shoulder is outside the thin-layer model. The
pressure-length warning does not catch every case in which the solution near a
shoulder is poorly resolved, so results close to a shoulder need separate
judgement even without a warning.

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
* No transition prediction, roughness, suction, blowing, transverse or
  longitudinal curvature terms, or chemical reactions, and no interaction with
  the outer flow. A transition position can be specified but is not predicted.
* Wall temperature changes over distances much longer than the boundary-layer
  thickness. Unlike a pressure change, a shorter change does not raise a warning.

SA-noft2 develops its own laminar start and transition near the leading edge.
This is not a transition prediction. Accuracy against experiment is not
established. Agreement with CFD that uses the same model does not establish it.
Wall-temperature changes as short as delta99 can be numerically resolved, but
neglect of streamwise diffusion then has no established accuracy.

Turbulent flow evaluates the viscosity model slightly above and below each
edge temperature. If an edge temperature lies at the very limit of a viscosity
model's range, that model can therefore raise or warn.

In strong expansions, or with a small ``edge_sa_ratio``, the outer part of
``nu_tilde_ratio`` can oscillate from one grid point to the next. In the cases
examined, wall quantities and delta99 were not affected: Cf changed by less
than 0.01% and delta99 by less than 0.06% when the outer grid was refined.
Treat the outer part of the ``nu_tilde_ratio`` profile itself as approximate.

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
