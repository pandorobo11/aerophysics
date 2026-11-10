.. _flat-plate-compressibility-corrections:

High-speed flat-plate models
============================

This page gives the models used with ``compressible=True``. The symbols are
defined in :doc:`flat_plate_boundary_layer`.

Recovery temperature
--------------------

When fast-moving gas is slowed by friction near an insulated wall, the wall
heats to the **recovery temperature** :math:`T_r`. It is a little below the
total temperature, because heat conduction carries part of the energy away.
The recovery factor :math:`r` expresses this:

.. math::

   T_r=T_e\left[1+r\frac{\gamma-1}{2}M^2\right],
   \qquad
   r_t=Pr^{1/3}.

This recovery-factor expression applies to turbulent points. For laminar points,
the adiabatic wall temperature is solved directly with zero wall heat flux;
:math:`r_l=\sqrt{Pr}` is only an initial estimate. If ``wall_temperature`` is
omitted, the corresponding adiabatic wall is used. With ``compressible=True``,
laminar points use the compressible similarity solution.

Laminar flow: similarity solution
---------------------------------

With ``compressible=True``, laminar points solve the Illingworth/Levy--Lees equations
(:ref:`Illingworth, 1949 <ref-illingworth-1949>`) for laminar friction and
thicknesses, including points before and at transition:

.. math::

   (Cf'')'+ff''=0,\qquad
   (Cg'/Pr)'+fg'+(\gamma-1)M_e^2C(f'')^2=0,
   \qquad g=T/T_e,\quad C=\frac{\mu(gT_e)}{\mu_{\mathrm{model}}(T_e)g}.

The coordinate is :math:`z=\eta/\sqrt{2}` in the constant-property limit.
No slip, no penetration, and unit edge speed and temperature give
:math:`f(0)=f'(0)=0`, :math:`f'(\infty)=g(\infty)=1`.
An omitted wall temperature imposes the **true adiabatic** boundary condition
:math:`g'(0)=0`, with :math:`T_{aw}=g(0)T_e`. A specified wall imposes
:math:`g(0)=T_w/T_e`; a separate adiabatic solve supplies
``recovery_temperature`` even in this case. No square-root-Prandtl recovery
approximation is imposed as a boundary condition.

With :math:`u/U_e=f'`, physical lengths follow from
:math:`\mathrm{d}y=\sqrt{2}g x/\sqrt{Re_x}\,\mathrm{d}z`:

.. math::

   c_f\sqrt{Re_x}=\sqrt{2}C_w f''(0),\qquad
   \frac{\delta_{99}\sqrt{Re_x}}{x}=\sqrt{2}\int_0^{z_{99}}g\,\mathrm{d}z,

   \frac{\delta^*\sqrt{Re_x}}{x}=\sqrt{2}\int_0^\infty(g-f')\,\mathrm{d}z,
   \qquad
   \frac{\theta_m\sqrt{Re_x}}{x}=\sqrt{2}\int_0^\infty f'(1-f')\,\mathrm{d}z.

Here :math:`f'(z_{99})=0.99`. Friction and Reynolds number use the edge
state. These expressions initially use :math:`\mu_{\mathrm{model}}(T_e)`;
when supplied :math:`\mu_e` differs, all four coefficients multiply
:math:`\sqrt{\mu_{\mathrm{model}}(T_e)/\mu_e}`. The momentum integral is
checked independently against friction before enforcing the exact ZPG
identity :math:`\theta_m=xc_f`. The average coefficient is exactly
:math:`\bar C_f=2c_f`. The specified-transition start :math:`Re_0` uses this
same laminar momentum thickness at :math:`Re_{\mathrm{tr}}`. Thus
post-transition friction and drag depend on it; for the wall-temperature mode
its momentum-integral and profile thicknesses also respond.

``effective_reynolds_number`` returns the external :math:`Re_x` at laminar
points. At zero Mach and :math:`T_w=T_e`, the numerical coefficients approach
Blasius: 0.6641147, 4.90999, and 1.720788 for friction/momentum, delta99,
and displacement. ``compressible=False`` retains its rounded 0.664, 4.9100, and 1.7208.
Displacement thickness is signed: strong wall cooling can make its integral
zero or negative, and the result is never clipped.

Distance :math:`x` starts at the leading edge and :math:`y` points outward
from the wall. Thicknesses are in m and temperatures in K. Assumptions are a
perfect gas, constant cp and Pr, zero pressure gradient, and no induced
pressure. For :math:`T_{aw}\gtrsim2000\ \mathrm{K}`, changing specific heats
and dissociation make results indicative only; this is not a universal
chemical threshold or a real-gas validity claim.

The numerical solution agrees with an independent solver within 0.1% for Mach
0.3 to 15; see :doc:`../verification/flat_plate_boundary_layer`. States that
cannot be solved reliably raise
:class:`~aerophysics.exceptions.ModelRangeError`.

Viscous interaction warning
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In compressible modes the laminar points are checked using

.. math::

   \bar\chi=M_e^3\sqrt{\frac{C_w}{Re_{x,e}}},\qquad
   C_w=\frac{T_e}{T_w}\frac{\mu(T_w)}{\mu_e}.

Wall viscosity uses the selected transport model and the laminar wall
temperature (the solved adiabatic temperature when omitted). If any laminar
point has :math:`\bar\chi>1`, one
:class:`~aerophysics.exceptions.ApplicabilityWarning` reports the maximum
laminar value and the count of points exceeding one. The pressure induced
near the leading edge invalidates the assumed ZPG flat-plate estimate.
Turbulent points and ``compressible=False`` are excluded. This diagnostic does not correct
the edge pressure or solve the interaction.

As order-of-magnitude guidance from :ref:`Hayes and Probstein (1959)
<ref-hayes-probstein-1959>`, at :math:`\bar\chi\approx1` the weak-interaction
pressure rise is about 10% for a cold wall and 30% for an adiabatic wall;
weak and strong interaction meet around :math:`\bar\chi\approx3`.
These are not universal pressure corrections. For the strong-interaction
similarity problem see :ref:`Bush (1966) <ref-bush-1966>`.

Turbulent flow: Eckert reference temperature
--------------------------------------------

Eckert's idea is to keep the incompressible formulas but evaluate density and
viscosity at a single **reference temperature** :math:`T^*` that represents
conditions inside the boundary layer, somewhere between the edge and the wall.
For turbulent portions of ``TurbulentModel.ECKERT``,

.. math::

   T^*=0.22T_r+0.28T_e+0.50T_w,
   \qquad
   Re_x^*=Re_x\frac{T_e}{T^*}\frac{\mu_e}{\mu^*}.

Here :math:`\mu^*=\mu(T^*)` is evaluated with the selected dynamic-viscosity
model, while :math:`\mu_e` is the supplied ``edge_dynamic_viscosity`` used in
the definition of :math:`Re_x`. For consistent results, supply an edge
viscosity calculated with the same model as ``viscosity_model``.
The implemented reference-temperature expression is recorded by
:ref:`Glass and Hunt (1988) <ref-glass-hunt-1988>`.

The incompressible friction correlations evaluated at :math:`Re_x^*` give the
coefficient referred to the reference density
:math:`\rho^*=\rho_eT_e/T^*`, which follows from the perfect-gas law at the
constant edge pressure:

.. math::

   C_f^*=\frac{\tau_w}{\tfrac{1}{2}\rho^*U_e^2}=C_{f,\mathrm{inc}}(Re_x^*),
   \qquad
   C_{f,x}=\frac{\rho^*}{\rho_e}C_f^*=\frac{T_e}{T^*}C_{f,\mathrm{inc}}(Re_x^*).

The average coefficient is converted in the same way, because :math:`T^*` is
constant along a plate with constant edge and wall conditions. The returned
coefficients, ``wall_shear_stress``, and ``drag_per_unit_width`` therefore use
the edge dynamic-pressure definitions given above
(:ref:`White, 2006 <ref-white-2006>`). Laminar portions instead use the
similarity solution described below and return the edge Reynolds number.

All three turbulent thicknesses are ``NaN`` with this model, also after a
specified transition. Friction and drag are available.

Turbulent flow: Van Driest II
-----------------------------

Van Driest II maps the compressible turbulent boundary layer onto an
equivalent incompressible one. It is derived from a mixing-length model with
the density variation across the layer included, and it is generally more
accurate than a reference-temperature method for turbulent flow at high Mach
number. For turbulent ``TurbulentModel.VAN_DRIEST_II``, the
Hopkins--Inouye form first defines

.. math::

   r=Pr^{1/3},
   \qquad
   m=r\frac{\gamma-1}{2}M^2,
   \qquad
   T_r=T_e(1+m),
   \qquad
   F=\frac{T_w}{T_e}.

The angular factors are

.. math::

   \alpha=\frac{m-1+F}{\sqrt{(m+1+F)^2-4F}},
   \qquad
   \beta=\frac{m+1-F}{\sqrt{(m+1+F)^2-4F}}.

The transformation factors are

.. math::

   F_C=\frac{m}{[\sin^{-1}(\alpha)+\sin^{-1}(\beta)]^2},
   \qquad
   F_\theta=\frac{\mu_e}{\mu_w},
   \qquad
   F_x=\frac{F_\theta}{F_C},
   \qquad
   Re_{x,i}=F_xRe_x.

Here too :math:`\mu_e` is the supplied ``edge_dynamic_viscosity``;
:math:`\mu_w=\mu(T_w)` is evaluated by ``viscosity_model``.

These transformation factors follow the form stated by :ref:`Gnoffo, Berry,
and Van Norman (2011) <ref-gnoffo-berry-van-norman-2011>`. The method and its
flat-plate skin-friction application are documented by :ref:`Lee and Faget
(1956) <ref-lee-faget-1956>`, :ref:`Hopkins and Inouye (1971)
<ref-hopkins-inouye-1971>`, and :ref:`Hopkins (1972) <ref-hopkins-1972>`.

The equivalent incompressible local coefficient is the positive solution of

.. math::

   \frac{0.242}{\sqrt{C_{f,i}}}
   =0.41+\log_{10}(Re_{x,i}C_{f,i}),
   \qquad
   C_f=\frac{C_{f,i}}{F_C}.

The local implicit relation is also given by :ref:`Willems and Gülhan
(2013) <ref-willems-gulhan-2013>`.

The average coefficient uses the corresponding equation without the local
intercept:

.. math::

   \frac{0.242}{\sqrt{\bar C_{f,i}}}
   =\log_{10}(Re_{x,i}\bar C_{f,i}),
   \qquad
   \bar C_f=\frac{\bar C_{f,i}}{F_C}.

Both equations are solved numerically. If no solution can be found, the call
raises :class:`~aerophysics.exceptions.ModelRangeError`. At the
incompressible adiabatic limit :math:`M\to0`, all factors equal one.

Laminar portions use the similarity solution above. The turbulent
:math:`\delta_{99}`, :math:`\delta^*`, and :math:`\theta_m` are ``NaN``,
also after a specified transition.

Turbulent flow: Van Driest II with wall-temperature correction
--------------------------------------------------------------

Select ``TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE`` to obtain
friction and profile-based thicknesses for compressible turbulent flow.
This empirical extension combines the VD2 factors above with a wall-temperature
correction and a momentum-integral closure. The local edge state sets
:math:`Re_{x,e}=\rho_e U_e x/\mu_e`, and the same VD2 state gives
:math:`F_C`, :math:`F_x`, :math:`T_r` with :math:`r=Pr^{1/3}`, and
:math:`T_w` (equal to :math:`T_r` if omitted).

The dimensionless wall-temperature factor is

.. math::

   w=\left(\frac{T_w}{T_r}\right)^{0.064}, \qquad R=F_x Re_{x,e}.

The incompressible friction pair is the Schlichting relation and its exact
accumulated-drag derivative, from
:ref:`Schlichting and Gersten (2000) <ref-schlichting-gersten-2000>`:

.. math::

   \bar C_{F,i}=\frac{0.455}{(\log_{10}R)^{2.58}}, \qquad
   c_{f,i}=\bar C_{F,i}\left(1-\frac{2.58}{\ln R}\right),

   c_f=\frac{w c_{f,i}}{F_C}, \qquad
   \bar C_F=\frac{w\bar C_{F,i}}{F_C}.

Thus :math:`w=1` for an adiabatic wall. This model uses the Schlichting pair;
``TurbulentModel.VAN_DRIEST_II`` uses the Willems pair. At constant edge and wall
conditions, :math:`c_f=\mathrm{d}(x\bar C_F)/\mathrm{d}x`.
For ``TRANSITIONAL``, the laminar part uses the similarity solution and the
turbulent part the equivalent plate of :ref:`flat-plate-specified-transition`
with this friction pair.

The zero-pressure-gradient momentum integral gives the full-edge momentum
thickness and wall shear stress:

.. math::

   \theta=\frac{x\bar C_F}{2}, \qquad
   \tau_w=\frac{1}{2}\rho_e U_e^2 c_f.

Here :math:`\bar C_F` is the selected turbulent or specified-transition
average, on the
edge dynamic-pressure basis. Distances and thicknesses are in m, wall shear
is in Pa, and positive shear acts against the flow on the fluid. The
streamwise coordinate starts at the leading edge; the profile coordinate
:math:`y` points outward from the wall.

The 99-percent thickness is the value :math:`d` satisfying

.. math::

   \theta_{\mathrm{prof}}(d)
    =\int_0^d \frac{\rho u}{\rho_e U_e}
                  \left(1-\frac{u}{U_e}\right)\,\mathrm{d}y
    =0.978\theta.

The profile is
:func:`~aerophysics.boundary_layer_profile.compressible_turbulent_boundary_layer_profile`
with the Van Driest velocity transformation, Walz temperature relation,
:math:`\kappa=0.41`, :math:`B=5.2`, and a Coles wake parameter chosen so that
:math:`u(d)=0.99U_e`. See :doc:`compressible_velocity_transformations`,
:ref:`Coles (1956) <ref-coles-1956>`, :ref:`Spalding (1961)
<ref-spalding-1961>`, and the Walz relation reproduced by
:ref:`Zhang et al. (2014) <ref-zhang-bi-hussain-she-2014>`:

.. math::

   T=T_w+(T_r-T_w)\frac{u}{U_e}
       +(T_e-T_r)\left(\frac{u}{U_e}\right)^2, \qquad
   \rho=\rho_e\frac{T_e}{T}.

The returned ``boundary_layer_thickness`` is :math:`d=\delta_{99}`;
``displacement_thickness`` is
:math:`\int_0^d[1-\rho u/(\rho_eU_e)]\,\mathrm{d}y`, integrated only through
:math:`\delta_{99}`; ``momentum_thickness`` is the full :math:`\theta`, before
multiplication by 0.978. The wake must satisfy :math:`0\le\Pi\le1`.
If no admissible thickness exists, :class:`~aerophysics.exceptions.ModelRangeError`
reports the affected broadcast array index and which wake bound is exceeded.

Range and accuracy
~~~~~~~~~~~~~~~~~~

The model was compared with CFD, DNS, and measurements over
:math:`M_e=0.3`--13.3, :math:`T_w/T_{aw}=0.07`--1, and
:math:`Re_{x,e}\approx10^6`--:math:`10^9`, for a smooth plate without
pressure gradient. The main findings are:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Compared with
     - Typical difference
   * - CFD with the Spalart--Allmaras model
     - Friction within about 2% and momentum thickness within about 3.5% on
       average. :math:`\delta_{99}` within 2--3% RMS up to
       :math:`Re_{x,e}\approx10^8`; the difference then drifts to about −10%
       at :math:`10^9` on cold supersonic walls.
   * - Low-speed measurements and DNS
     - :math:`\delta_{99}` too thick by 19% at :math:`Re_x=10^6`, 15% at
       :math:`10^7`, and 10% at :math:`10^8`. Most of this comes from the
       Schlichting friction formula.
   * - High-speed DNS at Mach 2.5 to 14
     - Friction within 4%. :math:`\delta_{99}` within 5% up to Mach 6, but
       15% too thick at Mach 14.

These comparisons show where the model has been checked. They do not
guarantee every combination of Mach number, wall temperature, and Reynolds
number. On very cold walls at very high Reynolds number (for example Mach 10
with :math:`T_w/T_{aw}=0.07` near :math:`Re_{x,e}\approx1.5\times10^9`), no
admissible profile exists and the call raises ``ModelRangeError``.

The exponent 0.064 and the factor 0.978 were fitted to CFD results and are
not universal constants. The full comparisons are in
:doc:`../verification/flat_plate_boundary_layer`.

Supply **local edge conditions** rather than freestream conditions. Viscous
interaction at a sharp leading edge can raise the edge pressure 10--27%
above freestream pressure at Mach 10. Evaluating a correlation with local
states does not turn it into a streamwise variable-edge-state solver: the
momentum/drag formulas retain their constant-state, zero-pressure-gradient
assumption.

Why Eckert and Van Driest II give no thickness
----------------------------------------------

Earlier versions estimated turbulent thicknesses for ``ECKERT`` and
``VAN_DRIEST_II`` from simple incompressible ratios. Compared with CFD over
Mach 2 to 10, those estimates of :math:`\delta_{99}` were on average 58% too
thick for Eckert and 17% too thick for Van Driest II, and the displacement
and momentum thicknesses were wrong by factors of 2 to 10. These models now
return ``NaN`` for the turbulent thicknesses rather than misleading values.
Their friction and drag are unchanged.

Symbols
-------

.. list-table:: Compressibility symbols
   :header-rows: 1
   :widths: 17 29 38 16

   * - Symbol
     - API name
     - Meaning
     - Unit
   * - :math:`T_e`
     - ``edge_temperature``
     - Edge temperature
     - K
   * - :math:`T_w`
     - ``wall_temperature``
     - Specified or adiabatic wall temperature
     - K
   * - :math:`T_r`
     - ``recovery_temperature``
     - Recovery temperature
     - K
   * - :math:`T^*`
     - internal
     - Eckert reference temperature
     - K
   * - :math:`Re_{x,i}`
     - ``effective_reynolds_number``
     - Equivalent incompressible VD2 friction Reynolds number (from
       :math:`Re_\mathrm{eq}` after a specified transition)
     - dimensionless
   * - :math:`F_C,F_\theta,F_x`
     - internal
     - Van Driest II transformation factors
     - dimensionless
