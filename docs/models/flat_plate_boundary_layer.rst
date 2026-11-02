.. _flat-plate-boundary-layer-model:

Flat-plate boundary-layer model
===============================

This page describes the smooth flat-plate boundary-layer correlations,
including their equations, variables, applicability, and API mapping. For
task-oriented examples, see :doc:`../guides/boundary_layers`. Mean-velocity
profile transformations are documented separately in
:doc:`compressible_velocity_transformations`.

Physical background
-------------------

Air sticks to a solid surface: its velocity at the wall is zero (the no-slip
condition). The thin region in which the velocity rises from zero to the
outer-flow value :math:`U_e` is the **boundary layer**. Friction inside it
produces skin-friction drag, and it grows thicker with distance from the
leading edge.

A boundary layer starts **laminar**, with smooth, orderly layers. Farther
downstream, small disturbances grow and the flow becomes **turbulent**, with
strong mixing. Turbulent mixing brings fast fluid closer to the wall, so a
turbulent boundary layer has much higher skin friction and grows faster. The
change, called transition, typically happens between :math:`Re_x\approx
3\times10^5` and :math:`3\times10^6`, depending on surface finish, noise, and
free-stream turbulence. This package does not predict where it happens.

The flat plate with no pressure gradient is the simplest case and the
reference for many engineering estimates:

* The laminar case has an exact similarity solution by **Blasius**.
* Turbulent results are **empirical correlations** fitted to measurements.

At high Mach number, friction heats the gas near the wall. The hot gas has
lower density and higher viscosity, which changes skin friction. The
**compressibility corrections** account for this by evaluating the
incompressible formulas at suitably chosen reference conditions.

Correlation definition
----------------------

:func:`aerophysics.boundary_layer.flat_plate_boundary_layer` models one side
of a smooth flat plate with a sharp leading edge, zero pressure gradient, and
constant edge conditions. Distance :math:`x` is measured from the leading
edge, and drag is reported per unit spanwise width.

The standard incompressible definitions and correlations used below are
summarized by :ref:`Schlichting and Gersten (2000)
<ref-schlichting-gersten-2000>` and :ref:`White (2006) <ref-white-2006>`.

The caller selects
:class:`~aerophysics.boundary_layer.BoundaryLayerRegime` ``LAMINAR``,
``TURBULENT``, or ``TRANSITIONAL``. A transitional calculation requires a
``transition_reynolds`` value; the API does not infer natural transition.

Reynolds number and dynamic pressure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   Re_x=\frac{\rho_eU_ex}{\mu_e},
   \qquad
   q_e=\frac{1}{2}\rho_eU_e^2.

.. list-table:: Common flat-plate symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`x`
     - ``distance``
     - Distance from the leading edge
     - m
   * - :math:`U_e`
     - ``edge_velocity``
     - Boundary-layer edge velocity
     - m/s
   * - :math:`\rho_e`
     - ``edge_density``
     - Edge density
     - kg/m³
   * - :math:`\mu_e`
     - ``edge_dynamic_viscosity``
     - Edge dynamic viscosity
     - Pa s
   * - :math:`Re_x`
     - ``reynolds_number``
     - Reynolds number based on :math:`x`
     - dimensionless
   * - :math:`q_e`
     - computed internally
     - Edge dynamic pressure
     - Pa

Thickness correlations
^^^^^^^^^^^^^^^^^^^^^^

The laminar Blasius relations are

.. math::

   \frac{\delta_{99}}{x}=\frac{5}{\sqrt{Re_x}},
   \qquad
   \frac{\delta^*}{x}=\frac{1.7208}{\sqrt{Re_x}},
   \qquad
   \frac{\theta_m}{x}=\frac{0.664}{\sqrt{Re_x}}.

For a turbulent boundary layer with ``NONE``, the one-fifth-power thickness and
one-seventh-power-profile estimates are

.. math::

   \delta_{99}=0.37xRe_x^{-1/5},
   \qquad
   \delta^*=\frac{\delta_{99}}{8},
   \qquad
   \theta_m=\frac{7}{72}\delta_{99}.

.. list-table:: Thickness symbols
   :header-rows: 1
   :widths: 17 30 37 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`\delta_{99}`
     - ``boundary_layer_thickness``
     - Conventional 99-percent thickness
     - m
   * - :math:`\delta^*`
     - ``displacement_thickness``
     - Displacement thickness
     - m
   * - :math:`\theta_m`
     - ``momentum_thickness``
     - Momentum thickness
     - m

.. _flat-plate-skin-friction-coefficients:

Skin-friction coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^

The local coefficient relates wall shear stress to edge dynamic pressure:

.. math::

   C_{f,x}=\frac{\tau_w(x)}{q_e}.

The average coefficient from the leading edge to :math:`x` is

.. math::

   \bar C_f(x)
   =\frac{1}{x}\int_0^xC_{f,\xi}\,\mathrm{d}\xi
   =\frac{D'(x)}{q_ex},
   \qquad
   D'(x)=\frac{1}{2}\rho_eU_e^2x\bar C_f.

:math:`D'` is one-sided accumulated drag per unit width. These definitions
assume constant :math:`q_e` along the plate. They are flat-plate
*skin-friction coefficients*, not Darcy or Fanning pipe-friction factors.

.. list-table:: Skin-friction outputs
   :header-rows: 1
   :widths: 16 32 36 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`C_{f,x}`
     - ``local_skin_friction_coefficient``
     - Local skin-friction coefficient
     - dimensionless
   * - :math:`\bar C_f`
     - ``average_skin_friction_coefficient``
     - Leading-edge-averaged coefficient
     - dimensionless
   * - :math:`\tau_w`
     - ``wall_shear_stress``
     - Local wall shear stress
     - Pa
   * - :math:`D'`
     - ``drag_per_unit_width``
     - Accumulated one-sided drag per unit width
     - N/m

For ``BoundaryLayerRegime.LAMINAR``, the Blasius solution gives

.. math::

   C_{f,x}=\frac{0.664}{\sqrt{Re_x}},
   \qquad
   \bar C_f=\frac{1.328}{\sqrt{Re_x}}.

For ``TurbulentCorrelation.POWER_LAW``,

.. math::

   C_{f,x}=0.0592Re_x^{-1/5},
   \qquad
   \bar C_f=0.074Re_x^{-1/5}.

The default ``TurbulentCorrelation.SCHLICHTING`` model uses

.. math::

   \bar C_f=\frac{0.455}{(\log_{10}Re_x)^{2.58}}.

For constant edge conditions, differentiating accumulated drag gives the
local relation

.. math::

   C_{f,x}=\frac{\mathrm{d}(x\bar C_f)}{\mathrm{d}x}
   =\bar C_f\left(1-\frac{2.58}{\ln Re_x}\right).

Turbulent use outside :math:`5\times10^5\le Re_x\le10^9` emits
:class:`~aerophysics.exceptions.ApplicabilityWarning`.

Specified transition
^^^^^^^^^^^^^^^^^^^^

For ``BoundaryLayerRegime.TRANSITIONAL``, local quantities switch from the
laminar to selected turbulent correlation at the caller-supplied
:math:`Re_\mathrm{tr}`. Average friction preserves the laminar drag accumulated
before transition:

.. math::

   \bar C_{f,\mathrm{mixed}}(Re_x)
   =\bar C_{f,t}(Re_x)
   +\frac{Re_\mathrm{tr}}{Re_x}
    \left[
      \bar C_{f,l}(Re_\mathrm{tr})
      -\bar C_{f,t}(Re_\mathrm{tr})
    \right].

This is a prescribed transition location, not a transition-prediction model.

.. _flat-plate-compressibility-corrections:

Compressibility corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compressibility is opt-in. When fast-moving gas is slowed by friction near
an insulated wall, the wall heats to the **recovery temperature**
:math:`T_r`. It is a little below the total temperature, because heat
conduction carries part of the energy away. The recovery factor :math:`r`
expresses this:

.. math::

   T_r=T_e\left[1+r\frac{\gamma-1}{2}M^2\right],
   \qquad
   r_l=\sqrt{Pr},
   \qquad
   r_t=Pr^{1/3}.

If ``wall_temperature`` is omitted, the adiabatic value :math:`T_w=T_r` is
used.

Eckert reference-temperature method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Eckert's idea is to keep the incompressible formulas but evaluate density and
viscosity at a single **reference temperature** :math:`T^*` that represents
conditions inside the boundary layer, somewhere between the edge and the wall.
For ``CompressibilityCorrection.ECKERT``,

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
(:ref:`White, 2006 <ref-white-2006>`). For laminar flow this is the familiar
Chapman--Rubesin form :math:`C_{f,x}=0.664\sqrt{C^*/Re_x}` with
:math:`C^*=\rho^*\mu^*/(\rho_e\mu_e)`. At :math:`M\to0` with
:math:`T_w=T_e` and a supplied edge viscosity :math:`\mu_e=\mu(T_e)`,
:math:`T^*=T_e` and the friction and laminar thicknesses reduce to ``NONE``. Turbulent
thicknesses remain unavailable. If the supplied edge
viscosity differs from the selected model, the effective Reynolds number need
not equal :math:`Re_x`.

Laminar thickness correlations are evaluated at :math:`Re_x^*` without a
density conversion. For laminar Blasius flow,
:math:`2\,\mathrm{d}\theta_m/\mathrm{d}x=C_f^*` rather than :math:`C_{f,x}`.
For turbulent portions, all three thicknesses are ``NaN``: the former
reference-property estimates were inaccurate in compressible flow. This
also applies after the prescribed transition. Friction and drag are unchanged.

Van Driest II method
~~~~~~~~~~~~~~~~~~~~

Van Driest II maps the compressible turbulent boundary layer onto an
equivalent incompressible one. It is derived from a mixing-length model with
the density variation across the layer included, and it is generally more
accurate than a reference-temperature method for turbulent flow at high Mach
number. For turbulent ``CompressibilityCorrection.VAN_DRIEST_II``, the
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

Laminar portions of a VD2 calculation still use Eckert, including the
:math:`T_e/T^*` conversion to the edge basis and the laminar thicknesses.
Turbulent :math:`\delta_{99}`, :math:`\delta^*`, and :math:`\theta_m` are
``NaN``, including post-transition points. The former thickness estimates
at :math:`F_\theta Re_x` were inaccurate in compressible flow.
``turbulent_correlation`` is ignored for VD2 turbulent portions and remains
active for ``NONE`` and ``ECKERT``. The Willems friction equations are unchanged.

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
     - Equivalent incompressible VD2 friction Reynolds number
     - dimensionless
   * - :math:`F_C,F_\theta,F_x`
     - internal
     - Van Driest II transformation factors
     - dimensionless

Van Driest II with wall-temperature correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select ``CompressibilityCorrection.VAN_DRIEST_II_WALL_TEMPERATURE`` to obtain
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

Thus :math:`w=1` for an adiabatic wall. This mode ignores
``turbulent_correlation`` and uses the Schlichting pair, rather than the
Willems pair used by the original VD2 mode. At constant edge and wall
conditions, :math:`c_f=\mathrm{d}(x\bar C_F)/\mathrm{d}x`.
For ``TRANSITIONAL``, the laminar part uses Eckert and the average uses the
specified-transition mixing equation above with this new turbulent average.

The zero-pressure-gradient momentum integral gives the full-edge momentum
thickness and wall shear stress:

.. math::

   \theta=\frac{x\bar C_F}{2}, \qquad
   \tau_w=\frac{1}{2}\rho_e U_e^2 c_f.

Here :math:`\bar C_F` is the selected turbulent or mixed average, on the
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

Applicability and external comparisons
++++++++++++++++++++++++++++++++++++++

The tested envelope is :math:`M_e=0.3`--13.3, :math:`T_w/T_{aw}=0.07`--1,
:math:`T_\infty=47`--271 K, and :math:`Re_{x,e}\approx10^6`--:math:`10^9`
(:math:`M_e<2` only up to :math:`1.6\times10^8`), for a smooth
zero-pressure-gradient plate and a
calorically perfect gas. Flow is turbulent from the leading edge, or uses
the prescribed transition model. This envelope is evidence of comparison,
not a guarantee for every combination of these parameters; the wake closure
can fail, including near a prescribed transition.

Supply **local edge conditions** rather than freestream conditions. Viscous
interaction at a sharp leading edge can raise the edge pressure 10--27%
above freestream pressure at Mach 10. Evaluating a correlation with local
states does not turn it into a streamwise variable-edge-state solver: the
momentum/drag formulas retain their constant-state, zero-pressure-gradient
assumption.

The following summaries are externally supplied validation results dated
2026-11-01, not simulations rerun by this package's test suite. Differences
are model minus reference, divided by reference:

* Against SU2 8.5.0 standard SA, :math:`\delta_{99}` RMS differences are
  **2.5%** for the eight fitted conditions (343 sections), **2.3%** for six
  held-out conditions (264 sections), and **2.1%** for long plates through
  :math:`Re_{x,e}=1.7\times10^8`. Above :math:`10^8`, the RMS difference is
  **5.2%**, spanning **−7.6 to +0.6%**.
* Displacement thickness integrated through :math:`\delta_{99}` differs by
  **2.3--2.7% RMS**.
* Four SU2 standard-SA cases at DNS-matched conditions vary freestream
  temperature and unit Reynolds number: :math:`T_\infty=47`--271 K,
  :math:`Re'=4.8`--:math:`10.6\times10^6\ \mathrm{m}^{-1}`,
  :math:`M_e=2.5`--13.3, :math:`T_w/T_{aw}=0.19`--0.99, and
  :math:`Re_{x,e}\le9\times10^7`. Casewise :math:`\delta_{99}` RMS differences
  are **1.8--4.0%** (range **−8.9 to +7.5%**); displacement-thickness RMS
  differences are **1.2--4.0%**.
* The :math:`\delta_{99}` difference falls almost linearly with
  :math:`\log_{10}Re_{x,e}`, by **5--6 percentage points per decade**, with no
  levelling off up to :math:`10^9`. Seven further SU2 standard-SA cases
  (M2, M6 and M10 at :math:`Re'=1.2\times10^8\ \mathrm{m}^{-1}` and
  :math:`T_w/T_{aw}=0.07`--0.5; M1.5, M0.8 and M0.3 at
  :math:`Re'=2\times10^7\ \mathrm{m}^{-1}` with the wall at recovery
  temperature) give :math:`\delta_{99}` differences of **+0.2 to +1.1%** at
  :math:`Re_{x,e}=10^7`, **−6.3 to −3.5%** at :math:`10^8`, and **−12.1 to
  −8.6%** at :math:`10^9` for the supersonic cold walls. For the
  recovery-temperature walls at :math:`M_e\le1.5`, the differences are
  **+9 to +11%** at :math:`10^6`, **+4 to +6%** at :math:`10^7`, and **−0.5 to
  +1.1%** at :math:`10^8`. The slope is the same at M0.3, so it comes from the
  profile shape, not from the compressibility transformation. Whether the
  correlation or the SA model is closer to measured :math:`\delta_{99}` at
  high Reynolds number has not been established.
* Over the same cases, :math:`C_f` differs by **−1.3 to +2.0%** mean and
  :math:`\theta` by **−1.3 to +3.5%** mean (RMS at most **2.2%** and
  **3.8%**), up to :math:`Re_{x,e}=10^9` and down to M0.3.
* Against incompressible zero-pressure-gradient measurements and DNS
  (:ref:`Österlund et al. (2000) <ref-osterlund-2000>`,
  :ref:`Schlatter and Örlü (2010) <ref-schlatter-orlu-2010>`,
  :ref:`Vallikivi, Hultmark and Smits (2015) <ref-vallikivi-hultmark-smits-2015>`;
  :math:`Re_\theta=2\times10^3`--:math:`2.3\times10^5`, about
  :math:`Re_x=10^6`--:math:`2.7\times10^8`), :math:`\delta_{99}` is too thick
  by about **+19%** at :math:`Re_x=10^6`, **+15%** at :math:`10^7`, **+10%**
  at :math:`10^8`, and **+9%** at :math:`1.5\times10^8`. SU2 standard SA is
  **+7 to +9%** thick over the same range, so most of the SA-relative drift
  above comes from the correlation approaching the measurements as
  :math:`Re_x` rises. The difference comes mainly from the incompressible
  Schlichting friction pair: at a given :math:`Re_\theta` its :math:`C_f` is
  **6--6.5%** above the measurements, which makes :math:`\theta` about
  **5.5%** thick and lowers the solved wake parameter (0.33--0.45 against
  about 0.6 measured). Given the measured :math:`\theta` and :math:`C_f`, the
  profile closure reproduces :math:`\delta_{99}` within **−6 to +2%** on
  average. The widely used friction pair is kept unchanged.
* The solved wake parameter rises with Reynolds number on cold hypersonic
  walls: at M10 with :math:`T_w/T_{aw}=0.07` it reaches :math:`\Pi=0.976` at
  :math:`Re_{x,e}=10^9`. Extrapolating that trend, the wake bound
  :math:`\Pi\le1` would be exceeded, and
  :class:`~aerophysics.exceptions.ModelRangeError` raised, near
  :math:`Re_{x,e}\approx1.5\times10^9` for that condition.
* At matched :math:`Re_\theta`, comparison with the DNS of
  :ref:`Zhang, Duan and Choudhari (2018) <ref-zhang-duan-choudhari-2018>`
  gives :math:`C_f` differences **−0.3, −0.6, −3.8, −2.2%**, and
  :math:`\delta_{99}` differences **−4.5, −0.5, +0.6, +15.3%**, respectively
  for M2.5, M6 with :math:`T_w/T_r=0.76`, M6 with :math:`T_w/T_r=0.25`, and
  M14 with :math:`T_w/T_r=0.18`. The M14 result shows that the good M6
  agreement should not be extrapolated to higher Mach numbers.

The exponent **0.064** is fitted to SU2 solutions; DNS also exhibits lower
cold-wall friction than VD2. The factor **0.978** is the SU2 mean ratio of
momentum thickness through :math:`\delta_{99}` to full-edge momentum thickness
(the truncation difference varies from **−5 to −1%**). Neither constant is
universal. The nine-row independent-implementation regression fixture uses
exponent 0.06391; rounding to 0.064 changes the factor by less than 0.03%.
Its ``su2_delta99`` column is comparison data, not a unit-test acceptance target.

Why the original compressible thicknesses are unavailable
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Across the supplied 1059-section SU2 comparison (standard SA, M2--M10,
:math:`T_w/T_{aw}=0.07`--1), Eckert's former :math:`\delta_{99}` estimate had
mean error **+58%** (range **+5 to +92%**); VD2 had mean error **+17%**
(**18% RMS**, range **−1 to +37%**). Fixed ratios :math:`\delta^*=\delta/8`
and :math:`\theta=7\delta/72` missed displacement thickness by factors of
2--5 (underprediction) and momentum thickness by factors of 2--10
(overprediction). Consequently ``ECKERT`` and ``VAN_DRIEST_II`` now return
``NaN`` for all three turbulent thicknesses, including post-transition
points. Friction, drag, laminar portions, and all ``NONE`` outputs retain
their existing definitions.

Worked friction example
^^^^^^^^^^^^^^^^^^^^^^^

For :math:`x=1\ \mathrm{m}`, :math:`U_e=100\ \mathrm{m/s}`,
:math:`\rho_e=1.225\ \mathrm{kg/m^3}`, and
:math:`\mu_e=1.7894\times10^{-5}\ \mathrm{Pa\,s}`, the Schlichting model
gives :math:`Re_x=6.84587\times10^6`, :math:`C_{f,x}=0.002670314`,
:math:`\bar C_f=0.003193859`, :math:`\tau_w=16.3557\ \mathrm{Pa}`, and
:math:`D'=19.5624\ \mathrm{N/m}`.

>>> from aerophysics import BoundaryLayerRegime, flat_plate_boundary_layer
>>> layer = flat_plate_boundary_layer(
...     1.0,
...     edge_velocity=100.0,
...     edge_density=1.225,
...     edge_dynamic_viscosity=1.7894e-5,
...     regime=BoundaryLayerRegime.TURBULENT,
... )
>>> f"{layer.reynolds_number:.5e}"
'6.84587e+06'
>>> round(layer.local_skin_friction_coefficient, 9)
0.002670314
>>> round(layer.average_skin_friction_coefficient, 9)
0.003193859
>>> round(layer.wall_shear_stress, 4)
16.3557
>>> round(layer.drag_per_unit_width, 4)
19.5624

The result object is
:class:`aerophysics.boundary_layer.FlatPlateBoundaryLayerResult`.

Applicability
-------------

The flat-plate model excludes surface roughness, streamwise pressure
gradients, separation, suction, blowing, and natural-transition prediction.
The gas remains calorically perfect and compressible transport properties use
the supplied dynamic-viscosity model. The thickness models remain engineering estimates, including the coupled
mean-profile closure; comparison with SA is not independent physical validation.
