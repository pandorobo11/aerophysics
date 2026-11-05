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

Set ``compressible`` independently of the turbulent model. It accepts Python
``bool`` and NumPy ``bool_`` values. Thermal inputs (``mach``,
``edge_temperature``, and ``wall_temperature``) require ``compressible=True``;
Mach and edge temperature are mandatory in that case. Temperatures use K.

.. list-table:: Valid model selections
   :header-rows: 1
   :widths: 22 16 62

   * - Regime
     - ``compressible``
     - ``turbulent_model``
   * - ``LAMINAR``
     - False or True
     - Must be None: Blasius for False, compressible similarity for True
   * - ``TURBULENT`` / ``TRANSITIONAL``
     - False
     - None defaults to SCHLICHTING; SCHLICHTING or SCHLICHTING_LOG_WAKE
   * - ``TURBULENT`` / ``TRANSITIONAL``
     - True
     - Required: ECKERT, VAN_DRIEST_II,
       or VAN_DRIEST_II_WALL_TEMPERATURE

Other combinations raise ``ValueError``; an irrelevant selection is never
silently ignored. Each member of
:class:`~aerophysics.boundary_layer.TurbulentModel` has a ``compressible``
property indicating which setting it requires.

.. list-table:: TurbulentModel friction and thickness treatment
   :header-rows: 1
   :widths: 28 38 34

   * - Member
     - Skin friction (equations below)
     - Thicknesses
   * - ``SCHLICHTING``
     - Logarithmic mean friction and accumulated-drag derivative
     - One-fifth-power delta99, one-seventh-power profile integrals
   * - ``SCHLICHTING_LOG_WAKE``
     - Same as ``SCHLICHTING``
     - Momentum-integral theta, Spalding--Coles delta99 and displacement
   * - ``ECKERT``
     - Schlichting at Eckert reference state, converted to edge basis
     - All three NaN at turbulent points
   * - ``VAN_DRIEST_II``
     - Transformed implicit Willems equations
     - All three NaN at turbulent points
   * - ``VAN_DRIEST_II_WALL_TEMPERATURE``
     - Transformed Schlichting pair multiplied by (Tw/Tr)**0.064
     - Momentum-integral theta, Van Driest/Walz delta99 and displacement

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

For ``compressible=False``, the laminar Blasius relations are

.. math::

   \frac{\delta_{99}}{x}=\frac{4.9100}{\sqrt{Re_x}},
   \qquad
   \frac{\delta^*}{x}=\frac{1.7208}{\sqrt{Re_x}},
   \qquad
   \frac{\theta_m}{x}=\frac{0.664}{\sqrt{Re_x}}.

For a turbulent boundary layer with ``compressible=False``, the two
incompressible models share the same friction and drag and differ only in
their thicknesses. Both are for smooth, zero-pressure-gradient flow.

.. _flat-plate-power-law-thickness:

One-fifth-power thickness (``SCHLICHTING``, default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \delta_{99}=0.37xRe_x^{-1/5},
   \qquad
   \delta^*=\frac{\delta_{99}}{8},
   \qquad
   \theta_m=\frac{7}{72}\delta_{99}.

This is the classical estimate of :ref:`von Kármán (1921)
<ref-von-karman-1921>`, as given by :ref:`Schlichting and Gersten (2000)
<ref-schlichting-gersten-2000>`. It assumes the one-seventh-power profile
:math:`u/U_e=(y/\delta)^{1/7}`, which gives the two ratios above, and Blasius'
pipe-friction law applied to the plate,
:math:`\tau_w/(\rho_eU_e^2)=0.0225(\nu/U_e\delta)^{1/4}`. Integrating the
momentum relation :math:`\mathrm{d}\theta_m/\mathrm{d}x=\tau_w/(\rho_eU_e^2)`
from the leading edge gives the coefficient
:math:`(0.0225\cdot\tfrac{72}{7}\cdot\tfrac{5}{4})^{4/5}=0.3707`.
The pipe law holds only at low Reynolds number, and the fixed profile gives a
constant shape factor :math:`H=9/7` and :math:`\delta_{99}/\theta_m=72/7`.
Measured layers have :math:`H\approx1.4` and :math:`\delta_{99}/\theta_m\approx8`
near :math:`Re_x=10^6`; equilibrium log-wake profiles give :math:`H` falling to
about 1.2 and :math:`\delta_{99}/\theta_m` rising to about 11 by :math:`10^9`.
The correlation is therefore nominally limited to
:math:`5\times10^5\le Re_x\le10^7`. Turbulent points above
:math:`Re_x=10^7` (for ``TRANSITIONAL``, the equivalent :math:`Re_\mathrm{eq}`
described under :ref:`flat-plate-specified-transition`) emit
:class:`~aerophysics.exceptions.ApplicabilityWarning`.

The Schlichting friction is not derived from these thicknesses:
:math:`\theta_m` differs from :math:`x\bar C_f/2`, the value implied by the
returned drag, by +2% at :math:`Re_x=10^6`, −15% at :math:`10^8`, and −27% at
:math:`10^9`.

For ``TRANSITIONAL``, :math:`x` and :math:`Re_x` in these equations are
:math:`x_\mathrm{eq}` and :math:`Re_\mathrm{eq}` of the equivalent turbulent
plate described under :ref:`flat-plate-specified-transition`, the same origin
used for friction. The power-law :math:`\theta_m` is then 1.3--2.2% above the
Blasius value at transition for :math:`Re_\mathrm{tr}=5\times10^5`--
:math:`5\times10^6`, the power law's own mismatch with Schlichting friction
at :math:`Re_\mathrm{eq}\approx1.4\times10^5`--:math:`6\times10^5`.
:math:`\delta_{99}` jumps there, because the assumed profile changes from
Blasius to the one-seventh-power shape.

.. _flat-plate-log-wake-thickness:

Log-wake profile thickness (``SCHLICHTING_LOG_WAKE``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SCHLICHTING_LOG_WAKE`` uses the momentum integral and an equilibrium
mean-velocity profile instead of fixed ratios:

.. math::

   \theta_m=\frac{x\bar C_f}{2},
   \qquad
   \tau_w=\tfrac{1}{2}\rho_eU_e^2C_{f,x}.

:math:`\bar C_f` is the Schlichting average, or the specified-transition
average, so :math:`\theta_m` always agrees with ``drag_per_unit_width`` and is
continuous at transition. :math:`\delta_{99}` is the value :math:`d` for which
the Spalding inner law (:ref:`Spalding, 1961 <ref-spalding-1961>`) plus the
Coles wake (:ref:`Coles, 1956 <ref-coles-1956>`), with :math:`\kappa=0.41`,
:math:`B=5.2`, :math:`0\le\Pi\le1`, and :math:`u(d)=0.99U_e`, satisfies

.. math::

   \int_0^d\frac{u}{U_e}\left(1-\frac{u}{U_e}\right)\mathrm{d}y=0.978\,\theta_m.

``displacement_thickness`` is :math:`\int_0^d(1-u/U_e)\,\mathrm{d}y`. The wall
unit is :math:`\mu_e/(\rho_eu_\tau)` with
:math:`u_\tau=\sqrt{\tau_w/\rho_e}`, using the supplied
``edge_dynamic_viscosity``. This is the constant-property limit of the
``VAN_DRIEST_II_WALL_TEMPERATURE`` closure described below, with the same
factor 0.978. The two models agree within 0.1% at Mach 0.03.

Where no wake parameter :math:`0\le\Pi\le1` matches, :math:`\delta_{99}` and
:math:`\delta^*` are ``NaN`` and one
:class:`~aerophysics.exceptions.ApplicabilityWarning` reports the number of
such points; :math:`\theta_m`, friction, and drag are still returned. Fully
turbulent plates are admissible for :math:`10^4\le Re_x\le10^{10}`. After a
prescribed transition, friction and :math:`\theta_m` are those of the fully
turbulent plate at :math:`Re_\mathrm{eq}`, so the closure is admissible
whenever that plate is; :math:`Re_\mathrm{eq}\approx1.4\times10^5` just after
:math:`Re_\mathrm{tr}=5\times10^5`.

Thickness accuracy
~~~~~~~~~~~~~~~~~~

The log-wake closure with Schlichting friction was compared with
incompressible zero-pressure-gradient measurements and DNS (see the
``Van Driest II with wall-temperature correction`` section). Its
:math:`\delta_{99}` is about 19% thick at :math:`Re_x=10^6`, 15% at
:math:`10^7`, 10% at :math:`10^8`, and 9% at :math:`1.5\times10^8`. Most of
this offset comes from the Schlichting friction, which makes :math:`\theta_m`
about 5.5% thick. Multiplying by the ratio of the two models gives the
estimated one-fifth-power error:

.. list-table:: Fully turbulent thickness ratios, ``SCHLICHTING`` /
   ``SCHLICHTING_LOG_WAKE``
   :header-rows: 1
   :widths: 14 16 16 16 38

   * - :math:`Re_x`
     - :math:`\delta_{99}`
     - :math:`\theta_m`
     - :math:`\delta^*`
     - Estimated ``SCHLICHTING`` :math:`\delta_{99}` error against measurements
   * - :math:`10^6`
     - 1.18
     - 1.02
     - 0.96
     - +40%
   * - :math:`10^7`
     - 1.01
     - 0.95
     - 0.98
     - +16%
   * - :math:`10^8`
     - 0.83
     - 0.85
     - 0.90
     - −8%
   * - :math:`1.5\times10^8`
     - 0.80
     - 0.83
     - 0.88
     - −12%
   * - :math:`10^9`
     - 0.67
     - 0.73
     - 0.79
     - not measured

The one-fifth-power error changes sign near :math:`Re_x\approx5\times10^7`
and grows quickly beyond :math:`10^8`, while the log-wake error decreases
steadily with Reynolds number. ``SCHLICHTING_LOG_WAKE`` is recommended when
:math:`\delta_{99}`, :math:`\delta^*`, or a :math:`\theta_m` consistent with drag
matters. ``SCHLICHTING`` remains the default for backward compatibility.

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

For ``BoundaryLayerRegime.LAMINAR`` with ``compressible=False``, the Blasius solution gives

.. math::

   C_{f,x}=\frac{0.664}{\sqrt{Re_x}},
   \qquad
   \bar C_f=\frac{1.328}{\sqrt{Re_x}}.

The default ``TurbulentModel.SCHLICHTING`` model uses

.. math::

   \bar C_f=\frac{0.455}{(\log_{10}Re_x)^{2.58}}.

For constant edge conditions, differentiating accumulated drag gives the
local relation

.. math::

   C_{f,x}=\frac{\mathrm{d}(x\bar C_f)}{\mathrm{d}x}
   =\bar C_f\left(1-\frac{2.58}{\ln Re_x}\right).

Turbulent use outside :math:`5\times10^5\le Re_x\le10^9` emits
:class:`~aerophysics.exceptions.ApplicabilityWarning`.

.. _flat-plate-specified-transition:

Specified transition
^^^^^^^^^^^^^^^^^^^^

For ``BoundaryLayerRegime.TRANSITIONAL``, local quantities switch from the
laminar to the selected turbulent model at the caller-supplied
:math:`Re_\mathrm{tr}`. Downstream, the layer continues as the fully turbulent
plate of that model whose momentum thickness equals the laminar value at
transition. For constant edge and wall conditions and zero pressure gradient,
the momentum integral gives :math:`\theta_m=x\bar C_f/2` for both the laminar
and the fully turbulent plate, so the turbulent start :math:`Re_0` solves

.. math::

   Re_0\,\bar C_{f,t}(Re_0)=Re_\mathrm{tr}\,\bar C_{f,l}(Re_\mathrm{tr})
   =2Re_{\theta,\mathrm{tr}}.

The equivalent Reynolds number and length then give

.. math::

   Re_\mathrm{eq}=Re_0+Re_x-Re_\mathrm{tr},
   \qquad
   x_\mathrm{eq}=x\frac{Re_\mathrm{eq}}{Re_x},

   C_{f,x}=C_{f,t}(Re_\mathrm{eq}),
   \qquad
   \bar C_f=\frac{Re_\mathrm{eq}}{Re_x}\bar C_{f,t}(Re_\mathrm{eq}).

Here :math:`C_{f,t}` and :math:`\bar C_{f,t}` are the fully turbulent
edge-basis coefficients of the selected model, and :math:`\bar C_{f,l}` is
Blasius, or the compressible similarity solution with ``compressible=True``.
The compressible models write these coefficients as
:math:`s\,C_{f,i}(kRe_x)`, with :math:`(s,k)=(\rho^*/\rho_e,Re^*/Re_x)` for
``ECKERT`` and :math:`(w/F_C,F_x)` for the Van Driest II models (:math:`w=1`
for ``VAN_DRIEST_II``), so the root is found in incompressible variables:

.. math::

   kRe_0\,\bar C_{f,i}(kRe_0)=2Re_{\theta,\mathrm{tr}}\frac{k}{s}.

For Van Driest II with :math:`w=1`, :math:`k/s=F_\theta`, the momentum-thickness
Reynolds number transformation. The Schlichting root is found by Brent's
method in :math:`\ln R`; the Willems average equation is inverted explicitly as
:math:`R=S(\log_{10}S/0.242)^2` with :math:`S=R\bar C_{f,i}`. Thickness
correlations, including the ``SCHLICHTING`` power law, use
:math:`x_\mathrm{eq}` and :math:`Re_\mathrm{eq}`, and
``effective_reynolds_number`` returns :math:`kRe_\mathrm{eq}` at turbulent
points (:math:`k=1` for ``compressible=False``). :math:`Re_0` for
:math:`Re_\mathrm{tr}=5\times10^5` is :math:`1.417\times10^5` with
Schlichting friction.

Accumulated drag, and with it the drag-based momentum thickness
:math:`x\bar C_f/2`, is continuous at transition for every model. Whether the
returned quantities satisfy stronger relations depends on the model:

.. list-table:: Exact and approximate relations after a specified transition
   :header-rows: 1

   * - Model
     - Returned :math:`\theta_m` at transition
     - :math:`\mathrm{d}(x\bar C_f)/\mathrm{d}x=C_{f,x}`
   * - ``SCHLICHTING``
     - Power law, 1.3--2.2% above the laminar value for
       :math:`Re_\mathrm{tr}=5\times10^5`--:math:`5\times10^6`
     - Exact
   * - ``SCHLICHTING_LOG_WAKE``
     - :math:`x\bar C_f/2`, continuous
     - Exact
   * - ``ECKERT``
     - ``NaN``
     - Exact
   * - ``VAN_DRIEST_II``
     - ``NaN``
     - Approximate: the Willems local and average equations are separate
       fits, 3.9% apart at :math:`R=10^6`
   * - ``VAN_DRIEST_II_WALL_TEMPERATURE``
     - :math:`x\bar C_f/2`, continuous
     - Exact

The ``SCHLICHTING`` and ``VAN_DRIEST_II`` departures belong to their fully
turbulent correlations and are unchanged by the equivalent origin. Local
friction, :math:`\delta_{99}`, and :math:`\delta^*` jump, because transition is
instantaneous here. Turbulent friction just after transition is consistent
with the drag-based momentum thickness, and the ``SCHLICHTING_LOG_WAKE`` and
``VAN_DRIEST_II_WALL_TEMPERATURE`` closures see a fully turbulent equilibrium
state.

Earlier versions kept the leading-edge origin for turbulent friction and
used the Prandtl-type mixed average
:math:`\bar C_{f,t}(Re_x)+(Re_\mathrm{tr}/Re_x)[\bar C_{f,l}-\bar C_{f,t}]_{Re_\mathrm{tr}}`.
That average conserves drag, but its local friction belongs to a layer whose
momentum thickness is 2.8 times the laminar value at
:math:`Re_\mathrm{tr}=5\times10^5`. With Schlichting friction, the local
coefficient just after transition is now 26% higher for
:math:`Re_\mathrm{tr}=5\times10^5` and 42% higher for
:math:`5\times10^6`. Average friction is 6--9% higher at
:math:`Re_x/Re_\mathrm{tr}=1.2`--2 for :math:`Re_\mathrm{tr}=5\times10^5`, up
to 15% for :math:`5\times10^6`, and less than 1% higher beyond
:math:`100Re_\mathrm{tr}`. Fully turbulent and laminar results are unchanged.

Drag continuity, the exact drag derivatives, the ``SCHLICHTING`` theta offset,
and the fully turbulent limits are verified analytically in the test suite. The equivalent plate has not been compared
with transitional measurements; a real transition region has finite length and
non-equilibrium turbulence. If :math:`2Re_{\theta,\mathrm{tr}}` is below the
smallest value the turbulent correlation can produce (about 4.5 for
Schlichting, :math:`Re_\mathrm{tr}\lesssim11` with Blasius), the call raises
:class:`~aerophysics.exceptions.ModelRangeError`.

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
   r_t=Pr^{1/3}.

This recovery-factor expression applies to turbulent points. For laminar points,
the adiabatic wall temperature is solved directly with zero wall heat flux;
:math:`r_l=\sqrt{Pr}` is only an initial estimate. If ``wall_temperature`` is
omitted, the corresponding adiabatic wall is used. With ``compressible=True``,
laminar points use the compressible similarity solution.

Eckert reference-temperature method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

For turbulent portions, all three thicknesses are ``NaN``: the former
reference-property estimates were inaccurate in compressible flow. This
also applies after the prescribed transition. Friction and drag are unchanged.

Compressible laminar similarity solution
++++++++++++++++++++++++++++++++++++++++

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

The calculation groups exactly equal states (Mach, edge temperature, wall
condition, Pr, gamma) within each call, without rounding or a persistent
state cache. Conservative second-order finite differences and damped Newton
iteration operate in batches of at most 256 states. A stretched, fixed
161/321-node pair and Richardson extrapolation are followed by 321/641
nodes only for uncertain states. The solver checks residuals, boundaries,
positive temperatures and viscosities, grid differences, far-field decay,
and the independent momentum integral. Remaining states use
``scipy.integrate.solve_bvp(tol=1e-8)`` with an enlarged-domain check.
Unresolved states raise ``ModelRangeError`` with their original input indices
and state. Blasius and Walz profiles supply initial guesses only; they are
not the returned physical approximation. There are no interpolation tables.

Distance :math:`x` starts at the leading edge and :math:`y` points outward
from the wall. Thicknesses are in m and temperatures in K. Assumptions are a
perfect gas, constant cp and Pr, zero pressure gradient, and no induced
pressure. For :math:`T_{aw}\gtrsim2000\ \mathrm{K}`, changing specific heats
and dissociation make results indicative only; this is not a universal
chemical threshold or a real-gas validity claim.

Independent similarity comparison
+++++++++++++++++++++++++++++++++

The unchanged independent ``solve_bvp`` reference generator
``docs/scripts/build_laminar_similarity_reference.py`` produces the CSV and
provenance JSON under
``tests/reference_data/viscous_flow/laminar_compressible_similarity``.
It uses Sutherland air, :math:`T_e=240\ \mathrm{K}`, :math:`Pr=0.72`,
:math:`\gamma=1.4`, :math:`M_e=0.3,2,4,6,10,15`, and
:math:`T_w/T_{aw}=0.1,0.3,0.5,1` (24 points). Isothermal API inputs use the
reference's solved :math:`T_{aw}` times the prescribed ratio; the last wall
condition is truly adiabatic. Errors below compare the library solution with
this independent reference, with equal weight for each point.

.. include:: ../_generated/laminar_similarity_accuracy.rst

The acceptance bound is 0.1% relative for friction, momentum and delta99.
For displacement coefficients of magnitude at least 0.05 it is also 0.1%;
closer to zero the absolute coefficient error must be at most 5e-5.
The table separates these error measures and suppresses absolute errors below
5e-10 to avoid reporting insignificant CPU differences. These are numerical
checks of constant-cp perfect-gas equations, not physical validation.
The generated-assets suite also checks 200 seeded random states over
Mach 0--20, Te 30--2000 K, and Tw/Taw 0.05--1.5, with selected independent
collocation comparisons; timing is not a CI requirement.

Independent adiabatic examples are recovered within 0.5%: at M2,
:math:`c_f\sqrt{Re_x}=0.6356105`, :math:`\delta_{99}\sqrt{Re_x}/x=6.414710`,
:math:`\delta^*\sqrt{Re_x}/x=3.303259`, :math:`T_{aw}/T_e=1.676431`;
at M10 the corresponding values are 0.4258457, 29.674188, 27.668758, and
17.628489. Tests also check the Blasius limit and reference convergence
with domain enlargement (14 to 20) and tolerance tightening (1e-8 to 1e-10).

Viscous interaction warning
+++++++++++++++++++++++++++

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

Van Driest II method
~~~~~~~~~~~~~~~~~~~~

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

With ``compressible=True``, laminar portions use the similarity solution,
with friction on the edge basis.
Turbulent :math:`\delta_{99}`, :math:`\delta^*`, and :math:`\theta_m` are
``NaN``, including post-transition points. The former thickness estimates
at :math:`F_\theta Re_x` were inaccurate in compressible flow.
``TurbulentModel.VAN_DRIEST_II`` selects the Willems friction equations.

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

Van Driest II with wall-temperature correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Applicability and external comparisons
++++++++++++++++++++++++++++++++++++++

The tested envelope is :math:`M_e=0.3`--13.3, :math:`T_w/T_{aw}=0.07`--1,
:math:`T_\infty=47`--271 K, and :math:`Re_{x,e}\approx10^6`--:math:`10^9`
(:math:`M_e<2` only up to :math:`1.6\times10^8`), for a smooth
zero-pressure-gradient plate and a
calorically perfect gas. Flow is turbulent from the leading edge, or uses
the prescribed transition model. This envelope is evidence of comparison,
not a guarantee for every combination of these parameters; the wake closure
can fail.

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

In-repository SA boundary-layer reference
+++++++++++++++++++++++++++++++++++++++++

``docs/scripts/build_turbulent_sa_reference.py`` is an independent
verification solver, not part of the package. It marches the compressible
zero-pressure-gradient thin-layer equations from the leading edge with the
:ref:`Spalart--Allmaras <ref-spalart-allmaras-1994>` model without the trip
term (SA-noft2), in Levy--Lees
variables, for a calorically perfect gas with constant :math:`Pr=0.72`,
turbulent Prandtl number 0.9, :math:`\gamma=1.4` and the Sutherland law
(:math:`S=110.4` K). The SA transport is the standard form, in which density
does not appear (:ref:`Allmaras, Johnson and Spalart, 2012
<ref-allmaras-johnson-spalart-2012>`, section 2.2), its edge value is
:math:`\tilde\nu=3\nu_e`, and the modified vorticity uses the smooth limiter of
the same paper. It uses no correlation, fitted
coefficient or package model. With :math:`\tilde\nu=0` it reproduces the
compressible laminar similarity solution within 0.01%. Across all tabulated
rows, the largest change in friction or a thickness is 0.03% with twice the
streamwise stations, about 0.1% with :math:`\eta` grid ratio 1.02 instead of
1.04, 0.07% with half the first spacing, and below 0.001% with twice the outer
boundary :math:`\eta_{\max}=150`. Each condition takes about 1--1.5 seconds to march with the coupled Newton
iteration and occasional Picard fallback (measured locally).

Its stored results are
``tests/reference_data/viscous_flow/turbulent_sa_boundary_layer.csv``. The
table compares them with this model at the same edge state and wall
temperature. Displacement thickness is integrated through
:math:`\delta_{99}` in both; momentum thickness is integrated to the edge.

.. include:: ../_generated/turbulent_sa_comparison.rst

Friction agrees within 3% at :math:`Re_{x,e}\ge10^7`, as in the SU2 standard-SA
summaries above. At M0.3 the :math:`\delta_{99}` differences also fall within
those summaries at :math:`Re_{x,e}=10^7` and :math:`10^8` and just above them
at :math:`10^6`. On cold supersonic walls they are lower than the summaries at
:math:`Re_{x,e}=10^7`: :math:`-3.0` to :math:`-3.5\%` here, against
:math:`+0.2` to :math:`+1.1\%` for SU2 at different walls and unit Reynolds
numbers. This gap has not been explained, so this reference does not replace
the SU2 comparison.

At :math:`Re_{x,e}=10^6` this model, which is fully turbulent from the leading
edge, gives 5.5--7.4% higher friction than SA, whose turbulence develops from
the edge value of :math:`\tilde\nu`, except at M10. There SA stays nearly
laminar up to about :math:`Re_{x,e}=10^5` and keeps that thicker layer
downstream: friction differs by +0.8% while :math:`\delta_{99}` differs by
:math:`-14.6\%`, narrowing to :math:`-3.2\%` at :math:`10^7`. The reference is a RANS model result, not a
measurement: SA itself gives incompressible :math:`\delta_{99}` about 7--9%
above the measurements cited above.

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
points. Friction and drag retain their existing definitions. Laminar thicknesses
use the closure above; ``compressible=False`` uses the Blasius constants stated earlier.

Wall heat transfer
^^^^^^^^^^^^^^^^^^

With ``compressible=True`` the result also holds the local Stanton number,
heat-transfer coefficient and wall heat flux, from the Reynolds analogy on the
local edge-basis friction:

.. math::

   St=s\,\frac{C_{f,x}}{2},
   \qquad
   h=St\,\rho_eU_ec_p,
   \qquad
   q_w=h\,(T_r-T_w).

:math:`q_w` is in W/m², :math:`h` in W/(m²·K), and :math:`q_w>0` means heat
flows into the wall. :math:`T_r` is the returned ``recovery_temperature``:
the solved adiabatic wall temperature at laminar points and
:math:`r_t=Pr^{1/3}` at turbulent points. An adiabatic wall has :math:`q_w=0`
and a defined :math:`St`. These fields are ``None`` for ``compressible=False``.

**Laminar points** use :math:`s=Pr^{-2/3}`. Against the wall flux of the
thin-layer equations with :math:`\tilde\nu=0` (see the
:ref:`SA reference below <sa-reference-heat-transfer>`), the difference is
within 1% for :math:`M_e=0.3`--10 and :math:`T_w/T_{aw}=0.2`--1.2.

**Turbulent points** use ``reynolds_analogy_factor``, with default
:math:`s=1.14` from the :ref:`generalized Reynolds analogy
<ref-zhang-bi-hussain-she-2014>`, so that the default temperature profile has
a matching wall gradient. Other common values are 1.16 from :ref:`Chi and
Spalding <ref-chi-spalding-1966>` and Colburn's :math:`Pr^{-2/3}\approx1.24`.
:ref:`Cary <ref-cary-1970>` finds the 1.16 fit a valid mean for :math:`M_e<4`--5
near adiabatic conditions, finds too little data and too much scatter above
that to define :math:`s`, and notes that it appears to fall as the wall is
cooled. In a Mach 6 cold-wall DNS, :ref:`Fu et al. <ref-fu-et-al-2021>` find
values within 5% of 1.16. Treat :math:`s` as an uncertainty of several percent
rather than a constant of the model.

.. _sa-reference-heat-transfer:

The :math:`q_w` column of the table in the SA reference section compares
``VAN_DRIEST_II_WALL_TEMPERATURE`` with :math:`s=1.14`. The difference is
:math:`-5.5` to :math:`+3.8\%` at :math:`Re_{x,e}\ge10^7` on the walls tested
(:math:`M_e=5`--10, :math:`T_w/T_r=0.20`--0.76). It reaches :math:`+18\%` at
:math:`Re_{x,e}=10^6`, :math:`T_w/T_r=0.76`, where SA is still developing and
:math:`T_r-T_w` is sensitive to the recovery factor. SA's own recovery factor
falls from about 0.91 at Mach 0.3 to 0.85 at Mach 10 instead of staying at
:math:`Pr^{1/3}=0.896`, so ``St`` should not be compared near :math:`T_w=T_r`.
``ECKERT`` and ``VAN_DRIEST_II`` friction also produce heat-transfer values, but
they carry their own friction differences and have not been compared in this
way. All of this assumes a perfect gas, constant :math:`c_p` and :math:`Pr`, an
isothermal wall and zero pressure gradient.

Velocity and temperature profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``flat_plate_boundary_layer_profile`` returns the velocity, temperature,
density, viscosity and local Mach number at one station as functions of
:math:`y/\delta_{99}`, together with the ``flat_plate_boundary_layer`` result.
The edge state must be physically consistent, because the edge Mach number and
viscosity follow from :math:`U_e`, :math:`T_e`, ``gas`` and ``viscosity_model``.

Laminar stations use the compressible similarity solution, so the wall
gradients of the profile reproduce :math:`\tau_w` and, within 1%, the analogy
:math:`q_w`. Turbulent stations pass the ``VAN_DRIEST_II_WALL_TEMPERATURE``
:math:`\tau_w` and :math:`\delta_{99}` to
``compressible_turbulent_boundary_layer_profile``; see :doc:`Compressible
boundary-layer profiles <../api/viscous_flow>` for that model. The thickness
closure uses the Walz relation, while the profile defaults to the generalized
Reynolds analogy so that its wall gradient matches :math:`s`. Measured locally
against the SA reference on the same :math:`y/\delta_{99}` coordinate, the
largest differences at :math:`Re_{x,e}\ge10^7` are 3.5--9% of :math:`U_e` in
velocity and 2.5--8% of the peak temperature. They grow to about 11% at
:math:`Re_{x,e}=10^6` and :math:`M_e=10`.

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
