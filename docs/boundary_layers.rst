Boundary-layer engineering models
=================================

This page describes the flat-plate boundary-layer and immersed-protrusion
models, including their equations, variables, applicability, and API mapping.
The separate :doc:`compressible_velocity_transformations` page covers
mean-velocity profile transformations in greater depth.

Flat-plate model
----------------

:func:`aerophysics.boundary_layer.flat_plate_boundary_layer` models one side
of a smooth flat plate with a sharp leading edge, zero pressure gradient, and
constant edge conditions. Distance :math:`x` is measured from the leading
edge, and drag is reported per unit spanwise width.

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

For a turbulent boundary layer, the one-fifth-power thickness and
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

Compressibility is opt-in. The recovery temperature is

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

For ``CompressibilityCorrection.ECKERT``,

.. math::

   T^*=0.22T_r+0.28T_e+0.50T_w,
   \qquad
   Re_x^*=Re_x\frac{T_e}{T^*}\frac{\mu_e}{\mu^*}.

Here :math:`\mu^*=\mu(T^*)` is evaluated with the selected dynamic-viscosity
model.
The incompressible thickness and friction correlations are evaluated at
:math:`Re_x^*`.

Van Driest II method
~~~~~~~~~~~~~~~~~~~~

For turbulent ``CompressibilityCorrection.VAN_DRIEST_II``, the
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

The inverse-sine arguments are clipped to :math:`[-1,1]` only to suppress
floating-point excursions at the endpoints. The transformation factors are

.. math::

   F_C=\frac{m}{[\sin^{-1}(\alpha)+\sin^{-1}(\beta)]^2},
   \qquad
   F_\theta=\frac{\mu_e}{\mu_w},
   \qquad
   F_x=\frac{F_\theta}{F_C},
   \qquad
   Re_{x,i}=F_xRe_x.

The equivalent incompressible local coefficient is the positive solution of

.. math::

   \frac{0.242}{\sqrt{C_{f,i}}}
   =0.41+\log_{10}(Re_{x,i}C_{f,i}),
   \qquad
   C_f=\frac{C_{f,i}}{F_C}.

The average coefficient uses the corresponding equation without the local
intercept:

.. math::

   \frac{0.242}{\sqrt{\bar C_{f,i}}}
   =\log_{10}(Re_{x,i}\bar C_{f,i}),
   \qquad
   \bar C_f=\frac{\bar C_{f,i}}{F_C}.

Both monotonic equations are solved on a verified positive bracket with
Brent's method. A missing bracket or failed solve raises
:class:`~aerophysics.exceptions.ModelRangeError`. At the incompressible
adiabatic limit :math:`M\to0`, the implementation assigns
:math:`F_C=F_\theta=F_x=1` directly.

Laminar portions of a VD2 calculation still use Eckert. Compressible
thicknesses are engineering estimates obtained from the incompressible
relations at :math:`F_\theta Re_x`, not the friction Reynolds number
:math:`Re_{x,i}`. ``turbulent_correlation`` is ignored for VD2 turbulent
portions and remains active for ``NONE`` and ``ECKERT``.

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

Mean-property profiles
----------------------

:func:`aerophysics.boundary_layer_profile.compressible_turbulent_boundary_layer_profile`
predicts a smooth-wall, zero-pressure-gradient mean profile from supplied edge
conditions, :math:`\delta_{99}`, and wall shear stress. It combines Spalding's
wall law and the Coles wake with either a Van Driest or Volpiani inverse
transformation. GRA is the default temperature--velocity closure; Walz is
selectable. The result includes velocity and thermodynamic arrays plus
profile-integrated displacement and momentum thicknesses.

This profile API is independent of ``flat_plate_boundary_layer``. The latter
can estimate :math:`\delta_{99}` and wall shear, but the caller passes them
explicitly. The model is limited to fully turbulent, smooth, approximately
zero-pressure-gradient flow. See
:ref:`compressible-velocity-transformations` for equations and examples.

Boundary-layer-immersed protrusions
-----------------------------------

:func:`aerophysics.protrusion.protrusion_drag` estimates direct drag on one
isolated protrusion by integrating the undisturbed local dynamic pressure over
its frontal area. For height :math:`h` and projected width :math:`b(y)`,

.. math::

   A_f=\int_0^hb(y)\,\mathrm{d}y,
   \qquad
   q_\mathrm{eff}
   =\frac{1}{A_f}\int_0^h
    \frac{1}{2}\rho(y)U(y)^2b(y)\,\mathrm{d}y.

The direct drag and shielding factor are

.. math::

   D=C_Dq_\mathrm{eff}A_f,
   \qquad
   \eta_q=\frac{q_\mathrm{eff}}{q_e},
   \qquad
   q_e=\frac{1}{2}\rho_eU_e^2.

The default turbulent velocity approximation is

.. math::

   \frac{U}{U_e}
   =\min\left[\left(\frac{y}{\delta}\right)^{1/7},1\right].

For constant width and density with :math:`h\le\delta`,

.. math::

   \eta_q=\frac{7}{9}\left(\frac{h}{\delta}\right)^{2/7}.

Measured or computed velocity and density profiles may be supplied instead.
For the optional compressible approximation, the Walz relation is

.. math::

   T(y)=T_w+(T_r-T_w)\frac{U}{U_e}
        +(T_e-T_r)\left(\frac{U}{U_e}\right)^2,

.. math::

   T_r=T_e\left[1+Pr^{1/3}\frac{\gamma-1}{2}M_e^2\right],
   \qquad
   \frac{\rho(y)}{\rho_e}=\frac{T_e}{T(y)}.

.. list-table:: Protrusion symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`h`
     - ``height``
     - Protrusion height
     - m
   * - :math:`b(y)`
     - ``frontal_width`` or profile input
     - Projected width at wall-normal position
     - m
   * - :math:`A_f`
     - ``frontal_area``
     - Projected frontal area
     - m²
   * - :math:`C_D`
     - ``drag_coefficient``
     - Supplied free-stream drag coefficient
     - dimensionless
   * - :math:`D`
     - ``direct_drag``
     - Estimated protrusion direct drag
     - N
   * - :math:`\eta_q`
     - ``shielding_factor``
     - Effective-to-edge dynamic-pressure ratio
     - dimensionless

This approximation does not solve the flow around the protrusion. It excludes
wall interference, horseshoe vortices, roughness-induced transition,
downstream skin-friction changes, multiple-element interference, and
shock/protrusion interaction. A transonic case with a single supplied drag
coefficient emits :class:`~aerophysics.exceptions.ApplicabilityWarning`.

>>> from aerophysics import protrusion_drag
>>> drag = protrusion_drag(
...     1.1,
...     height=0.02,
...     frontal_width=0.01,
...     edge_velocity=100.0,
...     edge_density=1.225,
...     boundary_layer_thickness=0.1,
... )
>>> round(drag.shielding_factor, 6)
0.491073
>>> round(drag.direct_drag, 6)
0.661721

Overall applicability
---------------------

The flat-plate model excludes surface roughness, streamwise pressure
gradients, separation, suction, blowing, and natural-transition prediction.
The gas remains calorically perfect and compressible transport properties use
the supplied dynamic-viscosity model. The thickness correlations and protrusion
corrections are engineering estimates rather than resolved profile or
three-dimensional flow solutions.
