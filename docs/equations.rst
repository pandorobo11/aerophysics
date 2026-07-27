Equations and variable definitions
==================================

This page collects the equations implemented by ``aerophysics`` and defines
their symbols, units, assumptions, and API mappings.

Perfect-gas and transport properties
------------------------------------

The :class:`~aerophysics.gas.PerfectGas` model is calorically perfect:
:math:`R` and :math:`\gamma` are constant. The ideal-gas law and the
implemented heat capacities and speed of sound are

.. math::

   p = \rho R T,
   \qquad
   c_p = \frac{\gamma R}{\gamma-1},
   \qquad
   c_v = \frac{R}{\gamma-1},
   \qquad
   a = \sqrt{\gamma R T}.

The :class:`~aerophysics.gas.SutherlandModel` computes dynamic viscosity from

.. math::

   \mu(T)
   = \mu_\mathrm{ref}
     \left(\frac{T}{T_\mathrm{ref}}\right)^{3/2}
     \frac{T_\mathrm{ref}+S}{T+S}.

The U.S. Standard Atmosphere air-conductivity correlation implemented by
:class:`~aerophysics.gas.USSAConductivityModel` is

.. math::

   k(T) = \frac{c_k T^{3/2}}{T + A_k 10^{-B_k/T}}.

``AIR`` uses :math:`R=8314.32/28.9644\ \mathrm{J/(kg\,K)}` and
:math:`\gamma=1.4`. ``AIR_VISCOSITY`` uses
:math:`\mu_\mathrm{ref}=1.7894\times10^{-5}\ \mathrm{Pa\,s}`,
:math:`T_\mathrm{ref}=288.15\ \mathrm{K}`, and :math:`S=110.4\ \mathrm{K}`.
``AIR_CONDUCTIVITY`` uses :math:`c_k=2.64638\times10^{-3}`,
:math:`A_k=245.4\ \mathrm{K}`, and :math:`B_k=12\ \mathrm{K}`.

.. list-table::
   :header-rows: 1
   :widths: 16 26 42 16

   * - Symbol
     - API name
     - Definition
     - SI unit
   * - :math:`p`
     - ``pressure``
     - Static pressure
     - Pa
   * - :math:`\rho`
     - ``density``
     - Gas density
     - kg/m³
   * - :math:`T`
     - ``temperature``
     - Absolute temperature
     - K
   * - :math:`R`
     - ``specific_gas_constant``
     - Specific gas constant
     - J/(kg K)
   * - :math:`\gamma`
     - ``heat_capacity_ratio``
     - Ratio :math:`c_p/c_v`
     - dimensionless
   * - :math:`c_p,\ c_v`
     - ``cp``, ``cv``
     - Constant-pressure and constant-volume specific heats
     - J/(kg K)
   * - :math:`a`
     - ``speed_of_sound``
     - Ideal-gas speed of sound
     - m/s
   * - :math:`\mu`
     - ``dynamic_viscosity``
     - Dynamic viscosity
     - Pa s
   * - :math:`k`
     - ``thermal_conductivity``
     - Thermal conductivity
     - W/(m K)

Standard atmosphere
-------------------

:func:`aerophysics.atmosphere.standard_atmosphere` accepts geometric altitude
:math:`h` and converts it to geopotential altitude :math:`H` using the
effective Earth radius :math:`r_0`:

.. math::

   H = \frac{r_0 h}{r_0+h},
   \qquad
   h = \frac{r_0 H}{r_0-H}.

Within each U.S. Standard Atmosphere layer, temperature is linear in
geopotential altitude:

.. math::

   T = T_b + L_b(H-H_b).

Hydrostatic pressure is evaluated with one of two expressions. For an
isothermal layer, :math:`L_b=0`,

.. math::

   p = p_b
       \exp\left[-\frac{g_0(H-H_b)}{R T_b}\right].

For a nonzero lapse rate,

.. math::

   p = p_b
       \left(\frac{T_b}{T}\right)^{g_0/(R L_b)}.

The remaining atmosphere properties follow from

.. math::

   \rho = \frac{p}{RT},
   \qquad
   g(h) = g_0\left(\frac{r_0}{r_0+h}\right)^2,
   \qquad
   \nu = \frac{\mu}{\rho},
   \qquad
   Pr = \frac{\mu c_p}{k}.

Here :math:`r_0=6\,356\,766\ \mathrm{m}` and
:math:`g_0=9.80665\ \mathrm{m/s^2}`. The temperature, pressure, and lapse-rate
base tables come from U.S. Standard Atmosphere 1976. The implemented geometric
altitude range is :math:`-5\,000 \le h \le 86\,000\ \mathrm{m}`.

.. list-table::
   :header-rows: 1
   :widths: 15 29 40 16

   * - Symbol
     - API name
     - Definition
     - SI unit
   * - :math:`h`
     - ``geometric_altitude``
     - Geometric altitude above the reference ellipsoid
     - m
   * - :math:`H`
     - ``geopotential_altitude``
     - Geopotential altitude used by layer hydrostatics
     - m
   * - :math:`H_b`
     - internal layer table
     - Base geopotential altitude of the active layer
     - m
   * - :math:`T_b,\ p_b`
     - internal layer tables
     - Base temperature and pressure of the active layer
     - K, Pa
   * - :math:`L_b`
     - internal layer table
     - Temperature lapse rate of the active layer
     - K/m
   * - :math:`g`
     - ``gravity``
     - Local gravitational acceleration
     - m/s²
   * - :math:`\nu`
     - ``kinematic_viscosity``
     - Kinematic viscosity
     - m²/s
   * - :math:`Pr`
     - ``prandtl_number``
     - Prandtl number
     - dimensionless

The transport properties :math:`\mu` and :math:`k` and the speed of sound
:math:`a` use the models in `Perfect-gas and transport properties`_.

Isentropic perfect-gas flow
---------------------------

Define the temperature factor

.. math::

   F(M) = 1 + \frac{\gamma-1}{2}M^2.

The total-to-static ratios returned by
:func:`aerophysics.isentropic.isentropic_ratios` are

.. math::

   \frac{T_0}{T} = F(M),
   \qquad
   \frac{p_0}{p} = F(M)^{\gamma/(\gamma-1)},
   \qquad
   \frac{\rho_0}{\rho} = F(M)^{1/(\gamma-1)}.

The inverse ratio functions solve these expressions algebraically for
:math:`M`. At :math:`M=1`, they also define the total-to-critical ratios
:math:`T_0/T^*`, :math:`p_0/p^*`, and :math:`\rho_0/\rho^*`.

The quasi-one-dimensional area--Mach relation is

.. math::

   \frac{A}{A^*}
   = \frac{1}{M}
     \left[
       \frac{2}{\gamma+1}F(M)
     \right]^{(\gamma+1)/[2(\gamma-1)]}.

For :math:`A/A^*>1`, the relation has both a subsonic and a supersonic root.
:func:`aerophysics.isentropic.mach_from_area_ratio` therefore requires an
explicit :class:`~aerophysics.isentropic.MachBranch` and solves the selected
root numerically.

The dimensionless mass-flow parameter and mass flux are

.. math::

   \operatorname{MFP}(M)
   = \frac{\sqrt{\gamma}\,M}
          {F(M)^{(\gamma+1)/[2(\gamma-1)]}},
   \qquad
   \frac{\dot m}{A}
   = \frac{p_0}{\sqrt{R T_0}}\operatorname{MFP}(M).

The choked mass flux is the value at :math:`M=1`.

.. list-table::
   :header-rows: 1
   :widths: 15 27 42 16

   * - Symbol
     - API name
     - Definition
     - SI unit
   * - :math:`M`
     - ``mach``
     - Mach number
     - dimensionless
   * - :math:`T,\ p,\ \rho`
     - static state
     - Static temperature, pressure, and density
     - K, Pa, kg/m³
   * - :math:`T_0,\ p_0,\ \rho_0`
     - ``total_*``
     - Isentropic total state
     - K, Pa, kg/m³
   * - :math:`A`
     - area
     - Local flow area
     - m²
   * - :math:`A^*`
     - critical area
     - Area at which :math:`M=1` for the same mass flow
     - m²
   * - :math:`\dot m/A`
     - ``mass_flux``
     - Mass flow per unit area
     - kg/(m² s)

Normal and oblique shocks
-------------------------

Shock ratios use downstream-over-upstream static quantities. For a normal
shock with upstream Mach number :math:`M_1`,

.. math::

   M_2^2
   = \frac{1+\frac{\gamma-1}{2}M_1^2}
          {\gamma M_1^2-\frac{\gamma-1}{2}},

.. math::

   \frac{p_2}{p_1}
   = 1+\frac{2\gamma}{\gamma+1}(M_1^2-1),
   \qquad
   \frac{\rho_2}{\rho_1}
   = \frac{(\gamma+1)M_1^2}
          {(\gamma-1)M_1^2+2},
   \qquad
   \frac{T_2}{T_1}
   = \frac{p_2/p_1}{\rho_2/\rho_1}.

The total-pressure loss is expressed as

.. math::

   \frac{p_{02}}{p_{01}}
   =
   \left[
     \frac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}
   \right]^{\gamma/(\gamma-1)}
   \left[
     \frac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}
   \right]^{1/(\gamma-1)}.

:func:`aerophysics.shocks.supersonic_pitot_pressure_ratio` returns the
Rayleigh--Pitot ratio

.. math::

   \frac{p_{02}}{p_1}
   =
   \left[\frac{(\gamma+1)M_1^2}{2}\right]^{\gamma/(\gamma-1)}
   \left[
     \frac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}
   \right]^{1/(\gamma-1)}.

For an oblique shock, the theta--beta--Mach relation is

.. math::

   \tan\theta
   =
   2\cot\beta\,
   \frac{M_1^2\sin^2\beta-1}
        {M_1^2(\gamma+\cos 2\beta)+2}.

The shock angle lies between the Mach angle
:math:`\mu_M=\sin^{-1}(1/M_1)` and :math:`\pi/2`. The selected weak or strong
root is found numerically. If :math:`\theta` exceeds the maximum attached-shock
deflection, :class:`~aerophysics.exceptions.NoAttachedShockError` is raised.

Oblique-shock state ratios are obtained from the normal component:

.. math::

   M_{n1}=M_1\sin\beta,
   \qquad
   M_2=\frac{M_{n2}}{\sin(\beta-\theta)}.

The normal-shock equations above are evaluated at :math:`M_{n1}`.

.. list-table::
   :header-rows: 1
   :widths: 16 27 41 16

   * - Symbol
     - API name
     - Definition
     - Unit
   * - :math:`M_1,\ M_2`
     - ``upstream_mach``, ``downstream_mach``
     - Upstream and downstream Mach numbers
     - dimensionless
   * - :math:`M_{n1},\ M_{n2}`
     - ``upstream_normal_mach``, ``downstream_normal_mach``
     - Mach components normal to an oblique shock
     - dimensionless
   * - :math:`\theta`
     - ``deflection_angle``
     - Flow-deflection angle
     - rad
   * - :math:`\beta`
     - ``shock_angle``
     - Shock angle measured from the upstream velocity
     - rad
   * - :math:`p_{01},\ p_{02}`
     - total pressure
     - Upstream and downstream total pressures
     - Pa

The normal and oblique calculations are exposed by
:func:`aerophysics.shocks.normal_shock` and
:func:`aerophysics.shocks.oblique_shock`.

Prandtl--Meyer expansion
------------------------

For a calorically perfect gas, the Prandtl--Meyer function is

.. math::

   \nu(M)
   =
   \sqrt{\frac{\gamma+1}{\gamma-1}}
   \tan^{-1}\left[
     \sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)}
   \right]
   - \tan^{-1}\left(\sqrt{M^2-1}\right).

Its limiting value as :math:`M\rightarrow\infty` is

.. math::

   \nu_\max
   = \frac{\pi}{2}
     \left(\sqrt{\frac{\gamma+1}{\gamma-1}}-1\right).

For a centered expansion through turn angle :math:`\delta`,

.. math::

   \nu(M_2)=\nu(M_1)+\delta.

:func:`aerophysics.expansion.mach_from_prandtl_meyer` solves this equation
numerically for :math:`M_2`. Total pressure and total temperature remain
constant, and downstream-over-upstream static ratios are

.. math::

   \frac{T_2}{T_1} = \frac{F(M_1)}{F(M_2)},
   \qquad
   \frac{p_2}{p_1}
   = \left(\frac{T_2}{T_1}\right)^{\gamma/(\gamma-1)},
   \qquad
   \frac{\rho_2}{\rho_1}
   = \left(\frac{T_2}{T_1}\right)^{1/(\gamma-1)}.

Here :math:`F(M)` is the `Isentropic perfect-gas flow`_ temperature factor.
Angles are in radians, :math:`M_1\ge1`, :math:`\delta\ge0`, and the downstream
Prandtl--Meyer angle must remain below :math:`\nu_\max`. The complete state
change is returned by
:func:`aerophysics.expansion.prandtl_meyer_expansion`.

.. list-table::
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Definition
     - Unit
   * - :math:`\nu`
     - ``prandtl_meyer_angle``
     - Prandtl--Meyer angle
     - rad
   * - :math:`\delta`
     - ``turn_angle``
     - Flow turning angle through the expansion
     - rad
   * - :math:`M_1,\ M_2`
     - ``upstream_mach``, ``downstream_mach``
     - Upstream and downstream Mach numbers
     - dimensionless

Flat-plate boundary-layer thicknesses
-------------------------------------

For the same smooth, zero-pressure-gradient flat plate used by the
skin-friction model, the laminar Blasius thickness relations are

.. math::

   \frac{\delta_{99}}{x} = \frac{5}{\sqrt{Re_x}},
   \qquad
   \frac{\delta^*}{x} = \frac{1.7208}{\sqrt{Re_x}},
   \qquad
   \frac{\theta_m}{x} = \frac{0.664}{\sqrt{Re_x}}.

The turbulent one-fifth-power thickness and the one-seventh-power-profile
estimates are

.. math::

   \delta_{99} = 0.37xRe_x^{-1/5},
   \qquad
   \delta^* = \frac{\delta_{99}}{8},
   \qquad
   \theta_m = \frac{7}{72}\delta_{99}.

.. list-table::
   :header-rows: 1
   :widths: 17 30 37 16

   * - Symbol
     - API name
     - Definition
     - SI unit
   * - :math:`\delta_{99}`
     - ``boundary_layer_thickness``
     - Conventional 99-percent boundary-layer thickness
     - m
   * - :math:`\delta^*`
     - ``displacement_thickness``
     - Displacement thickness
     - m
   * - :math:`\theta_m`
     - ``momentum_thickness``
     - Momentum thickness
     - m

For ``BoundaryLayerRegime.TRANSITIONAL``, local quantities switch from the
laminar to the selected turbulent correlation when :math:`Re_x` exceeds the
caller-supplied :math:`Re_\mathrm{tr}`. Average skin friction preserves the
laminar drag accumulated before transition:

.. math::

   \bar C_{f,\mathrm{mixed}}(Re_x)
   =
   \bar C_{f,t}(Re_x)
   + \frac{Re_\mathrm{tr}}{Re_x}
     \left[
       \bar C_{f,l}(Re_\mathrm{tr})
       - \bar C_{f,t}(Re_\mathrm{tr})
     \right].

No natural-transition criterion is inferred.

.. _flat-plate-compressibility-corrections:

Compressibility corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compressibility is opt-in. The recovery temperature is

.. math::

   T_r
   = T_e\left[
       1+r\frac{\gamma-1}{2}M^2
     \right],
   \qquad
   r_l=\sqrt{Pr},
   \qquad
   r_t=Pr^{1/3}.

If ``wall_temperature`` is omitted, the model uses the adiabatic value
:math:`T_w=T_r`.

For :class:`~aerophysics.boundary_layer.CompressibilityCorrection`\ ``.ECKERT``,
the reference temperature and effective Reynolds number are

.. math::

   T^* = 0.22T_r + 0.28T_e + 0.50T_w,
   \qquad
   Re_x^*
   = Re_x\frac{T_e}{T^*}\frac{\mu_e}{\mu^*}.

Here :math:`\mu^*=\mu(T^*)` is evaluated with the selected Sutherland model.
The incompressible thickness and skin-friction correlations are then evaluated
at :math:`Re_x^*`.

For turbulent
:class:`~aerophysics.boundary_layer.CompressibilityCorrection`\ ``.VAN_DRIEST_II``,
the Reynolds-number transform is

.. math::

   F_\theta = \frac{\mu_e}{\mu_w},
   \qquad
   Re_{x,\mathrm{eff}} = F_\theta Re_x.

The implemented friction transform defines

.. math::

   a_v
   =
   \sqrt{
     \frac{1}{2}Pr^{1/3}(\gamma-1)M^2\frac{T_e}{T_w}
   },
   \qquad
   b_v = \frac{T_r}{T_w}-1,

.. math::

   \alpha_v
   = \frac{2a_v^2-b_v}{\sqrt{b_v^2+4a_v^2}},
   \qquad
   \beta_v
   = \frac{b_v}{\sqrt{b_v^2+4a_v^2}},

.. math::

   F_c
   =
   \frac{T_r/T_e-1}
        {\left[
          \sin^{-1}(\alpha_v)+\sin^{-1}(\beta_v)
        \right]^2},
   \qquad
   C_f
   =
   \frac{C_{f,\mathrm{inc}}(Re_{x,\mathrm{eff}})}{F_c}.

Laminar portions still use the Eckert method. Thicknesses are engineering
estimates obtained by evaluating the incompressible thickness relations at
the effective Reynolds number; they are not transformed velocity-profile
solutions.

.. list-table::
   :header-rows: 1
   :widths: 17 29 38 16

   * - Symbol
     - API name
     - Definition
     - Unit
   * - :math:`T_e`
     - ``edge_temperature``
     - Boundary-layer edge temperature
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
   * - :math:`Re_{x,\mathrm{eff}}`
     - ``effective_reynolds_number``
     - Reynolds number supplied to the selected correlation
     - dimensionless
   * - :math:`F_c`
     - internal
     - Van Driest II skin-friction factor
     - dimensionless

.. _flat-plate-skin-friction-coefficients:

Flat-plate skin-friction coefficients
-------------------------------------

The flat-plate correlations describe the skin friction on one side of a smooth
plate with a sharp leading edge, zero pressure gradient, and constant
boundary-layer edge conditions. Distance is measured from the leading edge,
and drag is reported per unit spanwise width.

The Reynolds number and edge dynamic pressure at distance :math:`x` are

.. math::

   Re_x = \frac{\rho_e U_e x}{\mu_e},
   \qquad
   q_e = \frac{1}{2}\rho_e U_e^2.

The local skin-friction coefficient relates wall shear stress to edge dynamic
pressure:

.. math::

   C_{f,x} = \frac{\tau_w(x)}{q_e}.

The average coefficient from the leading edge to :math:`x` is

.. math::

   \bar C_f(x)
   = \frac{1}{x}\int_0^x C_{f,\xi}\,\mathrm{d}\xi
   = \frac{D'(x)}{q_e x},

where :math:`D'` is the accumulated one-sided drag per unit width. These
definitions assume that :math:`q_e` is constant along the plate.

The coefficients on this page are flat-plate *skin-friction coefficients*.
They are not the Darcy or Fanning friction factors used for internal pipe
flow.

Variables
^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 12 29 43 16

   * - Symbol
     - API name
     - Definition
     - SI unit
   * - :math:`x`
     - ``distance``
     - Distance from the sharp leading edge
     - m
   * - :math:`U_e`
     - ``edge_velocity``
     - Boundary-layer edge velocity
     - m/s
   * - :math:`\rho_e`
     - ``edge_density``
     - Boundary-layer edge density
     - kg/m³
   * - :math:`\mu_e`
     - ``edge_dynamic_viscosity``
     - Boundary-layer edge dynamic viscosity
     - Pa s
   * - :math:`Re_x`
     - ``reynolds_number``
     - Reynolds number based on :math:`x`
     - dimensionless
   * - :math:`q_e`
     - computed internally
     - Boundary-layer edge dynamic pressure
     - Pa
   * - :math:`C_{f,x}`
     - ``local_skin_friction_coefficient``
     - Local skin-friction coefficient at :math:`x`
     - dimensionless
   * - :math:`\bar C_f`
     - ``average_skin_friction_coefficient``
     - Skin-friction coefficient averaged from 0 to :math:`x`
     - dimensionless
   * - :math:`\tau_w`
     - ``wall_shear_stress``
     - Local wall shear stress
     - Pa
   * - :math:`D'`
     - ``drag_per_unit_width``
     - One-sided drag accumulated from 0 to :math:`x`, per unit width
     - N/m

Laminar correlation
^^^^^^^^^^^^^^^^^^^

The Blasius flat-plate solution gives

.. math::

   C_{f,x} = \frac{0.664}{\sqrt{Re_x}},
   \qquad
   \bar C_f = \frac{1.328}{\sqrt{Re_x}}.

Select this correlation with ``BoundaryLayerRegime.LAMINAR`` from
:class:`~aerophysics.boundary_layer.BoundaryLayerRegime`.

Turbulent power-law correlation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The one-fifth-power correlation is

.. math::

   C_{f,x} = 0.0592 Re_x^{-1/5},
   \qquad
   \bar C_f = 0.074 Re_x^{-1/5}.

Select it with ``TurbulentCorrelation.POWER_LAW`` from
:class:`~aerophysics.boundary_layer.TurbulentCorrelation`.

Schlichting turbulent correlation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default turbulent model uses the Schlichting average correlation

.. math::

   \bar C_f = \frac{0.455}{\left(\log_{10} Re_x\right)^{2.58}}.

For constant edge conditions, the local coefficient follows from
:math:`C_{f,x}=\mathrm{d}(x\bar C_f)/\mathrm{d}x`. Because
:math:`Re_x` is proportional to :math:`x`, the implemented expression is

.. math::

   C_{f,x}
   = \bar C_f\left(1 - \frac{2.58}{\ln Re_x}\right).

Select it with ``TurbulentCorrelation.SCHLICHTING`` from
:class:`~aerophysics.boundary_layer.TurbulentCorrelation`.

Applicability
^^^^^^^^^^^^^

The caller prescribes whether the boundary layer is laminar, turbulent, or
transitional; the API does not predict natural transition. Turbulent use
outside :math:`5\times10^5 \le Re_x \le 10^9` emits
:class:`~aerophysics.exceptions.ApplicabilityWarning`.

This section states the base incompressible correlations. The specified
transition and compressible transforms are defined in
:ref:`flat-plate-compressibility-corrections`.

Worked example
^^^^^^^^^^^^^^

Consider a turbulent plate at :math:`x=1\,\mathrm{m}` with
:math:`U_e=100\,\mathrm{m/s}`, :math:`\rho_e=1.225\,\mathrm{kg/m^3}`, and
:math:`\mu_e=1.7894\times10^{-5}\,\mathrm{Pa\,s}`. The Schlichting
correlation gives :math:`Re_x=6.84587\times10^6`,
:math:`C_{f,x}=0.002670314`, :math:`\bar C_f=0.003193859`,
:math:`\tau_w=16.3557\,\mathrm{Pa}`, and
:math:`D'=19.5624\,\mathrm{N/m}`.

The same calculation through the public API is

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

See :func:`aerophysics.boundary_layer.flat_plate_boundary_layer` for the full
input contract and returned
:class:`aerophysics.boundary_layer.FlatPlateBoundaryLayerResult`.

Integrated flight condition
---------------------------

:class:`aerophysics.flight.FlightCondition` combines the standard atmosphere
with either a supplied Mach number or velocity. At the selected geometric
altitude,

.. math::

   V = Ma,
   \qquad
   M = \frac{V}{a},
   \qquad
   q = \frac{1}{2}\rho V^2.

The Reynolds number per unit length and, when a positive characteristic
length :math:`L` is supplied, the dimensionless Reynolds number are

.. math::

   \frac{Re_L}{L} = \frac{\rho V}{\mu},
   \qquad
   Re_L = \frac{\rho V L}{\mu}.

Total conditions use the isentropic ratios:

.. math::

   T_0 = T F(M),
   \qquad
   p_0 = p F(M)^{\gamma/(\gamma-1)},
   \qquad
   \rho_0 = \rho F(M)^{1/(\gamma-1)}.

The atmosphere supplies :math:`T`, :math:`p`, :math:`\rho`, :math:`a`, and
:math:`\mu`; :math:`F(M)` is defined in `Isentropic perfect-gas flow`_.

.. list-table::
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Definition
     - SI unit
   * - :math:`V`
     - ``velocity``
     - Flight speed
     - m/s
   * - :math:`q`
     - ``dynamic_pressure``
     - Dynamic pressure
     - Pa
   * - :math:`L`
     - ``characteristic_length``
     - Optional Reynolds-number reference length
     - m
   * - :math:`Re_L/L`
     - ``reynolds_number_per_length``
     - Reynolds number per unit length
     - 1/m
   * - :math:`Re_L`
     - ``reynolds_number``
     - Reynolds number based on :math:`L`
     - dimensionless
   * - :math:`T_0,\ p_0,\ \rho_0`
     - ``total_temperature``, ``total_pressure``, ``total_density``
     - Isentropic total state
     - K, Pa, kg/m³

Unit conversions
----------------

The conversion functions in :mod:`aerophysics.units` apply explicit scale or
affine transformations. Forward and inverse functions use the same exact
constants.

.. list-table::
   :header-rows: 1
   :widths: 35 43 22

   * - Conversion
     - Equation
     - API
   * - International feet to metres
     - :math:`x_\mathrm{m}=0.3048x_\mathrm{ft}`
     - ``feet_to_meters``
   * - Knots to metres per second
     - :math:`V_\mathrm{m/s}=(1852/3600)V_\mathrm{kt}`
     - ``knots_to_meters_per_second``
   * - Degrees Fahrenheit to kelvin
     - :math:`T_\mathrm{K}=(T_\mathrm{^\circ F}-32)(5/9)+273.15`
     - ``fahrenheit_to_kelvin``
   * - Pounds-force per square inch to pascals
     - :math:`p_\mathrm{Pa}=6894.757293168\,p_\mathrm{psi}`
     - ``psi_to_pascals``
   * - Pounds-force per square foot to pascals
     - :math:`p_\mathrm{Pa}=47.8802589803\,p_\mathrm{psf}`
     - ``psf_to_pascals``
   * - Pounds mass to kilograms
     - :math:`m_\mathrm{kg}=0.45359237\,m_\mathrm{lbm}`
     - ``pounds_mass_to_kilograms``
   * - Slugs to kilograms
     - :math:`m_\mathrm{kg}=14.5939029372\,m_\mathrm{slug}`
     - ``slugs_to_kilograms``
   * - Degrees to radians
     - :math:`\theta_\mathrm{rad}=(\pi/180)\theta_\mathrm{deg}`
     - ``degrees_to_radians``

The inverse APIs divide by the listed scale factor. Temperature uses
:math:`T_\mathrm{^\circ F}=(T_\mathrm{K}-273.15)(9/5)+32`. Values below
absolute zero are rejected by the temperature converters.
