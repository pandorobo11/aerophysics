.. _aerodynamic-heating-model:

Aerodynamic heating
===================

These engineering estimates use SI units. Heat flux :math:`q_w` is
in W/m² and is positive into the wall; a sufficiently hot wall produces a
negative value. All enthalpies are in J/kg, and angles are in radians.

Equations
---------

For a calorically perfect gas, :math:`h_\infty=c_pT_\infty`,
:math:`h_w=c_pT_w`, :math:`H_s=h_\infty+V_\infty^2/2=c_pT_0`, and
:math:`\rho_\infty=p_\infty/(RT_\infty)`. The sphere stagnation pressure
:math:`p_s` is the total pressure downstream of a normal shock.

The Fay--Riddell stagnation-point correlation, including the Lewis-number and
dissociation-enthalpy term, can be written as

.. math::

   q_w=0.763 Pr_w^{-0.6}(\rho_w\mu_w)^{0.1}(\rho_e\mu_e)^{0.4}
       \left[1+(Le^{0.52}-1)\frac{h_D}{h_e}\right](h_e-h_w)
       \sqrt{\frac{1}{r_n}\sqrt{\frac{2(p_e-p_\infty)}{\rho_e}}}.

Here :math:`Le` is the Lewis number and :math:`h_D` is the dissociation
enthalpy per unit mass at the boundary-layer edge. When :math:`Le=1`, the
bracket is one. With perfect-gas enthalpy and the Newtonian velocity gradient,
this gives the implemented simplified form without a dissociation term:

.. math::

   q_w=0.763 Pr^{-0.6}(\rho_e\mu_e)^{0.4}(\rho_w\mu_w)^{0.1}
       \sqrt{(du_e/dx)_s}\,c_p(T_0-T_w),
   \quad (du_e/dx)_s=\frac{1}{r_n}\sqrt{\frac{2(p_s-p_\infty)}{\rho_e}}.

Here :math:`\rho_e=p_s/(RT_0)`, :math:`\rho_w=p_s/(RT_w)`, and the
viscosities are evaluated at :math:`T_0` and :math:`T_w`, respectively.
The original correlation uses the wall Prandtl number :math:`Pr_w`. This API
uses a constant ``prandtl_number``, assuming :math:`Pr_w=Pr`.

The :ref:`TFAWS lecture <ref-tfaws-2012-aerothermodynamics>` presents the
:ref:`Fay and Riddell (1958) <ref-fay-riddell-1958>` correlation with
coefficient 0.763, wall Prandtl number :math:`Pr_w`, and velocity-gradient
factor :math:`[(du_e/dx)_s]^{0.5}`. :ref:`West and Brandis
<ref-west-brandis-ntrs>` give the same property exponents and dissociation
term, rounding the coefficient to 0.76. Their equation (1) omits the outer
square root on the velocity gradient in print; this is treated as a printing
omission because dimensional consistency requires it. Specifically, the
combined property factor has units kg m⁻² s⁻⁰·⁵ and the enthalpy difference
has units m² s⁻², so a factor with units s⁻⁰·⁵ is required to obtain W/m².

Tauber's air correlations, equations (40) and (41), are

.. math::

   q_{w,s}=1.83\times10^{-4}\sqrt{\rho_\infty/r_n}\,V_\infty^3(1-h_w/H_s),

   h_{aw}=h_\infty+\frac{V_\infty^2}{2}(1-0.18\sin^2\Lambda),

   q_{w,c}=1.29\times10^{-4}\sqrt{\rho_\infty/r}\,
   (1-0.18\sin^2\Lambda)V_\infty^3(1-h_w/h_{aw})\cos\Lambda.

.. list-table:: Symbols and API mapping
   :header-rows: 1

   * - Symbol
     - API name and SI unit
   * - :math:`M, T_\infty, p_\infty`
     - ``mach``, ``temperature`` [K], ``pressure`` [Pa]
   * - :math:`r_n, r, \Lambda`
     - ``nose_radius``, ``radius`` [m], ``sweep_angle`` [rad]
   * - :math:`T_w, Pr`
     - ``wall_temperature`` [K], ``prandtl_number`` (wall value, assumed constant)
   * - :math:`q_w, V_\infty, \rho_\infty`
     - ``heat_flux`` [W/m²], ``freestream_velocity`` [m/s], ``freestream_density`` [kg/m³]
   * - :math:`T_0, p_s, H_s, h_{aw}`
     - ``stagnation_temperature`` [K], ``stagnation_pressure`` [Pa], ``stagnation_enthalpy`` and ``recovery_enthalpy`` [J/kg]
   * - :math:`(du_e/dx)_s`
     - ``velocity_gradient`` [1/s], Fay--Riddell only

Assumptions and applicability
-----------------------------

Fay--Riddell uses a perfect gas, Lewis number one and no dissociation.
Its Newtonian velocity gradient
is intended for hypersonic flow; :ref:`Tauber <ref-tauber-1989>` reports it
is accurate for a cylinder but 7.7% low for a sphere. Above 2000 K stagnation
temperature the model warns because perfect-gas edge properties omit
potential dissociation. The source gives no Mach threshold for the Newtonian
gradient, so the API does not issue a low-speed warning for this approximation.

The Tauber correlations are for air at high speeds where
:math:`h_\infty\ll V_\infty^2/2`. Figure 3 of TP-2914 shows four correlations
within the approximately ±25% data spread up to 14 km/s. The derivation
assumes wall temperatures below about 4000 K; the API warns above 14 km/s or
4000 K. It also warns when the freestream enthalpy ratio
:math:`h_\infty/(V_\infty^2/2)=2/((\gamma-1)M^2)` exceeds 0.1, corresponding
to :math:`M<7.07` approximately for :math:`\gamma=1.4`. This 0.1 threshold
is a practical interpretation of :math:`\ll 1`, approximately the ratio at
the lowest speed of about 2 km/s in figure 3 for air at a representative
freestream temperature of 200 K. Tauber notes that using the velocity gradient
from his reference 8
would change the sphere coefficient from 1.83 to 1.90.

Equation (41) assumes equal sphere and cylinder velocity gradients. At zero
sweep its ratio to equation (40) is :math:`1.29/1.83\approx1/\sqrt{2}`;
as :math:`\Lambda\to\pi/2`, the cosine drives heat flux to zero. Only the
laminar stagnation line is represented: attachment-line transition and
flat-plate-like heating on highly swept real wing leading edges (TP-2914
equation 49) are not included. All models require a supersonic freestream.

Turbulent swept-cylinder attachment line
----------------------------------------

:func:`aerophysics.heating.swept_cylinder_turbulent_heat_flux` with
:class:`aerophysics.heating.AttachmentLineHeatingModel` ``POLL`` estimates
local convective heat flux at the front attachment line of a long circular
swept cylinder, with negligible end effects. The user explicitly selects a
smooth-wall, **fully turbulent** assumption. There is no Reynolds-number
transition test, automatic switching, or intermittency correction.
Circumferential distributions, area averages, finite-length ends, shock
interactions, separation, dissociation, recombination, catalysis, radiation,
ablation, and turbulent sphere heating are outside this model's scope.

.. list-table:: Attachment-line symbols and API mapping
   :header-rows: 1
   :widths: 30 70

   * - Symbol
     - Meaning, API name and SI unit
   * - :math:`W_e=V_A=V_e`
     - Axial boundary-layer edge velocity along the attachment line;
       ``edge_velocity`` [m/s]. The normal component alone is stopped.
   * - :math:`A=(du_e/ds)_{AL}`
     - Circumferential edge-velocity gradient, not a wall-normal gradient;
       ``velocity_gradient`` [1/s].
   * - :math:`R_{cyl}, D=2R_{cyl}`
     - Cylinder radius ``radius`` [m]; literature diameter :math:`D` [m].
   * - :math:`\Lambda`
     - ``sweep_angle`` [rad], zero when the cylinder axis is perpendicular
       to the freestream; require :math:`0<\Lambda<\pi/2`.
   * - :math:`T_e,p_e,\rho_e,M_e`
     - ``edge_temperature`` [K], ``edge_pressure`` [Pa], ``edge_density``
       [kg/m³], ``edge_mach``.
   * - :math:`T_{aw},T^*,T_w`
     - ``recovery_temperature``, ``reference_temperature``,
       ``wall_temperature`` [K].
   * - :math:`\overline R^*=W_e/\sqrt{\nu^*A}`
     - ``reynolds_number``; distinct from diameter-based :math:`Re_D`
       and flat-plate :math:`Re_x`.
   * - :math:`St^*=q_w/[\rho^*W_ec_p(T_{aw}-T_w)]`
     - ``stanton_number``; :math:`\rho^*=p_e/(RT^*)` uses the reference
       temperature and attachment-line edge pressure.
   * - :math:`q_w`
     - ``heat_flux`` [W/m²], positive from fluid into wall.

The external edge state uses a modified-Newtonian engineering approximation:
a normal shock and isentropic deceleration stop only the crossflow component,
while conserving the axial component. With :math:`\Pi_0(M_\perp)` the
normal-shock total-pressure ratio,

.. math::

   V_\infty=M\sqrt{\gamma RT_\infty},\qquad
   W_e=V_\infty\sin\Lambda,\qquad M_\perp=M\cos\Lambda,

   F=1+\frac{\gamma-1}{2}M_\perp^2,\qquad T_e=T_\infty F,

   p_e=p_\infty F^{\gamma/(\gamma-1)}\Pi_0(M_\perp),\qquad
   \rho_e=\frac{p_e}{RT_e},\qquad
   A=\frac{1}{R_{cyl}}\sqrt{\frac{2(p_e-p_\infty)}{\rho_e}}.

This requires :math:`M\cos\Lambda>1`; satisfying that inequality alone does
not ensure an accurate external-flow approximation. The full freestream
Mach number is not used in the normal shock, and :math:`W_e` is not stopped.
The edge state conserves :math:`c_pT_e+W_e^2/2=c_pT_\infty+V_\infty^2/2`.
``gas`` supplies perfect-gas properties for these conversions, but the
correlation is based on air data and has no validation for arbitrary gases.

The implemented **Poll reference-temperature correlation as reproduced in
NASA CR-181893** uses equations (5.7a,b), (5.9), and (5.10a,b) of
:ref:`the report <ref-nasa-cr-181893>`:

.. math::

   T_{aw}=T_e+\frac{rW_e^2}{2c_p},\qquad
   T^*=T_e+0.1(T_w-T_e)+0.6(T_{aw}-T_e),

   \rho^*=\frac{p_e}{RT^*},\qquad
   \mu^*=\mu(T^*),\qquad \nu^*=\frac{\mu^*}{\rho^*},\qquad
   \overline R^*=\frac{W_e}{\sqrt{\nu^*A}},

   St^*=\frac{0.0689}{2Pr^{2/3}}(\overline R^*)^{-0.42},\qquad
   q_w=St^*\rho^*W_ec_p(T_{aw}-T_w),\qquad
   M_e=\frac{W_e}{\sqrt{\gamma RT_e}}.

The default ``recovery_factor`` :math:`r=0.89` comes from the description
of Poll's method in :ref:`Bénard et al. (2006) <ref-benard-2006>`, equation
(8); CR-181893 does not specify this value. ``prandtl_number`` defaults to
0.72. Neither an additional factor 1.24 nor any further Prandtl correction
is applied.

The reference-density normalization follows algebraically from equivalence
with CR-181893 (5.12), the freestream-normalized form. Set
:math:`U_{1,0}=AD/(V_\infty\cos\Lambda)` and
:math:`Re_D=\rho_\infty V_\infty D/\mu_\infty`. Then

.. math::

   St_\infty=\frac{q_w}{\rho_\infty V_\infty c_p(T_{aw}-T_w)}
   =\frac{0.0689}{2Pr^{2/3}}
   \left(\frac{\rho_eT_e}{\rho_\infty T^*}\right)^{0.79}
   \left(\frac{\mu^*U_{1,0}}{\mu_\infty Re_D}\right)^{0.21}
   \sin^{0.58}\Lambda\cos^{0.21}\Lambda.

Here :math:`\rho^*=\rho_eT_e/T^*` follows from constant attachment-line
pressure; the density exponent is :math:`1-0.21=0.79` and the sine exponent
is :math:`1-0.42=0.58`. The cosine factor appears when the actual gradient
coefficient :math:`AD/V_\infty` is written as :math:`U_{1,0}\cos\Lambda`.
The factor :math:`1/2` is supported by CR-181893 (5.11): its laminar
coefficient :math:`1.141=2\times0.5705` matches the axial wall gradient of
swept Hiemenz flow. With :math:`C_f=\tau_w/(\rho W_e^2/2)`, the turbulent
analogy is :math:`St^*/(C_f^*/2)=Pr^{-2/3}`, with
:math:`C_f^*/2\simeq0.0345(\overline R^*)^{-0.42}`. At :math:`Pr=0.72`
the ratio is about 1.245, consistent with Poll's 1.24 ratio reported by
Bénard et al.

``ApplicabilityWarning`` is emitted once per call for each exceeded range
(including array calls). The literature-based checks use Bénard et al.'s
272-point turbulent heat-transfer evaluation database:
:math:`0.35\le M_e\le6.89` and :math:`0.14\le T_w/T_0\le0.96`, where
:math:`T_0=T_\infty[1+(\gamma-1)M^2/2]=T_e+W_e^2/(2c_p)`.
These are not confirmed derivation ranges of Poll's original correlation
or of CR-181893. On this database, Bénard et al. report a standard deviation
of 26.7% for Poll's reference-temperature method. Being within the ranges
is not an accuracy guarantee.

A separate warning above :math:`T_0=2000` K is an implementation-policy
threshold, shared with the existing Fay--Riddell implementation, rather
than a literature database limit. Perfect-gas properties omit dissociation;
a numerical prediction at high Mach number is not validation of
high-enthalpy reacting flow. All quantities use SI units and radians.
A wall hotter than :math:`T_{aw}` returns negative heat flux without clipping
or an additional hot-wall warning. Applicability in that cooling direction
is not established by the literature; the database-range warnings still apply.

Comparison with :ref:`Poll (1981) <ref-poll-1981>` itself remains incomplete;
identity with the original paper is not assumed. Bénard et al.'s separate
correlation, equations (11)--(12), is **not implemented**: the ``e0`` subscript
in its Stanton density :math:`\rho_{e0}` is undefined in the nomenclature,
and the high-Reynolds-number Stanton plateau in figure 2, titled "Typical
transition case", reads approximately :math:`St=2.0\times10^{-3}` at
:math:`\overline R^*\approx330`--380. The figure conditions are
:math:`\Lambda=60^\circ`, :math:`M_\infty=7.14`,
:math:`T_w/T_0\approx0.38`, and :math:`M_e\approx3.28`.
Evaluating equation (11) at :math:`\overline R^*=350`, assuming
:math:`T_0=750` K within the paper's 650--800 K wind-tunnel range, and using
this library's air viscosity model to convert from the :math:`T^*` to the
:math:`\overline T` Reynolds-number basis, gives
:math:`St\approx4.6\times10^{-3}` (about 2.3 times
the figure reading). This is inconsistent with the reported 11.7% standard
deviation. No inferred coefficient correction is used.

Tests verify analytical relations, algebraic equivalence with CR-181893
(5.12), direct edge-state calculations, and numerical regression values.
Independent comparison against experimental data has not yet been performed.
The existing laminar APIs and their angle conventions are unchanged.
The following example has :math:`T_0=1857.3` K and lies inside both database
ranges, so it emits no applicability warning:

.. doctest::

   >>> from math import pi
   >>> from aerophysics import swept_cylinder_turbulent_heat_flux
   >>> result = swept_cylinder_turbulent_heat_flux(6, 226.5, 1197, 0.1, pi/6, wall_temperature=300)
   >>> round(result.heat_flux)
   450252

API and verification
--------------------

:func:`aerophysics.heating.sphere_stagnation_heat_flux` selects Fay--Riddell
or Tauber; :func:`aerophysics.heating.swept_cylinder_heat_flux` evaluates the
swept cylinder. The Tauber equations have been checked against TP-2914. The
original :ref:`Fay and Riddell paper <ref-fay-riddell-1958>` is available only
as a paid publisher version and has not been obtained. The correlation was
compared with two public NASA sources, the TFAWS lecture and West and Brandis;
the latter's printed omission of the outer square root is discussed above.
Tests check analytical relations, model comparisons, and regression values,
not physical validation.

.. doctest::

   >>> import warnings
   >>> from aerophysics import sphere_stagnation_heat_flux
   >>> from aerophysics.exceptions import ApplicabilityWarning
   >>> with warnings.catch_warnings():
   ...     warnings.simplefilter("ignore", ApplicabilityWarning)
   ...     result = sphere_stagnation_heat_flux(10, 226.5, 1197, 0.1, wall_temperature=300)
   >>> round(result.heat_flux)
   1986435
