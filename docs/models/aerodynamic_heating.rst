.. _aerodynamic-heating-model:

Aerodynamic heating
===================

These engineering estimates give the convective heat flux at the stagnation
point of a sphere and along the stagnation line (attachment line) of a swept
cylinder. The Fay--Riddell and Tauber correlations assume a laminar boundary
layer. The Poll correlation for the swept-cylinder attachment line assumes a
fully turbulent boundary layer; see
:ref:`aerodynamic-heating-turbulent-attachment-line`. All use SI units. Heat flux :math:`q_w` is in W/m² and is positive
into the wall; a wall hotter than the gas gives a negative value. Enthalpies
are in J/kg, and angles are in radians. For worked examples, see :doc:`../guides/aerodynamic_heating`.

Physical background
-------------------

At high speed, the air that meets a body is slowed almost to rest in a thin
layer near the surface. Its kinetic energy turns into heat, so the gas next to
the wall approaches the total temperature :math:`T_0`. At Mach 10 and 220 K
this is several thousand kelvin. The boundary layer then conducts heat into
the cooler wall.

The heating is usually highest at the **stagnation point** at the front of a
blunt nose, where the boundary layer is thinnest. Two simple scaling results
explain the main trends of the correlations below:

* Heat flux grows roughly with :math:`\sqrt{\rho_\infty}\,V_\infty^3`: the
  energy flux of the oncoming air rises with the cube of speed.
* Heat flux falls with :math:`1/\sqrt{r_n}`: a larger nose radius thickens
  the boundary layer and lowers the heating. This is why re-entry vehicles
  use blunt noses.

The driving potential is the enthalpy difference between the gas near the
surface and the wall, :math:`h_e-h_w` or :math:`H_s-h_w`. A hot wall therefore
receives less heat.

For a swept cylinder, such as a wing leading edge, only the velocity component
normal to the cylinder is brought to rest. Sweeping the leading edge back by
:math:`\Lambda` reduces the heating, roughly in proportion to
:math:`\cos\Lambda`. The flow along the attachment line of a strongly swept
cylinder can become turbulent, which raises the heating above the laminar
estimate. This library does not predict that transition; the user chooses
the laminar or the fully turbulent model.

Equations
---------

The freestream state is taken from Mach number, static temperature, and
static pressure with a calorically perfect gas: :math:`h_\infty=c_pT_\infty`,
:math:`h_w=c_pT_w`, :math:`H_s=h_\infty+V_\infty^2/2=c_pT_0`, and
:math:`\rho_\infty=p_\infty/(RT_\infty)`. The sphere stagnation pressure
:math:`p_s` is the total pressure behind a normal shock.

Fay--Riddell sphere stagnation point
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The general Fay--Riddell correlation, including the Lewis-number and
dissociation-enthalpy term, is

.. math::

   q_w=0.763 Pr_w^{-0.6}(\rho_w\mu_w)^{0.1}(\rho_e\mu_e)^{0.4}
       \left[1+(Le^{0.52}-1)\frac{h_D}{h_e}\right](h_e-h_w)
       \sqrt{\frac{1}{r_n}\sqrt{\frac{2(p_e-p_\infty)}{\rho_e}}}.

Here :math:`Le` is the Lewis number (the ratio of mass diffusion to heat
conduction) and :math:`h_D` is the dissociation enthalpy per unit mass at the
boundary-layer edge. The last factor is the square root of the velocity
gradient at the stagnation point, which describes how fast the flow
accelerates away from that point. It is estimated with Newtonian theory, a
hypersonic approximation that obtains the surface pressure from the angle at
which the flow meets the surface.

For a gas that does not dissociate, :math:`h_D=0` and the bracket equals one.
With perfect-gas enthalpy, this gives the form implemented here:

.. math::

   q_w=0.763 Pr^{-0.6}(\rho_e\mu_e)^{0.4}(\rho_w\mu_w)^{0.1}
       \sqrt{(du_e/dx)_s}\,c_p(T_0-T_w),
   \quad (du_e/dx)_s=\frac{1}{r_n}\sqrt{\frac{2(p_s-p_\infty)}{\rho_e}}.

Here :math:`\rho_e=p_s/(RT_0)`, :math:`\rho_w=p_s/(RT_w)`, and the
viscosities are evaluated at :math:`T_0` and :math:`T_w` with
``viscosity_model``. The original correlation uses the wall Prandtl number
:math:`Pr_w`; this API uses a constant ``prandtl_number`` and assumes
:math:`Pr_w=Pr`.

The :ref:`TFAWS lecture <ref-tfaws-2012-aerothermodynamics>` presents the
:ref:`Fay and Riddell (1958) <ref-fay-riddell-1958>` correlation with
coefficient 0.763, wall Prandtl number :math:`Pr_w`, and the square root of
the velocity gradient. :ref:`West and Brandis <ref-west-brandis-ntrs>` give
the same exponents and dissociation term with the coefficient rounded to
0.76. Their printed equation (1) omits the outer square root on the velocity
gradient. This is treated as a printing error, because the square root is
required for the result to have units of W/m².

Tauber sphere and swept cylinder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`Tauber (1989) <ref-tauber-1989>` gives simpler air correlations,
equations (40) and (41), in which the heat flux depends directly on
freestream density and velocity:

.. math::

   q_{w,s}=1.83\times10^{-4}\sqrt{\rho_\infty/r_n}\,V_\infty^3(1-h_w/H_s),

   h_{aw}=h_\infty+\frac{V_\infty^2}{2}(1-0.18\sin^2\Lambda),

   q_{w,c}=1.29\times10^{-4}\sqrt{\rho_\infty/r}\,
   (1-0.18\sin^2\Lambda)V_\infty^3(1-h_w/h_{aw})\cos\Lambda.

The coefficients carry the units needed for SI inputs and are calibrated for
air. The ``gas`` argument supplies perfect-gas properties for velocity,
density, and enthalpy conversions, including the wall-enthalpy corrections.
It does not change the fitted coefficients or establish validity for other
gases. The ``prandtl_number`` and ``viscosity_model`` arguments are not used
by the Tauber correlations.

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

This section covers the laminar Fay--Riddell and Tauber models. They require
a supersonic freestream and represent only the stagnation point or stagnation
line. The fully turbulent Poll model has its own assumptions and limits,
described in :ref:`aerodynamic-heating-turbulent-attachment-line`.

**Fay--Riddell** uses a perfect gas, Lewis number one, and no dissociation.
The Newtonian velocity gradient is intended for hypersonic flow; Tauber
reports that it is accurate for a cylinder but 7.7% low for a sphere. The
source gives no Mach-number threshold for this approximation, so no low-speed
warning is issued. Above 2000 K stagnation temperature the call emits
:class:`~aerophysics.exceptions.ApplicabilityWarning`, because real air
begins to dissociate and the perfect-gas edge properties no longer describe
it.

**Tauber** correlations are for air at high speed, where the freestream
enthalpy is small compared with the kinetic energy,
:math:`h_\infty\ll V_\infty^2/2`. Figure 3 of TP-2914 shows four correlations
within an approximately ±25% data spread up to 14 km/s, and the derivation
assumes wall temperatures below about 4000 K. The API warns:

* above 14 km/s;
* above a 4000 K wall temperature;
* when :math:`h_\infty/(V_\infty^2/2)=2/((\gamma-1)M^2)` exceeds 0.1, that is,
  below about Mach 7.07 for :math:`\gamma=1.4`.

The 0.1 threshold is a practical reading of ":math:`\ll`". It is close to the
ratio at the lowest speed in figure 3 (about 2 km/s) for air at a typical
freestream temperature of 200 K. Tauber notes that a different
velocity-gradient estimate would change the sphere coefficient from 1.83 to
1.90.

**Tauber swept cylinder:** ``sweep_angle`` must lie between 0 and
:math:`\pi/2`, endpoints included. (The turbulent Poll model is stricter; it
excludes both endpoints.)
Equation (41) assumes equal sphere and cylinder velocity gradients. At zero
sweep its ratio to equation (40) is :math:`1.29/1.83\approx1/\sqrt{2}`, and as
:math:`\Lambda\to\pi/2` the heat flux falls to zero. Transition along the
attachment line and the flat-plate-like heating on highly swept real wing
leading edges (TP-2914 equation 49) are not included. For a fully turbulent
attachment line, select the Poll model below.

.. _aerodynamic-heating-turbulent-attachment-line:

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

A separate warning above :math:`T_0=2000` K is this package's own threshold,
shared with the Fay--Riddell model, rather than a limit of the database. Perfect-gas properties omit dissociation;
a numerical prediction at high Mach number is not validation of
high-enthalpy reacting flow. All quantities use SI units and radians.
A wall hotter than :math:`T_{aw}` returns negative heat flux without clipping
or an additional hot-wall warning. Applicability in that cooling direction
is not established by the literature; the database-range warnings still apply.

The original :ref:`Poll (1981) <ref-poll-1981>` paper has not been compared
directly, and the model has not been compared independently with experimental
data. The separate correlation proposed by Bénard et al. is not implemented,
because its published form could not be reconciled with the figure in the same
paper.

The following example has :math:`T_0=1857.3` K and lies inside both database
ranges, so it emits no applicability warning:

.. doctest::

   >>> from math import pi
   >>> from aerophysics import swept_cylinder_turbulent_heat_flux
   >>> result = swept_cylinder_turbulent_heat_flux(6, 226.5, 1197, 0.1, pi/6, wall_temperature=300)
   >>> round(result.heat_flux)
   450252

API and examples
----------------

:func:`aerophysics.heating.sphere_stagnation_heat_flux` evaluates the sphere
stagnation point. Choose the correlation with
:class:`~aerophysics.heating.StagnationHeatingModel` ``FAY_RIDDELL`` (default)
or ``TAUBER``. :func:`aerophysics.heating.swept_cylinder_heat_flux` evaluates
the laminar swept cylinder with Tauber's equation (41).
:func:`aerophysics.heating.swept_cylinder_turbulent_heat_flux` evaluates the
fully turbulent attachment line with the Poll model; its example is given in
:ref:`aerodynamic-heating-turbulent-attachment-line`.

.. doctest::

   >>> import math
   >>> import warnings
   >>> from aerophysics import sphere_stagnation_heat_flux, swept_cylinder_heat_flux
   >>> from aerophysics.exceptions import ApplicabilityWarning
   >>> from aerophysics.heating import StagnationHeatingModel
   >>> with warnings.catch_warnings():
   ...     warnings.simplefilter("ignore", ApplicabilityWarning)
   ...     result = sphere_stagnation_heat_flux(10, 226.5, 1197, 0.1, wall_temperature=300)
   >>> round(result.heat_flux)
   1986435
   >>> tauber = sphere_stagnation_heat_flux(
   ...     10, 226.5, 1197, 0.1, wall_temperature=300,
   ...     model=StagnationHeatingModel.TAUBER,
   ... )
   >>> round(tauber.heat_flux)
   2020352
   >>> cylinder = swept_cylinder_heat_flux(
   ...     10, 226.5, 1197, 0.1, math.radians(30.0), wall_temperature=300
   ... )
   >>> round(cylinder.heat_flux)
   1174326

Verification
------------

The Tauber equations have been checked against TP-2914. The original
:ref:`Fay and Riddell paper <ref-fay-riddell-1958>` is available only as a
paid publisher version and has not been obtained; the correlation was checked
against two public NASA sources, the TFAWS lecture and West and Brandis. The
tests check analytical relations, agreement between the models, and
regression values. They are not a physical validation against flight or
wind-tunnel data. The verification status of the turbulent Poll model is
described at the end of :ref:`aerodynamic-heating-turbulent-attachment-line`.
