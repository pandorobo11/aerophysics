.. _aerodynamic-heating-model:

Aerodynamic heating
===================

These laminar engineering estimates use SI units. Heat flux :math:`q_w` is
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
