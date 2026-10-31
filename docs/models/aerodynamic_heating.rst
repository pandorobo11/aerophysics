.. _aerodynamic-heating-model:

Aerodynamic heating
===================

These laminar engineering estimates give the convective heat flux at the
stagnation point of a sphere and along the stagnation line of a swept
cylinder. They use SI units. Heat flux :math:`q_w` is in W/m² and is positive
into the wall; a wall hotter than the gas gives a negative value. Enthalpies
are in J/kg, and angles are in radians.

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
:math:`\cos\Lambda`.

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

The coefficients carry the units needed for SI inputs. The ``gas`` argument
only converts the freestream Mach number, temperature, and pressure into
velocity and density; ``prandtl_number`` and ``viscosity_model`` are not used
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

All models are laminar, require a supersonic freestream, and represent only
the stagnation point or stagnation line.

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

**Swept cylinder:** ``sweep_angle`` must lie between 0 and :math:`\pi/2`.
Equation (41) assumes equal sphere and cylinder velocity gradients. At zero
sweep its ratio to equation (40) is :math:`1.29/1.83\approx1/\sqrt{2}`, and as
:math:`\Lambda\to\pi/2` the heat flux falls to zero. Transition along the
attachment line and the flat-plate-like heating on highly swept real wing
leading edges (TP-2914 equation 49) are not included.

API and examples
----------------

:func:`aerophysics.heating.sphere_stagnation_heat_flux` evaluates the sphere
stagnation point. Choose the correlation with
:class:`~aerophysics.heating.StagnationHeatingModel` ``FAY_RIDDELL`` (default)
or ``TAUBER``. :func:`aerophysics.heating.swept_cylinder_heat_flux` evaluates
the swept cylinder with Tauber's equation (41).

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
wind-tunnel data.
