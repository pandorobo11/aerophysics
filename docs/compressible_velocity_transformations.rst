.. _compressible-velocity-transformations:

Compressible turbulent mean-velocity transformations
=====================================================

Purpose and scope
-----------------

The transformations on this page map a known compressible turbulent mean
velocity profile to coordinates in which it can be compared with an
incompressible law of the wall.  They are useful for comparing experiments and
simulations, assessing wall models, and constructing reduced-order models.

They are not, by themselves, closed models that predict a dimensional velocity
profile from edge Mach number and wall temperature.  A forward transformation
requires profiles of velocity and thermodynamic properties.  An inverse
prediction also requires an incompressible wall law or wall model, a
temperature--velocity closure, property models, outer matching data, and a
method for determining wall shear stress.

The discussion is limited to the inner layer of turbulent wall-bounded flows.
It does not provide a laminar or transition model, an outer-layer wake model,
or a universal treatment of roughness, three-dimensionality, separation, or
thermochemical nonequilibrium.

Common notation
---------------

Let :math:`\widetilde U` denote the Favre-averaged streamwise velocity,
:math:`y` the distance from the wall, and :math:`\tau_w`, :math:`\rho_w`, and
:math:`\mu_w` the wall shear stress, density, and dynamic viscosity.  Wall
units are

.. math::

   u_\tau = \sqrt{\frac{\tau_w}{\rho_w}},
   \qquad
   y^+ = \frac{\rho_w u_\tau y}{\mu_w},
   \qquad
   \widetilde U^+ = \frac{\widetilde U}{u_\tau}.

Define the local mean-property ratios

.. math::

   r_\rho = \frac{\overline\rho}{\rho_w},
   \qquad
   r_\mu = \frac{\overline\mu}{\mu_w},

and the semi-local wall coordinate

.. math::

   y^*
   = y^+\frac{\sqrt{r_\rho}}{r_\mu}
   = \frac{y\sqrt{\overline\rho\,\tau_w}}{\overline\mu}.

The overbar denotes a Reynolds average.  In a mean-profile data set,
:math:`\overline\mu` is normally evaluated from the mean thermodynamic state.
The edge Mach number is written :math:`M_e`; the symbol :math:`M` is not used
for the viscosity ratio on this page.

Three different uses of the Van Driest name
--------------------------------------------

The following models are distinct and should not be interchanged:

* The **Van Driest velocity transformation** from 1951 is the density-weighted
  mean-velocity mapping described below.
* The **Van Driest II transformation** maps compressible skin-friction and
  momentum-thickness correlations.  The
  :class:`~aerophysics.boundary_layer.CompressibilityCorrection`\ ``.VAN_DRIEST_II``
  option implements an engineering form of that skin-friction correction; it
  does not return a transformed mean velocity profile.
* The **Van Driest damping function** from 1956 modifies a near-wall
  mixing-length or eddy-viscosity model.  It is neither of the two
  transformations above.

Van Driest transformation
-------------------------

Van Driest assumed an approximately constant total-stress layer, negligible
viscous stress in the logarithmic region, and a Prandtl mixing length
:math:`\ell_m=\kappa y`.  The transformed velocity increment is

.. math::

   \mathrm{d}U_{VD}^+
   = \sqrt{r_\rho}\,\mathrm{d}\widetilde U^+,
   \qquad
   U_{VD}^+(y)
   = \int_0^{\widetilde U^+(y)}
     \sqrt{r_\rho}\,\mathrm{d}\widetilde U^+.

The wall coordinate remains :math:`y^+`.  When the assumptions hold,
:math:`U_{VD}^+` follows the incompressible logarithmic law,

.. math::

   U_{VD}^+ \simeq \frac{1}{\kappa}\ln y^+ + B.

The method requires :math:`\widetilde U`, :math:`\overline\rho`, and
:math:`\tau_w`; it does not use the viscosity profile explicitly.  Its
simplicity and long experimental history make it a useful baseline for
adiabatic or nearly adiabatic, approximately zero-pressure-gradient boundary
layers.  Its density-only scaling is generally inaccurate for strong wall
heating or cooling and is not a reliable general transformation for internal
flows.  See :ref:`Van Driest (1951) <ref-van-driest-1951>`.

Trettel--Larsson transformation
-------------------------------

Trettel and Larsson combined semi-local scaling with near-wall momentum
conservation.  The transformed wall distance is :math:`y^*`, and the velocity
increment is

.. math::

   \mathrm{d}U_{TL}^+
   = r_\mu\frac{\mathrm{d}y^*}{\mathrm{d}y^+}
     \,\mathrm{d}\widetilde U^+.

Expanding the coordinate derivative gives

.. math::

   \mathrm{d}U_{TL}^+
   =
   \sqrt{r_\rho}
   \left[
     1
     + \frac{y^+}{2r_\rho}
       \frac{\mathrm{d}r_\rho}{\mathrm{d}y^+}
     - \frac{y^+}{r_\mu}
       \frac{\mathrm{d}r_\mu}{\mathrm{d}y^+}
   \right]
   \mathrm{d}\widetilde U^+.

The transformation therefore needs velocity, density, and viscosity profiles,
as well as :math:`\tau_w`.  It recovers the asymptotically correct viscous
sublayer behavior and has performed particularly well for compressible
channel flows with substantial property variation.  Subsequent assessments
show that, in boundary layers, small sustained errors in the transformed
shear can accumulate into a shifted logarithmic-law intercept, especially
with wall heat transfer.  See
:ref:`Trettel and Larsson (2016) <ref-trettel-larsson-2016>`.

Volpiani data-driven transformation
-----------------------------------

Volpiani, Iyer, Pirozzoli, and Larsson considered coordinate and velocity
mappings with power-law dependence on :math:`r_\rho` and :math:`r_\mu`.
Near-wall stress similarity constrains the relationship between the powers;
the remaining exponents were calibrated against direct numerical simulation
(DNS) data.  The resulting mapping is

.. math::

   \mathrm{d}y_V^+
   = r_\rho^{1/2}r_\mu^{-3/2}\,\mathrm{d}y^+,
   \qquad
   \mathrm{d}U_V^+
   = r_\rho^{1/2}r_\mu^{-1/2}\,
     \mathrm{d}\widetilde U^+.

Unlike the Trettel--Larsson transformation, it does not require numerical
derivatives of the property profiles.  It was calibrated on four boundary
layers and tested on heated, cooled, and adiabatic boundary layers from
:math:`M_e=2.3` through :math:`13.64`, plus one channel-flow case.  The
original study reported good inner-layer collapse under strong wall cooling
and heating, where earlier transformations had difficulty.

The exponents are empirical DNS-calibration results rather than a unique
consequence of the governing equations.  The original study's largest
conventional friction Reynolds number was approximately 650; it called for
additional high-Reynolds-number evaluation.  Later comparisons found degraded
performance in some channel, pipe, and non-canonical flows.  See
:ref:`Volpiani et al. (2020) <ref-volpiani-2020>`.

Griffin--Fu--Moin total-stress transformation
---------------------------------------------

Griffin, Fu, and Moin (GFM) argued that distinct physics should be used for
the viscous and turbulent parts of the inner layer.  Define the
quasi-equilibrium and Trettel--Larsson-scaled mean shears as

.. math::

   S_{eq}^+
   = \frac{1}{r_\mu}
     \frac{\mathrm{d}\widetilde U^+}{\mathrm{d}y^*},
   \qquad
   S_{TL}^+
   = r_\mu
     \frac{\mathrm{d}\widetilde U^+}{\mathrm{d}y^+}.

Here :math:`S_{eq}^+` represents the logarithmic-region scaling motivated by
quasi-equilibrium between turbulence production and dissipation, whereas
:math:`S_{TL}^+` supplies the viscous-sublayer scaling.  If
:math:`\tau^+=\tau/\tau_w` is the local total shear stress, the composite
transformed shear is

.. math::

   S_t^+
   = \frac{\tau^+ S_{eq}^+}
          {\tau^+ + S_{eq}^+ - S_{TL}^+},
   \qquad
   U_t^+(y^*) = \int_0^{y^*} S_t^+\,\mathrm{d}\eta.

When a total-stress profile is unavailable, the constant-stress-layer
approximation :math:`\tau^+\simeq1` gives

.. math::

   S_t^+
   \simeq
   \frac{S_{eq}^+}{1 + S_{eq}^+ - S_{TL}^+}.

The exact form requires velocity, density, viscosity, wall shear, velocity
gradients, and the total-stress profile.  The simplified form omits the last
input but retains the constant-stress approximation.  The original study
reported collapse across the viscous, buffer, and logarithmic regions for
adiabatic, heated, and cooled boundary layers, channel and pipe flows, and a
boundary layer downstream of shock impingement.  Its test set covered edge
Mach numbers from 0 through 15 and friction Reynolds numbers from about 200
through 2000.  See
:ref:`Griffin, Fu, and Moin (2021) <ref-griffin-fu-moin-2021>`.

GFM has the broadest validation set among the four methods for the canonical
flows considered in its original paper, but this does not make it universally
applicable.  Later work found mixed results for high-enthalpy and supercritical
flows and noted limitations of semi-local wall-coordinate mappings.  Strong
pressure gradients and other non-equilibrium conditions can also invalidate
the stress and quasi-equilibrium assumptions.

Comparison and selection
------------------------

.. list-table::
   :header-rows: 1
   :widths: 11 17 17 18 18 19

   * - Method
     - Recommended starting point
     - Principal forward inputs
     - Validation evidence
     - Advantages
     - Limitations
   * - Van Driest
     - Adiabatic or nearly adiabatic, approximately ZPG boundary layers
     - :math:`\widetilde U`, :math:`\overline\rho`, :math:`\tau_w`
     - Long experimental and DNS history for adiabatic boundary layers
     - Simplest method; no viscosity profile or property derivatives
     - Density-only scaling; unreliable with strong wall heat transfer
   * - Trettel--Larsson
     - Variable-property channel flows and viscous-sublayer scaling
     - :math:`\widetilde U`, :math:`\overline\rho`,
       :math:`\overline\mu`, :math:`\tau_w`
     - Supersonic channels and boundary layers with several thermal conditions
     - Physics-based semi-local coordinate; correct near-wall asymptote
     - Property derivatives amplify noise; boundary-layer log intercept can
       shift
   * - Volpiani et al.
     - Heated or cooled ZPG boundary layers when full property profiles exist
     - :math:`\widetilde U`, :math:`\overline\rho`,
       :math:`\overline\mu`, :math:`\tau_w`
     - Original DNS set: :math:`M_e=2.3`--13.64 and varied wall temperature
     - No profile derivatives; good original high-Mach thermal-wall results
     - DNS-fitted powers; limited original :math:`Re_\tau`; weaker in some
       internal and non-canonical flows
   * - GFM
     - Canonical flows spanning adiabatic and diabatic walls when gradient or
       stress information is available
     - Other methods' inputs plus velocity gradients; exact form also needs
       :math:`\tau(y)`
     - Original set: boundary layers, channels, pipes, and post-shock flow at
       :math:`M_e=0`--15
     - Treats viscous and turbulent stresses separately; broad canonical-flow
       performance
     - Most complex; relies on stress and quasi-equilibrium assumptions;
       non-canonical performance is not uniform

Here ZPG means zero pressure gradient.  A practical first choice is Van Driest
for an adiabatic ZPG baseline, Trettel--Larsson for an internal flow or when
near-wall semi-local scaling is the main interest, Volpiani for a thermally
loaded ZPG boundary layer when derivative-free processing is valuable, and
GFM when the broader canonical-flow behavior justifies its additional inputs.
When the required inputs are available, applying more than one transformation
is preferable to treating this selection as a universal ranking.

Forward numerical workflow
--------------------------

A reproducible forward comparison can use the following procedure:

1. Assemble one-dimensional profiles of :math:`y`, :math:`\widetilde U`,
   :math:`\overline\rho`, and :math:`\overline\mu`, together with
   :math:`\tau_w`.  Require strictly increasing :math:`y`, positive density
   and viscosity, and a wall point or a justified near-wall extrapolation.
2. Compute :math:`u_\tau`, :math:`y^+`, :math:`\widetilde U^+`,
   :math:`r_\rho`, :math:`r_\mu`, and :math:`y^*` on one common grid.  At a
   resolved no-slip wall, :math:`\widetilde U^+=0` and
   :math:`r_\rho=r_\mu=1`.
3. Evaluate each differential mapping and use cumulative trapezoidal
   integration from the wall.  For Trettel--Larsson, evaluating
   :math:`r_\mu\,\mathrm{d}y^*/\mathrm{d}y^+` directly is algebraically
   equivalent to differentiating density and viscosity separately and is
   often easier to implement consistently.
4. Use one-sided derivatives at the wall and a derivative method appropriate
   for the grid spacing.  Experimental profiles may require uncertainty-aware
   smoothing before applying Trettel--Larsson or GFM; smoothing must not alter
   the wall value or create non-monotone coordinates.
5. Plot the transformed velocity against its corresponding transformed wall
   coordinate and compare with an incompressible reference at a comparable
   Reynolds number.  Check both :math:`U^+\simeq y^+` in the viscous sublayer
   and the logarithmic slope and intercept.  Do not interpret wake-region
   scatter as an inner transformation error without an outer-layer analysis.

The integrations are path integrals along the sampled profile: factors such
as :math:`r_\rho` and :math:`r_\mu` must be evaluated at the same local
position as each velocity or wall-distance increment.  Substituting one edge
or wall value throughout generally changes the transformation.

Inverse use for profile prediction
----------------------------------

An inverse calculation starts from an incompressible law of the wall or an
incompressible wall-model equation and maps its strain rate back to the
compressible variables.  The following loop is required in a practical model:

#. Guess or solve for :math:`\tau_w` and, when needed, wall heat flux.
#. Use an algebraic temperature--velocity relation or an energy equation to
   obtain :math:`T(y)` from the current velocity estimate.
#. Evaluate :math:`\overline\rho(y)` and :math:`\overline\mu(y)` from an
   equation of state and transport model.
#. Integrate the selected inverse velocity mapping and match the result to a
   supplied inner-layer state and boundary-layer-edge state.
#. Iterate the coupled velocity, temperature, properties, and wall quantities
   to convergence.

Consequently, the inverse problem cannot be reproduced from :math:`M_e`,
:math:`T_e`, and :math:`T_w` alone.  Griffin, Fu, and Moin coupled the inverse
GFM transformation to a generalized Reynolds-analogy temperature relation and
an incompressible wall model; their tested implementation predicts velocity,
temperature, wall shear, and heat flux from matching data.  See
:ref:`Griffin, Fu, and Moin (2023) <ref-griffin-fu-moin-2023>`.

Applicability and evidence
--------------------------

The numerical ranges quoted above are validation ranges, not hard mathematical
bounds.  Extrapolation needs independent evidence.  In particular:

* The original papers primarily address smooth, statistically
  two-dimensional, fully turbulent wall flows and inner-layer similarity.
* Bai, Griffin, and Fu found that no one transformation delivered uniform
  logarithmic-region performance across high-enthalpy, supercritical, and
  pressure-gradient databases.  All methods performed similarly for the weak,
  adiabatic pressure-gradient boundary layers they considered, whereas none
  was satisfactory for the supercritical boundary layers.  See
  :ref:`Bai, Griffin, and Fu (2022) <ref-bai-griffin-fu-2022>`.
* Danis and Durbin showed that the accuracy of a transformation below the
  logarithmic layer is closely related to eddy-viscosity equivalence.  Their
  assessment identified shortcomings in semi-local wall-coordinate mappings
  for strongly cooled hypersonic ZPG data, including mappings used by
  Trettel--Larsson and GFM.  See
  :ref:`Danis and Durbin (2024) <ref-danis-durbin-2024>`.
* Outer-layer wake behavior, rough walls, three-dimensional boundary layers,
  strong adverse pressure gradients, separation, shock/boundary-layer
  interaction away from the tested cases, high-enthalpy chemistry,
  supercritical properties, and strong thermodynamic nonequilibrium require
  separate validation.

The original papers' favorable results and later assessments' limitations are
therefore complementary evidence: the former establish useful operating
regimes, while the latter show why those regimes should not be generalized
without qualification.
