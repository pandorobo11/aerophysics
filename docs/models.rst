Models and conventions
======================

Perfect gas and transport properties
------------------------------------

``PerfectGas`` is calorically perfect: the specific gas constant and heat
capacity ratio are constant. ``AIR`` uses the molar gas constant and sea-level
molecular weight from U.S. Standard Atmosphere 1976. Sutherland viscosity and
the standard's air thermal-conductivity correlation are separate transport
models.

All temperatures must be positive. High-temperature dissociation, chemical
reactions, vibrational excitation, and thermodynamic nonequilibrium are outside
version 0.3.

Standard atmosphere
-------------------

``standard_atmosphere`` accepts geometric altitude from -5,000 m through
86,000 m. It converts to geopotential altitude before applying the hydrostatic
layer equations. An input outside this implemented range raises
``ModelRangeError``; no silent extrapolation is performed.

The returned state includes temperature, pressure, density, speed of sound,
gravity, viscosity, thermal conductivity, and Prandtl number. The model is a
reference atmosphere, not a weather forecast.

Isentropic flow
---------------

The implementation assumes steady adiabatic flow of a calorically perfect gas
with no entropy production. State ratios are total over static, so they are at
least one for non-negative Mach number.

The area-Mach inverse is double-valued for ``A/A* > 1``. Callers must select
``MachBranch.SUBSONIC`` or ``MachBranch.SUPERSONIC``. Mass flux uses total
pressure in pascals and total temperature in kelvin.

Shock and expansion waves
-------------------------

The shock and Prandtl–Meyer models assume steady flow of a calorically perfect
gas. Shock state ratios are downstream over upstream static quantities;
``total_pressure_ratio`` is downstream over upstream total pressure. Expansion
state ratios are likewise downstream over upstream, while total temperature
and total pressure remain constant.

All angles are in radians. Convert degrees explicitly with
``aerophysics.units.degrees_to_radians`` and
``aerophysics.units.radians_to_degrees``.

The theta–beta–Mach inverse requires an explicit ``ShockBranch.WEAK`` or
``ShockBranch.STRONG`` selection. ``oblique_shock`` defaults to the weak branch.
If the requested turning angle exceeds the maximum attached-shock angle,
``NoAttachedShockError`` is raised; the calculation does not silently replace
the requested solution with a detached normal shock.

``mach_from_prandtl_meyer`` accepts angles from zero up to, but not including,
the finite limiting Prandtl–Meyer angle. The numerical inverses use fixed
high-accuracy tolerances. The governing equations and reference values follow
NACA Report 1135.

Flat-plate boundary layers
--------------------------

``flat_plate_boundary_layer`` models one side of a smooth flat plate with a
sharp leading edge, zero pressure gradient, and constant edge conditions.
Distance is measured from the leading edge. Surface roughness, pressure
gradients, separation, suction, blowing, and natural-transition prediction are
outside version 0.3.

The caller explicitly selects ``BoundaryLayerRegime.LAMINAR``,
``TURBULENT``, or ``TRANSITIONAL``. A transitional calculation requires
``transition_reynolds``; no transition Reynolds number is inferred. The
laminar model uses the Blasius relations

.. math::

   \delta_{99}/x = 5/\sqrt{Re_x},\quad
   \delta^*/x = 1.7208/\sqrt{Re_x},\quad
   \theta/x = 0.664/\sqrt{Re_x}.

The local and average laminar skin-friction coefficients are
``0.664/sqrt(Re_x)`` and ``1.328/sqrt(Re_x)``. Turbulent thickness uses
``delta_99/x = 0.37 Re_x**(-1/5)`` with a one-seventh-power profile estimate
for displacement and momentum thickness.

``TurbulentCorrelation.POWER_LAW`` uses local and average coefficients
``0.0592 Re_x**(-1/5)`` and ``0.074 Re_x**(-1/5)``.
``TurbulentCorrelation.SCHLICHTING``, the default, uses

.. math::

   \bar C_f = 0.455/(\log_{10} Re_x)^{2.58}

and obtains the local coefficient by differentiating accumulated drag.
Turbulent use outside ``5e5 <= Re_x <= 1e9`` emits
``ApplicabilityWarning``.

For a specified transition, local quantities switch at the requested
Reynolds number. Average skin friction preserves the laminar drag accumulated
before transition by applying a leading-edge offset to the selected turbulent
correlation.

Compressibility
^^^^^^^^^^^^^^^

Compressibility is opt-in. ``CompressibilityCorrection.ECKERT`` evaluates
the incompressible correlation at reference properties using

.. math::

   T^* = 0.22 T_r + 0.28 T_e + 0.50 T_w.

The effective Reynolds number accounts for ideal-gas density and Sutherland
viscosity ratios. ``VAN_DRIEST_II`` applies Eckert to laminar portions and a
Reynolds-number-based Van Driest II engineering transform to turbulent
portions: ``F_theta = mu_e/mu_w``, the selected incompressible correlation is
evaluated at ``F_theta Re_x``, and skin friction is divided by ``F_c``.

The original Van Driest II momentum-thickness transformation does not uniquely
define ``delta_99``. Version 0.3 estimates compressible thicknesses by applying
the corresponding incompressible thickness relation at the effective Reynolds
number. These thicknesses are engineering estimates, not profile solutions.

If wall temperature is omitted, the wall is adiabatic. Recovery factors are
``sqrt(Pr)`` for laminar flow and ``Pr**(1/3)`` for turbulent flow. Results
expose the recovery and wall temperatures used. The gas remains calorically
perfect and transport properties follow the supplied Sutherland model.

Boundary-layer-immersed protrusions
-----------------------------------

``protrusion_drag`` estimates the direct drag of one isolated protrusion by
integrating the undisturbed local dynamic pressure over its frontal area. The
caller supplies a free-stream drag coefficient. The default undisturbed
turbulent velocity profile is the one-seventh-power approximation; measured
or computed velocity and density profiles can be supplied instead.

The optional compressible approximation uses the Walz temperature relation,
constant static pressure normal to the wall, and ideal-gas density scaling.
This is an effective-dynamic-pressure correction, not a solution of the flow
around the protrusion. Wall interference, horseshoe vortices, roughness-induced
transition, downstream skin-friction changes, multiple-element interference,
and shock/protrusion interactions are excluded. A transonic calculation with a
single supplied drag coefficient emits ``ApplicabilityWarning``.

Flight condition
----------------

``FlightCondition`` combines the standard atmosphere with either a Mach number
or velocity. Dynamic pressure is ``rho V² / 2``. Reynolds number per length is
always returned; the dimensionless Reynolds number is returned only when a
positive characteristic length is supplied.
