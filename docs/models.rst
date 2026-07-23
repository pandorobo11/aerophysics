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
version 0.2.

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

Flight condition
----------------

``FlightCondition`` combines the standard atmosphere with either a Mach number
or velocity. Dynamic pressure is ``rho V² / 2``. Reynolds number per length is
always returned; the dimensionless Reynolds number is returned only when a
positive characteristic length is supplied.
