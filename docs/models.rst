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
version 0.1.

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

Flight condition
----------------

``FlightCondition`` combines the standard atmosphere with either a Mach number
or velocity. Dynamic pressure is ``rho V² / 2``. Reynolds number per length is
always returned; the dimensionless Reynolds number is returned only when a
positive characteristic length is supplied.

