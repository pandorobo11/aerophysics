Aerodynamic heating
===================

At supersonic and hypersonic speed, the air that is brought to rest at the
front of a body becomes very hot and heats the surface. The heating is
strongest at the **stagnation point** of a blunt nose and along the leading
edge of a wing or fin. This guide shows how to estimate the heat flux there.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Shape
     - Function
   * - Nose of a sphere
     - :func:`~aerophysics.heating.sphere_stagnation_heat_flux`
   * - Leading edge of a swept cylinder, laminar
     - :func:`~aerophysics.heating.swept_cylinder_heat_flux`
   * - Leading edge of a swept cylinder, turbulent
     - :func:`~aerophysics.heating.swept_cylinder_turbulent_heat_flux`

All three take the free-stream Mach number, static temperature in K, and static
pressure in Pa, followed by the geometry and the wall temperature in K. The
free stream must be supersonic. Heat flux is returned in W/m². For heating on a
flat surface away from the nose, see :doc:`boundary_layers`.

The equations and their ranges are in :doc:`../models/aerodynamic_heating`.

Nose of a sphere
----------------

This example is a nose of 0.1 m radius at Mach 6 and 30 km altitude, with the
wall at 300 K:

>>> from aerophysics import sphere_stagnation_heat_flux, standard_atmosphere
>>> air = standard_atmosphere(30_000.0)
>>> heating = sphere_stagnation_heat_flux(
...     6.0,
...     air.temperature,
...     air.pressure,
...     0.1,
...     wall_temperature=300.0,
... )
>>> round(heating.heat_flux / 1e3, 1)
391.8
>>> round(heating.stagnation_temperature, 1)
1857.4

The heat flux is about 392 kW/m². The result also holds the stagnation pressure
and enthalpy and the free-stream speed and density.

A blunter nose is heated less. The heat flux is proportional to one over the
square root of the nose radius, so a radius four times larger halves it:

>>> blunt = sphere_stagnation_heat_flux(
...     6.0, air.temperature, air.pressure, 0.4, wall_temperature=300.0
... )
>>> round(blunt.heat_flux / heating.heat_flux, 2)
0.5

Two methods are available through the ``model`` argument:

* ``StagnationHeatingModel.FAY_RIDDELL`` (the default) calculates the heat flux
  from the properties of the air at the stagnation point.
* ``StagnationHeatingModel.TAUBER`` is a simpler formula fitted for hypersonic
  flight. It is intended for speeds above about Mach 7 and warns below that.

Both assume a laminar boundary layer.

Leading edge of a swept cylinder
--------------------------------

A wing or fin leading edge can be treated as a cylinder swept back from the
flow. The sweep angle is in radians, with zero meaning the cylinder lies square
across the flow:

>>> from math import radians
>>> from aerophysics import swept_cylinder_heat_flux
>>> edge = swept_cylinder_heat_flux(
...     6.0,
...     air.temperature,
...     air.pressure,
...     0.05,
...     radians(45.0),
...     wall_temperature=300.0,
... )

Sweeping the leading edge reduces the laminar heating:

>>> unswept = swept_cylinder_heat_flux(
...     6.0, air.temperature, air.pressure, 0.05, 0.0, wall_temperature=300.0
... )
>>> round(edge.heat_flux / unswept.heat_flux, 2)
0.63

This function uses the Tauber formula, so at Mach 6 both calls warn that the
speed is below the range the formula is meant for. See :doc:`vectorization_errors`
for how to handle warnings.

Turbulent leading edge
~~~~~~~~~~~~~~~~~~~~~~

On a swept leading edge, air flows along the edge, and this flow can become
turbulent. Turbulent heating is higher than laminar heating:

>>> from aerophysics import swept_cylinder_turbulent_heat_flux
>>> turbulent = swept_cylinder_turbulent_heat_flux(
...     6.0,
...     air.temperature,
...     air.pressure,
...     0.05,
...     radians(45.0),
...     wall_temperature=300.0,
... )
>>> round(turbulent.heat_flux / edge.heat_flux, 2)
1.89

The package does not decide whether the leading edge is laminar or turbulent.
You choose by calling one function or the other.

The turbulent function needs a sweep angle strictly between 0° and 90°, and the
part of the flow that meets the cylinder head-on must be supersonic:
:math:`M_\infty\cos\Lambda > 1`. Otherwise it raises ``ValueError``.

Limits to keep in mind
----------------------

* These are estimates for the stagnation point or leading edge only. They do
  not give the heating over the rest of the body.
* Air is treated as a gas with constant heat capacity. Above a stagnation
  temperature of about 2000 K real air starts to break up into atoms, which
  these methods do not model. The functions warn when this happens.
* Radiation from the hot gas, which matters at very high entry speeds, is not
  included.
