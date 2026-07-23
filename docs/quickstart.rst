Quick start
===========

Atmosphere and flight condition
-------------------------------

Evaluate sea-level standard atmosphere:

>>> from aerophysics import FlightCondition, standard_atmosphere
>>> sea_level = standard_atmosphere(0.0)
>>> sea_level.temperature
288.15
>>> round(sea_level.speed_of_sound, 3)
340.294

Create a Mach-defined flight condition. Altitude and characteristic length are
in metres:

>>> condition = FlightCondition.from_mach(0.0, 0.8, characteristic_length=1.5)
>>> round(condition.dynamic_pressure, 1)
45393.6
>>> round(condition.total_temperature, 4)
325.0332

Isentropic flow
---------------

State ratios use the total-to-static convention:

>>> from aerophysics.isentropic import isentropic_ratios
>>> ratios = isentropic_ratios(2.0)
>>> round(ratios.total_temperature_ratio, 6)
1.8
>>> round(ratios.total_pressure_ratio, 6)
7.824449

Shock and expansion waves
-------------------------

Angles are supplied in radians and branch selection is explicit:

>>> from aerophysics import ShockBranch, oblique_shock, prandtl_meyer_expansion
>>> from aerophysics.units import degrees_to_radians, radians_to_degrees
>>> shock = oblique_shock(
...     2.0, degrees_to_radians(10.0), branch=ShockBranch.WEAK
... )
>>> round(radians_to_degrees(shock.shock_angle), 3)
39.314
>>> expansion = prandtl_meyer_expansion(2.0, degrees_to_radians(10.0))
>>> round(expansion.downstream_mach, 6)
2.384887

Flat-plate boundary layer
-------------------------

Select the boundary-layer state explicitly. This laminar example uses SI
edge conditions and returns the one-sided drag per unit width:

>>> from aerophysics import BoundaryLayerRegime, flat_plate_boundary_layer
>>> layer = flat_plate_boundary_layer(
...     1.0,
...     edge_velocity=10.0,
...     edge_density=1.0,
...     edge_dynamic_viscosity=1e-5,
...     regime=BoundaryLayerRegime.LAMINAR,
... )
>>> round(layer.reynolds_number, 1)
1000000.0
>>> round(layer.boundary_layer_thickness, 6)
0.005
>>> round(layer.drag_per_unit_width, 6)
0.0664

Unit conversion
---------------

Customary units are converted explicitly rather than inferred:

>>> from aerophysics.units import feet_to_meters, knots_to_meters_per_second
>>> feet_to_meters(10_000.0)
3048.0
>>> round(knots_to_meters_per_second(100.0), 6)
51.444444
