Viscous-flow and property pages
===============================

These notes cover pages in the **粘性流れ** (viscous flow) and **熱物性**
(thermophysical properties) groups. The **平板境界層** and
**圧縮性境界層プロファイル** pages follow their theory pages directly; how
they pass results to each other and to **突起抗力** is described in
:doc:`index`.

Aerodynamic heating
-------------------

The **空力加熱** (aerodynamic heating) page takes either the standard
atmosphere at a given altitude or manually entered free-stream conditions,
together with the geometry, wall temperature, and radius. For a sphere,
choose Fay--Riddell, Tauber, or a comparison of both. For a swept cylinder,
the Tauber correlation is used. Sweeps over Mach number, altitude, or the
cylinder sweep angle are available. Changing the geometry selects a model
that suits it and clears the previous result. Heat flux is calculated in W/m²
and can be displayed in kW/m² or MW/m². See :ref:`aerodynamic-heating-model`
for the assumptions and ranges.

Protrusion CSV inputs
---------------------

The **突起抗力** page also accepts measured or externally generated CSV
inputs. Download its templates first. A profile file needs
``wall_distance,velocity,density`` columns, and a projected-shape file needs
``height,width`` columns. Values are read in the currently selected display
units. Wall distance or height must start at zero and increase strictly, and
at least two rows of finite numbers are required.

Thermophysical property pages
-----------------------------

The **粘性係数** (viscosity) page evaluates the Sutherland, Keyes, and
Blottner/Wilke dry-air viscosity models at one temperature or over a
temperature sweep, or compares all three.

The **熱化学** (thermochemistry) page evaluates NASA7 or NASA9 frozen dry air,
or compares them, at one temperature or over a temperature sweep. It shows
``cp``, ``gamma``, sound speed, and the sensible enthalpy measured from a
reference temperature that you choose. Temperatures outside 200--6000 K are
rejected unless extrapolation is explicitly allowed.
