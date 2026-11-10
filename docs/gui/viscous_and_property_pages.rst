Viscous-flow and property pages
===============================

These notes cover pages in the **粘性流れ** (viscous flow) and **熱物性**
(thermophysical properties) groups. How the boundary-layer pages pass results
to each other and to **突起抗力** (protrusion drag) is described in
:doc:`index`.

Flat-plate boundary layer
-------------------------

The **平板境界層** page has a **圧縮性を考慮** (include compressibility)
checkbox and, for turbulent or transitional flow, a **乱流モデル** (turbulent
model) selector. The selector lists only the models that match the checkbox:

.. list-table::
   :header-rows: 1
   :widths: 22 44 34

   * - Compressibility
     - Model
     - Thickness results
   * - Off
     - **Schlichting（1/5乗則の厚さ）** (default)
     - Simple one-fifth-power formula
   * - Off
     - **Schlichting（対数則＋ウェイク則の厚さ）**
     - More accurate at high Reynolds number; same friction and drag
   * - On
     - **Van Driest II＋壁温補正（Schlichting）** (default)
     - Calculated from the velocity profile
   * - On
     - **Eckert基準温度法（Schlichting）**, **Van Driest II（Willems の式）**
     - Not available

Changing the checkbox resets a model that no longer fits to the default.

When a model gives no thickness, the page shows **—** and leaves the thickness
columns empty in the table, plots, and CSV file. Friction and drag are still
shown. Such a result cannot be saved as the boundary-layer case for the profile
page, because that page needs the thickness; **現在の境界層ケースとして保存**
shows a warning instead.

For laminar flow with compressibility switched on, leaving the wall temperature
empty means an insulated wall. A warning appears when the boundary layer near
the leading edge is thick enough to change the pressure on the plate, which
the flat-plate model does not include.

Settings files saved by older versions of the GUI are converted when loaded. A
few old model combinations no longer exist; the page then asks you to select a
model again.

The models are explained in :doc:`../models/flat_plate_boundary_layer`, and
Python examples are in :doc:`../guides/boundary_layers`.

Flat-plate heat transfer and laminar profiles
---------------------------------------------

With **圧縮性を考慮** selected, the **平板境界層** page also reports the Stanton
number, the heat-transfer coefficient, and the wall heat flux, with heat flux
into the wall positive. The page uses the default Reynolds-analogy factors
described in :doc:`../models/flat_plate_heat_transfer`; they cannot be changed
in the GUI.

The **層流境界層プロファイル** page gives the velocity, temperature, density,
viscosity, and local Mach number through :math:`\delta_{99}` at one station, with
the same wall quantities. It takes only the edge velocity, density, and
temperature: the Mach number and edge viscosity follow from them, using air and
the Sutherland law. The **平板境界層** page takes those two values as separate
inputs, so it matches this page only if they agree. The profile page assumes a
laminar layer at the given position and does not test for transition. It has no
sweep mode.

Aerodynamic heating
-------------------

The **空力加熱** (aerodynamic heating) page takes either the standard
atmosphere at a given altitude or manually entered free-stream conditions,
together with the geometry, wall temperature, and radius. For a sphere,
choose Fay--Riddell, Tauber, or a comparison of both; these are laminar. For
a swept cylinder, choose **Tauber（層流）** (laminar), **Poll（完全乱流・付着線）**
(fully turbulent attachment line), or **Tauber 層流 / Poll 乱流 比較**, which
shows both and their ratio. Sweeps over Mach number, altitude, or the cylinder
sweep angle are available. Changing the geometry selects a model that suits it
and clears the previous result. Heat flux is calculated in W/m² and can be
displayed in kW/m² or MW/m².

A fully turbulent boundary layer is an assumption you make by choosing Poll.
The page does not test for transition and never switches between laminar and
turbulent models on its own. The Poll model needs a sweep angle strictly
between 0 and 90° and :math:`M_\infty\cos\Lambda>1`, so that the flow normal
to the cylinder is supersonic. Points outside this range have no turbulent
value; the table shows them as missing and the page reports how many were
skipped. In a comparison, the laminar Tauber values at those points are still
shown. Tauber itself accepts sweep angles of 0 and 90°.

See :ref:`aerodynamic-heating-model` for the assumptions and ranges, and
:ref:`aerodynamic-heating-turbulent-attachment-line` for the Poll equations,
the Stanton-number definition, and its validity limits.

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
