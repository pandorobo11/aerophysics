.. _gui-guide:

GUI manual
==========

The optional GUI runs the package's calculators in a local browser. It calls
the same calculation APIs and applies the same model limits as the Python
interface; it is not a separate physical model. The GUI labels are in
Japanese. This manual gives each Japanese page name with an English
translation.


Installation and launch
-----------------------

Install the ``gui`` extra from the current GitHub Release wheel and start the
launcher:

.. code-block:: console

   $ python -m pip install "aerophysics[gui] @ https://github.com/pandorobo11/aerophysics/releases/download/v0.7.0/aerophysics-0.7.0-py3-none-any.whl"
   $ aerophysics-gui

The launcher opens the GUI in your browser at a local address. It accepts
connections only from your own computer. Stop it with :kbd:`Ctrl-C` in the
terminal where you started it. Calculations and this manual are served
locally, so no internet connection is needed once the package is installed.

Pages and related documentation
-------------------------------

The sidebar groups the pages as follows. The last column links to the pages
that explain the equations and assumptions. Each calculator page also has a
**📖 理論・モデルの説明** (theory and model) button that opens the matching
page of this manual directly.

.. list-table:: GUI pages
   :header-rows: 1
   :widths: 18 30 52

   * - Group
     - Page
     - Theory and model details
   * - 飛行状態 (flight)
     - 大気・飛行条件 (atmosphere and flight conditions)
     - :doc:`../models/standard_atmosphere`, :doc:`../models/flight_conditions`,
       :doc:`../guides/atmosphere_flight`
   * - 圧縮性流れ (compressible flow)
     - 連続流れ計算 (sequential flow)
     - :doc:`../guides/flow_sequence`
   * -
     - 等エントロピー流れ (isentropic flow)
     - :doc:`../models/isentropic_flow`
   * -
     - 垂直衝撃波 (normal shock), 斜め衝撃波 (oblique shock)
     - :doc:`../models/normal_oblique_shocks`
   * -
     - 円錐衝撃波 (conical shock)
     - :doc:`../models/conical_shocks`
   * -
     - 離脱衝撃波 (detached shock)
     - :doc:`../models/detached_shocks`
   * -
     - Prandtl–Meyer膨張 (Prandtl--Meyer expansion)
     - :doc:`../models/expansion_waves`
   * - 粘性流れ (viscous flow)
     - 平板境界層 (flat-plate boundary layer)
     - :doc:`../models/flat_plate_boundary_layer`
   * -
     - 空力加熱 (aerodynamic heating)
     - :doc:`../models/aerodynamic_heating`
   * -
     - 圧縮性境界層プロファイル (compressible boundary-layer profile)
     - :doc:`../models/compressible_velocity_transformations`
   * -
     - 突起抗力 (protrusion drag)
     - :doc:`../models/protrusion_drag`
   * - 熱物性 (thermophysical properties)
     - 粘性係数 (viscosity)
     - :doc:`../models/transport_properties`
   * -
     - 熱化学 (thermochemistry)
     - :doc:`../models/gas_models`
   * - ヘルプ (help)
     - ドキュメント (documentation)
     - This manual

Most pages offer a **単点** (single point) mode and a **スイープ** (sweep)
mode that varies one input over a range. Most calculator pages also have a
**モデルの前提・適用範囲** (model assumptions and applicability) panel that
summarizes the limits.

Display units and the SI core
-----------------------------

The Python calculation APIs use SI units and radians. The GUI sidebar can
display other units for length, area, speed, pressure, temperature, density,
force, inverse length, angle, and heat flux (W/m², kW/m², or MW/m²). Choices
include millimetres and inches, feet per second, kilopascals, hectopascals,
pounds-force per square foot, Celsius, Rankine, and pounds mass per cubic
foot. Length and inverse length share one selector, so choosing feet also
shows inverse-length results per foot.

GUI inputs are converted to SI before the calculation, and results are
converted back only for display and CSV export. Changing a sidebar unit keeps
the physical value of each input; for example, 1000 m becomes 3280.84 ft
rather than 1000 ft.

Quantities without a sidebar selector, such as mass flux, drag per unit width,
and dynamic viscosity in ``Pa s``, keep their SI label. Always use the unit
shown next to an input or table column.

Passing a case between pages
----------------------------

Related pages can pass a result to the next page in the same browser session,
at full precision rather than the rounded displayed values:

1. Run a single-point calculation on **大気・飛行条件** and save the flight
   case. The **平板境界層** (flat-plate boundary layer) page can then use it
   as the edge condition.
2. Run a single-point, fully turbulent flat-plate case and save the
   boundary-layer case. The **圧縮性境界層プロファイル** page can then use its
   edge state, thickness, and wall shear stress.
3. Save one generated boundary-layer profile. The **突起抗力** (protrusion
   drag) page can then use its velocity and density profile.

Saved cases last only until the GUI is closed. Download a settings JSON file
to reproduce a calculation later.

CSV and settings JSON
---------------------

Every completed calculation offers two downloads:

* **Result CSV** contains the displayed result table. The column headings
  include the display units. The file is UTF-8 with a byte-order mark so that
  spreadsheet programs read it correctly.
* **Settings JSON** stores the page's model choices, inputs (in SI), sweep
  range, and display units. Load it on the same page to repeat the
  calculation. A file from a different page or an unsupported format version
  is rejected, and the current inputs are left unchanged.

Loading a settings file checks its structure and number ranges, but not
whether the physics will succeed; a calculation can still fail with the usual
errors after loading.

Offline documentation
---------------------

The installed package includes this manual. When ``aerophysics-gui`` starts,
it also serves the manual locally, and the **ドキュメント** page shows it
without contacting an external site.

In a source checkout, build the manual before launching the GUI:

.. code-block:: console

   $ uv run sphinx-build -W --keep-going -b html docs docs/_build/html
   $ aerophysics-gui

The launcher finds a valid ``docs/_build/html`` directory automatically. An
extracted release documentation ZIP can also be read directly by opening its
top-level ``index.html`` in a browser.

Page notes
----------

.. toctree::
   :maxdepth: 1

   compressible_flow_pages
   viscous_and_property_pages

The **大気・飛行条件** page and the **連続流れ計算** page are described in
:doc:`../guides/atmosphere_flight` and :doc:`../guides/flow_sequence`.
