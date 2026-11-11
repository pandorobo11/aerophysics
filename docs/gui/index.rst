.. _gui-guide:

GUI manual
==========

The optional GUI lets you run the package's calculations in a web browser on
your own computer, without writing Python. It uses the same models and the
same limits as the Python functions. The labels in the GUI are in Japanese;
this manual gives each Japanese name with an English translation.


Installation and launch
-----------------------

Install the ``gui`` extra from the current GitHub Release wheel and start the
launcher:

.. code-block:: console

   $ python -m pip install "aerophysics[gui] @ https://github.com/pandorobo11/aerophysics/releases/download/v0.8.0/aerophysics-0.8.0-py3-none-any.whl"
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
     - :doc:`../models/standard_atmosphere`, :doc:`../models/flight_conditions`
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
     - 軸対称物体流れ（特性曲線法） (axisymmetric body flow)
     - :doc:`../models/axisymmetric_characteristics`
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
     - 境界層マーチング解析 (marching boundary layer)
     - :doc:`../models/marching_boundary_layer`
   * -
     - 空力加熱 (aerodynamic heating)
     - :doc:`../models/aerodynamic_heating`
   * -
     - 圧縮性境界層プロファイル (compressible boundary-layer profile)
     - :doc:`../models/compressible_velocity_transformations`
   * -
     - 層流境界層プロファイル (laminar boundary-layer profile)
     - :doc:`../models/flat_plate_boundary_layer`
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
   case. The **平板境界層** (flat-plate boundary layer) and
   **層流境界層プロファイル** (laminar boundary-layer profile) pages can then
   use it as the edge condition.
2. Run a single-point, fully turbulent flat-plate case and save the
   boundary-layer case. The **圧縮性境界層プロファイル** page can then use its
   edge state, thickness, and wall shear stress.
3. Save one generated boundary-layer profile. The **突起抗力** (protrusion
   drag) page can then use its velocity and density profile.

4. Run **軸対称物体流れ（特性曲線法）** and save the body-flow case. The
   **境界層マーチング解析** (marching boundary layer) page can then use its
   surface pressure, radius, and post-shock total conditions as the edge.

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
  calculation. A loaded file restores its saved inputs as manual input, even if
  they came from a saved case and that page has a different case saved now;
  select the case again to use it. A file from a different page or an
  unsupported format version is rejected, and the current inputs are left
  unchanged.

Loading a settings file checks that the file is well formed and that its
numbers are within the limits of the input boxes. It does not check whether
the calculation will succeed, so a loaded case can still report an error when
you run it.

Offline documentation
---------------------

The installed package includes this manual. When ``aerophysics-gui`` starts,
it also serves the manual locally, and the **ドキュメント** page shows it
without contacting an external site.

The manual can also be read without the GUI: download the documentation ZIP
from the release page, extract it, and open its top-level ``index.html`` in a
browser. If you run the GUI from a source checkout, build the manual first as
described in the repository's ``DEVELOPMENT.md``.

Page notes
----------

.. toctree::
   :maxdepth: 1

   flight_and_sequence_pages
   compressible_flow_pages
   viscous_and_property_pages
