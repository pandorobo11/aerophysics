.. _gui-guide:

Local GUI
=========

The optional GUI exposes the package's atmospheric, compressible-flow,
boundary-layer, and thermophysical calculators in a local browser. It uses the
same calculation APIs and model limits as the Python interface; the GUI is not
a separate physical model.

Installation and launch
-----------------------

Install the ``gui`` extra and start the launcher:

.. code-block:: console

   $ python -m pip install "aerophysics[gui]"
   $ aerophysics-gui

The launcher starts Streamlit and opens its local URL. Stop it with
:kbd:`Ctrl-C` in the launching terminal. Calculations and the bundled manual
are served locally; an internet connection is not required after the wheel and
its dependencies are installed.

Display units and the SI core
-----------------------------

The public Python calculation APIs use SI units and radians. The GUI sidebar
can instead display selected aviation units for length, speed, pressure,
temperature, density, and angle. GUI inputs are converted to SI before calling
the calculation core, and results are converted from SI only for presentation
and CSV export. Changing a sidebar unit preserves the represented physical
input rather than reinterpreting its numeric value.

Model parameters that are intrinsically SI, such as dynamic viscosity in
``Pa s``, remain labelled as such. Always use the unit shown beside an input or
table column.

Case handoff
------------

Related calculators can pass a case through the current GUI session without
rounding through displayed values:

1. Run a single-point **Atmosphere and flight conditions** calculation and
   save the current flight case. The **Flat-plate boundary layer** page can
   select it as the edge-condition source.
2. Run a single-point, fully turbulent flat-plate case and save the current
   boundary-layer case. The **Compressible boundary-layer profile** page can
   use its edge state, thickness, and wall shear stress.
3. Save one generated boundary-layer profile. The **Protrusion drag** page can
   use the saved SI velocity and density profile.

These handoffs live only in the current Streamlit session. Download a settings
JSON file when a calculation must be reproduced after restarting the GUI.

CSV and settings JSON
---------------------

Every completed calculator provides two reproducibility downloads:

* **Result CSV** contains the displayed result table. Its headings include the
  active display units, and the file is UTF-8 with a byte-order mark for
  spreadsheet compatibility.
* **Settings JSON** is a versioned calculation configuration. It stores the
  calculator and model selections, canonical SI inputs, sweep definition when
  present, and the display-unit preferences. Load it from the same calculator
  page; a configuration for a different calculator or schema is rejected.

The protrusion calculator also accepts measured or externally generated CSV
inputs. Download its templates before preparing data. A profile file requires
``wall_distance,velocity,density`` columns, and a projected-shape file requires
``height,width`` columns. Values use the currently selected display units;
wall distance or height must start at zero and increase strictly. At least two
finite data rows are required. The detached-shock calculator can additionally
export the computed shock-shape coordinates as CSV.

Offline documentation
---------------------

The distributed wheel contains the rendered Sphinx manual. On launch,
``aerophysics-gui`` starts a loopback-only documentation server and the
**Documentation** page opens the bundled topics without contacting an external
site.

In a source checkout, build the manual before launching the GUI:

.. code-block:: console

   $ uv run sphinx-build -W --keep-going -b html docs docs/_build/html
   $ aerophysics-gui

The launcher automatically detects a valid ``docs/_build/html`` directory. An
extracted release documentation ZIP can also be read directly by opening its
top-level ``index.html`` in a browser.
