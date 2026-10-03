.. _gui-guide:

Local GUI
=========

The optional GUI exposes the package's atmospheric, compressible-flow,
boundary-layer, and thermophysical calculators in a local browser. It uses the
same calculation APIs and model limits as the Python interface; the GUI is not
a separate physical model.

Installation and launch
-----------------------

Install the ``gui`` extra from the current GitHub Release wheel and start the
launcher:

.. code-block:: console

   $ python -m pip install "aerophysics[gui] @ https://github.com/pandorobo11/aerophysics/releases/download/v0.6.0/aerophysics-0.6.0-py3-none-any.whl"
   $ aerophysics-gui

The launcher binds Streamlit to the IPv4 loopback address, opens its local URL,
and does not expose the application to other hosts by default. Stop it with
:kbd:`Ctrl-C` in the launching terminal. Calculations and the bundled manual
are served locally; an internet connection is not required after the wheel and
its dependencies are installed.

Display units and the SI core
-----------------------------

The public Python calculation APIs use SI units and radians. The GUI sidebar
can instead display selected aviation units for length, area, speed, pressure,
temperature, density, force, inverse length, and angle. Length and inverse
length share one selector, so choosing feet also displays inverse-length
results per foot. Choices include
millimetres and inches, feet per second, kilopascals, hectopascals, pounds-force
per square foot, Celsius, Rankine, and pounds mass per cubic foot. GUI inputs
are converted to SI before calling the calculation core, and results are
converted from SI only for presentation and CSV export. Changing a sidebar
unit preserves the represented physical input rather than reinterpreting its
numeric value.

Quantities without a sidebar selector, such as mass flux, drag per unit width,
and dynamic viscosity in ``Pa s``, remain labelled as such. Always use the unit
shown beside an input or table column.

Flow schematics
---------------

Single-point normal-shock results show a red shock perpendicular to the
left-to-right flow, with solved upstream and downstream Mach numbers.
Velocity directions are collinear; arrow lengths do not encode velocity
magnitude. At upstream Mach 1, the caption identifies the zero-strength limit.

Single-point Prandtl--Meyer results show a grey convex wall turning clockwise
through the positive expansion angle ``theta``, with flow above it. A green
fan spans directions ``mu_1`` and ``mu_2 - theta``, where each Mach angle
``mu = asin(1/M)`` is measured from its local flow direction. The downstream
arrow follows the wall at ``-theta``; large supported turns retain that actual
direction. The fan's interior rays are representative Mach lines. At zero
turn, only one Mach line is shown and labelled as no expansion. Angle arcs
use the selected degree/radian unit; the Prandtl--Meyer function ``nu`` is
not drawn as a geometric Mach angle. Sweep results retain the trend plots.

Single-point oblique and conical results include an equally scaled section:
the grey filled body is bounded by its wall, and the red ray is the attached
shock. The dashed horizontal line is the upstream direction and, for a cone,
its axis. Angle arcs measure both the wall angle (``theta`` or cone half-angle
``theta_c``) and shock angle ``beta`` from that horizontal reference, using
the selected display angle unit.

Blue arrows show flow directions, rather than Mach-label leader lines. The
upstream arrow is horizontal; the oblique downstream arrow is parallel to the
wedge. Labels include the solved Mach numbers to three significant figures,
matching the normal-shock and expansion schematics. The conical downstream
arrow denotes only the surface direction and
surface Mach ``M_s``: the intervening Taylor--Maccoll flow turns continuously
from the shock toward the surface. Arrow lengths are schematic and do not
encode velocity magnitude. The cone view is a meridional section with its
upper shock shown.

The detached-shock view uses the same grey body, red shock and blue freestream
arrow. Its origin ``O`` is the nose-curvature centre, with positive ``x``
pointing upstream and the flow arrow toward negative ``x``. The horizontal
axis is reversed so upstream appears on the left and the flow runs left to
right, matching the attached-shock views. Plot coordinates, hover values and
exported CSV coordinates retain the upstream-positive convention. The radius
dimension ``Rn`` runs from that centre to the circular nose. The standoff
bracket ``Delta`` spans from the body vertex at ``x = Rn`` to the shock vertex
at ``x = Rn + Delta``; both lengths use the selected display unit. The shock
curve is Billig's shape with Ambrosio--Wortman standoff, including when the
table also compares a Seiff estimate. See :doc:`../models/shock_waves` for
the correlation assumptions and coordinate conventions.

Thermally perfect shocks
------------------------

The **Oblique shock** page defaults to constant-gamma ``AIR``. Select
``NASA7``, ``NASA9``, or ``HARMONIC_OSCILLATOR`` and enter the upstream
**static** temperature to use temperature-dependent heat capacity. Both
branches and Mach/deflection sweeps are supported. Results include the gas
selection and upstream/downstream static temperatures in the selected units;
settings JSON preserves the gas model and static temperature in SI.
Older settings without those fields retain the ``AIR`` model.

An unavailable thermodynamic state is marked ``out_of_range`` in a sweep;
physical detachment is marked ``no_attached_shock``. A weak shock can be valid
even when its polar maximum exceeds the database range; the unavailable
maximum is then left blank. See :doc:`../models/shock_waves` for frozen-gas
assumptions and temperature ranges.

The **Normal shock** page also supports these gas models with upstream static
temperature, a Mach sweep, and settings replay. Results include upstream and
downstream temperatures. A range failure is retained as an ``out_of_range``
sweep row. The pitot ratio ``p02/p1`` additionally requires an in-range
stagnation temperature: if only that state is unavailable, the pitot value is
left blank with a warning while the valid shock ratios remain available.
Legacy settings continue to select ``AIR``.

The **Conical shock** page offers the same gas models and static-temperature
input, with Mach/cone-half-angle sweeps and settings replay. Thermal results
include upstream and cone-surface temperatures. The entire trajectory from
the shock to the surface must stay in the model range. An unavailable physical
attached limit is left blank while valid weak solutions remain usable.
The angle sweep starts with an editable 0--30 degree range. Switching modes or
editing inputs does not solve the attached limit; numerical calculations start
when the calculation button is pressed.

Thermally perfect expansions
----------------------------

The **Prandtl--Meyer expansion** page supports ``AIR``, ``NASA7``, ``NASA9``,
and ``HARMONIC_OSCILLATOR`` with upstream static temperature for thermal gases.
Results include downstream Mach, static temperature/pressure/density ratios,
upstream/downstream static temperatures, and the maximum turn allowed by the
temperature lower bound. This limit is labelled separately from the perfect
gas's physical limiting turn.

Mach and turn sweeps retain range failures as ``out_of_range`` rows; numerical
convergence failures are ``error`` rows and do not discard successful points.
Single-point failures are displayed on the page. If the sonic reference is
outside the gas range, the valid expansion still appears while only
``nu1/nu2`` are blank, with an explanation. Settings JSON stores gas selection
and SI static temperature; older settings retain ``AIR``. Display-unit changes
preserve the input temperature. See :doc:`../models/expansion_waves` for the
governing equations and model ranges.

Local flow properties and unit Reynolds number
----------------------------------------------

The normal, oblique and conical shock tables include local heat-capacity
ratios, the full velocity-magnitude ratio, dynamic-pressure ratio and
dimensionless entropy rise. Oblique-shock ``V2/V1`` uses the full speed,
including the tangential component; it is different from the normal velocity
ratio. Conical outputs compare the freestream with the cone surface, rather
than the state immediately behind the shock. For frozen ideal gases,
``Delta s / R = -ln(p02/p01)``. Local sound speeds and velocities use
``a = sqrt(gamma(T) R T)`` and ``V = M a``. These require static temperature
but no pressure. Constant-gamma ``AIR`` offers an optional temperature input.
The isentropic page likewise displays local gamma, sound speed and velocity
whenever temperature is available. Select the heat-capacity details checkbox
to include local ``cp`` and ``cv`` in the table and CSV.

Enable absolute-state outputs on a shock page and enter the upstream
**static pressure** (freestream pressure for cones). The isentropic page uses
its existing **total pressure** input. Pressure-dependent columns include
local static pressure, density, dynamic pressure, dynamic viscosity and
unit Reynolds number:

.. math::

   Re' = \frac{\rho V}{\mu(T)},\qquad Re_L = Re' L.

Unit Reynolds number has units of inverse length. Its display follows the
length selector (1/m, 1/mm, 1/ft or 1/in). A separate optional representative
length adds dimensionless ``Re_L``. Ideal-gas shock density is ``p/(R T)``;
the Beattie--Bridgeman isentropic page retains its actual equation-of-state
density and sound speed.

Select a viscosity model independently of the heat-capacity model. These are
the existing dry-air Sutherland, Keyes (79--1845 K) and frozen Blottner/Wilke
(1000--30000 K) correlations from :doc:`../models/transport_properties`.
They evaluate **local static temperature**, never total temperature. Keyes
and Blottner/Wilke are not extrapolated: an unavailable local viscosity leaves
only that state's viscosity and Reynolds columns blank, with a warning. Valid
flow states and other sweep rows remain available. Sutherland has no registered
upper bound; this does not establish accuracy throughout a thermal gas model's
temperature range. These temperature-only transport models describe dilute
air; a high-density thermodynamic model does not add pressure-dependent
transport corrections. The definition follows the `NASA Reynolds-number
description <https://www.grc.nasa.gov/www/BGH/viscosity.html>`_; use of local
static temperature follows the `NASA Wind-US transport documentation
<https://www.grc.nasa.gov/www/winddocs/user/keywords/viscosity.html>`_.

All optional inputs and output selections survive settings JSON replay.
Older settings retain the pressure-free defaults; optional columns appear
when requested and are shared by the displayed table and CSV.

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

Settings are checked for required and unsupported fields, value types, model
choices, finite SI numbers, and numeric input-widget limits before they are
applied. JSON ``NaN`` and infinity values are not accepted. Sweep grids use
the same bound-ordering and spacing validation as calculations. The settings
loader does not predict physical validity or reject individual sweep points:
calculators retain their existing errors and missing-point behavior.
A structural configuration error does not alter the current inputs or queue
values for the next rerun.

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
