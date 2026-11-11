# GUI internals

User-facing GUI behavior is described in the `docs/gui/` section of the
manual. This file records settings-file fields, documentation links, and
failure handling that contributors need when changing a GUI page.

## Altitude coordinate in settings JSON

- The atmosphere/flight page stores the resolved geometric altitude in
  `inputs_si.geometric_altitude` and the selected coordinate in
  `models.altitude_basis` (`"geometric"` or `"geopotential"`).
- Sweep bounds are stored in the selected coordinate, in metres.
- Older settings without `altitude_basis` use geometric altitude.
- The flow-sequence page stores `inputs_si.initial.altitude` in the
  coordinate named by `inputs_si.initial.altitude_basis` (geometric when the
  field is omitted). The Python `FlowState.from_atmosphere` API always uses
  geometric altitude.
- Altitude limits are checked after conversion to geometric altitude:
  -5,000 to 86,000 m geometric, approximately -5,003.936 to 84,852.046 m
  geopotential.

## Settings validation

Settings are checked for required and unsupported fields, value types, model
choices, finite SI numbers, and numeric input-widget limits before they are
applied. JSON `NaN` and infinity values are rejected. Sweep grids use the same
bound-ordering and spacing validation as calculations.

The loader does not predict physical validity or reject individual sweep
points; calculators keep their usual errors and missing-point behavior. A
structural configuration error leaves the current inputs unchanged and does
not queue values for the next Streamlit rerun.

Older settings files without the gas-model fields added for thermal shocks
and expansions load with the constant-gamma `AIR` model. Older settings
without the optional pressure, viscosity, and representative-length inputs
load with those outputs disabled.

## Legacy flat-plate settings

Saved schema-version-1 flat-plate configurations that use the former separate
correlation and correction keys are converted when loaded, except for
turbulent or transitional one-fifth-power friction with NONE or ECKERT: those
settings raise an error asking the user to select Schlichting. Legacy
power-law selections for laminar or VD2 flow still load because the friction
selection was unused in those cases. Exports use only `compressible` and
`turbulent_model` (null for laminar flow).

## Documentation links

`render_documentation_link` in `gui/documentation.py` adds a theory-page
button below each calculator title. It renders nothing when the local manual
is not available. Every linked path must also be a value of
`DOCUMENTATION_TOPICS`; `tests/gui/test_documentation.py` checks this. When a
manual page is renamed or split, update both the topic list and the page
links.

## Solver failures on GUI pages

Sweep rows use these status values:

| Status | Meaning |
| --- | --- |
| `out_of_range` | A required thermodynamic state is outside the gas model's range (`ModelRangeError`). |
| `no_attached_shock` | The requested deflection or cone angle exceeds the physical attached limit (`NoAttachedShockError`). |
| `error` | Numerical convergence failure (`ShockConvergenceError` or `ExpansionConvergenceError`). |

- A single-point failure is displayed on the page.
- A failed sweep point is kept with its status, and the remaining points are
  still evaluated.
- If only an auxiliary value fails after a valid state was found (for
  example the attached-limit angle of a conical shock, or the absolute
  Prandtl--Meyer angles when the sonic reference is outside the gas range),
  the valid result is kept and only that value is left blank.
- Unrelated runtime errors are not suppressed.

The conical-shock page does not solve the attached-limit angle while the user
switches modes or edits inputs; numerical work starts only when the
calculation button is pressed.

## Axisymmetric body-flow page

`axisymmetric_body` is single-point only, with `AIR` fixed and no sweep.
`body_flow_adapter.py` owns pure geometry construction, atmosphere resolution,
SI rows and Cp; it does not change the public solver. `body_flow_page.py`
uses the `body_flow_` widget prefix and `body_flow_payload`. The shared
`calculation_button` accepts an optional explicit key so this page's calculate
button also follows that prefix. Unit changes convert existing input state.

Settings `inputs_si` contains `mach`, resolved `pressure` and `temperature`,
`altitude` (null in static mode; otherwise metres in `models.altitude_basis`),
`nose_length`, `base_radius`, `ogive_segments` (2--200), `cylinder_length`,
nullable `boattail_length` / `boattail_radius`, nullable number lists
`meridian_x` / `meridian_r`, `characteristic_count` (4--1000), and
`fan_angle_step` (0.05--5 degrees stored as radians). Geometry lengths are
metres. Models are `geometry_source` (`parametric`/`csv`), `nose`
(`cone`/`tangent_ogive`), boolean `boattail`, `freestream_basis`
(`static`/`atmosphere`), `altitude_basis` (`geometric`/`geopotential`) and
`pressure_quantity` (`cp`/`pressure_rise`) and `net_coloring`
(`none`/`mach`/`cp`, default `none` on legacy imports). CSV settings embed SI vertices
and can replay without an upload. An upload is parsed once per Streamlit
`file_id` and its SI vertices are cached in `body_flow_meridian_upload_si`, so
a display-unit change never reinterprets the same file in the new unit. Geometric/physical validity is checked at
calculation time, including concavity in the solver. Tangent ogives require
`base_radius <= nose_length` for the specified circle parameterization.

The payload holds surface `CalculationResult`, shock rows, the raw
`AxisymmetricBodyFlowResult` and the calculated configuration. Display changes
reuse it. Export updates both display selectors and display units but keeps
the calculated inputs. Input changes show a previous-result warning.
`ValueError` (including CSV errors, `ModelRangeError` and
`NoAttachedShockError`) and `ShockConvergenceError` display their original
message via `st.error`; a new failed calculation clears the old payload.
Other exception types propagate. `current_axisymmetric_body_case` holds the
raw result only after an explicit save; the marching boundary-layer page
consumes it.


The adapter requests `include_characteristic_net=True` and retains it in
`raw`; no solver invocation occurs on pressure, coloring or unit changes.
The fourth figure tab uses one `None`-separated trace per characteristic
family (C+ blue, C- green; body and shock follow the cone-shock styling), capped at 150 rows and 80 C- polylines by evenly
spaced indices including the endpoints. `layout.meta.downsampled` drives the
page caption. Optional `Scattergl` markers sample at most 20,000 nodes from
the displayed rows. Cp uses the saved calculation's freestream, even when
current inputs are stale. The long-format `axisymmetric_body_net` table schema
converts all nodes to display units for CSV; display sampling never affects
export. Row/node indices are zero-based integer strings, preserving the
shared `CellValue` contract without changing other calculators; `ends_at_shock` describes the whole
row, not each node. GUI regressions cover four tabs, caps, CSV contents,
legacy settings, and selector changes without recalculation or warnings.


`body_flow_net_csv_cache` retains `(calculation, display_units, csv_text)` in
session state. Calculation identity is compared with `is`, avoiding reused
numeric IDs. CSV generation runs once for the current calculation and unit
selection; coloring, pressure quantity and unsaved input edits reuse the
text. Unit changes replace the entry. A new calculation attempt or missing
payload clears it, including failed calculations; reset and settings import
clear it through the `body_flow_` prefix. The download button receives the
cached string directly, without depending on callable data support.

## Marching boundary-layer page

`marching_boundary_layer` is single-case only, with `AIR` and Sutherland
viscosity fixed and no sweep. `marching_adapter.py` owns the pure SI edge
construction (`uniform_edge`, `body_case_edge`), the solver call, station rows
and profile rows; it does not change the public solver. `marching_page.py`
uses the `march_bl_` widget prefix and `march_bl_payload`.

Settings `inputs_si` contains `mach` (at least 0.01), `pressure`,
`temperature`, `length`, `station_count` (2--1000), nullable `cone_half_angle`
(0.1--89 degrees stored as radians; null unless the uniform cone is selected),
resolved `total_pressure` / `total_temperature`, nullable number lists
`edge_distance` / `edge_pressure` / `edge_radius` (null for the uniform edge;
`edge_radius` is also null for planar CSV input), nullable `wall_temperature`
(null for an adiabatic wall), nullable `transition_distance` (null unless
transitional), `prandtl_number`, `turbulent_prandtl_number` and
`edge_sa_ratio`. Models are `edge_source` (`uniform`/`csv`/`body`), `geometry`
(`planar`/`cone`, used by the uniform edge), `regime`
(`laminar`/`turbulent`/`transitional`) and `wall` (`adiabatic`/`isothermal`).
The uniform-edge fields keep their last values when another source is
selected.

CSV and body-case edges are embedded as SI arrays, so settings replay without
an upload or a saved case. An upload is parsed once per Streamlit `file_id`
and cached in `march_bl_edge_upload_si`; removing the upload drops that cache
and the arrays it supplied, while arrays from an imported settings file (which
leave no cache) stay in use. `edge_csv_template` writes the template in the
current display units. Selecting another edge source drops the embedded arrays. For the body source the session's
`current_axisymmetric_body_case` is converted on every rerun, except after a
settings import: `march_bl_keep_embedded` then keeps the loaded arrays until
the user presses the switch button, which is rendered only when both exist.

The payload holds `MarchingCalculation` (station `CalculationResult` with
captured `ApplicabilityWarning` messages, the raw result, the edge and an
optional separation distance) and the calculated configuration. The profile
station selector (`march_bl_profile_station`, reset by each calculation) and
display units reuse the payload. `BoundaryLayerSeparationError` with at least
one solved station is returned by the adapter as a truncated result plus
`separation_distance`; the page shows an error banner above that result.
A separation without solved stations, other `ValueError`s and
`BoundaryLayerConvergenceError` display the message via `st.error` and clear
the old payload. Other exception types propagate.
