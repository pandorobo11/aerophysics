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
