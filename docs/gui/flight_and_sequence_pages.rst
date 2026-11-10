Flight-condition and sequential-flow pages
==========================================

These notes cover the **大気・飛行条件** (atmosphere and flight conditions)
page and the **連続流れ計算** (sequential flow) page. The calculations behind
them are explained in :doc:`../guides/atmosphere_flight` and
:doc:`../guides/flow_sequence`.

Atmosphere and flight conditions
--------------------------------

Enter an altitude and either a Mach number or a speed. The page shows the
standard-atmosphere properties at that altitude together with the dynamic
pressure, total conditions, and Reynolds number of the flight. In sweep mode
it varies one input and plots the results. A single-point result can be saved
as the flight case for the boundary-layer pages; see :doc:`index`.

Choosing the kind of altitude
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The page accepts either of two altitudes:

* **Geometric altitude** ``h`` (the default) is the actual height above sea
  level.
* **Geopotential altitude** ``H`` is an adjusted height that allows for gravity
  becoming weaker with height. Standard-atmosphere tables are often given in
  it.

The two are related by ``H = R h / (R + h)`` with ``R = 6,356,766 m``. They
differ by about 16 m at 10 km and about 1.1 km at 86 km.

Your choice applies to the altitude input, the sweep range, and the plot axes.
If you change it, the numbers you entered stay the same and are read as the
other kind of altitude. Both altitudes always appear in the result table and
the CSV file. The allowed range is -5,000 m to 86,000 m geometric, which is
about -5,004 m to 84,852 m geopotential.

The Python functions always use geometric altitude.

Sequential flow
---------------

This page calculates a flow that passes through several shocks and expansions
in a row.

1. Select the gas model.
2. Enter the starting flow. It can be given as static conditions, total
   conditions, or an altitude in the standard atmosphere, with either a Mach
   number or a speed. With the standard atmosphere, the altitude can be
   geometric or geopotential, as on the flight-condition page.
3. Add the stages in order. Each stage can be named, duplicated, deleted, and
   moved up or down.
4. Press **計算** (calculate).

The results show the starting flow and the flow leaving each stage, with a
plot along the stages. If a stage cannot be calculated, the page shows the
reason and leaves the later stages empty; the stages before it are still
shown.

If you change an input after calculating, the page marks the old result as out
of date until you calculate again.

Inputs and the result table follow the display units selected in the sidebar.
The plot and the CSV file on this page are always in SI units. The settings
JSON file stores the whole case, including the stages, so that it can be
loaded again later.
