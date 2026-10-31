Compressible-flow pages
=======================

These notes cover the **等エントロピー流れ**, **垂直衝撃波**, **斜め衝撃波**,
**円錐衝撃波**, **離脱衝撃波**, and **Prandtl–Meyer膨張** pages. For the
theory, use each page's **📖 理論・モデルの説明** button or the table in
:doc:`index`.

Gas models on flow pages
------------------------

The **斜め衝撃波** (oblique shock), **垂直衝撃波** (normal shock),
**円錐衝撃波** (conical shock), and **Prandtl–Meyer膨張** pages use the
constant-gamma ``AIR`` model by default. Select ``NASA7``, ``NASA9``, or
``HARMONIC_OSCILLATOR`` and enter the upstream **static** temperature to use a
temperature-dependent heat capacity. Mach-number and angle sweeps work with
all models, and the results show the upstream and downstream temperatures.
Older settings files without a gas selection load with ``AIR``.

Some states fall outside a gas model's temperature range. The GUI then keeps
every valid result and marks only the affected item:

* In a sweep, a point outside the gas temperature range is marked
  ``out_of_range``. A deflection or cone angle beyond the attached-shock limit
  is marked ``no_attached_shock``. A numerical convergence failure is marked
  ``error``. The other points are still calculated.
* A weak oblique or conical shock can be valid even when the attached-shock
  limit is outside the gas data. The limit is then left blank.
* On the normal-shock page, the pitot ratio ``p02/p1`` also needs the
  stagnation temperature to be in range. If only that is out of range, the
  pitot value is left blank with a warning.
* On the expansion page, the maximum turn allowed by the gas data's lower
  temperature limit is shown separately from the perfect-gas physical limit.
  If the sonic reference state is out of range, the absolute angles
  ``nu1/nu2`` are left blank with an explanation.

On the conical-shock page, the cone-angle sweep starts with an editable
0--30 degree range. The attached-limit angle is only calculated when you
press the calculation button.

Local flow properties and unit Reynolds number
----------------------------------------------

The normal, oblique, and conical shock tables include the local heat-capacity
ratio, the ratio of flow speeds ``V2/V1``, the dynamic-pressure ratio, and the
dimensionless entropy rise. For an oblique shock, ``V2/V1`` uses the full
speed, including the component along the shock, so it differs from the
normal-velocity ratio. Conical outputs compare the free stream with the cone
surface, not with the state just behind the shock. For a frozen ideal gas,
``Delta s / R = -ln(p02/p01)``.

Local sound speed and velocity use ``a = sqrt(gamma(T) R T)`` and ``V = M a``.
They need the static temperature but not the pressure. The constant-gamma
``AIR`` model offers an optional temperature input for this purpose. The
**等エントロピー流れ** page also shows local gamma, sound speed, and velocity
when a temperature is available. Tick the heat-capacity details checkbox to
add local ``cp`` and ``cv`` to the table and CSV.

To add absolute-state outputs, enable them on a shock page and enter the
upstream **static pressure** (free-stream pressure for cones). The
isentropic page uses its **total pressure** input instead. The extra columns
are static pressure, density, dynamic pressure, dynamic viscosity, and unit
Reynolds number:

.. math::

   Re' = \frac{\rho V}{\mu(T)},\qquad Re_L = Re' L.

Unit Reynolds number is the Reynolds number per unit length, so its unit
follows the length selector (1/m, 1/mm, 1/ft, or 1/in). Enter an optional
representative length to add the dimensionless ``Re_L``. Shock pages use the
ideal-gas density ``p/(R T)``; the Beattie--Bridgeman option on the isentropic
page uses its own equation-of-state density and sound speed.

The viscosity model is chosen separately from the gas model. The options are
the dry-air Sutherland, Keyes (79--1845 K), and Blottner/Wilke (1000--30000 K)
correlations from :doc:`../models/transport_properties`. They are evaluated at
the **local static temperature**, never the total temperature. Keyes and
Blottner/Wilke are not extrapolated: outside their ranges only the viscosity
and Reynolds-number columns of that state are left blank, with a warning.
Sutherland has no upper limit in the code, but this does not mean it is
accurate over a thermal gas model's whole range. These viscosity models
describe low-density air and include no high-pressure corrections. The
definition follows the `NASA Reynolds-number description
<https://www.grc.nasa.gov/www/BGH/viscosity.html>`_; the use of local static
temperature follows the `NASA Wind-US transport documentation
<https://www.grc.nasa.gov/www/winddocs/user/keywords/viscosity.html>`_.

Flow schematics
---------------

Single-point results on the shock and expansion pages include a schematic of
the flow. Flow runs from left to right. Arrows show flow **direction** only;
their lengths do not represent speed. Angles are drawn in the selected display
unit.

* **Normal shock:** a red shock perpendicular to the flow, labelled with the
  upstream and downstream Mach numbers. At upstream Mach 1 the caption notes
  that the shock has zero strength.
* **Prandtl--Meyer expansion:** a grey convex wall turning clockwise through
  the expansion angle ``theta``, with the flow above it. A green fan spans from
  the upstream Mach line (at ``mu_1`` from the upstream flow) to the
  downstream Mach line (at ``mu_2`` from the downstream flow), where
  ``mu = asin(1/M)``. The downstream arrow follows the wall. The lines inside
  the fan are representative Mach lines. At zero turn, a single Mach line is
  shown and labelled as no expansion. The Prandtl--Meyer function ``nu`` is
  not drawn, because it is not a geometric angle.
* **Oblique and conical shocks:** a section drawn to scale. The grey body is
  the wedge or cone, the red line is the attached shock, and the dashed
  horizontal line is the upstream direction (and the cone axis). Arcs show the
  wall angle (``theta`` or the cone half-angle ``theta_c``) and the shock angle
  ``beta``. The oblique downstream arrow is parallel to the wedge. For the
  cone, the downstream arrow shows only the surface direction and surface Mach
  ``M_s``; between the shock and the surface the flow turns gradually. Only the
  upper half of the cone section is shown.
* **Detached shock:** the grey body, red bow shock, and blue free-stream arrow.
  The origin ``O`` is the centre of nose curvature. ``Rn`` marks the nose
  radius, and the bracket ``Delta`` marks the standoff distance between the
  nose and the shock. The shock curve is Billig's shape with the
  Ambrosio--Wortman standoff, even when the table also shows a Seiff estimate.
  The model uses positive ``x`` pointing upstream; the plot reverses the
  horizontal axis so that the flow runs left to right. Hover values and
  exported CSV coordinates keep the upstream-positive convention. See
  :doc:`../models/detached_shocks`.

Sweep results show trend plots instead of schematics.

Shock-shape export
------------------

The **離脱衝撃波** (detached shock) page can export the calculated
shock-shape coordinates as CSV, in the upstream-positive coordinate system
described above.
