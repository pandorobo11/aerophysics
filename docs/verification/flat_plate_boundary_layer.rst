Flat-plate boundary-layer comparisons
=====================================

This record collects the comparisons behind the accuracy statements in
:doc:`../models/flat_plate_boundary_layer` and its subpages.

Compressible laminar similarity solution
----------------------------------------

The unchanged independent ``solve_bvp`` reference generator
``docs/scripts/build_laminar_similarity_reference.py`` produces the CSV and
provenance JSON under
``tests/reference_data/viscous_flow/laminar_compressible_similarity``.
It uses Sutherland air, :math:`T_e=240\ \mathrm{K}`, :math:`Pr=0.72`,
:math:`\gamma=1.4`, :math:`M_e=0.3,2,4,6,10,15`, and
:math:`T_w/T_{aw}=0.1,0.3,0.5,1` (24 points). Isothermal API inputs use the
reference's solved :math:`T_{aw}` times the prescribed ratio; the last wall
condition is truly adiabatic. Errors below compare the library solution with
this independent reference, with equal weight for each point.

.. include:: ../_generated/laminar_similarity_accuracy.rst

The acceptance bound is 0.1% relative for friction, momentum and delta99.
For displacement coefficients of magnitude at least 0.05 it is also 0.1%;
closer to zero the absolute coefficient error must be at most 5e-5.
The table separates these error measures and suppresses absolute errors below
5e-10 to avoid reporting insignificant CPU differences. These are numerical
checks of constant-cp perfect-gas equations, not physical validation.
The generated-assets suite also checks 200 seeded random states over
Mach 0--20, Te 30--2000 K, and Tw/Taw 0.05--1.5, with selected independent
collocation comparisons; timing is not a CI requirement.

Independent adiabatic examples are recovered within 0.5%: at M2,
:math:`c_f\sqrt{Re_x}=0.6356105`, :math:`\delta_{99}\sqrt{Re_x}/x=6.414710`,
:math:`\delta^*\sqrt{Re_x}/x=3.303259`, :math:`T_{aw}/T_e=1.676431`;
at M10 the corresponding values are 0.4258457, 29.674188, 27.668758, and
17.628489. Tests also check the Blasius limit and reference convergence
with domain enlargement (14 to 20) and tolerance tightening (1e-8 to 1e-10).

Van Driest II with wall-temperature correction
----------------------------------------------

External comparisons
~~~~~~~~~~~~~~~~~~~~

The tested envelope is :math:`M_e=0.3`--13.3, :math:`T_w/T_{aw}=0.07`--1,
:math:`T_\infty=47`--271 K, and :math:`Re_{x,e}\approx10^6`--:math:`10^9`
(:math:`M_e<2` only up to :math:`1.6\times10^8`), for a smooth
zero-pressure-gradient plate and a
calorically perfect gas. Flow is turbulent from the leading edge, or uses
the prescribed transition model. This envelope is evidence of comparison,
not a guarantee for every combination of these parameters; the wake closure
can fail.

The following summaries are externally supplied validation results dated
2026-11-01, not simulations rerun by this package's test suite. Differences
are model minus reference, divided by reference. "Standard SA" below means
``SA_OPTIONS=NONE`` in SU2 8.5.0, whose convection differs from the published
model at high speed; see :ref:`the note on SU2 convection <su2-sa-convection>`.

* Against SU2 8.5.0 standard SA, :math:`\delta_{99}` RMS differences are
  **2.5%** for the eight fitted conditions (343 sections), **2.3%** for six
  held-out conditions (264 sections), and **2.1%** for long plates through
  :math:`Re_{x,e}=1.7\times10^8`. Above :math:`10^8`, the RMS difference is
  **5.2%**, spanning **−7.6 to +0.6%**.
* Displacement thickness integrated through :math:`\delta_{99}` differs by
  **2.3--2.7% RMS**.
* Four SU2 standard-SA cases at DNS-matched conditions vary freestream
  temperature and unit Reynolds number: :math:`T_\infty=47`--271 K,
  :math:`Re'=4.8`--:math:`10.6\times10^6\ \mathrm{m}^{-1}`,
  :math:`M_e=2.5`--13.3, :math:`T_w/T_{aw}=0.19`--0.99, and
  :math:`Re_{x,e}\le9\times10^7`. Casewise :math:`\delta_{99}` RMS differences
  are **1.8--4.0%** (range **−8.9 to +7.5%**); displacement-thickness RMS
  differences are **1.2--4.0%**.
* The :math:`\delta_{99}` difference falls almost linearly with
  :math:`\log_{10}Re_{x,e}`, by **5--6 percentage points per decade**, with no
  levelling off up to :math:`10^9`. Seven further SU2 standard-SA cases
  (M2, M6 and M10 at :math:`Re'=1.2\times10^8\ \mathrm{m}^{-1}` and
  :math:`T_w/T_{aw}=0.07`--0.5; M1.5, M0.8 and M0.3 at
  :math:`Re'=2\times10^7\ \mathrm{m}^{-1}` with the wall at recovery
  temperature) give :math:`\delta_{99}` differences of **+0.2 to +1.1%** at
  :math:`Re_{x,e}=10^7`, **−6.3 to −3.5%** at :math:`10^8`, and **−12.1 to
  −8.6%** at :math:`10^9` for the supersonic cold walls. For the
  recovery-temperature walls at :math:`M_e\le1.5`, the differences are
  **+9 to +11%** at :math:`10^6`, **+4 to +6%** at :math:`10^7`, and **−0.5 to
  +1.1%** at :math:`10^8`. The slope is the same at M0.3, so it comes from the
  profile shape, not from the compressibility transformation. Whether the
  correlation or the SA model is closer to measured :math:`\delta_{99}` at
  high Reynolds number has not been established.
* Over the same cases, :math:`C_f` differs by **−1.3 to +2.0%** mean and
  :math:`\theta` by **−1.3 to +3.5%** mean (RMS at most **2.2%** and
  **3.8%**), up to :math:`Re_{x,e}=10^9` and down to M0.3.
* Against incompressible zero-pressure-gradient measurements and DNS
  (:ref:`Österlund et al. (2000) <ref-osterlund-2000>`,
  :ref:`Schlatter and Örlü (2010) <ref-schlatter-orlu-2010>`,
  :ref:`Vallikivi, Hultmark and Smits (2015) <ref-vallikivi-hultmark-smits-2015>`;
  :math:`Re_\theta=2\times10^3`--:math:`2.3\times10^5`, about
  :math:`Re_x=10^6`--:math:`2.7\times10^8`), :math:`\delta_{99}` is too thick
  by about **+19%** at :math:`Re_x=10^6`, **+15%** at :math:`10^7`, **+10%**
  at :math:`10^8`, and **+9%** at :math:`1.5\times10^8`. SU2 standard SA is
  **+7 to +9%** thick over the same range, so most of the SA-relative drift
  above comes from the correlation approaching the measurements as
  :math:`Re_x` rises. The difference comes mainly from the incompressible
  Schlichting friction pair: at a given :math:`Re_\theta` its :math:`C_f` is
  **6--6.5%** above the measurements, which makes :math:`\theta` about
  **5.5%** thick and lowers the solved wake parameter (0.33--0.45 against
  about 0.6 measured). Given the measured :math:`\theta` and :math:`C_f`, the
  profile closure reproduces :math:`\delta_{99}` within **−6 to +2%** on
  average. The widely used friction pair is kept unchanged.
* The solved wake parameter rises with Reynolds number on cold hypersonic
  walls: at M10 with :math:`T_w/T_{aw}=0.07` it reaches :math:`\Pi=0.976` at
  :math:`Re_{x,e}=10^9`. Extrapolating that trend, the wake bound
  :math:`\Pi\le1` would be exceeded, and
  :class:`~aerophysics.exceptions.ModelRangeError` raised, near
  :math:`Re_{x,e}\approx1.5\times10^9` for that condition.
* At matched :math:`Re_\theta`, comparison with the DNS of
  :ref:`Zhang, Duan and Choudhari (2018) <ref-zhang-duan-choudhari-2018>`
  gives :math:`C_f` differences **−0.3, −0.6, −3.8, −2.2%**, and
  :math:`\delta_{99}` differences **−4.5, −0.5, +0.6, +15.3%**, respectively
  for M2.5, M6 with :math:`T_w/T_r=0.76`, M6 with :math:`T_w/T_r=0.25`, and
  M14 with :math:`T_w/T_r=0.18`. The M14 result shows that the good M6
  agreement should not be extrapolated to higher Mach numbers.

In-repository SA boundary-layer reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``docs/scripts/build_turbulent_sa_reference.py`` is an independent
verification solver, not part of the package. It marches the compressible
zero-pressure-gradient thin-layer equations from the leading edge with the
:ref:`Spalart--Allmaras <ref-spalart-allmaras-1994>` model without the trip
term (SA-noft2), in Levy--Lees
variables, for a calorically perfect gas with constant :math:`Pr=0.72`,
turbulent Prandtl number 0.9, :math:`\gamma=1.4` and the Sutherland law
(:math:`S=110.4` K). The SA transport is the standard form, in which density
does not appear (:ref:`Allmaras, Johnson and Spalart, 2012
<ref-allmaras-johnson-spalart-2012>`, section 2.2), its edge value is
:math:`\tilde\nu=3\nu_e`, and the modified vorticity uses the smooth limiter of
the same paper. It uses no correlation, fitted
coefficient or package model. With :math:`\tilde\nu=0` it reproduces the
compressible laminar similarity solution within 0.01%. Across all tabulated
rows, the largest change in friction or a thickness is 0.03% with twice the
streamwise stations, about 0.1% with :math:`\eta` grid ratio 1.02 instead of
1.04, 0.07% with half the first spacing, and below 0.001% with twice the outer
boundary :math:`\eta_{\max}=150`. Each condition takes about 1--1.5 seconds to march with the coupled Newton
iteration and occasional Picard fallback (measured locally).

Its stored results are
``tests/reference_data/viscous_flow/turbulent_sa_boundary_layer.csv``. The
table compares them with this model at the same edge state and wall
temperature. Displacement thickness is integrated through
:math:`\delta_{99}` in both; momentum thickness is integrated to the edge.

.. include:: ../_generated/turbulent_sa_comparison.rst

Friction agrees within 3% at :math:`Re_{x,e}\ge10^7`, as in the SU2 summaries
above. At M0.3 the :math:`\delta_{99}` differences also fall within those
summaries at :math:`Re_{x,e}=10^7` and :math:`10^8` and just above them at
:math:`10^6`. On cold supersonic walls they are lower than the summaries at
:math:`Re_{x,e}=10^7`: :math:`-3.0` to :math:`-3.5\%` here, against
:math:`+0.2` to :math:`+1.1\%` for SU2, because the two solvers do not solve the
same SA convection, as described next.

.. _su2-sa-convection:

**SU2 convection.** :ref:`SU2 8.5.0 <ref-su2-sa-convection>` discretizes SA
convection as :math:`\nabla\cdot(\mathbf{u}\tilde\nu)` with the volume flux,
whereas the published model has :math:`\mathbf{u}\cdot\nabla\tilde\nu`. The two
differ by :math:`\tilde\nu\,\nabla\cdot\mathbf{u}`, which is not negligible in a
high-speed boundary layer. Adding this term to a scratch copy of this reference
reproduces the stored SU2 fields at matched local edge states: at
:math:`M_e\approx9.8`, :math:`T_w=300` K, :math:`Re_{x,e}\approx10^7`, the peak
:math:`\tilde\nu/\nu_e` is 343 against 345 in SU2 and 380 in the reference, and
:math:`u/U_e` agrees within 0.001 through the layer. Over the eight SU2
conditions at this :math:`Re_{x,e}`, the difference of the reference from SU2
falls from :math:`+0.5` to :math:`+7.4\%` in friction to :math:`-0.1` to
:math:`+2.0\%`, and from :math:`+0.5` to :math:`+5.0\%` in :math:`\delta_{99}` to
:math:`-2.8` to :math:`-0.1\%`.

SU2 itself confirms the cause. The same :math:`M_\infty=10`, :math:`T_w=300` K
case was rerun on the same grid with development commit ``f4b67d2``, starting
from the converged 8.5.0 solution. With the default ``SCALAR_UPWIND`` that commit
holds the 8.5.0 solution. ``CONV_NUM_METHOD_TURB= BOUNDED_SCALAR`` subtracts the
face-flux divergence from the convection, leaving
:math:`\mathbf{u}\cdot\nabla\tilde\nu`, and also takes its upwind weights from
the flow solver's mass flux; these two changes were not run separately. At
:math:`x=0.5` m, where :math:`Re_{x,e}\approx10^7`, friction rises 3.4%,
:math:`\delta_{99}` 6.8%, momentum thickness 4.2% and wall heat flux 3.6%, and the
peak :math:`\tilde\nu/\nu_e` rises from 346 to 385. The scratch-copy comparison
above predicts 3.2%, 6.1% and 3.8% for the first three, and 343 to 380 for the
peak.

Over :math:`Re_{x,e}\approx10^6` to :math:`2\times10^7`, the reference then
differs from this rerun by :math:`+1.7` to :math:`+4.0\%` in friction and
:math:`-4.6` to :math:`-1.3\%` in :math:`\delta_{99}`, close to the
:math:`+1.7` to :math:`+4.3\%` and :math:`-3.6` to :math:`-0.8\%` between the
scratch copy and SU2 8.5.0. Convection therefore does not explain this remainder,
which has not been separated from the viscous-interaction pressure history, the
grid or the leading-edge development. Only this one condition was rerun.

At the :ref:`DNS <ref-zhang-duan-choudhari-2018>` momentum-thickness Reynolds
numbers, friction from the published form differs from DNS by :math:`-5.9`,
:math:`-8.3`, :math:`-1.5` and :math:`+1.2\%` for M2.5, M6 with
:math:`T_w/T_r=0.25`, M6 with :math:`T_w/T_r=0.76` and M14 with
:math:`T_w/T_r=0.18`; SU2 differs by :math:`-7.2`, :math:`-11.4`, :math:`-4.9`
and :math:`-4.3\%`. These reference runs approximate the edge temperature by
the freestream value. The constants below were fitted to the SU2 results and
therefore include this convection difference; they have not been re-evaluated
against the published form.

At :math:`Re_{x,e}=10^6` this model, which is fully turbulent from the leading
edge, gives 5.5--7.4% higher friction than SA, whose turbulence develops from
the edge value of :math:`\tilde\nu`, except at M10. There SA stays nearly
laminar up to about :math:`Re_{x,e}=10^5` and keeps that thicker layer
downstream: friction differs by +0.8% while :math:`\delta_{99}` differs by
:math:`-14.6\%`, narrowing to :math:`-3.2\%` at :math:`10^7`. The reference is a RANS model result, not a
measurement: SA itself gives incompressible :math:`\delta_{99}` about 7--9%
above the measurements cited above.

The exponent **0.064** is an empirical fit to SU2 solutions; DNS also exhibits
lower cold-wall friction than VD2. At the four
:ref:`DNS <ref-zhang-duan-choudhari-2018>` conditions listed in the SU2 note above,
with the same edge-temperature approximation and matched momentum-thickness
Reynolds number, friction differs from DNS by :math:`-0.4`, :math:`-3.2`,
:math:`+0.8` and :math:`+0.5\%` with this factor, and by :math:`-0.4`,
:math:`+7.8`, :math:`+2.8` and :math:`+14.1\%` without it. These four conditions
are the only DNS check of the exponent, and it was not fitted to them. The factor
**0.978** is the SU2 mean ratio of
momentum thickness through :math:`\delta_{99}` to full-edge momentum thickness
(the truncation difference varies from **−5 to −1%**). Neither constant is
universal. The nine-row independent-implementation regression fixture uses
exponent 0.06391; rounding to 0.064 changes the factor by less than 0.03%.
Its ``su2_delta99`` column is comparison data, not a unit-test acceptance target.

Thicknesses formerly returned by Eckert and Van Driest II
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Across the supplied 1059-section SU2 comparison (standard SA, M2--M10,
:math:`T_w/T_{aw}=0.07`--1), Eckert's former :math:`\delta_{99}` estimate had
mean error **+58%** (range **+5 to +92%**); VD2 had mean error **+17%**
(**18% RMS**, range **−1 to +37%**). Fixed ratios :math:`\delta^*=\delta/8`
and :math:`\theta=7\delta/72` missed displacement thickness by factors of
2--5 (underprediction) and momentum thickness by factors of 2--10
(overprediction). Consequently ``ECKERT`` and ``VAN_DRIEST_II`` now return
``NaN`` for all three turbulent thicknesses, including post-transition
points. Friction and drag retain their existing definitions. Laminar thicknesses
use the closure above; ``compressible=False`` uses the Blasius constants stated earlier.

Wall heat transfer
------------------

.. _sa-reference-heat-transfer:

The :math:`q_w` column of the table in the SA reference section compares
``VAN_DRIEST_II_WALL_TEMPERATURE`` with :math:`s=1.14`. The difference is
:math:`-5.4` to :math:`+5.4\%` at :math:`Re_{x,e}\ge10^7` on the walls tested
(:math:`M_e=5`--10, :math:`T_w/T_r=0.20`--0.76). It reaches :math:`+20\%` at
:math:`Re_{x,e}=10^6`, :math:`T_w/T_r=0.76`, where SA is still developing and
:math:`T_r-T_w` is sensitive to the recovery factor. SA's own recovery factor
falls from about 0.91 at Mach 0.3 to 0.84--0.87 at Mach 10 for
:math:`Re_{x,e}\ge10^7` instead of staying at :math:`Pr^{1/3}=0.896`, so ``St`` should not be compared near :math:`T_w=T_r`.
``ECKERT`` and ``VAN_DRIEST_II`` friction also produce heat-transfer values, but
they carry their own friction differences and have not been compared in this
way. All of this assumes a perfect gas, constant :math:`c_p` and :math:`Pr`, an
isothermal wall and zero pressure gradient.
