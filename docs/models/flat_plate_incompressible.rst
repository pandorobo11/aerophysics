Low-speed flat-plate correlations
=================================

This page gives the formulas used with ``compressible=False``. The symbols
are defined in :doc:`flat_plate_boundary_layer`.

Laminar flow (Blasius)
----------------------

The laminar flat plate has an exact solution by Blasius. Its skin friction is

.. math::

   C_{f,x}=\frac{0.664}{\sqrt{Re_x}},
   \qquad
   \bar C_f=\frac{1.328}{\sqrt{Re_x}}.

and its thicknesses are

.. math::

   \frac{\delta_{99}}{x}=\frac{4.9100}{\sqrt{Re_x}},
   \qquad
   \frac{\delta^*}{x}=\frac{1.7208}{\sqrt{Re_x}},
   \qquad
   \frac{\theta_m}{x}=\frac{0.664}{\sqrt{Re_x}}.

Turbulent skin friction (Schlichting)
-------------------------------------

Both low-speed turbulent models use Schlichting's formula for the average
coefficient:

.. math::

   \bar C_f=\frac{0.455}{(\log_{10}Re_x)^{2.58}}.

For constant edge conditions, differentiating accumulated drag gives the
local relation

.. math::

   C_{f,x}=\frac{\mathrm{d}(x\bar C_f)}{\mathrm{d}x}
   =\bar C_f\left(1-\frac{2.58}{\ln Re_x}\right).

Turbulent use outside :math:`5\times10^5\le Re_x\le10^9` emits
:class:`~aerophysics.exceptions.ApplicabilityWarning`.

Turbulent thicknesses
---------------------

The two low-speed turbulent models share the friction and drag above and
differ only in their thicknesses.

.. _flat-plate-power-law-thickness:

One-fifth-power thickness (``SCHLICHTING``, default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \delta_{99}=0.37xRe_x^{-1/5},
   \qquad
   \delta^*=\frac{\delta_{99}}{8},
   \qquad
   \theta_m=\frac{7}{72}\delta_{99}.

This is the classical estimate of :ref:`von Kármán (1921)
<ref-von-karman-1921>`, as given by :ref:`Schlichting and Gersten (2000)
<ref-schlichting-gersten-2000>`. It assumes the one-seventh-power profile
:math:`u/U_e=(y/\delta)^{1/7}`, which gives the two ratios above, and Blasius'
pipe-friction law applied to the plate,
:math:`\tau_w/(\rho_eU_e^2)=0.0225(\nu/U_e\delta)^{1/4}`. Integrating the
momentum relation :math:`\mathrm{d}\theta_m/\mathrm{d}x=\tau_w/(\rho_eU_e^2)`
from the leading edge gives the coefficient
:math:`(0.0225\cdot\tfrac{72}{7}\cdot\tfrac{5}{4})^{4/5}=0.3707`.
The pipe law holds only at low Reynolds number, and the fixed profile gives a
constant shape factor :math:`H=9/7` and :math:`\delta_{99}/\theta_m=72/7`.
Measured layers have :math:`H\approx1.4` and :math:`\delta_{99}/\theta_m\approx8`
near :math:`Re_x=10^6`; equilibrium log-wake profiles give :math:`H` falling to
about 1.2 and :math:`\delta_{99}/\theta_m` rising to about 11 by :math:`10^9`.
The correlation is therefore nominally limited to
:math:`5\times10^5\le Re_x\le10^7`. Turbulent points above
:math:`Re_x=10^7` (for ``TRANSITIONAL``, the equivalent :math:`Re_\mathrm{eq}`
described under :ref:`flat-plate-specified-transition`) emit
:class:`~aerophysics.exceptions.ApplicabilityWarning`.

The Schlichting friction is not derived from these thicknesses. The momentum
thickness therefore differs from :math:`x\bar C_f/2`, the value implied by
the returned drag, by +2% at :math:`Re_x=10^6`, −15% at :math:`10^8`, and −27%
at :math:`10^9`.

.. _flat-plate-log-wake-thickness:

Log-wake profile thickness (``SCHLICHTING_LOG_WAKE``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SCHLICHTING_LOG_WAKE`` uses the momentum integral and an equilibrium
mean-velocity profile instead of fixed ratios:

.. math::

   \theta_m=\frac{x\bar C_f}{2},
   \qquad
   \tau_w=\tfrac{1}{2}\rho_eU_e^2C_{f,x}.

:math:`\bar C_f` is the Schlichting average, or the specified-transition
average, so :math:`\theta_m` always agrees with ``drag_per_unit_width`` and is
continuous at transition. :math:`\delta_{99}` is the value :math:`d` for which
the Spalding inner law (:ref:`Spalding, 1961 <ref-spalding-1961>`) plus the
Coles wake (:ref:`Coles, 1956 <ref-coles-1956>`), with :math:`\kappa=0.41`,
:math:`B=5.2`, :math:`0\le\Pi\le1`, and :math:`u(d)=0.99U_e`, satisfies

.. math::

   \int_0^d\frac{u}{U_e}\left(1-\frac{u}{U_e}\right)\mathrm{d}y=0.978\,\theta_m.

``displacement_thickness`` is :math:`\int_0^d(1-u/U_e)\,\mathrm{d}y`. The wall
unit is :math:`\mu_e/(\rho_eu_\tau)` with
:math:`u_\tau=\sqrt{\tau_w/\rho_e}`, using the supplied
``edge_dynamic_viscosity``. This is the constant-property limit of the
``VAN_DRIEST_II_WALL_TEMPERATURE`` closure in
:doc:`flat_plate_compressible`, with the same
factor 0.978. The two models agree within 0.1% at Mach 0.03.

Where no wake parameter :math:`0\le\Pi\le1` matches, :math:`\delta_{99}` and
:math:`\delta^*` are ``NaN`` and one
:class:`~aerophysics.exceptions.ApplicabilityWarning` reports the number of
such points; :math:`\theta_m`, friction, and drag are still returned. Fully
turbulent plates can be calculated for :math:`10^4\le Re_x\le10^{10}`.

Thickness accuracy
~~~~~~~~~~~~~~~~~~

The log-wake thickness was compared with low-speed measurements and DNS (see
:doc:`../verification/flat_plate_boundary_layer`). Its
:math:`\delta_{99}` is about 19% thick at :math:`Re_x=10^6`, 15% at
:math:`10^7`, 10% at :math:`10^8`, and 9% at :math:`1.5\times10^8`. Most of
this offset comes from the Schlichting friction, which makes :math:`\theta_m`
about 5.5% thick. Multiplying by the ratio of the two models gives the
estimated one-fifth-power error:

.. list-table:: Fully turbulent thickness ratios, ``SCHLICHTING`` /
   ``SCHLICHTING_LOG_WAKE``
   :header-rows: 1
   :widths: 14 16 16 16 38

   * - :math:`Re_x`
     - :math:`\delta_{99}`
     - :math:`\theta_m`
     - :math:`\delta^*`
     - Estimated ``SCHLICHTING`` :math:`\delta_{99}` error against measurements
   * - :math:`10^6`
     - 1.18
     - 1.02
     - 0.96
     - +40%
   * - :math:`10^7`
     - 1.01
     - 0.95
     - 0.98
     - +16%
   * - :math:`10^8`
     - 0.83
     - 0.85
     - 0.90
     - −8%
   * - :math:`1.5\times10^8`
     - 0.80
     - 0.83
     - 0.88
     - −12%
   * - :math:`10^9`
     - 0.67
     - 0.73
     - 0.79
     - not measured

The one-fifth-power error changes sign near :math:`Re_x\approx5\times10^7`
and grows quickly beyond :math:`10^8`, while the log-wake error decreases
steadily with Reynolds number. ``SCHLICHTING_LOG_WAKE`` is recommended when
:math:`\delta_{99}`, :math:`\delta^*`, or a :math:`\theta_m` consistent with drag
matters. ``SCHLICHTING`` remains the default for backward compatibility.
