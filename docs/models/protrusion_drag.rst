.. _protrusion-drag-model:

Boundary-layer protrusion drag
==============================

:func:`aerophysics.protrusion.protrusion_drag` estimates direct drag on one
isolated protrusion by integrating the undisturbed local dynamic pressure over
its frontal area. For a task-oriented example and guidance on supplying a
computed boundary-layer profile, see :doc:`../guides/boundary_layers`.
Background on aircraft excrescence drag and engineering drag coefficients is
provided by :ref:`Young and Paterson (1981) <ref-young-paterson-1981>` and
:ref:`Hoerner (1965) <ref-hoerner-1965>`.

Model equations
---------------

For height :math:`h` and projected width :math:`b(y)`,

.. math::

   A_f=\int_0^hb(y)\,\mathrm{d}y,
   \qquad
   q_\mathrm{eff}
   =\frac{1}{A_f}\int_0^h
    \frac{1}{2}\rho(y)U(y)^2b(y)\,\mathrm{d}y.

The direct drag and shielding factor are

.. math::

   D=C_Dq_\mathrm{eff}A_f,
   \qquad
   \eta_q=\frac{q_\mathrm{eff}}{q_e},
   \qquad
   q_e=\frac{1}{2}\rho_eU_e^2.

The default turbulent velocity approximation is

.. math::

   \frac{U}{U_e}
   =\min\left[\left(\frac{y}{\delta}\right)^{1/7},1\right].

For constant width and density with :math:`h\le\delta`,

.. math::

   \eta_q=\frac{7}{9}\left(\frac{h}{\delta}\right)^{2/7}.

Measured or computed velocity and density profiles may be supplied instead.
For the optional compressible approximation, the Walz relation is

.. math::

   T(y)=T_w+(T_r-T_w)\frac{U}{U_e}
        +(T_e-T_r)\left(\frac{U}{U_e}\right)^2,

.. math::

   T_r=T_e\left[1+Pr^{1/3}\frac{\gamma-1}{2}M_e^2\right],
   \qquad
   \frac{\rho(y)}{\rho_e}=\frac{T_e}{T(y)}.

Experiments with artificially thickened boundary layers and immersed circular
cylinders illustrate why case-specific profiles and drag coefficients can be
important; see :ref:`Johnson and Mitchell (1971)
<ref-johnson-mitchell-1971>` and :ref:`Stallings, Lamb, and Howell (1973)
<ref-stallings-lamb-howell-1973>`.

Symbols
-------

.. list-table:: Protrusion symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`h`
     - ``height``
     - Protrusion height
     - m
   * - :math:`b(y)`
     - ``frontal_width`` or profile input
     - Projected width at wall-normal position
     - m
   * - :math:`A_f`
     - ``frontal_area``
     - Projected frontal area
     - m²
   * - :math:`C_D`
     - ``drag_coefficient``
     - Supplied free-stream drag coefficient
     - dimensionless
   * - :math:`D`
     - ``direct_drag``
     - Estimated protrusion direct drag
     - N
   * - :math:`\eta_q`
     - ``shielding_factor``
     - Effective-to-edge dynamic-pressure ratio
     - dimensionless

Worked example
--------------

>>> from aerophysics import protrusion_drag
>>> drag = protrusion_drag(
...     1.1,
...     height=0.02,
...     frontal_width=0.01,
...     edge_velocity=100.0,
...     edge_density=1.225,
...     boundary_layer_thickness=0.1,
... )
>>> round(drag.shielding_factor, 6)
0.491073
>>> round(drag.direct_drag, 6)
0.661721

Applicability
-------------

This approximation does not solve the flow around the protrusion. It excludes
wall interference, horseshoe vortices, roughness-induced transition,
downstream skin-friction changes, multiple-element interference, and
shock/protrusion interaction. A transonic case with a single supplied drag
coefficient emits :class:`~aerophysics.exceptions.ApplicabilityWarning`.
The correction is an engineering estimate rather than a three-dimensional
flow solution.
