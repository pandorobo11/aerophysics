.. _marching-boundary-layer-model:

Marching boundary-layer model
=============================

This page describes the laminar boundary-layer solver for edge conditions that
change along a wall, including its equations, assumptions, and API mapping. For
a task-oriented example, see :doc:`../guides/boundary_layers`.

Physical background
-------------------

The flat-plate models in :doc:`flat_plate_boundary_layer` assume that the flow
outside the boundary layer does not change along the wall. On most real shapes
it does: the outer flow speeds up through an expansion and slows down through a
compression.

A boundary layer remembers what happened upstream. Accelerating flow thins it
and raises the wall shear, while decelerating flow thickens it and lowers the
wall shear. If the deceleration is strong enough, the wall shear falls to zero
and the flow **separates**. Two stations with the same local outer flow can
therefore have different skin friction, because their upstream histories differ.

This solver follows that history. It starts at a sharp leading edge and marches
downstream, using each solved station as the upstream state of the next one.

The outer flow is an input. The solver does not calculate it, does not represent
a body shape, and does not feed the boundary layer back into the outer flow.

Governing equations
-------------------

Distance :math:`s` is measured along the wall from the leading edge. The thin
boundary-layer equations are written in Levy--Lees variables
(:ref:`Levy, 1954 <ref-levy-1954>`; :ref:`Lees, 1956 <ref-lees-1956>`):

.. math::

   \xi(s) = \int_0^s \rho_e \mu_e U_e \, ds, \qquad
   \eta = \frac{U_e}{\sqrt{2\xi}} \int_0^y \rho \, dy, \qquad
   f' = \frac{u}{U_e}, \qquad g = \frac{T}{T_e}

With :math:`C = \rho\mu/(\rho_e\mu_e)`, momentum and energy become

.. math::

   \left(C f''\right)' + f f'' + \beta\left(g - f'^2\right)
   = 2\xi\left(f' f'_\xi - f'' f_\xi\right)

.. math::

   \left(\frac{C}{Pr} g'\right)' + f g' + (\gamma - 1) M_e^2 C f''^2
   = 2\xi\left(f' g_\xi - g' f_\xi\right)

where primes denote :math:`\partial/\partial\eta` and

.. math::

   \beta = 2\xi \frac{d \ln U_e}{d\xi}

is positive for accelerating flow. The energy equation has no explicit
:math:`\beta` term. For an isentropic edge the pressure work cancels the change
in edge enthalpy, in the same way that a fluid particle expanding isentropically
keeps a constant :math:`T/T_e`.

At the wall, :math:`f = f' = 0`, and :math:`g` is either the specified wall
temperature ratio or follows from a zero wall heat flux. At the outer boundary,
:math:`f' = g = 1`.

Results
~~~~~~~

.. math::

   C_f = \frac{\tau_w}{\tfrac12 \rho_e U_e^2}
       = \frac{2 C_w \mu_e f''_w}{\sqrt{2\xi}}, \qquad
   q_w = \frac{C_w \rho_e \mu_e U_e c_p T_e}{Pr \sqrt{2\xi}} \, g'_w

Heat flux into the wall is positive. :math:`q_w` comes from the computed wall
temperature gradient, not from the Reynolds analogy that
``flat_plate_boundary_layer`` uses. With :math:`Y = \int_0^\eta g \, d\eta`,

.. math::

   y = \frac{\sqrt{2\xi}}{\rho_e U_e} Y, \qquad
   \delta^* = \frac{\sqrt{2\xi}}{\rho_e U_e} \int_0^{\eta_{max}} (g - f') \, d\eta,
   \qquad
   \theta = \frac{\sqrt{2\xi}}{\rho_e U_e} \int_0^{\eta_{max}} f' (1 - f') \, d\eta

:math:`\delta_{99}` is the first wall distance where :math:`u/U_e = 0.99`.

Edge state
----------

The equations require :math:`dp_e = -\rho_e U_e \, dU_e` and a constant total
enthalpy along the edge, so pressure, temperature, and velocity are not
independent. ``BoundaryLayerEdge.from_pressure`` takes the edge pressure and the
edge total pressure and temperature, then calculates the rest from isentropic
relations. Downstream of a shock, the total pressure is the post-shock value.

Between the supplied stations, :math:`\ln p_e` follows a monotone piecewise cubic
(PCHIP) interpolation. Supply a smooth, well-resolved pressure distribution:
:math:`\beta` depends on the pressure gradient, so noise in the pressure appears
directly in the solution.

Assumptions and range
---------------------

* Steady, planar, laminar, attached flow from a sharp leading edge.
* Calorically perfect gas with constant :math:`Pr`, and pressure constant across
  the layer.
* Isentropic edge flow with constant total enthalpy. Entropy layers from blunt
  noses are not represented.
* No transition, roughness, suction, blowing, wall curvature, or chemical
  reactions, and no interaction with the outer flow.

``BoundaryLayerRegime.TURBULENT`` and ``TRANSITIONAL`` raise
``ModelRangeError``.

.. list-table:: Reported conditions
   :header-rows: 1
   :widths: 40 60

   * - Condition
     - Report
   * - Invalid input or an edge state that is not isentropic
     - ``ValueError``
   * - The wall shear collapses in an adverse pressure gradient
     - ``BoundaryLayerSeparationError``, holding ``separation_distance`` and the
       ``result`` at the stations upstream. The march cannot continue because
       the thin-layer equations are not valid past separation.
   * - The solution cannot be resolved
     - ``BoundaryLayerConvergenceError``
   * - :math:`\delta_{99}\,|d \ln p_e/ds| > 1`
     - ``ApplicabilityWarning``. The pressure changes too quickly for the
       thin-layer approximation, as near a sharp corner, where viscous
       interaction with the outer flow matters.
   * - :math:`\bar\chi = M_e^3 \sqrt{C_w/Re_s} > 1`
     - ``ApplicabilityWarning``. Leading-edge viscous interaction changes the
       pressure, which the supplied edge state may not include.

Verification
------------

:doc:`../verification/viscous_flow` records the comparisons. They use the
Falkner--Skan wedge flows (:ref:`Falkner and Skan, 1931
<ref-falkner-skan-1931>`), similar solutions with wall heat transfer in the form
of :ref:`Cohen and Reshotko, 1956 <ref-cohen-reshotko-1956>`, Howarth's linearly
retarded flow (:ref:`Howarth, 1938 <ref-howarth-1938>`), and its compressible
image under the transformation of :ref:`Illingworth, 1949
<ref-illingworth-1949>` and :ref:`Stewartson, 1949 <ref-stewartson-1949>`.

API mapping
-----------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Symbol
     - Field of ``MarchingBoundaryLayerResult``
   * - :math:`s`
     - ``distance``
   * - :math:`\rho_e U_e s/\mu_e`
     - ``reynolds_number``
   * - :math:`\delta_{99},\ \delta^*,\ \theta`
     - ``boundary_layer_thickness``, ``displacement_thickness``,
       ``momentum_thickness``
   * - :math:`C_f,\ \tau_w`
     - ``local_skin_friction_coefficient``, ``wall_shear_stress``
   * - :math:`T_w,\ q_w`
     - ``wall_temperature``, ``wall_heat_flux``
   * - :math:`\beta`
     - ``pressure_gradient_parameter``
   * - :math:`\eta,\ y,\ u/U_e,\ T/T_e`
     - ``eta``, ``wall_distance``, ``velocity_ratio``, ``temperature_ratio``

The leading edge is excluded from the result, because :math:`\tau_w` and
:math:`q_w` are singular there.
