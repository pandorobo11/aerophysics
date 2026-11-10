.. _flat-plate-specified-transition:

Specified transition on a flat plate
====================================

With ``BoundaryLayerRegime.TRANSITIONAL`` the boundary layer is laminar up to
the Reynolds number :math:`Re_\mathrm{tr}` that you supply, and turbulent
after it. This is a transition position that you choose, not a prediction.

The equivalent turbulent plate
------------------------------

Downstream of transition, the layer continues as a fully turbulent plate
that has the same momentum thickness as the laminar layer at the transition
point. For constant edge and wall conditions and zero pressure gradient,
the momentum integral gives :math:`\theta_m=x\bar C_f/2` for both the laminar
and the fully turbulent plate, so the turbulent start :math:`Re_0` solves

.. math::

   Re_0\,\bar C_{f,t}(Re_0)=Re_\mathrm{tr}\,\bar C_{f,l}(Re_\mathrm{tr})
   =2Re_{\theta,\mathrm{tr}}.

The equivalent Reynolds number and length then give

.. math::

   Re_\mathrm{eq}=Re_0+Re_x-Re_\mathrm{tr},
   \qquad
   x_\mathrm{eq}=x\frac{Re_\mathrm{eq}}{Re_x},

   C_{f,x}=C_{f,t}(Re_\mathrm{eq}),
   \qquad
   \bar C_f=\frac{Re_\mathrm{eq}}{Re_x}\bar C_{f,t}(Re_\mathrm{eq}).

Here :math:`C_{f,t}` and :math:`\bar C_{f,t}` are the fully turbulent
edge-basis coefficients of the selected model, and :math:`\bar C_{f,l}` is
Blasius, or the compressible similarity solution with ``compressible=True``.
The compressible models write these coefficients as
:math:`s\,C_{f,i}(kRe_x)`, with :math:`(s,k)=(\rho^*/\rho_e,Re^*/Re_x)` for
``ECKERT`` and :math:`(w/F_C,F_x)` for the Van Driest II models (:math:`w=1`
for ``VAN_DRIEST_II``), so the root is found in incompressible variables:

.. math::

   kRe_0\,\bar C_{f,i}(kRe_0)=2Re_{\theta,\mathrm{tr}}\frac{k}{s}.

For Van Driest II with :math:`w=1`, :math:`k/s=F_\theta`, the
momentum-thickness Reynolds number transformation. Thickness formulas,
including the ``SCHLICHTING`` power law, use :math:`x_\mathrm{eq}` and
:math:`Re_\mathrm{eq}`, and ``effective_reynolds_number`` returns
:math:`kRe_\mathrm{eq}` at turbulent points (:math:`k=1` for
``compressible=False``). As an example, :math:`Re_\mathrm{tr}=5\times10^5`
gives :math:`Re_0=1.417\times10^5` with Schlichting friction.

What is continuous at transition
--------------------------------

Accumulated drag, and with it the drag-based momentum thickness
:math:`x\bar C_f/2`, is continuous at transition for every model. Whether the
returned quantities satisfy stronger relations depends on the model:

.. list-table:: Exact and approximate relations after a specified transition
   :header-rows: 1

   * - Model
     - Returned :math:`\theta_m` at transition
     - :math:`\mathrm{d}(x\bar C_f)/\mathrm{d}x=C_{f,x}`
   * - ``SCHLICHTING``
     - Power law, 1.3--2.2% above the laminar value for
       :math:`Re_\mathrm{tr}=5\times10^5`--:math:`5\times10^6`
     - Exact
   * - ``SCHLICHTING_LOG_WAKE``
     - :math:`x\bar C_f/2`, continuous
     - Exact
   * - ``ECKERT``
     - ``NaN``
     - Exact
   * - ``VAN_DRIEST_II``
     - ``NaN``
     - Approximate: the Willems local and average equations are separate
       fits, 3.9% apart at :math:`R=10^6`
   * - ``VAN_DRIEST_II_WALL_TEMPERATURE``
     - :math:`x\bar C_f/2`, continuous
     - Exact

The ``SCHLICHTING`` and ``VAN_DRIEST_II`` departures belong to their fully
turbulent correlations and are unchanged by the equivalent origin. Local
friction, :math:`\delta_{99}`, and :math:`\delta^*` jump, because transition is
instantaneous here. Turbulent friction just after transition is consistent
with the drag-based momentum thickness, and the ``SCHLICHTING_LOG_WAKE`` and
``VAN_DRIEST_II_WALL_TEMPERATURE`` closures see a fully turbulent equilibrium
state.

With the ``SCHLICHTING`` power law, the returned :math:`\delta_{99}` also
jumps at transition, because the assumed profile changes from Blasius to the
one-seventh-power shape.

Limits
------

The equivalent plate has not been compared with measurements of transitional
flow. A real transition region has a finite length, and the turbulence just
after it is not yet in equilibrium. If the laminar momentum thickness at
transition is smaller than any turbulent plate can have (about
:math:`Re_\mathrm{tr}\lesssim11` with Blasius and Schlichting), the call
raises :class:`~aerophysics.exceptions.ModelRangeError`.
