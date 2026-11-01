Shock waves
===========

This page explains what a shock wave is and the conservation laws behind all
shock models. The equations for each kind of shock are on the following
pages. State ratios are downstream over upstream, and angles are in radians.
Convert explicitly with :func:`aerophysics.units.degrees_to_radians` and
:func:`aerophysics.units.radians_to_degrees`. For task-oriented examples, see
:doc:`../guides/compressible_flow`.

Physical background
-------------------

A shock wave is a very thin region in which a supersonic flow slows down
abruptly. Pressure, temperature, and density rise across it almost
instantly. Viscosity and heat conduction act only inside this thin layer, so
the flow on either side can be treated as inviscid. Because the change is
sudden and irreversible, entropy increases and some of the flow's ability to
do work is lost. That loss appears as a drop in **total pressure**. No heat
is added or removed, so **total temperature** (more generally, total
enthalpy) stays the same.

Apply conservation of mass, momentum, and energy to a thin control volume
that encloses the shock. With :math:`u` the velocity component normal to the
shock and :math:`h` the specific enthalpy,

.. math::

   \rho_1 u_1=\rho_2 u_2,\qquad
   p_1+\rho_1u_1^2=p_2+\rho_2u_2^2,\qquad
   h_1+\frac{u_1^2}{2}=h_2+\frac{u_2^2}{2}.

These are the Rankine--Hugoniot conditions. Together with an equation of
state, they determine the downstream state from the upstream state. The
second law of thermodynamics allows only the compressive solution, so a
shock always needs a supersonic normal Mach number upstream and always
produces a subsonic normal Mach number downstream.

* A **normal shock** stands perpendicular to the flow. An **oblique shock**
  forms when a wall turns into a supersonic stream; only the velocity
  component normal to the shock changes, so the normal-shock relations apply
  to that component. See :doc:`normal_oblique_shocks`.
* A **conical shock** forms on a sharp cone. The flow keeps compressing
  smoothly between the shock and the cone surface, so it is weaker than the
  oblique shock on a wedge with the same half-angle. See
  :doc:`conical_shocks`.
* A **detached (bow) shock** stands in front of a blunt body, or in front of
  a wedge or cone that turns the flow more than an attached shock allows. See
  :doc:`detached_shocks`.

.. toctree::
   :maxdepth: 1

   normal_oblique_shocks
   conical_shocks
   detached_shocks
