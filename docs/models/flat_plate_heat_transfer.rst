Flat-plate heat transfer and profiles
=====================================

Wall heat transfer
------------------

With ``compressible=True`` the result also holds the local Stanton number,
heat-transfer coefficient and wall heat flux, from the Reynolds analogy on the
local edge-basis friction:

.. math::

   St=s\,\frac{C_{f,x}}{2},
   \qquad
   h=St\,\rho_eU_ec_p,
   \qquad
   q_w=h\,(T_r-T_w).

:math:`q_w` is in W/m², :math:`h` in W/(m²·K), and :math:`q_w>0` means heat
flows into the wall. :math:`T_r` is the returned ``recovery_temperature``:
the solved adiabatic wall temperature at laminar points and
:math:`r_t=Pr^{1/3}` at turbulent points. An adiabatic wall has :math:`q_w=0`
and a defined :math:`St`. These fields are ``None`` for ``compressible=False``.

**Laminar points** use :math:`s=Pr^{-2/3}`. Against the wall flux of the
thin-layer equations with :math:`\tilde\nu=0` (see the
:ref:`SA reference below <sa-reference-heat-transfer>`), the difference is
within 1% for :math:`M_e=0.3`--10 and :math:`T_w/T_{aw}=0.2`--1.2.

**Turbulent points** use ``reynolds_analogy_factor``, with default
:math:`s=1.14` from the :ref:`generalized Reynolds analogy
<ref-zhang-bi-hussain-she-2014>`, so that the default temperature profile has
a matching wall gradient. Other common values are 1.16 from :ref:`Chi and
Spalding <ref-chi-spalding-1966>` and Colburn's :math:`Pr^{-2/3}\approx1.24`.
:ref:`Cary <ref-cary-1970>` finds the 1.16 fit a valid mean for :math:`M_e<4`--5
near adiabatic conditions, finds too little data and too much scatter above
that to define :math:`s`, and notes that it appears to fall as the wall is
cooled. In a Mach 6 cold-wall DNS, :ref:`Fu et al. <ref-fu-et-al-2021>` find
values within 5% of 1.16. Treat :math:`s` as an uncertainty of several percent
rather than a constant of the model.

Compared with a boundary-layer solution using the Spalart--Allmaras model,
the heat flux of ``VAN_DRIEST_II_WALL_TEMPERATURE`` with :math:`s=1.14`
differs by about ±5% at :math:`Re_{x,e}\ge10^7` on the cooled walls tested
(:math:`M_e=5`--10, :math:`T_w/T_r=0.20`--0.76). The difference is larger at
lower Reynolds number and when the wall temperature is close to the recovery
temperature, where :math:`T_r-T_w` is sensitive to the recovery factor.
``ECKERT`` and ``VAN_DRIEST_II`` also return heat-transfer values, but these
have not been compared in the same way. All of this assumes constant
:math:`c_p` and :math:`Pr`, an isothermal wall, and zero pressure gradient.
Details are in :doc:`../verification/flat_plate_boundary_layer`.

Velocity and temperature profiles
---------------------------------

``flat_plate_boundary_layer_profile`` returns the velocity, temperature,
density, viscosity and local Mach number at one station as functions of
:math:`y/\delta_{99}`, together with the ``flat_plate_boundary_layer`` result.
The edge state must be physically consistent, because the edge Mach number and
viscosity follow from :math:`U_e`, :math:`T_e`, ``gas`` and ``viscosity_model``.

Laminar stations use the compressible similarity solution, so the wall
gradients of the profile reproduce :math:`\tau_w` and, within 1%, the analogy
:math:`q_w`. Turbulent stations pass the ``VAN_DRIEST_II_WALL_TEMPERATURE``
:math:`\tau_w` and :math:`\delta_{99}` to
``compressible_turbulent_boundary_layer_profile``; see :doc:`Compressible
boundary-layer profiles <../api/viscous_flow>` for that model. The thickness
closure uses the Walz relation, while the profile defaults to the generalized
Reynolds analogy so that its wall gradient matches :math:`s`. Measured locally
against the SA reference on the same :math:`y/\delta_{99}` coordinate, the
largest differences at :math:`Re_{x,e}\ge10^7` are 4--8% of :math:`U_e` in
velocity and 1--7% of the peak temperature. At :math:`Re_{x,e}=10^6` they are
7--9% and 5--9%.
