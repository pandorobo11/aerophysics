Verification
============

These records provide reproducible checks of equations, reference tables,
cross-implementation results, and physical or mathematical invariants. A
verified result establishes agreement with the stated acceptance evidence.
Non-gating physical-accuracy assessments are reported separately and do not
change that status.

Current status
--------------

.. list-table:: Verification records
   :header-rows: 1
   :widths: 25 22 31 22

   * - Domain
     - Status
     - Principal evidence
     - Record
   * - U.S. Standard Atmosphere 1976
     - Verified with observations
     - Official tables, ``fluids`` snapshot, dense-grid invariants
     - :doc:`standard_atmosphere`
   * - Compressible flow
     - Verified
     - NACA Report 1135, NASA SP-3004, source equations and invariants
     - :doc:`compressible_flow`
   * - Thermophysical properties
     - Verified
     - USSA values, published equations, Cantera, NIST assessment, invariants
     - :doc:`thermophysical`
   * - Viscous flow
     - Verified
     - Published correlations, NASA TN D-6945, profile and integration checks
     - :doc:`viscous_flow`

``Verified with observations`` means every acceptance criterion passes but a
stricter diagnostic or non-acceptance comparison is retained for review.

Reproduction
------------

Regenerate all committed verification fragments and figures with:

.. code-block:: console

   uv run python docs/scripts/generate_verification.py

Check them without changing files with:

.. code-block:: console

   uv run python docs/scripts/generate_verification.py --check

Each record documents its source-specific snapshots, tolerances, limitations,
and individual regeneration commands.

.. toctree::
   :maxdepth: 1

   standard_atmosphere
   compressible_flow
   thermophysical
   viscous_flow
