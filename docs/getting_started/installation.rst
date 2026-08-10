Installation
============

``aerophysics`` requires Python 3.12 or newer. Install the calculation package
with pip:

.. code-block:: console

   python -m pip install aerophysics

Confirm that the package can evaluate the reference atmosphere:

>>> from aerophysics import standard_atmosphere
>>> standard_atmosphere(0.0).temperature
288.15

Optional local GUI
------------------

Install the GUI extra and launch the local Streamlit application:

.. code-block:: console

   python -m pip install "aerophysics[gui]"
   aerophysics-gui

The wheel includes an offline copy of this manual. See :doc:`../guides/gui`
for display-unit handling, case handoffs, CSV and JSON files, and the embedded
documentation page.

Development installation
------------------------

Repository contributors should use the locked environment described in the
top-level ``DEVELOPMENT.md`` rather than installing editable dependencies by
hand.
