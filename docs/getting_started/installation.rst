Installation
============

``aerophysics`` requires Python 3.12 or newer. Until a PyPI distribution is
available, install the current wheel from the `latest GitHub Release`_:

.. code-block:: console

   python -m pip install "https://github.com/pandorobo11/aerophysics/releases/download/v0.6.0/aerophysics-0.6.0-py3-none-any.whl"

Check the release page for a newer version before copying a versioned URL.

Confirm that the package can evaluate the reference atmosphere:

>>> from aerophysics import standard_atmosphere
>>> standard_atmosphere(0.0).temperature
288.15

Optional local GUI
------------------

Install the GUI extra and launch the local Streamlit application:

.. code-block:: console

   python -m pip install "aerophysics[gui] @ https://github.com/pandorobo11/aerophysics/releases/download/v0.6.0/aerophysics-0.6.0-py3-none-any.whl"
   aerophysics-gui

The wheel includes an offline copy of this manual. See :doc:`../guides/gui`
for display-unit handling, case handoffs, CSV and JSON files, and the embedded
documentation page.

Development installation
------------------------

Repository contributors should use the locked environment described in the
top-level ``DEVELOPMENT.md`` rather than installing editable dependencies by
hand.

.. _latest GitHub Release: https://github.com/pandorobo11/aerophysics/releases/latest
