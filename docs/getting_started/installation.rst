Installation
============

``aerophysics`` requires Python 3.12 or newer. Until a PyPI distribution is
available, install the current wheel from the `latest GitHub Release`_:

.. code-block:: console

   python -m pip install "https://github.com/pandorobo11/aerophysics/releases/download/v0.8.0/aerophysics-0.8.0-py3-none-any.whl"

Check the release page for a newer version before copying a versioned URL.

Check that the installation works by asking for the sea-level temperature of
the standard atmosphere:

>>> from aerophysics import standard_atmosphere
>>> standard_atmosphere(0.0).temperature
288.15

Optional browser interface
--------------------------

The GUI is an optional extra. Install it and start it with:

.. code-block:: console

   python -m pip install "aerophysics[gui] @ https://github.com/pandorobo11/aerophysics/releases/download/v0.8.0/aerophysics-0.8.0-py3-none-any.whl"
   aerophysics-gui

The GUI opens in your web browser and runs entirely on your own computer. The
package includes a copy of this manual, which the GUI shows on its
**ドキュメント** (documentation) page. See :doc:`../gui/index`.

Development installation
------------------------

Repository contributors should use the locked environment described in the
top-level ``DEVELOPMENT.md`` rather than installing editable dependencies by
hand.

.. _latest GitHub Release: https://github.com/pandorobo11/aerophysics/releases/latest
