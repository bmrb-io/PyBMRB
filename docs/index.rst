Welcome to PyBMRB!
======================================

A Python module for visualizing Nuclear Magnetic Resonance(NMR)  chemical shift data from Biological Magnetic
Resonance data Bank(`BMRB <http://bmrb.ip>`_) and  from NMR-STAR :footcite:`Ulrich2019` format files

|BuildStatus| |License| |Wheel| |PythonVersions|

Python versions supported: 3.6, 3.7, 3.8, and 3.9

Previous python versions (back to 2.6) are supported by the v2 branch
(version 2.x releases on PyPI).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage/nmrstar-introduction
   usage/quick-start
   usage/examples
   usage/full
   release-notes
Module documentation
======================================

Spectra class
~~~~~~~~~~~

.. autoclass:: pybmrb2.Spectra
   :special-members:
   :members:


.. autoclass:: pybmrb2.ChemicalShift
   :special-members:
   :members:

.. |PythonVersions| image:: https://img.shields.io/pypi/pyversions/pynmrstar.svg
   :target: https://pypi.org/project/PyBMRB

.. |License| image::  https://img.shields.io/github/license/kumar-physics/PyBMRB
   :target: https://pypi.org/project/PyBMRB

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pynmrstar.svg
   :target: https://pypi.org/project/PyNMRSTAR


.. |BuildStatus| image:: https://img.shields.io/github/workflow/status/kumar-physics/PyBMRB/CI/dev
   :target: https://pypi.org/project/PyNMRSTAR

.. footbibliography::