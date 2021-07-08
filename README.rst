Welcome to PyBMRB!
======================================

A Python module for visualizing Nuclear Magnetic Resonance(NMR)  chemical shift data from Biological Magnetic
Resonance data Bank (`BMRB <http://bmrb.ip>`_) and  from NMR-STAR  format files. PyBMRB helps the
users to view the one dimensional chemical shift list as multi-dimensional NMR spectra. Chemical shift distributions
can be studied using this library by plotting the histograms using differnt filtering criteria.

This package uses `PyBMRB <https://github.com/uwbmrb/PyNMRSTAR>`_ to parse the NMR-STAR files
and `BMRB-API <https://github.com/uwbmrb/BMRB-API>`_ to fetch the data directly from BMRB, which avoids the hustle of
downloading and parsing the data from the BMRB for visualizations.

`Read the documents here <https://pybmrb.readthedocs.org>`__.


.. |PythonVersions| image:: https://img.shields.io/pypi/pyversions/pynmrstar.svg
   :target: https://github.com/uwbmrb/PyBMRB

.. |License| image::  https://img.shields.io/github/license/kumar-physics/PyBMRB
   :target: https://github.com/uwbmrb/PyBMRB


.. |BuildStatus| image:: https://img.shields.io/github/workflow/status/kumar-physics/PyBMRB/CI/dev
   :target: https://github.com/uwbmrb/PyBMRB


|BuildStatus| |License| |PythonVersions|