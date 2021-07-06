Welcome to PyBMRB!
======================================

A Python module for visualizing Nuclear Magnetic Resonance(NMR)  chemical shift data from Biological Magnetic
Resonance data Bank (`BMRB <http://bmrb.ip>`_) and  from NMR-STAR :footcite:`Ulrich2019` format files. PyBMRB helps the
users to view the one dimensional chemical shift list as multi-dimensional NMR spectra. Chemical shift distributions
can be studied using this library by plotting the histograms using differnt filtering criteria.

This package uses `PyBMRB <https://github.com/uwbmrb/PyNMRSTAR>`_ to parse the NMR-STAR files
and `BMRB-API <https://github.com/uwbmrb/BMRB-API>`_ to fetch the data directly from BMRB, which avoids the hustle of
downloading and parsing the data from the BMRB for visualizations.



|BuildStatus| |License| |Wheel| |PythonVersions|