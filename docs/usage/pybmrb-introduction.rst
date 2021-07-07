Introduction to PyBMRB
------------------------

Biological Magnetic Resonance data Bank (BMRB)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`BMRB <http://bmrb.ip>`_ is the global archive of NMR spectroscopic derived from biological
molecules like proteins, nucleic acids and metabolites. BMRB collects, annotates, archives,
and disseminates (worldwide in the public domain) the important spectral and quantitative data
derived from NMR spectroscopic investigations of biologically relevant molecules. The goal is
to empower scientists in their analysis of the structure, dynamics, and chemistry of
biological systems and to support further development of the field of biomolecular
NMR spectroscopy.

NMR-STAR
~~~~~~~~~

BMRB uses NMR-STAR :footcite:`Ulrich2019` data
model, which is an object oriented data model driven by
`NMR-STAR dictionary <https://github.com/uwbmrb/nmr-star-dictionary>`_ . Each BMRB entry includes
the meta data such as the information about the biological sample and solvent, experimental conditions,
instrument details, author information and so on along the observed NMR chemical shift data. NMR chemical
shift  is the core information available in NMR-STAR formatted files provided by BMRB. These files may also
contain additional derived data such as restraints, Residual Dipolar couplings (RDC), J-couplings, Chemical Shfit
Aniostrophy(CSA) data and so on. BMRB provides python parser (`PyNMRSTAR <https://github.com/uwbmrb/PyNMRSTAR>`_ )
to read, write and edit NMR-STAR files. BMRB data is also available through `BMRB-API <https://github.com/uwbmrb/BMRB-API>`_
for machine to machine communication and programmatic access.

Why do we need PyBMRB?
~~~~~~~~~~~~~~~~~~~~~~~~
BMRB is mostly used by  biologists and biochemists who have very little or no programing experience.Though NMR-STAR files are
resourceful for the research community, one needs to understand the NMR-STAR data model and to know little bit of programing
to extract the desired information from NMR-STAR files. NMR spectroscopists may want to view the data as NMR spectra and
compare with their spectra measured using their sample. PyBMRB can do this in a single command, aviods the hustle
of downloading and parsing the data for visualization. This would greatly benefit the research community by allowing them
to quickly and easily compare their data with any BMRB entry and  visualizing the BMRB data in different types of 2D spectra.

How it works?
~~~~~~~~~~~~~~
PyBMRB extracts the assigned chemical shift list from NMR-STAR files and combines them with certain rules
defined by the NMR experiment to generate the peak positions. This peak list is displayed on a 2D plane using
interactive data visualization tool `plotly <https://plotly.com>`_ . It can also generate chemical shift histograms
by fetching the data directly from BMRB through BMRB-API.

.. footbibliography::