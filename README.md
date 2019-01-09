# PyBMRB

## Biological Magnetic Resonance data Bank (BMRB)
[BMRB](http://www.bmrb.wisc.edu/) collects, annotates, archives, and disseminates (worldwide in the public domain)
 the important spectral and quantitative data derived from NMR spectroscopic 
 investigations of biological macromolecules and metabolites. The goal is to empower 
 scientists in their analysis of the structure, dynamics, and chemistry of biological 
 systems and to support further development of the field of 
 biomolecular NMR spectroscopy. 
 
 BMRB maintains its data in a relational data base and also as flat files in [NMR-STAR](https://doi.org/10.1007/s10858-018-0220-3)
  format. It also provides data access via API and provide software tools to handle NMR-STAR files. 
 [PyNMRSTAR](https://github.com/uwbmrb/PyNMRSTAR) is one such tool to help read, write and parse
 NMR-STAR files. [PyBMRB](https://github.com/uwbmrb/PyBMRB) has been developed to facilitate easy data access and visualization
 of chemical shift data from BMRB. 
 
 ## Data visualization
 
PyBMRB uses modern python visualization tool [plotly](https://plot.ly/python/) for 
its visualization. Current version (version 1.1) features chemical shift histograms of 
standard amino acids from BMRB data base and N<sup>15</sup>-HSQC peak position simulation
for any given BMRB entry (or) list of BMRB entries (or) from a local assigned chemical 
shift list in NMR-STAR format. Random coil N<sup>15</sup>-HSQC using BMRB statistics has 
been implemented as a pilot project. 

This library can produce stand alone interactive visualizations as an html file or it can 
be used in a [Jupyter Notebook](https://jupyter.org/). Example notebooks can be found 
[here](https://github.com/uwbmrb/PyBMRB/tree/master/pybmrb/examples)