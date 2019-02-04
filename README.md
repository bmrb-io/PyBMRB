# PyBMRB
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/uwbmrb/PyBMRB/master?filepath=pybmrb%2Fexamples) Click the launcher to try PyBMRB on a live Jupyter Notebook.
## Biological Magnetic Resonance data Bank (BMRB)
[BMRB](http://www.bmrb.wisc.edu/) collects, annotates, archives, and disseminates (worldwide in the public domain)
 the important spectral and quantitative data derived from NMR spectroscopic 
 investigations of biological macromolecules and metabolites. The goal is to empower 
 scientists in their analysis of the structure, dynamics, and chemistry of biological 
 systems and to support further development of the field of 
 biomolecular NMR spectroscopy. 
 
 BMRB maintains its data in a relational data base and also as flat files in [NMR-STAR](https://doi.org/10.1007/s10858-018-0220-3)
  format. It also provides data access via API and provide software tools to handle NMR-STAR files. 
 [PyNMRSTAR](https://github.com/uwbmrb/PyNMRSTAR) is one such tool to facilitate read, write and parse
 NMR-STAR files. [PyBMRB](https://github.com/uwbmrb/PyBMRB) enhances the usability of BMRB data 
 and data access and visualization of chemical shift data from BMRB. 
 
 ## Data visualization
 
PyBMRB uses modern python visualization tool [plotly](https://plot.ly/python/) for 
its visualization. It features chemical shift histograms of 
standard amino acids from BMRB data base and N<sup>15</sup>-HSQC peak position simulation
for any given BMRB entry (or) list of BMRB entries (or) from a local assigned chemical 
shift list in NMR-STAR format. Random coil N<sup>15</sup>-HSQC simulation based on BMRB 
database statistics is also included.  

This library can produce stand alone interactive visualizations as an html file or it can 
be used in a [Jupyter Notebook](https://jupyter.org/). Example notebooks can be found 
[here](https://github.com/uwbmrb/PyBMRB/tree/master/pybmrb/examples)

### Installation
PyBMRB is available in [Python Package Index (PIP)](https://pypi.org/project/pybmrb/). It can be easily 
installed with the following command
>pip install pybmrb

To install from the source, clone or download the source and go to the folder which contains 'setup.py' 
and run the following command
>pip install .

If necessary use 'sudo'. 
### Chemical shift Histogram
BMRB serves Nuclear Magnetic Resonance(NMR) community by providing 
high quality curated chemical shift data of various biologically important 
macro molecules like proteins and nucleic acids and small molecules like ligands , co-factors,
small peptides and metabolites. Chemical shift distribution of an atom in an amino acid or 
nucleic acid helps us to understand its diverse chemical environment.
Chemical shift histograms from a database like BMRB will help
us to understand the biophysical nature bio-molecules and provides a priory 
probabilities for resonance assignments.

#### Conditional Histogram
Chemical shift assignment in a multi-dimensional NMR experiment is 
a biggest challenge. PyBMRB could help the process by generating 
conditional histogram. If partial assignments of a particular 
amnio acid is available, then PyBMRB could estimate the chemical shift
distribution of an unassigned atom in that amino acid only from those entries
 whose chemical shifts are similar to the given list for that 
 particular amino acid. 


#### Examples in Notebook
Chemical shift histogram [examples notebook](https://github.com/uwbmrb/PyBMRB/tree/master/pybmrb/examples/histogram.ipynb) 
can be found [here](https://github.com/uwbmrb/PyBMRB/tree/master/pybmrb/examples/histogram.ipynb). 
The output might look like [this](https://htmlpreview.github.io/?https://github.com/uwbmrb/PyBMRB/blob/master/pybmrb/examples/histogram.html)  

### N<sup>15</sup>-HSQC Simulation.

N<sup>15</sup>-HSQC is one of the most common and useful NMR experiment.
It carries signatures of folded state, ligand interaction and experimental environments 
like pH and temperature. BMRB contains thousands of assigned chemical shifts
of various bio-molecules. PyBMRB make use of this data to simulate N<sup>15</sup>-HSQC 
peak positions for any BMRB entry. Its also possible to compare N<sup>15</sup>-HSQC 
spectra of different entries or BMRB entries with user generated
 chemical shift information in NMR-STAR format. 
 
 #### N<sup>15</sup>-HSQC spectrum from sequence. 
 As an experiment, BMRB developed tools to simulate random coil
 N<sup>15</sup>-HSQC spectrum for a given protein sequence using 
 database statistics. Assuming the fact that in BMRB, small segments 
 of given sequence might exists in all possible confirmations, 
 one could calculate the chemical shift of an atom in random 
 confirmation by averaging the chemical shifts of small segments
 over the entire database
 
 #### Examples in Notebook
 N<sup>15</sup>-HSQC [examples notebook](https://github.com/uwbmrb/PyBMRB/tree/master/pybmrb/examples/n15hsqc.ipynb) 
 can be found [here](https://github.com/uwbmrb/PyBMRB/tree/master/pybmrb/examples/n15hsqc.ipynb)
 The output of the notebook looks like [this](https://htmlpreview.github.io/?https://github.com/uwbmrb/PyBMRB/blob/master/pybmrb/examples/n15hsqc.html)
 
