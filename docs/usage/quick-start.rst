PyBMRB quick start
======================

First, pull up an interactive python session and import the module:

.. code:: python

    import pybmrb

you may also import specific modules as follows

.. code:: python

    from pybmrb import Spectra, Histogram

Use the following command if you want to simulate |n15| -HSQC for entry 15060. By default legends are
disabled. If you want residue names as legend use legend='residue'
It will open the visualization on your default web browser

View BMRB entry or NMR-STAR file as spectra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose you are working with a protein called arsenate reductase and you have your data in a NMR-STAR format.
You found out that there are already two arsenate reductase  entries (17076,17077) in the BMRB. You may now easily
compare your data with BMRB as overlying |n15| - HSQC spectra using the following command

.. code:: python

    peak_list = Spectra.n15hsqc(bmrb_ids=[17076,17077],
                file_names='tests/test_data/MyData.str',
                legend='dataset')

This will open the visualization on your default web browser. When you mouseover the tool-tip will show the informatoin
about each peak. You may turn on and off the data set using legend on the right.
`Click here to view the output1 <../_static/n15hsqc_compare.html>`_

If you want the output as an image and not to open the visualization on web browser then use the following option

.. code:: python

    peak_list = Spectra.n15hsqc(bmrb_ids=[17076,17077],
                file_names='tests/test_data/MyData.str',
                legend='dataset',
                output_format='jpg',
                output_file='../docs/_images/n15hsqc_compare',
                show_visualization = False)

.. figure:: ../_images/n15hsqc_compare.jpg
    :alt: n15hsqc
    :align: center

    Comparing of local data with BMRB entries

The command will output the peak list information in the variable peak_list.

If you want to trace the chemical shift changes, use the following command

.. code:: python

    peak_list = Spectra.n15hsqc(bmrb_ids=[17076,17077],
                file_names='tests/test_data/MyData.str',
                legend='dataset',
                draw_trace = True)

`Click here to view the output2 <../_static/n15hsqc_compare2.html>`_

If you don't have your data in NMR-STAR format, then no problem. You may extract the peak list from any NMR spectra as
a csv file. You may use the csv file to compare your peak list with any BMRB entry

.. code:: python

    peak_list = Spectra.n15hsqc(bmrb_ids=[17076,17077],
                file_names='tests/test_data/MyData.str',
                legend='dataset',
                draw_trace = True)

.. figure:: ../_images/n15_peaklist.jpg
    :alt: n15hsqc
    :align: center

    Comparing of peak list with BMRB entries

Chemical shift histograms

You may easily generate chemical shift histogram of any atom or list of atoms or any residue with single command. The
same set of above parameters can be used to write output as static image

.. code:: python

    cs_data = Histogram.hist(residue='TYR', atom='CB')

.. figure:: ../_images/tyr-cb.jpg
    :alt: tyr-cb
    :align: center

    Chemical shift distribution of TYR CB

Different plot types (box, violin) are also supported. Click the figure caption for html version. When you mouseover the
box and violin plots will show all the statistical properties of the distribution

.. code:: python

    cs_data = Histogram.hist(residue='CYS', atom='CB',plot_type='box')

.. figure:: ../_images/cys-cb-box.jpg
    :alt: tyr-cb
    :align: center

    `Box plot <../_static/cys-cb-box.html>`_

.. code:: python

    cs_data = Histogram.hist(residue='CYS', atom='CB',plot_type='violin')

.. figure:: ../_images/cys-cb-violin.jpg
    :alt: tyr-cb
    :align: center

    `Violin plot <../_static/cys-cb-violin.html>`_


You may also use the wildcard

.. code:: python

    cs_data = Histogram.hist(residue='TYR', atom='H*')

.. figure:: ../_images/tyr-h.jpg
    :alt: tyr-cb
    :align: center

    Chemical shift distribution of TYR protons


Leaving out the residue will plot CB chemical shift distribution of all 20 standard amino acids

.. code:: python

    cs_data = Histogram.hist( atom='CB')

.. figure:: ../_images/cb.jpg
    :alt: tyr-cb
    :align: center

    Chemical shift distribution of CB

You may also plot 2D chemical shift correlation plot for two atoms in the same residue

.. code:: python

    cs_data = Histogram.hist2d(residue='CYS',atom1='N', atom2='CB')

.. figure:: ../_images/cys-n-cb.jpg
    :alt: tyr-cb
    :align: center

    Chemical shift correlation


.. |n15| replace:: :sup:`1` H - :sup:`15` N
.. |c13| replace:: :sup:`1` H - :sup:`13` C
.. |hh| replace:: :sup:`1` H - :sup:`1` H

