Examples
========

Load the package into your python environment using the following command

.. code:: python

    import pybmrb

.. |n15| replace:: :sup:`1` H - :sup:`15` N
.. |c13| replace:: :sup:`1` H - :sup:`13` C
.. |hh| replace:: :sup:`1` H - :sup:`1` H


Spectra simulation
--------------------
PyBMRB combines the chemical shift information from BMRB entries or NMR-STAR files using certain rules defined by the
experiment type. In this method it can generate |n15| - HSQC, |c13| -HSQC and |HH|-TOCSY. You may define which atom to
be on X axis and which atom to be on Y axis to generate through bond correlation 2D spectrum.


|n15| - HSQC peak position simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Example 1: Single entry from BMRB

.. code:: python

    pybmrb.Spectra.n15hsqc(bmrb_ids=15060, legend='residue')

.. figure:: ../_images/15060_n15.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for n15-HSQC from BMRB entry 15060 <../_static/15060_n15.html>`_

- Example 2: Multiple entries from BMRB along with a NMR-STAR file.

For multiple data set use data set as legend, so that you may turn on and off different data set.
You may also use residue as legend to turn on and off different residue types

.. code:: python

    pybmrb.Spectra.n15hsqc(bmrb_ids=[17076,17077],
    input_file_names='test_data/MyData.str',
    legend='dataset')

.. figure:: ../_images/multi_n15.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for n15-HSQC from BMRB entries along with NMR-STAR file <../_static/multi_n15.html>`_

- Example 3: Multiple entries from BMRB along with a NMR-STAR file and a peak list in csv format

.. code:: python

    pybmrb.Spectra.n15hsqc(bmrb_ids=[17076,17077],
    input_file_names='test_data/MyData.str',
    peak_list='test_data/my_peak_list.csv',
    legend='dataset')

.. figure:: ../_images/multi_n152.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for n15-HSQC from BMRB entries along with NMR-STAR file and a peak list <../_static/multi_n152.html>`_


- Example 4: Multiple entries from BMRB with chemical shift tracking

.. code:: python

    pybmrb.Spectra.n15hsqc(bmrb_ids=[17076,17077],
    input_file_names='test_data/MyData.str',
    legend='dataset',
    draw_trace=True)

.. figure:: ../_images/multi2_n15.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for n15-HSQC from BMRB entries 17074,17076 and 17076 with trace <../_static/multi2_n15.html>`_


|c13| - HSQC peak position simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Example 5: Single entry from BMRB

.. code:: python

    pybmrb.Spectra.c13hsqc(bmrb_ids=15060, legend='residue')

.. figure:: ../_images/15060_c13.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for c13-HSQC from BMRB entry 15060 <../_static/15060_c13.html>`_

- Example 6: Multiple entries from BMRB

.. code:: python

    pybmrb.Spectra.c13hsqc(bmrb_ids=[17074,17076,17077], legend='dataset')

.. figure:: ../_images/multi_c13.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for c13-HSQC from BMRB entries 17074,17076 and 17076 <../_static/multi_c13.html>`_

- Example 7: Multiple entries from BMRB with chemical shift tracking

.. code:: python

    pybmrb.Spectra.c13hsqc(bmrb_ids=[17074,17076,17077], legend='dataset', draw_trace=True)

.. figure:: ../_images/multi2_c13.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for c13-HSQC from BMRB entries 17074,17076 and 17076 with trace <../_static/multi2_c13.html>`_


|hh| - TOCSY peak position simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Example 8: Single entry from BMRB

.. code:: python

    pybmrb.Spectra.tocsy(bmrb_ids=15060, legend='residue')

.. figure:: ../_images/15060_tocsy.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for TOCSY from BMRB entry 15060 <../_static/15060_tocsy.html>`_

- Example 9: Multiple entries from BMRB

.. code:: python

    pybmrb.Spectra.tocsy(bmrb_ids=[17074,17076,17077], legend='dataset')

.. figure:: ../_images/multi_tocsy.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for TOCSY from BMRB entries 17074,17076 and 17076 <../_static/multi_tocsy.html>`_

- Example 10: Multiple entries from BMRB with residues as legend

.. code:: python

    pybmrb.Spectra.tocsy(bmrb_ids=[17074,17076,17077], legend='residue')

.. figure:: ../_images/multi_tocsy2.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for TOCSY from BMRB entries 17074,17076 and 17076 with residues as legend <../_static/multi_tocsy2.html>`_

- Example 11 : Multiple entries from BMRB with chemical shift tracking

.. code:: python

    pybmrb.Spectra.tocsy(bmrb_ids=[17074,17076,17077], legend='dataset', draw_trace=True)

.. figure:: ../_images/multi2_tocsy.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for TOCSY from BMRB entries 17074,17076 and 17076 with trace <../_static/multi2_tocsy.html>`_

Please not the above TOCSY with chemical shift visualization will take some time to load, because of hundreds of traces


Generic 2D peak position simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may use any two atoms in a residue to generate a generic 2D spectrum. For the following examples, N chemical shifts
were used as  x axis and CB chemical shifts were was used a Y axis.

- Example 12: Single entry from BMRB

.. code:: python

    pybmrb.Spectra.c13hsqc(bmrb_ids=15060,
    atom_x='N',
    atom_y='CB',
    legend='residue')

.. figure:: ../_images/15060_2d.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for generic 2D spectrum from BMRB entry 15060 <../_static/15060_2d.html>`_

- Example 13: Multiple entries from BMRB

.. code:: python

    pybmrb.Spectra.c13hsqc(bmrb_ids=[17074,17076,17077],
    atom_x='N',
    atom_y='CB',
    legend='dataset')

.. figure:: ../_images/multi_2d.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for generic 2D spectrum from BMRB entries 17074,17076 and 17076 <../_static/multi_2d.html>`_

- Example 14: Multiple entries from BMRB with chemical shift tracking

.. code:: python

    pybmrb.Spectra.c13hsqc(bmrb_ids=[17074,17076,17077],
    atom_x='N',
    atom_y='CB',
    legend='dataset',
    draw_trace=True)

.. figure:: ../_images/multi2_2d.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for generic 2D spectrum from BMRB entries 17074,17076 and 17076 with trace <../_static/multi2_2d.html>`_

Chemical shift Histograms
---------------------------

PyBMRB is able to fetch database wide chemical shift data and plot the distribution in different ways.
The default distribution would be Chemical shift vs number of instances(count). You may also plot the 'percent' or 'probability'
or 'probability density' by providing desired value for 'histnorm'. These distributions
can be filtered using temperature range and PH range. Here are some of the examples.


Single distribution
^^^^^^^^^^^^^^^^^^^^

- Example 15: Chemical shift distribution of CYS-CB

.. code:: python

    pybmrb.Histogram.hist(residue='CYS', atom='CB')

.. figure:: ../_images/cys_cb_hist.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for CYS-CB histogram <../_static/cys_cb_hist.html>`_

- Example 16: Chemical shift distribution of CYS-CB with standard deviation filter

You may exclude extreme values by using the filter based on standard deviation. sd_limit=5 would  exclude
the values beyond 5 times standard deviation on moth sides of the mean

.. code:: python

    pybmrb.Histogram.hist(residue='CYS', atom='CB', sd_limt=5 )

.. figure:: ../_images/cys_cb_hist_sd5.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for CYS-CB histogram with in 5 std on both sides of the mean <../_static/cys_cb_hist_sd5.html>`_

- Example 17: Chemical shift distribution of CYS-CB with Ph filter

You may use experimental conditions like  Ph or temperature values as a filter

.. code:: python

    pybmrb.Histogram.hist(residue='CYS', atom='CB', sd_limt=5,
    ph_min=7.0, ph_max=8.2)

.. figure:: ../_images/cys_cb_hist_ph.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for CYS-CB histogram with Ph filter <../_static/cys_cb_hist_ph.html>`_

- Example 18: Chemical shift distribution of CYS-CB as box plot

Box plot and Violin plot will show all the statistical properties of the distribution, when you mouse over the distribution.

.. code:: python

    pybmrb.Histogram.hist(residue='CYS', atom='CB',
    plot_type='box')

.. figure:: ../_images/cys_cb_box_sd5.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for CYS-CB box plot <../_static/cys_cb_box_sd5.html>`_

- Example 19: Chemical shift distribution of CYS-CB as violin plot

Box plot and Violin plot will show all the statistical properties of the distribution, when you mouse over the distribution.

.. code:: python

    pybmrb.Histogram.hist(residue='CYS', atom='CB',
    plot_type='violin')

.. figure:: ../_images/cys_cb_violin_sd5.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for CYS-CB violin plot <../_static/cys_cb_violin_sd5.html>`_

Multiple distribution
^^^^^^^^^^^^^^^^^^^^^^

- Example 20: Histogram from list of atoms

You may also provide list of atoms as input

.. code:: python

    pybmrb.Histogram.hist(list_of_atoms=['GLN-CB','CYS-CB','TYR-CB'],
    histnorm='probability density')

.. figure:: ../_images/multi_hist.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for multiple distribution <../_static/multi_hist.html>`_

- Example 21: Violin plot for list of atoms



.. code:: python

    pybmrb.Histogram.hist(list_of_atoms=['GLN-CB','CYS-CB','TYR-CB'],
    plot_type='violin')

.. figure:: ../_images/multi_violin.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for violin plot <../_static/multi_violin.html>`_


- Example 22: Histogram method supports wildcard

If you want to see the chemical shift distribution of protons in GLN, then you may use the following command.
You may chose histnorm as 'probability density' to compare distributions

.. code:: python

    pybmrb.Histogram.hist(residue='GLN', atom='H*',
    hist_norm='probability density')

.. figure:: ../_images/gln_h_hist.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for GLN-H* distribution <../_static/gln_h_hist.html>`_

- Example 23: Distribution of all atoms from a residue

If you want to see the chemical shift distribution of all atoms from a residue you may use atom='*' or simply leave out atom.

.. code:: python

    pybmrb.Histogram.hist(residue='ASP', atom='*')

or

.. code:: python

    pybmrb.Histogram.hist(residue='ASP')


.. figure:: ../_images/asp_hist.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for ASP distribution <../_static/asp_hist.html>`_

- Example 24: Distribution of specific atom type from al residues

If you want to see the chemical shift distribution CG atoms from all 20 standard residues you may use residue=*' or simply leave out residue.

.. code:: python

    pybmrb.Histogram.hist(residue='*', atom='CG*',
    hist_norm='percent')

or

.. code:: python

    pybmrb.Histogram.hist(atom='CG*',
    hist_norm='percent')


.. figure:: ../_images/cg_hist.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for CG* distribution <../_static/cg_hist.html>`_

2D Histograms
^^^^^^^^^^^^^^^^^^^^

- Example 25: Chemical shift correlation as 2d heatmap


.. code:: python

    pybmrb.Histogram.hist2d(residue='CYS', atom1='CA', atom2='CB', sd_limut=5)

.. figure:: ../_images/cys-ca-cb.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive 2D heatmap  <../_static/cys-ca-cb.html>`_

- Example 26: Chemical shift correlation as contour heatmap


.. code:: python

    pybmrb.Histogram.hist2d(residue='GLN', atom1='HE21', atom2='HE22',
    sd_limut=5, plot_type='contour')

.. figure:: ../_images/gln-2d.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive contour plot  <../_static/gln-2d.html>`_

Conditional histogram
^^^^^^^^^^^^^^^^^^^^^^

- Example 27: Conditional histogram with chemical shift filtering

You may filter the chemical shift distribution of an atom in a residue based on the chemical shift values of one or
more atom in the same residue. In the following example CYS-CB values are filtered based on CYS-H=8.9. During the seach
0.1ppm tolerance for protons and 2.0 ppm tolerance for heavy atoms is used.

.. code:: python

    pybmrb.Histogram.conditional_hist(residue='CYS', atom='CB', histnorm='percent'
    filtering_rules=[('H',8.9)])

.. figure:: ../_images/filt1.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive conditional histogram with one rule  <../_static/filt1.html>`_

- Example 28: Conditional histogram with chemical shift list


.. code:: python

    pybmrb.Histogram.conditional_hist(residue='CYS', atom='CB', histnorm='percent'
    filtering_rules=[('H', 8.9), ('CA', 61)])

.. figure:: ../_images/filt2.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive conditional histogram with list of rules  <../_static/filt2.html>`_
