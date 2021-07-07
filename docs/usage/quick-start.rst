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

.. code:: python

    peak_list = Spectra.n15hsqc(bmrb_ids=15060,legend='residue')

.. |n15| replace:: :sup:`1` H - :sup:`15` N
.. |c13| replace:: :sup:`1` H - :sup:`13` C
.. |hh| replace:: :sup:`1` H - :sup:`1` H
