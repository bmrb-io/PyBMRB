Examples
========

Load the moduel
----------------

.. code:: python

    import pybmrb

.. |n15| replace:: :sup:`1` H - :sup:`15` N
.. |c13| replace:: :sup:`1` H - :sup:`13` C
.. |hh| replace:: :sup:`1` H - :sup:`1` H


|n15| - HSQC peak position simulation
---------------------------------------

Single entry from BMRB

.. code:: python

    pybmrb.Spectra.n15hsqc(bmrb_ids=15060, legend='residue')

.. figure:: ../_images/15060_n15.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for n15-HSQC from BMRB entry 15060 <../_static/15060_n15.html>`_

Multiple entries from BMRB

.. code:: python

    pybmrb.Spectra.n15hsqc(bmrb_ids=[17074,17076,17077], legend='dataset')

.. figure:: ../_images/multi_n15.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for n15-HSQC from BMRB entries 17074,17076 and 17076 <../_static/multi_n15.html>`_

Multiple entries from BMRB with chemical shift tracking

.. code:: python

    pybmrb.Spectra.n15hsqc(bmrb_ids=[17074,17076,17077], legend='dataset', draw_trace=True)

.. figure:: ../_images/multi2_n15.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for n15-HSQC from BMRB entries 17074,17076 and 17076 with trace <../_static/multi2_n15.html>`_


|c13| - HSQC peak position simulation
---------------------------------------

Single entry from BMRB

.. code:: python

    pybmrb.Spectra.c13hsqc(bmrb_ids=15060, legend='residue')

.. figure:: ../_images/15060_c13.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for c13-HSQC from BMRB entry 15060 <../_static/15060_c13.html>`_

Multiple entries from BMRB

.. code:: python

    pybmrb.Spectra.c13hsqc(bmrb_ids=[17074,17076,17077], legend='dataset')

.. figure:: ../_images/multi_c13.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for c13-HSQC from BMRB entries 17074,17076 and 17076 <../_static/multi_c13.html>`_

Multiple entries from BMRB with chemical shift tracking

.. code:: python

    pybmrb.Spectra.c13hsqc(bmrb_ids=[17074,17076,17077], legend='dataset', draw_trace=True)

.. figure:: ../_images/multi2_c13.jpg
    :alt: n15hsqc
    :align: center

    `Click here for interactive html for c13-HSQC from BMRB entries 17074,17076 and 17076 with trace <../_static/multi2_c13.html>`_

