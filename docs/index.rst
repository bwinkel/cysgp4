
:tocdepth: 3

#####################
cysgp4 Documentation
#####################

The `cysgp4` package is a `Cython <https://www.cython.org>`_-powered wrapper
of the `sgp4lib <https://www.danrw.com/sgp4/>`_ library (by Daniel Warner) to
compute satellite positions from two-line elements (TLE).

It provides similar functionality as the well-known `sgp4
<https://pypi.org/project/sgp4/>`_ Python package (by Brandon Rhodes), which
uses `Numba <http://numba.pydata.org/>`_ internally to speed-up the
calculations. In contrast to `sgp4`_, `cysgp4` can work well with arrays of
TLEs and/or times and make use of multi-core platforms (via `OpenMP
<https://www.openmp.org/>`_) to boost processing times a lot.

***************
Getting Started
***************

.. toctree::
   :maxdepth: 1

   install
   importing_cysgp4
   quick_tour

******************
User Documentation
******************

.. toctree::
   :maxdepth: 1

   user_manual
   Tutorials (Jupyter notebooks) <http://nbviewer.jupyter.org/github/bwinkel/cysgp4/blob/master/notebooks/>


***************
License
***************

.. toctree::
   :maxdepth: 1

   license


***************
Acknowledgments
***************

This code makes use of the excellent work provided by the
`Astropy <http://www.astropy.org/>`__ community. Cysgp4 uses the
`Astropy Package Template <https://github.com/astropy/package-template>`__
for the packaging.
