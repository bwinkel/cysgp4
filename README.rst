******
cysgp4
******

.. container::

    |License-GPL| |License-Apache| |License-BSD3|

.. container::

    |PyPI Badge| |PyPI Downloads|

.. container::

    |Conda-Version| |Conda-Platforms-Badge| |Conda-Downloads-Badge|

- *Version:* 0.3.6
- *Lead author/package maintainer:* Benjamin Winkel
- *Contributor(s):* Gyula I.G. Józsa
- *Repository:* `on GitHub <https://github.com/bwinkel/cysgp4>`__
- *Bug tracker:* `on GitHub <https://github.com/bwinkel/cysgp4/issues>`__
- *User manual:* `stable <https://bwinkel.github.io/cysgp4/>`__ |  `developer <https://bwinkel.github.io/cysgp4/latest/>`__

Project Status
==============

.. container::

    |Azure Status|

Purpose
=======

The `cysgp4` package is a `Cython <https://www.cython.org>`_-powered wrapper
of the `sgp4lib <https://www.danrw.com/sgp4/>`_ library (by Daniel Warner) to
compute satellite positions from two-line elements (TLE).

It provides similar functionality as the well-known `sgp4
<https://pypi.org/project/sgp4/>`_ Python package (by Brandon Rhodes), which
uses `Numba <http://numba.pydata.org/>`_ internally to speed-up the
calculations. In contrast to `sgp4`_, `cysgp4` can work well with arrays of
TLEs and/or times and make use of multi-core platforms (via OpenMP) to boost
processing times a lot.

Installation
============

We highly recommend to use `cysgp4` with the `Anaconda Python distribution <https://www.anaconda.com/>`_, in which
case installiation is as easy as ::

    conda install -c conda-forge cysgp4

Otherwise, you should install cysgp4 via `pip`::

    python -m pip install cysgp4

The installation is also possible from source. `Detailed installation
instructions <https://bwinkel.github.io/cysgp4/latest/install.html>`_
can be found in the user manual.

Dependencies
------------

We kept the dependencies as minimal as possible. The following packages are
required:

- `Python 3.8` or later
- `numpy 1.20` or later

If you want to run the notebooks yourself, you will also need the Jupyter
server and install matplotlib. To run the tests, you'll need `sgp4
<https://pypi.org/project/sgp4/>`_.

Note, for compiling the C-extension, OpenMP is used for parallelization. If you use gcc, for example, you should have at least version 4.8 otherwise the setup-script may fail. Again, see `Detailed installation instructions` for
more information.

Usage
=====

Using cysgp4 is possible via an object-oriented interface or with a
fast numpy-array functional approach. The former works like this:

.. code-block:: python

    import cysgp4

    # Define a date/time and an observer
    pydt = cysgp4.PyDateTime.from_mjd(58805.57)
    lon_deg, lat_deg = 6.88375, 50.525
    alt_km = 0.366
    obs = cysgp4.PyObserver(lon_deg, lat_deg, alt_km)

    # Define satellite properties/orbit via two-line element string (TLE)
    hst_tle = cysgp4.PyTle(
        'HST',
        '1 20580U 90037B   19321.38711875  .00000471  00000-0  17700-4 0  9991',
        '2 20580  28.4699 288.8102 0002495 321.7771 171.5855 15.09299865423838',
        )

    # Create a satellite object for querying coordinates
    sat = cysgp4.Satellite(hst_tle, obs, pydt)
    sat.eci_pos().loc  # ECI cartesian position, km
    (5879.5931344459295, 1545.7455647032068, 3287.4155452595)
    sat.eci_pos().vel  # ECI cartesian velocity, km/s
    (-1.8205895517672226, 7.374044252723081, -0.20697960810978586)
    sat.geo_pos()  # geographic (geodetic) position, lon/lat/alt
    <PyCoordGeodetic: 112.2146d, 28.5509d, 538.0186km>
    sat.topo_pos()  # topocentric position, az/el/dist/dist_rate
    <PyCoordTopocentric: 60.2453d, -35.6844d, 8314.5683km, 3.5087km/s>

    # One can change time to determine positions at another moment
    sat.mjd += 1 / 720.  # one minute later
    sat.topo_pos()
    <PyCoordTopocentric: 54.8446d, -38.2749d, 8734.9195km, 3.4885km/s>

In many cases, however, one probably wants to calculate coordinates for a
(large) number of satellites, observer locations, and/or observing times. For
this, the function `~cysgp4.propagate_many` is useful. This is an array
interface to the sgp4 calculations, which allows to perform calculations for
different satellite TLEs, observers and times in a parallelized manner.
`~numpy` broadcasting `rules
<https://docs.scipy.org/doc/numpy/user/basics.broadcasting.html>`_ apply:

.. code-block:: python

        import requests
        import numpy as np
        from cysgp4 import PyTle, PyObserver, propagate_many

        # Download many TLEs from a website
        url = 'http://celestrak.com/NORAD/elements/science.txt'
        ctrak_science = requests.get(url)
        all_lines = ctrak_science.text.split('\\r\\n')

        # Need to convert them to a list of tuples (each tuple consisting
        # of the three TLE strings)
        tle_list = list(zip(*tuple(
            all_lines[idx::3] for idx in range(3)
            )))
        # Create an array of PyTle and PyObserver objects, and MJDs
        tles = np.array([
            PyTle(*tle) for tle in tle_list
            ])[np.newaxis, np.newaxis, :20]  # use first 20 TLEs
        observers = np.array([
            PyObserver(6.88375, 50.525, 0.366),
            PyObserver(16.88375, 50.525, 0.366),
            ])[np.newaxis, :, np.newaxis]
        mjds = np.linspace(
            58805.5, 58806.5, 1000  # 1000 time steps
            )[:, np.newaxis, np.newaxis]

        # The result is a dictionary
        result = propagate_many(mjds, tles, observers)
        print(result.keys())
        dict_keys(['eci_pos', 'eci_vel', 'geo', 'topo'])

        # Returned array shapes are as follows; last array dimension
        # contains the coordinate pairs.
        print(np.broadcast(mjds, tles, observers).shape)
        (1000, 2, 20)
        print(result['eci_pos'].shape, result['topo'].shape)
        (1000, 2, 20, 3) (1000, 2, 20, 4)

        # One can also skip over coordinate frames.
        result = propagate_many(
            mjds, tles, observers,
            do_eci_pos=False, do_eci_vel=False, do_geo=False, do_topo=True
            )
        print(result.keys())
        dict_keys(['topo'])


More use-cases and tutorials
----------------------------

Check out the `user manual <https://bwinkel.github.io/cysgp4/latest/>`_ or the
`Jupyter tutorial notebooks <https://github.com/bwinkel/cysgp4/tree/master/notebooks>`_
in the repository for further examples of how to use `cysgp4`. Note that you
can only view the notebooks on GitHub, if you want to edit something
it is necessary to clone the repository or download a notebook to run it on
your machine.

Who do I talk to?
=================

If you encounter any problems or have questions, do not hesitate to raise an
issue or make a pull request. Moreover, you can contact the devs directly:

- <bwinkel@mpifr.de>

Licenses
========

`cysgp4` itself is published under `GPL v3 <https://www.github.com/bwinkel/cysgp4/blob/master/COPYING.GPLv3.txt>`_, an open-source license. The package
is a `Cython <https://www.cython.org>`_-powered wrapper of the `sgp4lib
<https://www.danrw.com/sgp4/>`_ library (by Daniel Warner) to compute
satellite positions from two-line elements (TLE). The sgp4lib source code is
licensed under `Apache-2.0 license
<https://www.github.com/bwinkel/cysgp4/blob/master/COPYING.Apache2>`_

The package is partly based on the `Astropy-affiliated package template <https://github.com/astropy/package-template>`_, which is under `BSD 3-clause <https://github.com/bwinkel/cysgp4/blob/master/TEMPLATE_LICENCE.rst>`_ license.




.. |PyPI Badge| image:: https://img.shields.io/pypi/v/cysgp4.svg
    :target: https://pypi.python.org/pypi/cysgp4
    :alt: PyPI tag

.. |PyPI Downloads| image:: https://img.shields.io/pypi/dm/cysgp4
   :target: https://pypi.python.org/pypi/cysgp4
   :alt: PyPI - Downloads

.. |License-GPL| image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: https://www.github.com/bwinkel/cysgp4/blob/master/COPYING.GPLv3.txt
    :alt: License-GPL3

.. |License-Apache| image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://www.github.com/bwinkel/cysgp4/blob/master/COPYING.Apache2
    :alt: License-Apache

.. |License-BSD3| image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
    :target: https://www.github.com/bwinkel/cysgp4/blob/master/TEMPLATE_LICENCE.rst
    :alt: License-BSD3

.. |Conda-Version| image:: https://anaconda.org/conda-forge/cysgp4/badges/version.svg
    :target: https://anaconda.org/conda-forge/cysgp4
    :alt: conda-forge platforms: Version on conda-forge

.. |Conda-Platforms-Badge| image:: https://anaconda.org/conda-forge/cysgp4/badges/platforms.svg
    :target: https://anaconda.org/conda-forge/cysgp4
    :alt: conda-forge platforms: linux-64, osx-64, win-64

.. |Conda-Downloads-Badge| image:: https://anaconda.org/conda-forge/cysgp4/badges/downloads.svg
    :target: https://anaconda.org/conda-forge/cysgp4
    :alt: conda-forge downloads


.. |Azure Status| image:: https://dev.azure.com/bwinkel78/Benjamin-Winkel-Projects/_apis/build/status/bwinkel.cysgp4?branchName=master
    :target: https://dev.azure.com/bwinkel78/Benjamin-Winkel-Projects/_build
    :alt: cysgp4's Azure Pipelines Status


