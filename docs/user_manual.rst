.. user-manual:

******************
Cysgp4 user manual
******************

.. currentmodule:: cysgp4



Introduction
============

The `~cysgp4` package can be used to calculate the orbits of satellites around
Earth, based on so-called `Two-line element strings (TLE)`_. They not only
encapsulate not only the basic `orbital elements`_ that describe Keplarian
orbits (based on six free parameters), but also include more complicated
physical processes, such as drag caused by Earth's atmosphere and other
pertubations owing to the ellipsoidal gravitational field, influence of other
solar bodies etc. One of the more famous `simplified perturbations models`_
is SGP4, which is directly related to the TLEs published by NORAD and NASA,
and many tools exist to do the calculation.

Here we provide a Python wrapper around the `C++ SGP4 Satellite Library`_ by
Daniel Warner. It provides similar functionality as the well-known `sgp4
<https://pypi.org/project/sgp4/>`_ Python package (by Brandon Rhodes), which
uses `Numba <http://numba.pydata.org/>`_ internally to speed-up the
calculations. In contrast to `sgp4`_, `cysgp4` can work well with arrays of
TLEs and observing times and make use of multi-core platforms (via `OpenMP
<https://www.openmp.org/>`_) to boost processing times a lot.

Using cysgp4
============

The Python package `~cysgp4` provides a relatively thin wrapper around the
`C++ SGP4 Satellite Library`_ on the one hand, but on the other hand also
contains some higher-level functions, which make it possible to bulk-process
satellite orbits with few lines of code.


.. _wrapped-approach-label:

Option 1: C++ SGP4 library wrapper
----------------------------------

The package provides the following Classes.

PyDateTime
^^^^^^^^^^
`~cysgp4.PyDateTime` specifies a datetime. Various conversions are provided
for convenience, from and to `Modified Julian Date (MJD)
<https://en.wikipedia.org/wiki/Julian_day#Variants>`_, Python's `~datetime`
class, so-called "ticks", which are the number of micro-seconds since 1.
January 0001, 00:00, but also the TLE-epoch format.

.. note::

    TLE epochs have a very strange format: first two digits are the year,
    next three digits are the day from beginning of year, then the fraction
    of a day is given, e.g. 20180.25 would be 2020, day 180, 6 hours (
    probably UT as no timezone is mentioned). See also `Two-line element
    strings (TLE)`_ or `Celestrak <https://celestrak.com/columns/v04n03/>`_.

There are several possibilities to create a `~cysgp4.PyDateTime` object,
e.g.::

    >>> from datetime import datetime
    >>> from cysgp4 import PyDateTime

    >>> # from Python's datetime
    >>> dt = datetime(2019, 1, 1, 12, 13, 14)
    >>> pydt = PyDateTime(dt)
    >>> pydt
    <PyDateTime: 2019-01-01 12:13:14.000000 UTC>

    >>> # from Modified Julian Date
    >>> mjd = 56458.123
    >>> pydt = PyDateTime.from_mjd(mjd)
    >>> pydt
    <PyDateTime: 2013-06-15 02:57:07.199999 UTC>

    >>> # from "Ticks"
    >>> ticks = 58628880000000000  # aka MJD zero point
    >>> pydt = PyDateTime.from_ticks(ticks)
    >>> pydt
    <PyDateTime: 1858-11-17 00:00:00.000000 UTC>

Accessing and modifying the date and time is done via the following
properties and methods::

    >>> print(pydt.datetime)
    1858-11-17 00:00:00
    >>> pydt.datetime = datetime(2010, 5, 1, 0, 10, 20)

    >>> print(pydt.mjd)
    55317.00717592592
    >>> pydt.mjd = 55489.01358

    >>> print(pydt.ticks)
    63423130773311999
    >>> pydt.ticks = 63681941594000000

    >>> observer_longitude = 6.67  # degrees
    >>> print('GMST: {:9.6f} hours'.format(pydt.gmst()))
    GMST:  4.959715 hours
    >>> print('LMST: {:9.6f} hours'.format(pydt.lmst(observer_longitude)))
    LMST:  5.076129 hours

One can also set the date and time manually with the `~cysgp4.PyDateTime.set`
method.

PyTle
^^^^^
This helper class holds the `Two-line element strings (TLE)`_ and provides
access to (some of) its parameters. Although the name suggests that two lines
of strings are used, in practice often three lines are defined, the first
containing a satellite name (and/or) ID. It is important to note that for
many satellites, the corresponding TLEs get outdated quickly. Especially for
low earth orbit satellites, one should query up-to-date information at least
once per day.

`~cysgp4.PyTle` is only used to parse the TLE strings and store the orbital
parameters for use by other SGP4 routines. TLEs can be obtained from many
sources, such as `Celestrak`_::

    >>> import requests
    >>> import cysgp4
    >>> url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=science&FORMAT=tle'

.. doctest-remote-data::
    >>> ctrak_science = requests.get(url).text
    >>> all_lines = ctrak_science.split('\r\n')
    >>> all_lines[:3]  # doctest: +SKIP
    ['AKEBONO (EXOS-D)        ',
     '1 19822U 89016A   19331.07700297  .00016037  93960-6  28103-3 0  9992',
     '2 19822  75.0289 316.1359 1764871 261.4338  78.3806 12.01292738 15108']
    >>> tle_list = list(zip(*tuple(
    ...     all_lines[idx::3]
    ...     for idx in range(3)
    ...     )))
    >>> len(tle_list)  # doctest: +SKIP
    65
    >>> print(*tle_list[1], sep='\n')  # doctest: +SKIP
    HST
    1 20580U 90037B   19321.38711875  .00000471  00000-0  17700-4 0  9991
    2 20580  28.4699 288.8102 0002495 321.7771 171.5855 15.09299865423838

.. note::
    As the above code example will download the latest (daily!) TLE from
    `Celestrak`_, you'll be seeing a different set of numbers.

Because re-structuring the TLE text from the download or file into a list of
tuples is often needed, there is a utility routine to do just that

.. doctest-remote-data::

    >>> tle_list = cysgp4.tle_tuples_from_text(ctrak_science)

And it will do the exact same thing. For testing, `~cysgp4` has one example
TLE file onboard::

    >>> ctrak_science = cysgp4.get_example_tles()
    >>> ctrak_science[:200]
    'AKEBONO (EXOS-D)        \r\n1 19822U 89016A   19321.49921565  .00016421  94291-6  28704-3 0  9992\r\n2 19822  75.0304 327.7460 1766579 276.4058  63.9734 12.00957719 13956\r\nHST                     \r\n1 2058'

    >>> tle_list = cysgp4.tle_tuples_from_text(ctrak_science)
    >>> tle_list  # doctest: +SKIP
    [('AKEBONO (EXOS-D)        ',
      '1 19822U 89016A   19321.49921565  .00016421  94291-6  28704-3 0  9992',
      '2 19822  75.0304 327.7460 1766579 276.4058  63.9734 12.00957719 13956'),
     ('HST                     ',
      '1 20580U 90037B   19321.38711875  .00000471  00000-0  17700-4 0  9991',
      '2 20580  28.4699 288.8102 0002495 321.7771 171.5855 15.09299865423838'),
     ...
     ('ZHANGZHENG-1 (CSES)     ',
      '1 43194U 18015C   19322.55091545  .00001006  00000-0  48427-4 0  9993',
      '2 43194  97.4135  85.6925 0016973 127.3937   1.4087 15.20935646 99453')]

Now, let's feed this into `~cysgp4.PyTle`::

    >>> hst_tle = cysgp4.PyTle(*tle_list[1])
    >>> hst_tle
    <PyTle: HST                     >
    >>> print(hst_tle)
    Norad Number:         20580
    Int. Designator:      90037B
    Epoch:                2019-11-17 09:17:27.060000 UTC
    Orbit Number:         42383
    Mean Motion Dt2:        0.00000471
    Mean Motion Ddt6:       0.00000000
    Eccentricity:           0.00024950
    BStar:                  0.00001770
    Inclination:           28.46990000
    Right Ascending Node: 288.81020000
    Argument Perigee:     321.77710000
    Mean Anomaly:         171.58550000
    Mean Motion:           15.09299865

PyCoordGeodetic
^^^^^^^^^^^^^^^
This is one of three coordinate frames, which are implemented. The `Geodetic
frame <https://en.wikipedia.org/wiki/Geodetic_datum>`_ is the same as the
well-known geographic longitude and latitude, plus an altitude above Earth's
surface. As such, it is fixed to the (rotating) Earth, in contrast to the
`ECI <https://en.wikipedia.org/wiki/Earth-centered_inertial>`_ system (see
below).

.. note::
    The `Geodetic frame`_ goes usually along with a specific definition of
    the Geoid (Earth's ellipsoid). For many application, the
    `WGS84 <https://en.wikipedia.org/wiki/World_Geodetic_System#A_new_World_Geodetic_System:_WGS_84>`_ is used.
    However, most TLEs, which can be downloaded, still use WGS72. The
    underlying `C++ SGP4 Satellite Library`_ unfortunately has the WGS72
    hardcoded.

Constructing and using a `~cysgp4.PyCoordGeodetic` object is straightforward::

    >>> from cysgp4 import PyCoordGeodetic

    >>> lon_deg, lat_deg = 6.88375, 50.525
    >>> alt_km = 0.366
    >>> geo = PyCoordGeodetic(lon_deg, lat_deg, alt_km)
    >>> geo
    <PyCoordGeodetic: 6.8838d, 50.5250d, 0.3660km>

    >>> # Access is also possible via properties, e.g.:
    >>> geo.lon
    6.88375
    >>> geo.alt = 0.4


PyEci
^^^^^
The `~cysgp4.PyEci` class holds an `ECI`_ location (latitude, longitude,
altitude) for a particular datetime. Note, internally, the coordinates (and
velocities) are stored in Cartesian form (read-only!). One can access these,
via the `~cysgp4.PyEci.loc` and `~cysgp4.PyEci.vel` properties. Setting
new parameters is only possible via a geographic location, i.e., via
`~cysgp4.PyCoordGeodetic`, which is then converted to Cartesian ECI based on
the provided `~cysgp4.PyDateTime`::

    >>> from cysgp4 import PyEci, PyCoordGeodetic, PyDateTime

    >>> pydt = PyDateTime.from_mjd(55555.)
    >>> lon_deg, lat_deg = 6.88375, 50.525
    >>> alt_km = 0.366
    >>> geo = PyCoordGeodetic(lon_deg, lat_deg, alt_km)
    >>> eci = PyEci(pydt, geo)
    >>> eci
    <PyEci: 6.8837d, 50.5250d, 0.3660km 2010-12-25 00:00:00.000000 UTC>

    >>> # Access is also possible via properties, e.g.:
    >>> eci.loc  # read-only!  # doctest: +FLOAT_CMP
    (-725.3304166274728, 3997.924210010933, 4900.402205553537)
    >>> eci.vel  # read-only!  # doctest: +FLOAT_CMP
    (-0.2915332651982093, -0.05289193431368386, 0.0)
    >>> eci.pydt = PyDateTime.from_mjd(55556.)

    >>> # or the update method:
    >>> eci.update(pydt, PyCoordGeodetic(0, 0, 0))


PyCoordTopocentric
^^^^^^^^^^^^^^^^^^
The third of the three coordinate frames is the `topocentric frame
<https://en.wikipedia.org/wiki/Horizontal_coordinate_system>`_ (or horizontal
frame), which is provided by `~cysgp4.PyCoordTopocentric`. It holds a
topocentric location (azimuth, elevation, range/distance and distance/range
rate).

.. note::
    The topocentric position of a satellite is always relative to
    an observer (in SGP4 the observer is defined in geographic coordinates).
    The distance and distance change (aka range rate) is thus the distance
    between the satellite and the observer. However,
    `~cysgp4.PyCoordTopocentric` only holds the azimuth, elevation, distance
    and distance rate parameters, but contains no information on the
    observer. It is only useful in conjuction with a `~cysgp4.PyObserver` (
    which holds a reference to a geographic location) and a datetime; see
    also the discussion of `~cysgp4.Satellite` here: :ref:`satellite-label`.

Constructing and using a `~cysgp4.PyCoordTopocentric` object is
straightforward::

    >>> from cysgp4 import PyCoordTopocentric

    >>> az_deg, el_deg = 130.1, 10.53
    >>> dist_km, dist_rate_kms = 1200., 0.03
    >>> topo = PyCoordTopocentric(az_deg, el_deg, dist_km, dist_rate_kms)
    >>> topo
    <PyCoordTopocentric: 130.1000d, 10.5300d, 1200.0000km, 0.0300km/s>

    >>> # Access is also possible via properties, e.g.:
    >>> topo.az
    130.1
    >>> topo.dist = 1000.

PyObserver
^^^^^^^^^^
The `~cysgp4.PyObserver` class holds the location (as `~cysgp4.PyEci`) of
an observer.

.. note::

    Usually, a `~cysgp4.PyDateTime` is attached to an ECI location in the
    underlying `C++ SGP4 Satellite Library`_.
    However, `~cysgp4.PyObserver` does not allow to specify the datetime
    and internally it is always set to zero.

Using `~cysgp4.PyObserver` is similar to working with
`~cysgp4.PyCoordGeodetic`::


    >>> from cysgp4 import PyObserver, PyCoordGeodetic

    >>> lon_deg, lat_deg = 6.88375, 50.525
    >>> alt_km = 0.366
    >>> obs = PyObserver(lon_deg, lat_deg, alt_km)
    >>> obs
    <PyObserver: 6.8838d, 50.5250d, 0.3660km>

    >>> # Access is also possible via location property:
    >>> obs.loc
    <PyCoordGeodetic: 6.8838d, 50.5250d, 0.3660km>
    >>> obs.loc = PyCoordGeodetic(1, 2, 3)

We will see in the next section, how the observer interplays with the
satellite orbit calculation.

.. _satellite-label:

Satellite interface
^^^^^^^^^^^^^^^^^^^

With the `~cysgp4.Satellite` class, one can calculate positions of a
satellite at given times. Here, the satellite is defined via an TLE (see
`~cysgp4.PyTle`). Furthermore, to be able to calculate apparent positions of
the satellite, an observer (see `~cysgp4.PyObserver`) needs to be defined.
The position calculations are lazy, i.e., only calculated if the newly
requested time differs by a certain amount from the last time at which a
calculation was performed. The granularity of this "cache" can be defined via
the `mjd_cache_resolution` parameter.

The following demonstrates how a typical use of the `~cysgp4.Satellite` class
would look like::

    >>> from cysgp4 import *

    >>> # Specify a datetime, for which the satellite position is desired
    >>> pydt = PyDateTime.from_mjd(58805.57)

    >>> # Define the observer...
    >>> lon_deg, lat_deg = 6.88375, 50.525
    >>> alt_km = 0.366
    >>> obs = PyObserver(lon_deg, lat_deg, alt_km)

    >>> # ... and TLE (i.e., satellite parameters)
    >>> hst_tle = PyTle(
    ... 'HST',
    ... '1 20580U 90037B   19321.38711875  .00000471  00000-0  17700-4 0  9991',
    ... '2 20580  28.4699 288.8102 0002495 321.7771 171.5855 15.09299865423838',
    ... )

    >>> # All of this goes into the Satellite constructor
    >>> sat = Satellite(hst_tle, obs, pydt)

    >>> # We can now query positions...
    >>> sat.eci_pos().loc  # ECI cartesian position  # doctest: +FLOAT_CMP
    (5879.5931344459295, 1545.7455647032068, 3287.4155452595)
    >>> sat.eci_pos().vel  # ECI cartesian velocity  # doctest: +FLOAT_CMP
    (-1.8205895517672226, 7.374044252723081, -0.20697960810978586)
    >>> sat.geo_pos()  # geographic position  # doctest: +FLOAT_CMP
    <PyCoordGeodetic: 112.2146d, 28.5509d, 538.0173km>
    >>> sat.topo_pos()  # topocentric position  # doctest: +FLOAT_CMP
    <PyCoordTopocentric: 60.2453d, -35.6845d, 8314.5681km, 3.5087km/s>

    >>> # ... also for different times, also for different times
    >>> sat.mjd += 1 / 720.  # one minute later
    >>> sat.topo_pos()  # doctest: +FLOAT_CMP
    <PyCoordTopocentric: 54.8446d, -38.2749d, 8734.9196km, 3.4885km/s>
    >>> sat.topo_pos().az, sat.topo_pos().el  # doctest: +FLOAT_CMP
    (54.84463503781068, -38.274852915850126)

    >>> # Change by less than cache resolution (1 ms) produces the same result
    >>> sat.mjd += 0.0005 / 86400.  # 0.5 ms
    >>> sat.topo_pos().az, sat.topo_pos().el  # doctest: +FLOAT_CMP
    (54.84463503781068, -38.274852915850126)

    >>> # Change by another 0.5 ms triggers re-calculation
    >>> sat.mjd += 0.00051 / 86400.
    >>> sat.topo_pos().az, sat.topo_pos().el  # doctest: +FLOAT_CMP
    (54.844568313870965, -38.274885794151324)

Option 2: Batch calculations
----------------------------
All of the above can of course be combined with for-loops to process many
different times, or several satellites. And this would be reasonably fast.
However, as the wrapper around `C++ SGP4 Satellite Library`_ is powered
with `Cython <https://www.cython.org>`_ it was possible to do implement the
for-loops directly in Cython to reach C-level speed. Furthermore, `OpenMP`_
can parallelize the computation such that it makes use of multi-core
platforms.

The batch processing is implemented with the `~cysgp4.propagate_many`
function. We start by defining datetimes (as MJDs), observers, and TLEs::

    >>> import numpy as np
    >>> from cysgp4 import PyTle, PyObserver, propagate_many
    >>> from cysgp4 import get_example_tles, tles_from_text

    >>> tle_text = get_example_tles()
    >>> tles = np.array(tles_from_text(tle_text))
    >>> observers = np.array([
    ...     PyObserver(6.88375, 50.525, 0.366),
    ...     PyObserver(16.88375, 50.525, 0.366),
    ...     ])
    >>> mjds = np.linspace(58805.5, 58806.5, 1000)  # 1000 time steps

The `~cysgp4.propagate_many` will apply numpy's `broadcasting rules
<https://docs.scipy.org/doc/numpy/user/basics.broadcasting.html>`_. As in
this case all input parameters are arrays, we need to add two (length-1) axes
to each of them to make them "compatible" (think of this as an "outer
product". Here we choose to make `mjds` vary on the 1st axis, `tles` on the
3rd axis, and observers on the 2nd axis.::

    >>> result = propagate_many(
    ...     mjds[:, np.newaxis, np.newaxis],
    ...     tles[np.newaxis, np.newaxis, :20],  # use first 20 TLEs,
    ...     observers[np.newaxis, :, np.newaxis],
    ...     )

Here, `result` is a Python dictionary, with the following entries:

    >>> print(sorted(result.keys()))
    ['eci_pos', 'eci_vel', 'geo', 'topo']

Each entry is an array containing the coordinates of the frame (on the last
axis). The broadcasted input arrays form the first axes of the output::

    >>> # Shapes are as follows
    >>> print(np.broadcast(
    ...     mjds[:, np.newaxis, np.newaxis],
    ...     tles[np.newaxis, np.newaxis, :20],  # use first 20 TLEs,
    ...     observers[np.newaxis, :, np.newaxis],
    ...     ).shape)
    (1000, 2, 20)
    >>> print(result['eci_pos'].shape, result['topo'].shape)
    (1000, 2, 20, 3) (1000, 2, 20, 4)
    >>> # First mjd, first sat, first observer:
    >>> print('{:.3f} {:.3f} {:.3f}'.format(*result['eci_pos'][0, 0, 0]))
    6611.457 -4125.652 762.758
    >>> print('{:.3f} {:.3f} {:.3f} {:.3f}'.format(*result['topo'][0, 0, 0]))
    90.982 -34.119 9361.523 -2.074

    >>> # One can also skip over coordinate frames.
    >>> result = propagate_many(
    ...     mjds[:, np.newaxis, np.newaxis],
    ...     tles[np.newaxis, np.newaxis, :20],  # use first 20 TLEs,
    ...     observers[np.newaxis, :, np.newaxis],
    ...     do_eci_pos=False, do_eci_vel=False, do_geo=False, do_topo=True
    ...     )
    >>> print(sorted(result.keys()))
    ['topo']



Benchmarks
==========

Because we make use of `Cython`_, `~cysgp4` turns out to be really fast. In
the following, a few benchmarks are run. First, it should be mentioned that
`~cysgp4.propagate_many` will utilize all CPU cores that are available on
your system. It is possible to change that by calling::

    >>> num_cores = 1
    >>> cysgp4.set_num_threads(num_cores)  # set desired number of cores

With the Python `~timeit` package, it is possible to do some speed testing.
We can also try out, how well the parallelization is working::

    >>> import numpy as np
    >>> import cysgp4
    >>> import timeit

    >>> def prepare_arrays():
    ...
    ...     tle_text = cysgp4.get_example_tles()
    ...     tles = np.array(cysgp4.tles_from_text(tle_text))
    ...     observers = np.array([
    ...         cysgp4.PyObserver(6.88375, 50.525, 0.366),
    ...         cysgp4.PyObserver(16.88375, 50.525, 0.366),
    ...         ])
    ...     mjds = np.linspace(58805.5, 58806.5, 1000)  # 1000 time steps
    ...
    ...     return (
    ...         mjds[:, np.newaxis, np.newaxis],
    ...         tles[np.newaxis, np.newaxis, :],
    ...         observers[np.newaxis, :, np.newaxis],
    ...         )

    >>> mjds, tles, observers = prepare_arrays()

    >>> for num_cores in [1, 2, 4, 8, 16]:  # doctest: +SKIP
    ...
    ...     cysgp4.set_num_threads(num_cores)
    ...     run_times = timeit.Timer(
    ...         'cysgp4.propagate_many(mjds, tles, observers)',
    ...         globals=globals()
    ...         ).repeat(3, number=10)  # do 3 repeats, 10 runs each
    ...
    ...     print('{:2d} cores: {:.2f} seconds (best of three)'.format(
    ...           num_cores, min(run_times)))
     1 cores: 7.16 seconds (best of three)
     2 cores: 3.63 seconds (best of three)
     4 cores: 1.86 seconds (best of three)
     8 cores: 0.99 seconds (best of three)
    16 cores: 0.58 seconds (best of three)

How does this relate to using Python for-loops with the wrapped classes (see
:ref:`wrapped-approach-label`)? There is a function,
`~cysgp4.propagate_many_slow` that does exactly that and is merely included
in the package for benchmarking. You can see its `source code here <_modules/
cysgp4/helpers.html#propagate_many_slow>`_. We can compare this to the faster
routine like this::

    >>> run_times = timeit.Timer(  # doctest: +SKIP
    ...     'cysgp4.propagate_many_slow(mjds, tles, observers)',
    ...     globals=globals()
    ...     ).repeat(3, number=10)  # do 3 repeats, 10 runs each
    >>> print('slow version: {:.2f} seconds (best of three)'.format(
    ...       min(run_times)))  # doctest: +SKIP
    slow version: 74.94 seconds (best of three)

Last but not least, we are interested in a comparison with the well-known
`sgp4`_ Python package (by Brandon Rhodes). Again, for the sole purpose of
benchmarking, a function was added, which uses `sgp4`_ but has the same
interface as `~cysgp4.propagate_many`. The source code of
`~cysgp4.propagate_many_sgp4` is available in the `manual <_modules/
cysgp4/helpers.html#propagate_many_sgp4>`_, as well. Here are the results::

    >>> run_times = timeit.Timer(  # doctest: +SKIP
    ...     'cysgp4.propagate_many_sgp4(mjds, tles, observers)',
    ...     globals=globals()
    ...     ).repeat(3, number=10)  # do 3 repeats, 10 runs each
    >>> print('sgp4 version: {:.2f} seconds (best of three)'.format(
    ...       min(run_times)))  # doctest: +SKIP
    sgp4 version: 236.46 seconds (best of three)

As you can see, there is some significant speed-up to be gained.


See Also
========

- `C++ SGP4 Satellite Library <https://www.danrw.com/sgp4/>`_ by Daniel Warner
- `sgp4 Python package <https://pypi.org/project/sgp4/>`_ by Brandon Rhodes
- `Two-line element strings (TLE)
  <https://en.wikipedia.org/wiki/Two-line_element_set>`_
- `Orbital elements <https://en.wikipedia.org/wiki/Orbital_elements>`_
- `Simplified perturbations models
  <https://en.wikipedia.org/wiki/Simplified_perturbations_models>`_


Reference/API
=============

.. automodapi:: cysgp4
    :inherited-members:
    :no-inheritance-diagram:
    :no-main-docstr:
..    :no-heading:

