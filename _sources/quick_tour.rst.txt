.. user-manual:

******************
Cysgp4 quick tour
******************

.. currentmodule:: cysgp4

Using cysgp4 is possible via an object-oriented interface or with a
fast `~numpy`-array functional approach. The former works like this::

    >>> import cysgp4

    >>> # Define a date/time and an observer
    >>> pydt = cysgp4.PyDateTime.from_mjd(58805.57)
    >>> lon_deg, lat_deg = 6.88375, 50.525
    >>> alt_km = 0.366
    >>> obs = cysgp4.PyObserver(lon_deg, lat_deg, alt_km)

    >>> # Define satellite properties/orbit via two-line element string (TLE)
    >>> hst_tle = cysgp4.PyTle(
    ...     'HST',
    ...     '1 20580U 90037B   19321.38711875  .00000471  00000-0  17700-4 0  9991',
    ...     '2 20580  28.4699 288.8102 0002495 321.7771 171.5855 15.09299865423838',
    ...     )

    >>> # Create a satellite object for querying coordinates
    >>> sat = cysgp4.Satellite(hst_tle, obs, pydt)
    >>> sat.eci_pos().loc  # ECI cartesian position, km  # doctest: +FLOAT_CMP
    (5879.5931344459295, 1545.7455647032068, 3287.4155452595)
    >>> sat.eci_pos().vel  # ECI cartesian velocity, km/s  # doctest: +FLOAT_CMP
    (-1.8205895517672226, 7.374044252723081, -0.20697960810978586)
    >>> sat.geo_pos()  # geographic (geodetic) position, lon/lat/alt  # doctest: +FLOAT_CMP
    <PyCoordGeodetic: 112.2146d, 28.5509d, 538.0173km>
    >>> sat.topo_pos()  # topocentric position, az/el/dist/dist_rate  # doctest: +FLOAT_CMP
    <PyCoordTopocentric: 60.2453d, -35.6845d, 8314.5681km, 3.5087km/s>

    >>> # One can change time to determine positions at another moment
    >>> sat.mjd += 1 / 720.  # one minute later
    >>> sat.topo_pos()  # doctest: +FLOAT_CMP
    <PyCoordTopocentric: 54.8446d, -38.2749d, 8734.9196km, 3.4885km/s>

In many cases, however, one probably wants to calculate coordinates for a
(large) number of satellites, observer locations, and/or observing times. For
this, the function `~cysgp4.propagate_many` is useful. This is an array
interface to the sgp4 calculations, which allows to perform calculations for
different satellite TLEs, observers and times in a parallelized manner.
`~numpy` broadcasting `rules
<https://docs.scipy.org/doc/numpy/user/basics.broadcasting.html>`_ apply::

        >>> import numpy as np
        >>> from cysgp4 import PyTle, PyObserver
        >>> from cysgp4 import get_example_tles, propagate_many

        >>> # Download many TLEs from a website
        >>> # import requests
        >>> # url = 'http://celestrak.com/NORAD/elements/science.txt'
        >>> # ctrak_science = requests.get(url).text

        >>> # or use built-in example TLE text file
        >>> ctrak_science = get_example_tles()

        >>> # Need to convert them to a list of tuples (each tuple consisting
        >>> # of the three TLE strings)
        >>> tle_list = cysgp4.tle_tuples_from_text(ctrak_science)
        >>> # Create an array of PyTle and PyObserver objects, and MJDs
        >>> tles = np.array([
        ...     PyTle(*tle) for tle in tle_list
        ...     ])[np.newaxis, np.newaxis, :20]  # use first 20 TLEs
        >>> observers = np.array([
        ...     PyObserver(6.88375, 50.525, 0.366),
        ...     PyObserver(16.88375, 50.525, 0.366),
        ...     ])[np.newaxis, :, np.newaxis]
        >>> mjds = np.linspace(
        ...     58805.5, 58806.5, 1000  # 1000 time steps
        ...     )[:, np.newaxis, np.newaxis]

        >>> # The result is a dictionary
        >>> result = propagate_many(mjds, tles, observers)
        >>> print(sorted(result.keys()))
        ['eci_pos', 'eci_vel', 'geo', 'topo']

        >>> # Returned array shapes are as follows; last array dimension
        >>> # contains the coordinate pairs.
        >>> print(np.broadcast(mjds, tles, observers).shape)
        (1000, 2, 20)
        >>> print(result['eci_pos'].shape, result['topo'].shape)
        (1000, 2, 20, 3) (1000, 2, 20, 4)

        >>> # One can also skip over coordinate frames.
        >>> result = propagate_many(
        ...     mjds, tles, observers,
        ...     do_eci_pos=False, do_eci_vel=False, do_geo=False, do_topo=True
        ...     )
        >>> print(result.keys())
        dict_keys(['topo'])

.. note::

    This "minimal example" is also contained in a fully-working `Jupyter
    tutorial notebook <https://github.com/bwinkel/cygrid/blob/master/notebooks/01_minimal_example.ipynb>`_.

