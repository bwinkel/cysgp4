#!/usr/bin/python
# -*- coding: utf-8 -*-
# cython: language_level=3
# cython: cdivision=True, boundscheck=False, wraparound=False
# cython: embedsignature=True

# ####################################################################
#
# title                  :cysgp4.pyx
# description            :Cython-powered wrapper of the sgp4lib library
# author                 :Benjamin Winkel
#
# ####################################################################
#  Copyright (C) 2014+ by Benjamin Winkel
#  bwinkel@mpifr.de
#  This file is part of cysgp4.
#
#  cysgp4 is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
#
# Note: cyaatm is a wrapper around sgp4lib library by Daniel Warner
#       (see package in cextern directory).
# ####################################################################

# import python3 compat modules
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

cimport cython
from cython.parallel import prange, parallel
cimport numpy as np
cimport openmp
from numpy cimport PyArray_MultiIter_DATA as Py_Iter_DATA
from cython.operator cimport dereference as deref
from cython.operator cimport address as addr
from cython.operator cimport preincrement as inc
from cpython cimport bool as python_bool
from libcpp cimport bool as cpp_bool
from libc.math cimport M_PI, floor, fabs, fmod, sqrt, sin, cos
from .cysgp4 cimport *

from datetime import datetime
import numpy as np

np.import_array()


cdef double NAN = np.nan
cdef double DEG2RAD = M_PI / 180.
cdef double RAD2DEG = 180. / M_PI
cdef double MJD_RESOLUTION = 0.001 / 24. / 3600.
cdef long long MJD0_TICKS = 58628880000000000


ctypedef SGP4* sgp4_ptr_t
ctypedef Tle* tle_ptr_t
ctypedef Observer* obs_ptr_t


__all__ = [
    'PyDateTime', 'PyTle', 'PyObserver',
    'PyCoordGeodetic', 'PyCoordTopocentric', 'PyEci',
    'Satellite', 'propagate_many', 'propagate_many_slow', 'set_num_threads',
    ]


def set_num_threads(int nthreads):
    '''
    Change maximum number of threads to use.

    Parameters
    ----------
    nthreads - int
        Number of threads to use.

    Notes
    -----
    - This can also be controlled by setting the environment variable
      `OMP_NUM_THREADS`.
    '''

    openmp.omp_set_num_threads(nthreads)


cdef inline long long ticks_from_mjd(double mjd) nogil:

    cdef:
        double days, fdays
        long long idays

    idays = <long long> mjd
    fdays = mjd - idays

    return (
        idays * 86400 * 1000000 +
        (<long long> (fdays * 86400 * 1000000)) +
        MJD0_TICKS
        )


cdef inline double mjd_from_ticks(long long ticks) nogil:

    cdef:
        double days, fdays
        long long idays

    ticks -= MJD0_TICKS
    return ticks / 8.64e10


cdef inline (double, double, double) ecef_from_geo(
        double lon_rad, double lat_rad, double alt_km
        ) nogil:
    '''
    Return ECEF (Cartesian) in km.

    https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
    '''

    cdef:

        double a = kXKMPER
        double b = kXKMPER * (1. - kF)

        double slam = sin(lon_rad)
        double clam = cos(lon_rad)
        double sphi = sin(lat_rad)
        double cphi = cos(lat_rad)

        double N = a ** 2 / sqrt(a ** 2 * cphi ** 2 + b ** 2 * sphi ** 2)

        double x = (N + alt_km) * cphi * clam
        double y = (N + alt_km) * cphi * slam
        double z = (b ** 2 / a ** 2 * N + alt_km) * sphi

    return x, y, z


cdef class PyDateTime(object):
    '''
    Thin wrapper around sgp4 (C++) DateTime class.

    Internally, DateTime objects store time on a linear scale, "ticks",
    which is the number of micro-seconds since 1. January 0001, 00:00.
    For convenience, one can convert to and from Python `datetime
    <https://docs.python.org/3.8/library/datetime.html>`_ objects and
    `Modified Julian Date (MJD)
    <https://en.wikipedia.org/wiki/Julian_day#Variants>`_, which is often
    used in astronomy.

    For calculating the position of a satellite for a given observer
    (i.e., in the horizontal frame, Azimuth and Elevation), one needs
    the `local mean sidereal time
    <https://en.wikipedia.org/wiki/Sidereal_time>`_  (derived from Greenwhich
    mean sidereal time), which can be queried using the
    `~cysgp4.PyDateTime.lmst()` and `~cysgp4.PyDateTime.gmst()` methods.

    Parameters
    ----------
    dt : `~datetime.datetime` object, or None (default: None)
        Construct the `~cysgp4.PyDateTime` instance with a given
        Python `~datetime.datetime` object, or - if None - use the
        current date and time (as given by `~datetime.datetime.now()`)
    init : Boolean (default: True)
        If set to true, the `~datetime.datetime` object will ignore the
        `dt` parameter and set the internal ticks variable to 0.
        This is only useful for computationally heavy tasks.

    Returns
    -------
    pydt : `~cysgp4.PyDateTime` object

    Examples
    --------
    There are several possibilities to create a `~cysgp4.PyDateTime` object::

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

    Accessing and Modifying the date and time is done via the following
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

    One can also set the date and time manually with the
    `~cysgp4.PyDateTime.set` method.
    '''

    # hold the C++ instance, which we're wrapping
    cdef DateTime _cobj

    def __init__(self, object dt=None, init=True):

        self._cobj = DateTime(0)
        if init:
            self._set_datetime(dt)

    @classmethod
    def from_ticks(cls, unsigned long long ticks):
        '''
        Creates a new `~cysgp4.PyDateTime` instance from "ticks".

        Internally, DateTime objects store time on a linear scale, "ticks",
        which is the number of micro-seconds since 1. January 0001, 00:00.

        Parameters
        ----------
        ticks : unsigned long long
            Number of micro-seconds since 1. January 0001, 00:00.

        Returns
        -------
        pydt : `~cysgp4.PyDateTime` object

        Examples
        --------
        "from_ticks" is a classmethod::

            >>> ticks = 58628880000000000  # aka MJD zero point
            >>> pydt = PyDateTime.from_ticks(ticks)
            >>> pydt
            <PyDateTime: 1858-11-17 00:00:00.000000 UTC>

        '''

        dt = cls(dt=None, init=False)
        dt.ticks = ticks

        return dt

    @classmethod
    def from_mjd(cls, double mjd):
        '''
        Creates a new `~cysgp4.PyDateTime` instance from MJD.

        `Modified Julian Date (MJD)
        <https://en.wikipedia.org/wiki/Julian_day#Variants>`_, is a
        timescale, which is often used in astronomy.

        Parameters
        ----------
        mjd : double
            Modified Julian Date.

        Returns
        -------
        pydt : `~cysgp4.PyDateTime` object

        Examples
        --------
        "from_mjd" is a classmethod::

            >>> mjd = 56458.123
            >>> pydt = PyDateTime.from_mjd(mjd)
            >>> pydt
            <PyDateTime: 2013-06-15 02:57:07.199999 UTC>

        '''

        dt = cls(dt=None, init=False)
        dt.mjd = mjd

        return dt

    def get_datetime_tuple(self):

        return (
            self._cobj.Year(),
            self._cobj.Month(),
            self._cobj.Day(),
            self._cobj.Hour(),
            self._cobj.Minute(),
            self._cobj.Second(),
            self._cobj.Microsecond(),
            )

    def _get_datetime(self):

        return datetime(*self.get_datetime_tuple())

    def _set_datetime(self, dt):
        '''
        Initialize PyDateTime from python datetime object
        '''

        if dt is None:
            dt = datetime.now()

        assert isinstance(dt, datetime)

        self._cobj.Initialise(
            <int> dt.year, <int> dt.month, <int> dt.day,
            <int> dt.hour, <int> dt.minute, <int> dt.second,
            <int> dt.microsecond
            )

    datetime = property(
        _get_datetime, _set_datetime, None,
        doc='datetime (see Class documentation).'
        )

    def set(
            self,
            int year, int month, int day,
            int hour, int minute, int second, int microsecond
            ):
        '''
        Creates a new `~cysgp4.PyDateTime` instance from given date and time.

        Parameters
        ----------
        year : int
            Year of desired date.
        month : int
            Month of desired date.
        day : int
            Day of desired date.
        hour : int
            Hours of desired time.
        minute : int
            Minutes of desired time.
        second : int
            Seconds of desired time.
        microsecond : int
            Microseconds of desired time.
        '''

        self._cobj.Initialise(
            <int> year, <int> month, <int> day,
            <int> hour, <int> minute, <int> second, <int> microsecond
            )

    def _get_ticks(self):
        return <long long> self._cobj.Ticks()

    def _set_ticks(self, unsigned long long ticks):

        # this is a bit ugly, but there is to setter method in the C++ code

        cdef:
            long long ticks_new = ticks
            long long ticks_old = self._get_ticks()

        # AddTicks returns a new instance...
        self._cobj = self._cobj.AddTicks(ticks_new - ticks_old)

    ticks = property(
        _get_ticks, _set_ticks, None, doc='Ticks (see Class documentation).'
        )

    def _get_mjd(self):
        return mjd_from_ticks(self._cobj.Ticks())

    def _set_mjd(self, double mjd):

        # this is a bit ugly, but there is to setter method in the C++ code

        cdef:
            long long ticks_new = ticks_from_mjd(mjd)
            long long ticks_old = self._get_ticks()

        # AddTicks returns a new instance...
        self._cobj = self._cobj.AddTicks(ticks_new - ticks_old)

    mjd = property(
        _get_mjd, _set_mjd, None, doc='MJD (see Class documentation).'
        )

    def __str__(self):

        return self._cobj.ToString().decode('UTF-8')

    def __repr__(self):

        return '<PyDateTime: ' + self.__str__() + '>'

    def gmst(self):
        '''
        Greenwhich mean sidereal time (GMST) of current date/time.

        Returns
        -------
        gmst : float
            Greenwhich mean sidereal time (GMST) in (fractional) hours.
        '''
        return self._cobj.ToGreenwichSiderealTime()

    def lmst(self, obslon_deg):
        '''
        Local mean sidereal time (LMST) of current date/time.

        Parameters
        ----------
        obslon_deg : float
            Geographic longitude of the observer.

        Returns
        -------
        lmst : float
            Local mean sidereal time (LMST) in (fractional) hours.
        '''

        return self._cobj.ToLocalMeanSiderealTime(DEG2RAD * obslon_deg)


cdef class PyTle(object):
    '''
    Thin wrapper around sgp4 (C++) Tle class.

    `Two-line element sets
    <https://en.wikipedia.org/wiki/Two-line_element_set>`_ contain the
    parameters necessary to describe a satellite's orbit (around Earth).
    Although the name suggests that two lines of strings are used,
    in practice often three lines are defined, the first containing a
    satellite name (and/or) ID. It is important to note that for many
    satellites, the corresponding TLEs get outdated quickly. Especially for
    low earth orbit satellites, one should query up-to-date information
    at least once per day.

    `~cysgp4.PyTle` is only used to parse the TLE strings and store the
    orbital parameters for use by other sgp4 routines.

    Parameters
    ----------
    name : str
        Name/ID of satellite.
    line_one : str
        First TLE line string.
    line_two : str
        Second TLE line string.

    Returns
    -------
    tle : `~cysgp4.PyTle` object

    Examples
    --------
    TLEs can be obtained from many sources, such as
    `<https://celestrak.com/>`_::

        >>> import requests
        >>> from cysgp4 import PyTle

        >>> url = 'http://celestrak.com/NORAD/elements/science.txt'
        >>> ctrak_science = requests.get(url)
        >>> all_lines = ctrak_science.text.split('\\r\\n')
        >>> tle_list = list(zip(*tuple(
        ...     all_lines[idx::3]
        ...     for idx in range(3)
        ...     )))
        >>> len(tle_list)
        >>> print(*tle_list[1], sep='\\n')
        HST
        1 20580U 90037B   19321.38711875  .00000471  00000-0  17700-4 0  9991
        2 20580  28.4699 288.8102 0002495 321.7771 171.5855 15.09299865423838

        >>> hst_tle = PyTle(*tle_list[1])
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

    '''

    # hold the C++ instance, which we're wrapping
    cdef Tle *thisptr

    def __init__(self, name, line_one, line_two):

        # Note: as there is no empty constructor ("Tle()") in the sgp4
        # library, we need to instantiate the object on the heap.
        # This is because Cython can not (yet?) do something like
        # my_tle = Tle(*args)
        # without the empty constructor being available, even though it
        # is not used.
        self.thisptr = new Tle(
            name.encode('UTF-8'),
            line_one.encode('UTF-8'),
            line_two.encode('UTF-8')
            )

    def __dealloc__(self):

        del self.thisptr

    def __str__(self):

        return self.thisptr.ToString().decode('UTF-8')

    def __repr__(self):

        return '<PyTle: ' + self.thisptr.Name().decode('UTF-8') + '>'

    def tle_strings(self):

        return (
            self.thisptr.Name().decode('UTF-8'),
            self.thisptr.Line1().decode('UTF-8'),
            self.thisptr.Line2().decode('UTF-8'),
            )


cdef class PyCoordGeodetic(object):
    '''
    Thin wrapper around sgp4 (C++) CoordGeodetic struct.

    The CoordGeodetic struct holds a geographic location (latitude,
    longitude, altitude).

    Parameters
    ----------
    lon_deg : float (default: 0.)
        Longitude of geographic location [deg].
    lat_deg : float (default: 0.)
        Latitude of geographic location [deg].
    alt_km : float (default: 0.)
        Altitude of geographic location [km].

    Returns
    -------
    geo : `~cysgp4.PyCoordGeodetic` object

    Examples
    --------
    Constructing and using a `~cysgp4.PyCoordGeodetic` object is
    straightforward::

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

    '''

    # hold the C++ instance, which we're wrapping
    cdef CoordGeodetic _cobj

    def __init__(
            self,
            double lon_deg=0,
            double lat_deg=0,
            double alt_km=0
            ):

        self._cobj = CoordGeodetic(
            lat_deg, lon_deg, alt_km
            )

    def __str__(self):

        return ', '.join([
            '{:.4f}d'.format(self.lon),
            '{:.4f}d'.format(self.lat),
            '{:.4f}km'.format(self.alt),
            ])

    def __repr__(self):

        return '<PyCoordGeodetic: ' + self.__str__() + '>'

    def _get_lon(self):

        return RAD2DEG * self._cobj.longitude

    def _set_lon(self, double lon_deg):

        self._cobj.longitude = DEG2RAD * lon_deg

    def _get_lat(self):

        return RAD2DEG * self._cobj.latitude

    def _set_lat(self, double lat_deg):

        self._cobj.latitude = DEG2RAD * lat_deg

    def _get_alt(self):

        return self._cobj.altitude

    def _set_alt(self, double alt_km):

        self._cobj.altitude = alt_km

    def _get_ecef(self):

        cdef double x, y, z

        (x, y, z) = ecef_from_geo(
            self._cobj.longitude,
            self._cobj.latitude,
            self._cobj.altitude,
            )

        return x, y, z

    lon = property(
        _get_lon, _set_lon, None,
        doc='Geographic longitude [deg] (see also Class documentation).'
        )
    lat = property(
        _get_lat, _set_lat, None,
        doc='Geographic latitude [deg] (see also Class documentation).'
        )
    alt = property(
        _get_alt, _set_alt, None,
        doc='Geographic altitude [km] (see also Class documentation).'
        )
    ecef = property(
        _get_ecef, None, None,
        doc='ECEF [km] (Earth-centered, Earth-fixed frame; x, y, z).'
        )


cdef class PyCoordTopocentric(object):
    '''
    Thin wrapper around sgp4 (C++) CoordTopocentric struct.

    The CoordTopocentric struct holds a topocentric location
    (azimuth, elevation, range/distance and distance/range rate).

    Note: the topocentric position of a satellite is always relative to
    an observer (in SGP4 the observer is defined in geographic coordinates).
    The distance and distance change (aka range rate) is thus the distance
    between the satellite and the observer. However, this struct/class
    only holds the azimuth, elevation, distance and distance rate parameters,
    but contains no information on the observer. It is only useful in
    conjuction with a `~cysgp4.PyObserver` (which holds a reference to a
    geographic location) and a datetime; see `~cysgp4.Satellite`.

    Parameters
    ----------
    az_deg : float (default: 0.)
        Azimuth of topocentric location [deg].
    el_deg : float (default: 0.)
        Elevation of topocentric location [deg].
    dist_km : float (default: 0.)
        Distance/range of topocentric location [km].
    dist_rate_km_per_s : float (default: 0.)
        Distance/range rate of topocentric location [km/s].

    Returns
    -------
    topo : `~cysgp4.PyCoordTopocentric` object

    Examples
    --------
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

    '''

    # hold the C++ instance, which we're wrapping
    cdef CoordTopocentric _cobj

    def __init__(
            self,
            double az_deg=0,
            double el_deg=0,
            double dist_km=0,
            double dist_rate_km_per_s=0,
            ):

        self._cobj = CoordTopocentric(
            az_deg * DEG2RAD, el_deg * DEG2RAD, dist_km, dist_rate_km_per_s
            )

    def __str__(self):

        return ', '.join([
            '{:.4f}d'.format(self.az),
            '{:.4f}d'.format(self.el),
            '{:.4f}km'.format(self.dist),
            '{:.4f}km/s'.format(self.dist_rate),
            ])

    def __repr__(self):

        return '<PyCoordTopocentric: ' + self.__str__() + '>'

    def _get_az(self):

        return RAD2DEG * self._cobj.azimuth

    def _set_az(self, double az_deg):

        self._cobj.azimuth = DEG2RAD * az_deg

    def _get_el(self):

        return RAD2DEG * self._cobj.elevation

    def _set_el(self, double el_deg):

        self._cobj.elevation = DEG2RAD * el_deg

    def _get_dist(self):

        return self._cobj.distance

    def _set_dist(self, double dist_km):

        self._cobj.distance = dist_km

    def _get_dist_rate(self):

        return self._cobj.distance_rate

    def _set_dist_rate(self, double dist_rate_km_per_s):

        self._cobj.distance_rate = dist_rate_km_per_s

    az = property(
        _get_az, _set_az, None,
        doc='Topocentric azimuth [deg] (see also Class documentation).'
        )
    el = property(
        _get_el, _set_el, None,
        doc='Topocentric elevation [deg] (see also Class documentation).'
        )
    dist = property(
        _get_dist, _set_dist, None,
        doc='Topocentric distance [km] (see also Class documentation).'
        )
    dist_rate = property(
        _get_dist_rate, _set_dist_rate, None,
        doc='Topocentric distance rate [km/s] (see also Class documentation).'
        )


cdef class PyObserver(object):
    '''
    Thin wrapper around sgp4 (C++) Observer class.

    The Observer class holds the location (as ECI, see `~cysgp4.PyEci`) of
    an observer.

    Note: Usually, a datetime is attached to an ECI location in sgp4.
    However, the Observer class has no interface to set the datetime
    and internally it is always set to zero. This is an odd design choice,
    because a geographic location is not fixed in the ECI system if
    time passes.

    Parameters
    ----------
    lon_deg : float (default: 0.)
        Longitude of observer (geographic location) [deg].
    lat_deg : float (default: 0.)
        Latitude of observer (geographic location) [deg].
    alt_km : float (default: 0.)
        Altitude of observer (geographic location) [km].

    Returns
    -------
    obs : `~cysgp4.PyObserver` object

    Examples
    --------
    Constructing and using a `~cysgp4.PyObserver` object is
    straightforward::

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

    '''

    # hold the C++ instance, which we're wrapping
    cdef:
        Observer _cobj

    def __init__(
            self,
            double lon_deg=0,
            double lat_deg=0,
            double alt_km=0,
            ):
        '''
        Constructor PyObserver(double lon_deg, double lat_deg, double alt_km)
        '''

        _obs_loc = PyCoordGeodetic(
            lon_deg=lon_deg,
            lat_deg=lat_deg,
            alt_km=alt_km
            )
        self._cobj = Observer(_obs_loc._cobj)

    def __str__(self):

        return self._get_location().__str__()

    def __repr__(self):

        return '<PyObserver: ' + self.__str__() + '>'

    def _get_location(self):

        _obs_loc = PyCoordGeodetic()
        _obs_loc._cobj = self._cobj.GetLocation()
        return _obs_loc

    def _set_location(self, PyCoordGeodetic loc):

        self._cobj.SetLocation(loc._cobj)

    loc = property(
        _get_location, _set_location, None,
        doc='Geographic location (see also Class documentation).'
        )


cdef class PyEci(object):
    '''
    Thin wrapper around sgp4 (C++) Eci class.

    The Eci class holds an `ECI location
    <https://en.wikipedia.org/wiki/Earth-centered_inertial>`_
    (latitude, longitude, altitude) for a particular datetime.

    Note, internally, the coordinates (and velocities) are stored in
    Cartesian form (read-only!). One can access these, via the
    `~cysgp4.PyEci.loc` and `~cysgp4.PyEci.vel` properties. Setting
    new parameters is only possible via a geographic location.

    Parameters
    ----------
    pydt : `~cysgp4.PyDateTime`
        Date and time.
    geo_loc : `~cysgp4.PyCoordGeodetic`
        Geographic location.

    Returns
    -------
    eci : `~cysgp4.PyEci` object

    Examples
    --------
    Constructing and using a `~cysgp4.PyEci` object is
    straightforward::

        >>> from cysgp4 import PyEci, PyCoordGeodetic, PyDateTime

        >>> pydt = PyDateTime.from_mjd(55555.)
        >>> lon_deg, lat_deg = 6.88375, 50.525
        >>> alt_km = 0.366
        >>> geo = PyCoordGeodetic(lon_deg, lat_deg, alt_km)
        >>> eci = PyEci(pydt, geo)
        >>> eci
        <PyEci: 6.8837d, 50.5250d, 0.3660km 2010-12-25 00:00:00.000000 UTC>

        >>> # Access is also possible via properties, e.g.:
        >>> eci.loc
        (-725.3304166274728, 3997.924210010933, 4900.402205553537)
        >>> eci.pydt = PyDateTime.from_mjd(55556.)

        >>> # or the update method:
        >>> eci.update(pydt, PyCoordGeodetic(0, 0, 0))

    '''
    cdef:
        # hold the C++ instance, which we're wrapping
        Eci _cobj

    def __init__(self, PyDateTime pydt=None, PyCoordGeodetic geo_loc=None):
        '''
        Constructor PyEci(PyDateTime pydt, PyCoordGeodetic geo_loc)
        '''
        self.update(pydt, geo_loc)

    def __str__(self):

        return self._get_geo_loc().__str__() + ' ' + self._get_dt().__str__()

    def __repr__(self):

        return '<PyEci: ' + self.__str__() + '>'

    def update(self, PyDateTime pydt=None, PyCoordGeodetic geo_loc=None):
        '''
        Update `~cysgp4.PyEci` object.

        Parameters
        ----------
        pydt : `~cysgp4.PyDateTime`
            Date and time.
        geo_loc : `~cysgp4.PyCoordGeodetic`
            Geographic location.
        '''

        if pydt is None:

            pydt = PyDateTime()

        if geo_loc is None:

            geo_loc = PyCoordGeodetic()

        self._cobj = Eci(pydt._cobj, geo_loc._cobj)

    def _get_loc(self):

        cdef:
            Vector _pos = self._cobj.Position()

        return _pos.x, _pos.y, _pos.z

    def _get_vel(self):

        cdef:
            Vector _vel = self._cobj.Velocity()

        return _vel.x, _vel.y, _vel.z

    def _get_geo_loc(self):

        geo_loc = PyCoordGeodetic()
        geo_loc._cobj = self._cobj.ToGeodetic()
        return geo_loc

    def _set_geo_loc(self, PyCoordGeodetic geo_loc):

        self.update(self.pydt, geo_loc)

    def _get_dt(self):

        pydt = PyDateTime()
        pydt._cobj = self._cobj.GetDateTime()
        return pydt

    def _set_dt(self, PyDateTime pydt):

        self.update(pydt, self.geo_loc)

    loc = property(
        _get_loc, None, None,
        doc='Cartesian location (readonly, see also Class documentation).'
        )
    vel = property(
        _get_vel, None, None,
        doc='Cartesian velocity (readonly, see also Class documentation).'
        )
    geo_loc = property(
        _get_geo_loc, _set_geo_loc, None,
        doc='Geographic location (see also Class documentation).'
        )
    pydt = property(
        _get_dt, _set_dt, None,
        doc='Datetime (see also Class documentation).'
        )


cdef class Satellite(object):
    '''
    Calculate position of a satellite at a given time.

    The satellite is defined via a TLE (see `~cysgp4.PyTle`). Furthermore,
    for apparent positions of the satellite, an observer (see
    `~cysgp4.PyObserver`) needs to be defined.

    The position calculations are lazy, i.e., only calculated if the newly
    requested time differs by a certain amount from the last time at which
    a calculation was performed. The granularity of this "cache" can be
    defined via the `mjd_cache_resolution` parameter.

    Parameters
    ----------
    tle : `~cysgp4.PyTle`
        TLE instance of the satellite of interest.
    obs : `~cysgp4.PyObserver` or None (default: None)
        Observer instance. If `None` then the observer location is set to
        (0 deg, 0 deg, 0 km).
    pydt : `~cysgp4.PyDateTime` or `None` (default: None)
        Date and time at which to calculate the position. If `None` then
        the current datetime is used, as given by a call to
        `~datetime.datetime.now()`.
    mjd_cache_resolution : double (default: 0.001)
        Granularity of the internal cache [s]. Satellite positions are
        calculated only if the newly requested datetime differs by more than
        this from the previous calculation.
    on_error : str, optional (either 'raise' or 'coerce_to_nan', default: 'raise')
        If the underlying SGP C++ library throws an error (which often
        happens if one works with times that are strongly deviating from
        the TLE epoch), a Python RuntimeError is usually thrown. For batch
        processing, this is not always desirable. If
        `on_error = 'coerce_to_nan'` then C++ errors will be suppressed and
        the resulting position vectors will be converted to NaN-values
        instead.

    Returns
    -------
    sat : `~cysgp4.Satellite` object

    Examples
    --------
    The following demonstrates how a typical use of the `~cysgp4.Satellite`
    class would look like::

        >>> from cysgp4 import *

        >>> pydt = PyDateTime.from_mjd(58805.57)
        >>> lon_deg, lat_deg = 6.88375, 50.525
        >>> alt_km = 0.366
        >>> obs = PyObserver(lon_deg, lat_deg, alt_km)

        >>> hst_tle = PyTle(
        ... 'HST',
        ... '1 20580U 90037B   19321.38711875  .00000471  00000-0  17700-4 0  9991',
        ... '2 20580  28.4699 288.8102 0002495 321.7771 171.5855 15.09299865423838',
        ... )

        >>> sat = Satellite(hst_tle, obs, pydt)
        >>> # can now query positions, also for different times
        >>> sat.eci_pos().loc  # ECI cartesian position
        (5879.5931344459295, 1545.7455647032068, 3287.4155452595)
        >>> sat.eci_pos().vel  # ECI cartesian velocity
        (-1.8205895517672226, 7.374044252723081, -0.20697960810978586)
        >>> sat.geo_pos()  # geographic position
        <PyCoordGeodetic: 112.2146d, 28.5509d, 538.0186km>
        >>> sat.topo_pos()  # topocentric position
        <PyCoordTopocentric: 60.2453d, -35.6844d, 8314.5683km, 3.5087km/s>

        >>> # change time
        >>> sat.mjd += 1 / 720.  # one minute later
        >>> sat.topo_pos()
        <PyCoordTopocentric: 54.8446d, -38.2749d, 8734.9195km, 3.4885km/s>

        >>> # change by less than cache resolution (1 ms)
        >>> sat.topo_pos().az, sat.topo_pos().el
        (54.84463503781068, -38.274852915850126)
        >>> sat.mjd += 0.0005 / 86400.  # 0.5 ms
        >>> sat.topo_pos().az, sat.topo_pos().el
        (54.84463503781068, -38.274852915850126)
        >>> # change by another 0.5 ms triggers re-calculation
        >>> sat.mjd += 0.00051 / 86400.
        >>> sat.topo_pos().az, sat.topo_pos().el
        (54.844568313870965, -38.274885794151324)

    '''

    cdef:
        # hold C++ instances, which we're wrapping
        SGP4 *thisptr

        PyTle _tle
        PyObserver _obs

        PyDateTime _pydt
        PyEci _eci
        PyCoordTopocentric _topo
        PyCoordGeodetic _geo

        # _cmjd holds the mjd for the last calculation time
        double _cmjd, _mjd_cache_resolution
        python_bool _pos_dirty, _tle_dirty, on_error_raise

    def __init__(
            self,
            PyTle tle,
            PyObserver obs=None,
            PyDateTime pydt=None,
            double mjd_cache_resolution=MJD_RESOLUTION,
            str on_error='raise',
            ):
        '''
        Constructs a new Satellite object from given TLE
        '''

        assert on_error in ['raise', 'coerce_to_nan']
        self.on_error_raise = True if on_error == 'raise' else False

        if obs is None:

            obs = PyObserver()

        if pydt is None:

            pydt = PyDateTime()

        self._mjd_cache_resolution = mjd_cache_resolution

        # want deep copies to avoid side effects!
        self._tle = PyTle(
            tle.thisptr.Name().decode('UTF-8'),
            tle.thisptr.Line1().decode('UTF-8'),
            tle.thisptr.Line2().decode('UTF-8'),
            )

        # important: the class members are not yet initialized at all???
        self._obs = PyObserver()
        self._obs._cobj.SetLocation(obs._cobj.GetLocation())
        self._pydt = PyDateTime(init=False)
        self._pydt._cobj = DateTime(pydt._cobj.Ticks())
        self._cmjd = self._pydt.mjd

        self.thisptr = new SGP4(deref(tle.thisptr))
        self._tle_dirty = <python_bool> False

        # the following is important, otherwise self._topo and self._geo will
        # just be None after _refresh_coords()

        # initialize workspaces, otherwise we cannot assign values in
        # _refresh_coords (would produce segfault)
        # note: it is not sufficient to define these as class members (above)
        self._eci = PyEci()
        self._topo = PyCoordTopocentric()
        self._geo = PyCoordGeodetic()

        self._pos_dirty = <python_bool> True

    def __dealloc__(self):

        # del self.eci_ptr
        del self.thisptr
        # del self._tle.thisptr

    def _get_mjd(self):

        return self._pydt.mjd

    def _set_mjd(self, double mjd):

        assert mjd < 1000000., 'warning, make sure to use mjd'

        self._pydt.mjd = mjd

        if fabs(self._pydt.mjd - self._cmjd) > self._mjd_cache_resolution:

            self._pos_dirty = <python_bool> True

    mjd = property(
        _get_mjd, _set_mjd, None,
        doc='Modified Julian Date (see also Class documentation).'
        )

    def _get_datetime(self):

        return self._pydt

    def _set_datetime(self, pydt):

        # must use the set_mjd method otherwise the caching would not work
        self._set_mjd(pydt.mjd)

    pydt = property(
        _get_datetime, _set_datetime, None,
        doc='Datetime (see also Class documentation).'
        )

    def topo_pos(self):
        '''
        Topocentric position of satellite w.r.t. given observer.

        Lazily calculated, i.e., only if difference to the timestamp
        for which the previous calculation was performed is larger
        than the MJD cache resolution (see class documentation).

        Returns
        ----------
        topo : `~cysgp4.PyCoordTopocentric`
            Topocentric position of satellite w.r.t. given observer.
        '''

        if self._pos_dirty:
            self._refresh_coords()

        if self._tle_dirty:
            return None

        return self._topo

    def geo_pos(self):
        '''
        Geographic (geodetic) position of satellite.

        Lazily calculated, i.e., only if difference to the timestamp
        for which the previous calculation was performed is larger
        than the MJD cache resolution (see class documentation).

        Returns
        ----------
        geo : `~cysgp4.PyCoordGeodetic`
            Geographic position of satellite.
        '''

        if self._pos_dirty:
            self._refresh_coords()

        if self._tle_dirty:
            return None

        return self._geo

    def eci_pos(self):
        '''
        ECI position of satellite.

        Lazily calculated, i.e., only if difference to the timestamp
        for which the previous calculation was performed is larger
        than the MJD cache resolution (see class documentation).

        Returns
        ----------
        eci : `~cysgp4.PyEci`
            ECI position of satellite.
        '''

        if self._pos_dirty:
            self._refresh_coords()

        if self._tle_dirty:
            return None

        return self._eci

    def _refresh_coords(self):

        # FindPosition doesn't update ECI time, need to do manually :-/
        self._eci = PyEci(pydt=self._pydt)
        if self.on_error_raise:
            self._eci._cobj = self.thisptr.FindPosition(self._pydt._cobj)
        else:
            self._eci._cobj = self.thisptr.FindPositionNaN(self._pydt._cobj)

        self._tle_dirty = <python_bool> False

        self._topo._cobj = self._obs._cobj.GetLookAngle(
            self._eci._cobj
            )
        self._geo._cobj = self._eci._cobj.ToGeodetic()

        self._cmjd = self._pydt.mjd
        self._pos_dirty = <python_bool> False


def propagate_many(
        mjds, tles, observers=None,
        bint do_eci_pos=True, bint do_eci_vel=True,
        bint do_geo=True, bint do_topo=True,
        str on_error='raise',
        ):
    '''
    Calculate positions of many satellites at a various times at once.

    This is an array interface to the sgp4 calculations, which allows to
    perform calculations for different satellite TLEs, observers and times
    in a parallelized manner. `~numpy` broadcasting rules apply.

    With the Boolean parameters, `do_eci_pos`, `do_eci_vel`, `do_geo`, and
    `do_topo` the user can decide, which position frames are returned
    in the output dictionary.

    Satellite are defined via TLEs (see `~cysgp4.PyTle`). The
    `~.cysgp4.propagate_many` function works with a single (scalar) PyTle
    object or a list/array of them. The same is true for the `mjds` and
    optional `observers` parameters, which must be a scalar (a double and
    `~cysgp4.PyObserver` instance, respectively) or a list/array of them.

    As in most use cases, a large number of positions is probably queried,
    the returned values do not use the `~cysgp4.PyEci` or
    `~cysgp4.PyCoordTopocentric` classes, but pack everything in arrays.
    This is not only faster, but makes it easier to further process the
    results.

    For parallelization, `OpenMP <https://www.openmp.org/>`_ is utilized.
    To change the number of CPU cores that are used, one can either set
    the environment variable `OMP_NUM_THREADS` or use
    `~cysgp4.set_num_threads`.

    Parameters
    ----------
    mjds : `~numpy.ndarray`, `~list`, or scalar of float
        Modified Julian Date.
    tles : `~numpy.ndarray`, `~list`, or scalar of `~cysgp4.PyTle`
        TLE instance of the satellite of interest.
    observers : `~numpy.ndarray`, `~list`, or scalar of `~cysgp4.PyObserver` or None (default: None)
        Observer instance. If `None` then the observer location is set to
        (0 deg, 0 deg, 0 km).
    do_eci_pos : Boolean, optional (default: True)
        Whether to include ECI cartesian positions in the result.
    do_eci_vel : Boolean, optional (default: True)
        Whether to include ECI cartesian velocity in the result.
    do_geo : Boolean, optional (default: True)
        Whether to include geographic/geodetic positions in the result.
    do_topo : Boolean, optional (default: True)
        Whether to include topocentric positions in the result.
    on_error : str, optional (either 'raise' or 'coerce_to_nan', default: 'raise')
        If the underlying SGP C++ library throws an error (which often
        happens if one works with times that are strongly deviating from
        the TLE epoch), a Python RuntimeError is usually thrown. For batch
        processing, this is not always desirable. If
        `on_error = 'coerce_to_nan'` then C++ errors will be suppressed and
        the resulting position vectors will be converted to NaN-values
        instead.

    Returns
    -------
    result : dict
        Resulting positions for each requested frame:

        - `eci_pos` : `~numpy.ndarray` of float
          Satellites ECI cartesian positions. Last dimension has length 3,
          one for each of `x`, `y`, and `z`. First dimensions are
          determined by the (broadcasted) shape of the inputs `mjd`,
          `tles`, and `observers`.

        - `eci_vel` : `~numpy.ndarray` of float
          Satellites ECI cartesian velicities. Last dimension has length 3,
          one for each of `v_x`, `v_y`, and `v_z`. First dimensions are
          determined by the (broadcasted) shape of the inputs `mjd`,
          `tles`, and `observers`.

        - `geo` : `~numpy.ndarray` of float
          Satellites Geodetic positions. Last dimension has length 3, one
          for each of `lon`, `lat`, and `alt`. First dimensions are
          determined by the (broadcasted) shape of the inputs `mjd`,
          `tles`, and `observers`.

        - `topo` : `~numpy.ndarray` of float
          Satellites Topocentric positions. Last dimension has length 4,
          one for each of `az`, `el`, `dist`, and `dist_rate`. First
          dimensions are determined by the (broadcasted) shape of the
          inputs `mjd`, `tles`, and `observers`.

    Examples
    --------
    The following demonstrates how to use the `~cysgp4.propagate_many`
    function::

        >>> import requests
        >>> import numpy as np
        >>> from cysgp4 import PyTle, PyObserver, propagate_many

        >>> url = 'http://celestrak.com/NORAD/elements/science.txt'
        >>> ctrak_science = requests.get(url)
        >>> all_lines = ctrak_science.text.split('\\r\\n')

        >>> tle_list = list(zip(*tuple(
        ...     all_lines[idx::3] for idx in range(3)
        ...     )))
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

        >>> result = propagate_many(mjds, tles, observers)
        >>> print(result.keys())
        dict_keys(['eci_pos', 'eci_vel', 'geo', 'topo'])

        >>> # shapes are as follows
        >>> print(np.broadcast(mjds, tles, observers).shape)
        (1000, 2, 20)
        >>> print(result['eci_pos'].shape, result['topo'].shape)
        (1000, 2, 20, 3) (1000, 2, 20, 4)

        >>> result = propagate_many(
        ...     mjds, tles, observers,
        ...     do_eci_pos=False, do_eci_vel=False, do_geo=False, do_topo=True
        ...     )
        >>> print(result.keys())
        dict_keys(['topo'])

    '''

    cdef:

        SGP4 *_sgp4_ptr
        Tle *_tle_ptr
        Observer _obs
        # PyEci py_nan_eci = PyEci(geo_loc=PyCoordGeodetic(NAN, NAN, NAN))
        Eci _eci
        # Eci _nan_eci = py_nan_eci._cobj
        DateTime _dt
        CoordGeodetic _geo
        CoordTopocentric _topo
        Vector _eci_pos, _eci_vel

        np.ndarray[double] mjd
        np.ndarray[object] tle, obs
        double[:, ::1] eci_pos_v
        double[:, ::1] eci_vel_v
        double[:, ::1] geo_v  # geodetic
        double[:, ::1] topo_v  # topocentric
        int i, size, n, pnum

        tle_ptr_t* _tle_ptr_array
        obs_ptr_t* _obs_ptr_array

        bint on_error_raise = True

    assert on_error in ['raise', 'coerce_to_nan']

    if on_error == 'coerce_to_nan':
        on_error_raise = False

    if observers is None:
        observers = PyObserver()

    # TODO: allocate only those arrays, which are actually requested
    # ("do_topo" etc.)

    pnum = 0
    out_dts = []
    if do_eci_pos:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    if do_eci_vel:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    if do_geo:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    if do_topo:
        out_dts.append(np.dtype(('float64', 4)))
        pnum += 1

    it = np.nditer(
        [tles, observers, mjds] + [None] * pnum,
        flags=['external_loop', 'buffered', 'delay_bufalloc', 'refs_ok'],
        op_flags=[['readonly']] * 3 + [['readwrite', 'allocate']] * pnum,
        op_dtypes=['object', 'object', 'float64'] + out_dts
        )

    # it would be better to use the context manager but
    # "with it:" requires numpy >= 1.14

    it.reset()

    for itup in it:

        tle = itup[0]
        obs = itup[1]
        mjd = itup[2]

        n = 3
        if do_eci_pos:
            eci_pos_v = itup[n]
            n += 1

        if do_eci_vel:
            eci_vel_v = itup[n]
            n += 1

        if do_geo:
            geo_v = itup[n]
            n += 1

        if do_topo:
            topo_v = itup[n]
            n += 1

        size = mjd.shape[0]
        _tle_ptr_array = array_new[tle_ptr_t](size)
        _obs_ptr_array = array_new[obs_ptr_t](size)

        for i in range(size):
            # unfortunately, it is not possible in nogil loop to access
            # the cdef'ed class members; therefore, we have to maintain
            # arrays of pointers to the sgp4 and observer objects
            _tle_ptr_array[i] = (<PyTle> tle[i]).thisptr
            _obs_ptr_array[i] = &(<PyObserver> obs[i])._cobj

        for i in prange(size, nogil=True):

            _tle_ptr = _tle_ptr_array[i]
            # _obs_ptr = _obs_ptr_array[i]

            # AddTicks returns a new instance...
            _dt = _dt.AddTicks(ticks_from_mjd(mjd[i]) - _dt.Ticks())

            # it is mandatory to create SGP4 instance here, as the
            # "FindPosition" method heavily changes class members
            # this could create problems with parallelization if instances
            # where re-used!
            _sgp4_ptr = new SGP4(deref(_tle_ptr))
            if on_error_raise:
                _eci = _sgp4_ptr.FindPosition(_dt)
            else:
                # "FindPositionE" will return NaNs if an error occured
                # at C++ level
                _eci = _sgp4_ptr.FindPositionNaN(_dt)

            del _sgp4_ptr

            if do_eci_pos:
                _eci_pos = _eci.Position()
                eci_pos_v[i, 0] = _eci_pos.x
                eci_pos_v[i, 1] = _eci_pos.y
                eci_pos_v[i, 2] = _eci_pos.z

            if do_eci_vel:
                _eci_vel = _eci.Velocity()
                eci_vel_v[i, 0] = _eci_vel.x
                eci_vel_v[i, 1] = _eci_vel.y
                eci_vel_v[i, 2] = _eci_vel.z

            if do_geo:
                _geo = _eci.ToGeodetic()
                geo_v[i, 0] = _geo.longitude * RAD2DEG
                geo_v[i, 1] = _geo.latitude * RAD2DEG
                geo_v[i, 2] = _geo.altitude

            if do_topo:
                _obs = Observer(_obs_ptr_array[i].GetLocation())
                _topo = _obs.GetLookAngle(_eci)
                topo_v[i, 0] = _topo.azimuth * RAD2DEG
                topo_v[i, 1] = _topo.elevation * RAD2DEG
                topo_v[i, 2] = _topo.distance
                topo_v[i, 3] = _topo.distance_rate

        array_delete(_tle_ptr_array)
        array_delete(_obs_ptr_array)

    result = {}

    n = 3
    if do_eci_pos:
        result['eci_pos'] = it.operands[n]
        n += 1

    if do_eci_vel:
        result['eci_vel'] = it.operands[n]
        n += 1

    if do_geo:
        result['geo'] = it.operands[n]
        n += 1

    if do_topo:
        result['topo'] = it.operands[n]
        n += 1

    return result


def propagate_many_slow(
        mjds, tles, observers=None,
        bint do_eci_pos=True, bint do_eci_vel=True,
        bint do_geo=True, bint do_topo=True,
        str on_error='raise',
        ):
    '''
    This is a slow (non-parallelized, Python-looping) version of
    `~cysgp4.propagate_many` that is meant for testing and benchmarking
    only. It has the same interface.
    '''

    assert on_error in ['raise', 'coerce_to_nan']

    pnum = 0
    out_dts = []
    if do_eci_pos:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    if do_eci_vel:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    if do_geo:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    if do_topo:
        out_dts.append(np.dtype(('float64', 4)))
        pnum += 1

    it = np.nditer(
        [tles, observers, mjds] + [None] * pnum,
        flags=['external_loop', 'buffered', 'delay_bufalloc', 'refs_ok'],
        op_flags=[['readonly']] * 3 + [['readwrite', 'allocate']] * pnum,
        op_dtypes=['object', 'object', 'float64'] + out_dts
        )

    it.reset()
    for itup in it:

        tle = itup[0]
        obs = itup[1]
        mjd = itup[2]

        size = mjd.shape[0]
        for i in range(size):

            sat = Satellite(
                tle[i], obs[i], PyDateTime.from_mjd(mjd[i]),
                on_error=on_error
                )
            eci = sat.eci_pos()

            eci_pos_x, eci_pos_y, eci_pos_z = eci.loc
            eci_vel_x, eci_vel_y, eci_vel_z = eci.vel

            n = 3
            if do_eci_pos:
                itup[n][i][0] = eci_pos_x
                itup[n][i][1] = eci_pos_y
                itup[n][i][2] = eci_pos_z
                n += 1

            if do_eci_vel:
                itup[n][i][0] = eci_vel_x
                itup[n][i][1] = eci_vel_y
                itup[n][i][2] = eci_vel_z
                n += 1

            if do_geo:
                geo_pos = sat.geo_pos()
                itup[n][i][0] = geo_pos.lon
                itup[n][i][1] = geo_pos.lat
                itup[n][i][2] = geo_pos.alt
                n += 1

            if do_topo:
                topo_pos = sat.topo_pos()
                itup[n][i][0] = topo_pos.az
                itup[n][i][1] = topo_pos.el
                itup[n][i][2] = topo_pos.dist
                itup[n][i][3] = topo_pos.dist_rate
                n += 1

    result = {}

    n = 3
    if do_eci_pos:
        result['eci_pos'] = it.operands[n]
        n += 1

    if do_eci_vel:
        result['eci_vel'] = it.operands[n]
        n += 1

    if do_geo:
        result['geo'] = it.operands[n]
        n += 1

    if do_topo:
        result['topo'] = it.operands[n]
        n += 1

    return result