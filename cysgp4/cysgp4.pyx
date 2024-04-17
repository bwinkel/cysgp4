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
from libc.math cimport M_PI, floor, fabs, fmod, sqrt, sin, cos, atan2, acos
from libc.stdint cimport uint32_t, int64_t
from .cysgp4 cimport *

import datetime
import numpy as np

np.import_array()


cdef double NAN = np.nan
cdef double DEG2RAD = M_PI / 180.
cdef double RAD2DEG = 180. / M_PI
cdef double MJD_RESOLUTION = 0.001 / 24. / 3600.
cdef int64_t MJD0_TICKS = 58628880000000000


ctypedef SGP4* sgp4_ptr_t
ctypedef Tle* tle_ptr_t
ctypedef Observer* obs_ptr_t


__all__ = [
    'PyDateTime', 'PyTle', 'PyObserver',
    'PyCoordGeodetic', 'PyCoordTopocentric', 'PyEci',
    'Satellite', 'set_num_threads',
    'eci_to_geo', 'geo_to_eci', 'lookangles',
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


cdef inline int64_t ticks_from_mjd(double mjd) noexcept nogil:

    cdef:
        double days, fdays
        int64_t idays

    idays = <int64_t> mjd
    fdays = mjd - idays

    return (
        idays * 86400 * 1000000 +
        (<int64_t> (fdays * 86400 * 1000000)) +
        MJD0_TICKS
        )


cdef inline double mjd_from_ticks(int64_t ticks) noexcept nogil:

    cdef:
        double days, fdays
        int64_t idays

    ticks -= MJD0_TICKS
    return ticks / 8.64e10


cdef inline DateTime datetime_from_mjd(double mjd) noexcept nogil:

    cdef:

        int64_t ticks = ticks_from_mjd(mjd)
        DateTime dt = DateTime(ticks)

    return dt


cdef inline (double, double, double) ecef_from_geo(
        double lon_rad, double lat_rad, double alt_km
        ) noexcept nogil:
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


cdef inline Vector normalize_vector(Vector a) noexcept nogil:

    cdef double norm = sqrt(a.x ** 2 + a.y ** 2 + a.z ** 2)
    a.x /= norm
    a.y /= norm
    a.z /= norm

    return a


cdef inline Vector cross_prod(Vector a, Vector b) noexcept nogil:

    return Vector(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
        )

cdef inline double dot_prod(Vector a, Vector b) noexcept nogil:

    return a.x * b.x + a.y * b.y + a.z * b.z


cdef inline (double, double, double, Vector, Vector, Vector) calc_sat_azel(
        Vector sat_pos, Vector sat_vel, Vector obs_pos
        ) noexcept nogil:

    cdef:

        double norm
        Vector e_x, e_y, e_z
        Vector sat_azel, diff_pos
        double b_x, b_y, b_z
        double sat_dist, sat_az, sat_el

    e_z = normalize_vector(sat_vel)
    e_x = normalize_vector(Vector(-sat_pos.x, -sat_pos.y, -sat_pos.z))
    e_y = normalize_vector(cross_prod(e_z, e_x))
    e_x = cross_prod(e_y, e_z)

    diff_pos.x = obs_pos.x - sat_pos.x
    diff_pos.y = obs_pos.y - sat_pos.y
    diff_pos.z = obs_pos.z - sat_pos.z
    b_x = dot_prod(diff_pos, e_x)
    b_y = dot_prod(diff_pos, e_y)
    b_z = dot_prod(diff_pos, e_z)

    sat_dist = sqrt(b_x ** 2 + b_y ** 2 + b_z ** 2)
    sat_el = M_PI / 2 - acos(b_z / sat_dist)
    sat_az = atan2(b_y, b_x)

    return sat_dist, sat_az, sat_el, e_x, e_y, e_z


cdef inline (double, double, double, Vector, Vector, Vector) calc_sat_iso(
        Vector sat_pos, Vector sat_vel, Vector obs_pos
        ) noexcept nogil:

    cdef:

        double norm
        Vector e_x, e_y, e_z
        Vector sat_iso, diff_pos
        double b_x, b_y, b_z
        double sat_dist, sat_theta, sat_phi

    e_x = normalize_vector(sat_vel)
    e_z = normalize_vector(Vector(-sat_pos.x, -sat_pos.y, -sat_pos.z))
    e_y = normalize_vector(cross_prod(e_z, e_x))
    e_z = cross_prod(e_x, e_y)

    diff_pos.x = obs_pos.x - sat_pos.x
    diff_pos.y = obs_pos.y - sat_pos.y
    diff_pos.z = obs_pos.z - sat_pos.z
    b_x = dot_prod(diff_pos, e_x)
    b_y = dot_prod(diff_pos, e_y)
    b_z = dot_prod(diff_pos, e_z)

    sat_dist = sqrt(b_x ** 2 + b_y ** 2 + b_z ** 2)
    sat_theta = acos(b_z / sat_dist)
    sat_phi = atan2(b_y, b_x)

    return sat_dist, sat_phi, sat_theta, e_x, e_y, e_z


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
        current date and time (as given by `~datetime.datetime.utcnow()`)
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
    def from_ticks(cls, int64_t ticks):
        '''
        Creates a new `~cysgp4.PyDateTime` instance from "ticks".

        Internally, DateTime objects store time on a linear scale, "ticks",
        which is the number of micro-seconds since 1. January 0001, 00:00.

        Parameters
        ----------
        ticks : int64_t
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

    @classmethod
    def from_tle_epoch(cls, double tle_epoch):
        '''
        Creates a new `~cysgp4.PyDateTime` instance from TLE epoch format.

        Note: TLE epochs have a very strange format: first two digits are the
        year, next three digits are the day from beginning of year, then
        the fraction of a day is given, e.g. 20180.25 would be 2020, day 180,
        6 hours (probably UT as no timezone is mentioned). See also
        `Wikipedia <https://en.wikipedia.org/wiki/Two-line_element_set>`_
        or `Celestrak <https://celestrak.com/columns/v04n03/>`_.

        Corresponding year must be between 1957 and 2056.

        Parameters
        ----------
        tle_epoch : double
            Datetime in TLE epoch format.

        Returns
        -------
        pydt : `~cysgp4.PyDateTime` object

        Examples
        --------
        "from_tle_epoch" is a classmethod::

            >>> tle_epoch = 19050.1
            >>> pydt = PyDateTime.from_tle_epoch(tle_epoch)
            >>> pydt
            <PyDateTime: 2019-02-19 02:24:00.000000 UTC>
        '''

        dt = cls(dt=None, init=False)
        dt.tle_epoch = tle_epoch

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

        return datetime.datetime(*self.get_datetime_tuple())

    def _set_datetime(self, dt):
        '''
        Initialize PyDateTime from python datetime object
        '''

        if dt is None:
            dt = datetime.datetime.utcnow()

        assert isinstance(dt, datetime.datetime)

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
        return <int64_t> self._cobj.Ticks()

    def _set_ticks(self, int64_t ticks):

        # this is a bit ugly, but there is to setter method in the C++ code

        cdef:
            int64_t ticks_new = ticks
            int64_t ticks_old = self._get_ticks()

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
            int64_t ticks_new = ticks_from_mjd(mjd)
            int64_t ticks_old = self._get_ticks()

        # AddTicks returns a new instance...
        self._cobj = self._cobj.AddTicks(ticks_new - ticks_old)

    mjd = property(
        _get_mjd, _set_mjd, None, doc='MJD (see Class documentation).'
        )

    def _get_tle_epoch(self):

        cdef:
            DateTime dt = self._cobj
            int year = dt.Year()
            double _tle_epoch = (
                (year % 100) * 1000. +
                dt.DayOfYear(year, dt.Month(), dt.Day()) +
                dt.Hour() / 24. +
                dt.Minute() / 1440. +
                (dt.Second() + dt.Microsecond() * 1.e-6) / 86400.
                )

        if year < 1957 or year > 2056:
            raise ValueError('Year must be between 1957 and 2056')

        return _tle_epoch

    def _set_tle_epoch(self, double _tle_epoch):

        cdef:

            int year = (<int> _tle_epoch) // 1000
            int idoy = (<int> _tle_epoch) % 1000
            double fdoy = _tle_epoch - <int> _tle_epoch
            double doy = idoy + fdoy

            DateTime dt

        year += 1900
        if year < 1957:
            year += 100

        dt = DateTime(year, doy)
        self._cobj = dt

    tle_epoch = property(
        _get_tle_epoch, _set_tle_epoch, None,
        doc='Datetime in TLE epoch format. See also `~cysgp4.from_tle_epoch`.'
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

    @property
    def epoch(self):
        '''
        Epoch of the satellites TLE.

        Note that TLEs get outdated quickly in which case, the orbit
        computations can have large errorbars.

        Returns
        ----------
        epoch : `~cysgp4.PyDateTime`
            Epoch of the satellites TLE.
        '''

        cdef PyDateTime pydt = PyDateTime(init=False)
        pydt._cobj = self.thisptr.Epoch()

        return pydt

    @property
    def int_designator(self):
        '''
        Satellite's international designator.

        First two digits is the (last two digits of the) launch year.
        Next three digits are the launch number (in that year).
        Last digit is the piece of the launch.
        (see `TLEs on Wikipedia
        <https://en.wikipedia.org/wiki/Two-line_element_set>`_)
        '''

        return self.thisptr.IntDesignator()

    @property
    def catalog_number(self):
        '''
        Satellite's (NORAD) catalog number.
        '''

        return self.thisptr.NoradNumber()

    @property
    def inclination(self):
        '''
        Satellite's inclination [deg].
        '''

        cdef bint in_degrees = True

        return self.thisptr.Inclination(in_degrees)

    @property
    def raan(self):
        '''
        Satellite's right ascension of the ascending node (RAAN) [deg].
        '''

        cdef bint in_degrees = True

        return self.thisptr.RightAscendingNode(in_degrees)

    @property
    def eccentricity(self):
        '''
        Satellite orbit's eccentricity [dimless].
        '''

        return self.thisptr.Eccentricity()

    @property
    def arg_perigee(self):
        '''
        Satellite's argument of perigee [deg].
        '''

        cdef bint in_degrees = True

        return self.thisptr.ArgumentPerigee(in_degrees)

    @property
    def mean_anomaly(self):
        '''
        Satellite's mean anomaly [deg].
        '''

        cdef bint in_degrees = True

        return self.thisptr.MeanAnomaly(in_degrees)

    @property
    def mean_motion(self):
        '''
        Satellite's mean motion [1 / day].
        '''

        return self.thisptr.MeanMotion()

    @property
    def rev_number(self):
        '''
        Satellite's revolution number at epoch (revolutions).
        '''

        return self.thisptr.OrbitNumber()

    @property
    def mean_motion_dt2(self):
        '''
        Satellite's first time derivative of the mean motion divided by two.
        '''

        return self.thisptr.MeanMotionDt2()

    @property
    def mean_motion_dt6(self):
        '''
        Satellite's second time derivative of mean motion divided by six.
        '''

        return self.thisptr.MeanMotionDdt6()

    @property
    def bstar(self):
        '''
        Satellite's BSTAR drag term.
        '''

        return self.thisptr.BStar()

    # need comparison operators for finding unique TLEs
    def __lt__(self, other):
        return self.tle_strings() < other.tle_strings()

    def __eq__(self, other):
        return self.tle_strings() == other.tle_strings()


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

    # need comparison operators for finding unique observers
    def __lt__(self, other):
        return self.loc.lon < other.loc.lon

    def __eq__(self, other):
        return (
            self.loc.lon == other.loc.lon and
            self.loc.lat == other.loc.lat and
            self.loc.alt == other.loc.alt
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
        `~datetime.datetime.utcnow()`.
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


def _propagate_many_cysgp4(
        mjds, tles, observers=None,
        bint do_eci_pos=True, bint do_eci_vel=True,
        bint do_geo=True, bint do_topo=True,
        bint do_obs_pos=False, bint do_sat_azel=False,
        bint do_sat_rotmat=False,
        str sat_frame='zxy', str on_error='raise',
        ):
    '''
    Calculate positions of many satellites at a various times at once.

    This is an array interface to the sgp4 calculations, which allows to
    perform calculations for different satellite TLEs, observers and times
    in a parallelized manner. `~numpy` broadcasting rules apply.

    With the Boolean parameters, `do_eci_pos`, `do_eci_vel`, `do_geo`,
    `do_topo`, `do_obs_pos`, and `do_sat_azel` the user can decide, which
    position frames are returned in the output dictionary.

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
    do_obs_pos : Boolean, optional (default: False)
        Whether to include the observer ECI position in the results.
    do_sat_azel : Boolean, optional (default: False)
        Whether to include the observer position as seen by the satellite
        (distance/azimuth/elevation) in the results.
    do_sat_rotmat : Boolean, optional (default: False)
        Whether to include the rotation matrix that converts the
        (moving and rotated) satellite frame (in cartesian) into
        cartesian ECI-aligned coordinates in the results. This can be useful
        for cases where the user needs to transform additional vectors
        between both frames (and is not only interested the observer
        position in the satellite frame as returned by `do_sat_azel').
    sat_frame : 'zxy' or 'xyz', optional (default: 'zxy')
        How the moving satellite frame is defined. Two options are
        implemented, 'zxy' and 'xyz'. If 'zxy' is chosen, the moving
        satellite frame is constructed such that the `z` axis is
        aligned with the satellite motion vector. The `y` axis is lies
        perpendicularly to the plane defined by the motion vector and
        the ECI zero point (aka the Earth centre). The resulting `x`
        axis, which is orthogonal to the `y` and `z` axes, is then
        approximately pointing towards nadir. Alternatively, if the
        frame is set as `xyz`, the `x` axis is the motion vector, `y`
        has the same meaning (but points into the opposite direction)
        and `z` is approximately pointing towards the nadir. The
        definition of the output polar angles is different for the two
        reference frames, see Returns.
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
          one for each of `x`, `y`, and `z`.

        - `eci_vel` : `~numpy.ndarray` of float

          Satellites ECI cartesian velicities. Last dimension has length 3,
          one for each of `v_x`, `v_y`, and `v_z`.

        - `geo` : `~numpy.ndarray` of float

          Satellites Geodetic positions. Last dimension has length 3, one
          for each of `lon`, `lat`, and `alt`.

        - `topo` : `~numpy.ndarray` of float

          Satellites Topocentric positions. Last dimension has length 4,
          one for each of `az`, `el`, `dist`, and `dist_rate`.

        - `obs_pos` : `~numpy.ndarray` of float

          Observer positions in ECI frame (Cartesian). Last dimension has
          length 3, one for each of `x`, `y`, and `z`.

        - `sat_azel` : `~numpy.ndarray` of float

          If `sat_frame` is 'zxy', `z` lies in the direction of motion,
          `y` perpendicular to the z-axis and the Earth center, `x` is
          pointing approximately towards nadir, see also `sat_frame`
          parameter description. The Observer positions in the
          (co-moving) satellite frame are given as azimuth, elevation
          in the specified reference frame, and distance (`az`, `el`,
          `dist`). `az` is the angle between the projection of the vector
          towards the Observer onto the xy-plane and the x-axis. -180
          deg < `az` < 180 deg. `el` is the angle between the normal
          vector and the xy-plane. -90 deg < `el` < 90 deg.

          If `sat_frame` is 'xyz', `x` lies in the direction of motion,
          `y` is perpendicular to `z` and the Earth center, `z` is pointing
          approximately towards nadir, see also `sat_frame` parameter
          description. The Observer positions in the (moving)
          satellite frame are given as azimuth and polar angle in the
          specified reference frame, and distance (`az`, `theta`, `dist`). `az`
          is the angle between the projection of the vector towards
          the observer onto the xy-plane and the x-axis. -180 deg < `az`
          < 180 deg. `theta` is the angle between the normal vector and
          the z-axis. -90 deg < `theta` < 90 deg.

        - `sat_rotmat` : `~numpy.ndarray` of float

          Rotation matrices which would transform a vector defined in the
          (moving and rotated) satellite frames (in cartesian) to the
          cartesian ECI-aligned basis frame. It is noted that the origin of
          this ECI-aligned frame is still at the satellite center.

          This can be useful for cases where the user needs to transform
          additional vectors between both frames (and is not only interested
          the observer position in the satellite frame as returned by
          `do_sat_azel').

          Likewise, the inverse of these rotation matrices (aka the
          transposed) can be used to rotate any vector from ECI-aligned
          satellite basis frame to the satellite frame.

        In all cases the first dimensions are determined by the
        (broadcasted) shape of the inputs `mjd`, `tles`, and `observers`.

    Examples
    --------
    The following demonstrates how to use the `~cysgp4.propagate_many`
    function::

        >>> import numpy as np
        >>> from cysgp4 import PyTle, PyObserver, propagate_many
        >>> from cysgp4 import get_example_tles, tles_from_text

        >>> tle_text = get_example_tles()
        >>> tles = np.array(
        ...     tles_from_text(tle_text)
        ...     )[np.newaxis, np.newaxis, :20]  # use first 20 TLEs
        >>> observers = np.array([
        ...     PyObserver(6.88375, 50.525, 0.366),
        ...     PyObserver(16.88375, 50.525, 0.366),
        ...     ])[np.newaxis, :, np.newaxis]
        >>> mjds = np.linspace(
        ...     58805.5, 58806.5, 1000  # 1000 time steps
        ...     )[:, np.newaxis, np.newaxis]

        >>> result = propagate_many(mjds, tles, observers)
        >>> print(sorted(result.keys()))
        ['eci_pos', 'eci_vel', 'geo', 'topo']

        >>> # shapes are as follows
        >>> print(np.broadcast(mjds, tles, observers).shape)
        (1000, 2, 20)
        >>> print(result['eci_pos'].shape, result['topo'].shape)
        (1000, 2, 20, 3) (1000, 2, 20, 4)

        >>> result = propagate_many(
        ...     mjds, tles, observers,
        ...     do_eci_pos=False, do_eci_vel=False, do_geo=False, do_topo=True
        ...     )
        >>> print(sorted(result.keys()))
        ['topo']

    '''

    cdef:

        SGP4 *_sgp4_ptr
        Tle *_tle_ptr
        Observer _obs
        # PyEci py_nan_eci = PyEci(geo_loc=PyCoordGeodetic(NAN, NAN, NAN))
        Eci _sat_eci, _obs_eci
        # Eci _nan_eci = py_nan_eci._cobj
        DateTime _dt
        CoordGeodetic _geo
        CoordTopocentric _topo
        Vector _sat_pos, _sat_vel, _obs_pos
        double _sat_dist, _sat_az, _sat_el  # observer pos in sat frame
        Vector _sat_e_x, _sat_e_y, _sat_e_z

        np.ndarray[double] mjd
        np.ndarray[object] tle, obs
        double[:, ::1] sat_pos_v
        double[:, ::1] sat_vel_v
        double[:, ::1] geo_v  # geodetic
        double[:, ::1] topo_v  # topocentric
        double[:, ::1] obs_pos_v  # observer eci position
        double[:, ::1] sat_azel_v  # observer pos in sat frame (as dist/az/el)
        double[:, ::1] sat_e_x_v, sat_e_y_v, sat_e_z_v  # sat frame basis
        int i, size, n, pnum

        tle_ptr_t* _tle_ptr_array
        obs_ptr_t* _obs_ptr_array

        bint on_error_raise = True

        int _sat_frame = 0  # 0: 'zxy', 1: 'xyz'

    if on_error not in ['raise', 'coerce_to_nan']:
        raise ValueError(
            'Argument "on_error" must be one of "raise" or "coerce_to_nan"'
            )
    if sat_frame not in ['zxy', 'xyz']:
        raise ValueError(
            'Argument "sat_frame" must be one of "zxy" or "xyz"'
            )

    if sat_frame == 'xyz':
        _sat_frame = 1

    if on_error == 'coerce_to_nan':
        on_error_raise = False

    if observers is None:
        observers = PyObserver()

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
    if do_obs_pos:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    if do_sat_azel:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    if do_sat_rotmat:
        out_dts.append(np.dtype(('float64', 3)))
        out_dts.append(np.dtype(('float64', 3)))
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 3

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
            sat_pos_v = itup[n]
            n += 1

        if do_eci_vel:
            sat_vel_v = itup[n]
            n += 1

        if do_geo:
            geo_v = itup[n]
            n += 1

        if do_topo:
            topo_v = itup[n]
            n += 1

        if do_obs_pos:
            obs_pos_v = itup[n]
            n += 1

        if do_sat_azel:
            sat_azel_v = itup[n]
            n += 1

        if do_sat_rotmat:
            sat_e_x_v = itup[n]
            n += 1
            sat_e_y_v = itup[n]
            n += 1
            sat_e_z_v = itup[n]
            n += 1

        size = mjd.shape[0]
        _tle_ptr_array = array_new[tle_ptr_t](size)
        _obs_ptr_array = array_new[obs_ptr_t](size)

        for i in range(size):
            # unfortunately, it is not possible in nogil loop to access
            # the cdef'ed class members; therefore, we have to maintain
            # arrays of pointers to the sgp4 and observer objects
            if not isinstance(tle[i], PyTle):
                print(type(tle[i]))
                raise TypeError(
                    'Argument "tles" must be of type "PyTle" '
                    '(or list/array of "PyTle")'
                    )
            if not isinstance(obs[i], PyObserver):
                print(type(obs[i]))
                raise TypeError(
                    'Argument "observers" must be of type "PyObserver" '
                    '(or list/array of "PyObserver")'
                    )

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
                _sat_eci = _sgp4_ptr.FindPosition(_dt)
            else:
                # "FindPositionE" will return NaNs if an error occured
                # at C++ level
                _sat_eci = _sgp4_ptr.FindPositionNaN(_dt)

            del _sgp4_ptr

            if do_eci_pos:
                _sat_pos = _sat_eci.Position()
                sat_pos_v[i, 0] = _sat_pos.x
                sat_pos_v[i, 1] = _sat_pos.y
                sat_pos_v[i, 2] = _sat_pos.z

            if do_eci_vel:
                _sat_vel = _sat_eci.Velocity()
                sat_vel_v[i, 0] = _sat_vel.x
                sat_vel_v[i, 1] = _sat_vel.y
                sat_vel_v[i, 2] = _sat_vel.z

            if do_geo:
                _geo = _sat_eci.ToGeodetic()
                geo_v[i, 0] = _geo.longitude * RAD2DEG
                geo_v[i, 1] = _geo.latitude * RAD2DEG
                geo_v[i, 2] = _geo.altitude

            if do_topo:
                _obs = Observer(_obs_ptr_array[i].GetLocation())
                _topo = _obs.GetLookAngle(_sat_eci)
                topo_v[i, 0] = _topo.azimuth * RAD2DEG
                topo_v[i, 1] = _topo.elevation * RAD2DEG
                topo_v[i, 2] = _topo.distance
                topo_v[i, 3] = _topo.distance_rate

            if do_obs_pos:
                _obs_eci = Eci(_dt, _obs_ptr_array[i].GetLocation())
                _obs_pos = _obs_eci.Position()
                obs_pos_v[i, 0] = _obs_pos.x
                obs_pos_v[i, 1] = _obs_pos.y
                obs_pos_v[i, 2] = _obs_pos.z

            if do_sat_azel or do_sat_rotmat:
                _obs_eci = Eci(_dt, _obs_ptr_array[i].GetLocation())
                _sat_pos = _sat_eci.Position()
                _sat_vel = _sat_eci.Velocity()
                _obs_pos = _obs_eci.Position()
                if _sat_frame == 0:
                    (
                        _sat_dist, _sat_az, _sat_el,
                        _sat_e_x, _sat_e_y, _sat_e_z
                        ) = calc_sat_azel(
                        _sat_pos, _sat_vel, _obs_pos
                        )
                elif _sat_frame == 1:
                    (
                        _sat_dist, _sat_az, _sat_el,
                        _sat_e_x, _sat_e_y, _sat_e_z
                        ) = calc_sat_iso(
                        _sat_pos, _sat_vel, _obs_pos
                        )

                if do_sat_azel:
                    sat_azel_v[i, 0] = _sat_az * RAD2DEG
                    sat_azel_v[i, 1] = _sat_el * RAD2DEG
                    sat_azel_v[i, 2] = _sat_dist

                if do_sat_rotmat:
                    sat_e_x_v[i, 0] = _sat_e_x.x
                    sat_e_x_v[i, 1] = _sat_e_x.y
                    sat_e_x_v[i, 2] = _sat_e_x.z
                    sat_e_y_v[i, 0] = _sat_e_y.x
                    sat_e_y_v[i, 1] = _sat_e_y.y
                    sat_e_y_v[i, 2] = _sat_e_y.z
                    sat_e_z_v[i, 0] = _sat_e_z.x
                    sat_e_z_v[i, 1] = _sat_e_z.y
                    sat_e_z_v[i, 2] = _sat_e_z.z

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

    if do_obs_pos:
        result['obs_pos'] = it.operands[n]
        n += 1

    if do_sat_azel:
        result['sat_azel'] = it.operands[n]
        n += 1

    if do_sat_rotmat:
        result['sat_rotmat'] = np.stack([
            it.operands[n + i] for i in range(3)
            ], axis=-1)
        n += 3


    return result


def eci_to_geo(eci_x, eci_y, eci_z, mjds):

    cdef:

        Eci _eci
        DateTime _dt
        CoordGeodetic _geo
        Vector _eci_pos

        np.ndarray[double] mjd
        const double[:] eci_x_v
        const double[:] eci_y_v
        const double[:] eci_z_v
        double[:] geo_lon_v
        double[:] geo_lat_v
        double[:] geo_alt_v
        int i, size

    it = np.nditer(
        [mjds, eci_x, eci_y, eci_z] + [None] * 3,
        flags=['external_loop', 'buffered', 'delay_bufalloc', 'refs_ok'],
        op_flags=[['readonly']] * 4 + [['readwrite', 'allocate']] * 3,
        op_dtypes=['float64'] * 7,
        )

    # it would be better to use the context manager but
    # "with it:" requires numpy >= 1.14

    it.reset()

    for itup in it:

        mjd = itup[0]
        eci_x_v = itup[1]
        eci_y_v = itup[2]
        eci_z_v = itup[3]
        geo_lon_v = itup[4]
        geo_lat_v = itup[5]
        geo_alt_v = itup[6]

        size = mjd.shape[0]

        for i in prange(size, nogil=True):

            # AddTicks returns a new instance...
            _dt = _dt.AddTicks(ticks_from_mjd(mjd[i]) - _dt.Ticks())
            _eci_pos = Vector(eci_x_v[i], eci_y_v[i], eci_z_v[i])
            _eci = Eci(_dt, _eci_pos)

            _geo = _eci.ToGeodetic()
            # attributes are in rad!
            geo_lon_v[i] = _geo.longitude * RAD2DEG
            geo_lat_v[i] = _geo.latitude * RAD2DEG
            geo_alt_v[i] = _geo.altitude

    return it.operands[4:7]


def geo_to_eci(lon, lat, alt, mjds):

    cdef:

        Eci _eci
        DateTime _dt
        CoordGeodetic _geo
        Vector _eci_pos

        np.ndarray[double] mjd
        double[:] eci_x_v
        double[:] eci_y_v
        double[:] eci_z_v
        const double[:] geo_lon_v
        const double[:] geo_lat_v
        const double[:] geo_alt_v
        int i, size

    it = np.nditer(
        [mjds, lon, lat, alt] + [None] * 3,
        flags=['external_loop', 'buffered', 'delay_bufalloc', 'refs_ok'],
        op_flags=[['readonly']] * 4 + [['readwrite', 'allocate']] * 3,
        op_dtypes=['float64'] * 7,
        )

    # it would be better to use the context manager but
    # "with it:" requires numpy >= 1.14

    it.reset()

    for itup in it:

        mjd = itup[0]
        geo_lon_v = itup[1]
        geo_lat_v = itup[2]
        geo_alt_v = itup[3]
        eci_x_v = itup[4]
        eci_y_v = itup[5]
        eci_z_v = itup[6]

        size = mjd.shape[0]

        for i in prange(size, nogil=True):

            # AddTicks returns a new instance...
            _dt = _dt.AddTicks(ticks_from_mjd(mjd[i]) - _dt.Ticks())
            # constructor uses radians by default
            _eci = Eci(_dt, geo_lat_v[i], geo_lon_v[i], geo_alt_v[i])
            _eci_pos = _eci.Position()

            eci_x_v[i] = _eci_pos.x
            eci_y_v[i] = _eci_pos.y
            eci_z_v[i] = _eci_pos.z

    return it.operands[4:7]


def lookangles(
        sat_pos_x, sat_pos_y, sat_pos_z,
        sat_vel_x, sat_vel_y, sat_vel_z,
        mjds, observers, str sat_frame='zxy',
        bint do_sat_rotmat=False,
        ):

    cdef:

        Eci _sat_eci, _obs_eci
        DateTime _dt
        Observer _obs
        # CoordGeodetic _geo
        CoordTopocentric _topo
        Vector _obs_pos, _sat_pos, _sat_vel
        double _sat_dist, _sat_az, _sat_el
        Vector _sat_e_x, _sat_e_y, _sat_e_z

        np.ndarray[double] mjd
        np.ndarray[object] obs

        const double[:] sat_pos_x_v
        const double[:] sat_pos_y_v
        const double[:] sat_pos_z_v
        const double[:] sat_vel_x_v
        const double[:] sat_vel_y_v
        const double[:] sat_vel_z_v
        double[:] obs_az_v, obs_el_v, sat_az_v, sat_el_v
        double[:] sat_e_x_x_v, sat_e_x_y_v, sat_e_x_z_v
        double[:] sat_e_y_x_v, sat_e_y_y_v, sat_e_y_z_v
        double[:] sat_e_z_x_v, sat_e_z_y_v, sat_e_z_z_v
        double[:] dist_v
        double[:] distrate_v

        obs_ptr_t* _obs_ptr_array

        int i, size
        int _sat_frame = 0  # 0: 'zxy', 1: 'xyz'
        int outpnum = 15 if do_sat_rotmat else 6

    assert sat_frame in ['zxy', 'xyz']

    if sat_frame == 'xyz':
        _sat_frame = 1

    it = np.nditer(
        [
            mjds, observers,
            sat_pos_x, sat_pos_y, sat_pos_z,
            sat_vel_x, sat_vel_y, sat_vel_z,
            ] + [None] * outpnum,
        flags=['external_loop', 'buffered', 'delay_bufalloc', 'refs_ok'],
        op_flags=[['readonly']] * 8 + [['readwrite', 'allocate']] * outpnum,
        op_dtypes=['float64', 'object'] + ['float64'] * (6 + outpnum),
        )

    # it would be better to use the context manager but
    # "with it:" requires numpy >= 1.14

    it.reset()

    for itup in it:

        mjd = itup[0]
        obs = itup[1]
        sat_pos_x_v = itup[2]
        sat_pos_y_v = itup[3]
        sat_pos_z_v = itup[4]
        sat_vel_x_v = itup[5]
        sat_vel_y_v = itup[6]
        sat_vel_z_v = itup[7]
        obs_az_v = itup[8]
        obs_el_v = itup[9]
        sat_az_v = itup[10]
        sat_el_v = itup[11]
        dist_v = itup[12]
        distrate_v = itup[13]
        if do_sat_rotmat:
            sat_e_x_x_v = itup[14]
            sat_e_x_y_v = itup[15]
            sat_e_x_z_v = itup[16]
            sat_e_y_x_v = itup[17]
            sat_e_y_y_v = itup[18]
            sat_e_y_z_v = itup[19]
            sat_e_z_x_v = itup[20]
            sat_e_z_y_v = itup[21]
            sat_e_z_z_v = itup[22]

        size = mjd.shape[0]
        _obs_ptr_array = array_new[obs_ptr_t](size)

        for i in range(size):
            # unfortunately, it is not possible in nogil loop to access
            # the cdef'ed class members; therefore, we have to maintain
            # arrays of pointers to the sgp4 and observer objects
            _obs_ptr_array[i] = &(<PyObserver> obs[i])._cobj

        for i in prange(size, nogil=True):

            # AddTicks returns a new instance...
            _dt = _dt.AddTicks(ticks_from_mjd(mjd[i]) - _dt.Ticks())
            _sat_pos = Vector(sat_pos_x_v[i], sat_pos_y_v[i], sat_pos_z_v[i])
            _sat_vel = Vector(sat_vel_x_v[i], sat_vel_y_v[i], sat_vel_z_v[i])
            _sat_eci = Eci(_dt, _sat_pos, _sat_vel)

            _obs = Observer(_obs_ptr_array[i].GetLocation())
            _topo = _obs.GetLookAngle(_sat_eci)
            obs_az_v[i] = _topo.azimuth * RAD2DEG
            obs_el_v[i] = _topo.elevation * RAD2DEG
            dist_v[i] = _topo.distance
            distrate_v[i] = _topo.distance_rate

            _obs_eci = Eci(_dt, _obs_ptr_array[i].GetLocation())
            _obs_pos = _obs_eci.Position()

            if _sat_frame == 0:
                (
                    _sat_dist, _sat_az, _sat_el,
                    _sat_e_x, _sat_e_y, _sat_e_z
                    ) = calc_sat_azel(
                    _sat_pos, _sat_vel, _obs_pos
                    )
            elif _sat_frame == 1:
                (
                    _sat_dist, _sat_az, _sat_el,
                    _sat_e_x, _sat_e_y, _sat_e_z
                    ) = calc_sat_iso(
                    _sat_pos, _sat_vel, _obs_pos
                    )
            sat_az_v[i] = _sat_az * RAD2DEG
            sat_el_v[i] = _sat_el * RAD2DEG
            if do_sat_rotmat:
                sat_e_x_x_v[i] = _sat_e_x.x
                sat_e_x_y_v[i] = _sat_e_x.y
                sat_e_x_z_v[i] = _sat_e_x.z
                sat_e_y_x_v[i] = _sat_e_y.x
                sat_e_y_y_v[i] = _sat_e_y.y
                sat_e_y_z_v[i] = _sat_e_y.z
                sat_e_z_x_v[i] = _sat_e_z.x
                sat_e_z_y_v[i] = _sat_e_z.y
                sat_e_z_z_v[i] = _sat_e_z.z

    return it.operands[8:23] if do_sat_rotmat else it.operands[8:14]
