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
from numpy cimport PyArray_MultiIter_DATA as Py_Iter_DATA
from cython.operator cimport dereference as deref
from cython.operator cimport address as addr
from cython.operator cimport preincrement as inc
from cpython cimport bool as python_bool
from libcpp cimport bool as cpp_bool
from libc.math cimport M_PI, floor, fabs, fmod
from .cysgp4 cimport *

from datetime import datetime
import numpy as np

np.import_array()


cdef double DEG2RAD = M_PI / 180.
cdef double RAD2DEG = 180. / M_PI
cdef double MJD_RESOLUTION = 0.001 / 24. / 3600.
cdef long long MJD0_TICKS = 58628880000000000


ctypedef SGP4* sgp4_ptr_t
ctypedef Observer* obs_ptr_t


__all__ = [
    'PyDateTime', 'PyTle', 'PyObserver',
    'PyCoordGeodetic', 'PyCoordTopocentric', 'PyEci',
    'Satellite', 'propagate_many',
    ]


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


cdef class PyDateTime(object):
    '''
    Wrapper around (C++) DateTime class
    '''

    # hold the C++ instance, which we're wrapping
    cdef DateTime _cobj

    def __init__(self, object dt=None, init=True):

        self._cobj = DateTime(0)
        if init:
            self._set_datetime(dt)

    @classmethod
    def from_ticks(cls, unsigned long long ticks):

        dt = cls(dt=None, init=False)
        dt.ticks = ticks

        return dt

    @classmethod
    def from_mjd(cls, double mjd):

        dt = cls(dt=None, init=False)
        dt.mjd = mjd

        return dt

    def _get_datetime(self):

        return datetime(
            self._cobj.Year(),
            self._cobj.Month(),
            self._cobj.Day(),
            self._cobj.Hour(),
            self._cobj.Minute(),
            self._cobj.Second(),
            self._cobj.Microsecond(),
            )

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

    datetime = property(_get_datetime, _set_datetime, None)

    def set(
            self,
            int year, int month, int day,
            int hour, int minute, int second, int microsecond
            ):
        '''
        Initialize PyDateTime
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

    ticks = property(_get_ticks, _set_ticks, None)

    def _get_mjd(self):
        return mjd_from_ticks(self._cobj.Ticks())

    def _set_mjd(self, double mjd):

        # this is a bit ugly, but there is to setter method in the C++ code

        cdef:
            long long ticks_new = ticks_from_mjd(mjd)
            long long ticks_old = self._get_ticks()

        # AddTicks returns a new instance...
        self._cobj = self._cobj.AddTicks(ticks_new - ticks_old)

    mjd = property(_get_mjd, _set_mjd, None)

    def __str__(self):

        return self._cobj.ToString().decode('UTF-8')

    def __repr__(self):

        return '<PyDateTime: ' + self.__str__() + '>'

    def gmst(self):

        return self._cobj.ToGreenwichSiderealTime()

    def lmst(self, obslon_deg):

        return self._cobj.ToLocalMeanSiderealTime(DEG2RAD * obslon_deg)


cdef class PyTle(object):
    '''
    Wrapper around (C++) Tle class
    '''

    # hold the C++ instance, which we're wrapping
    cdef Tle *thisptr

    def __init__(self, name, line_one, line_two):

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


cdef class PyCoordGeodetic(object):
    '''
    Wrapper around (C++) CoordGeodetic struct
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

    lon = property(_get_lon, _set_lon, None)
    lat = property(_get_lat, _set_lat, None)
    alt = property(_get_alt, _set_alt, None)


cdef class PyCoordTopocentric(object):
    '''
    Wrapper around (C++) CoordTopocentric struct
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

    az = property(_get_az, _set_az, None)
    el = property(_get_el, _set_el, None)
    dist = property(_get_dist, _set_dist, None)
    dist_rate = property(_get_dist_rate, _set_dist_rate, None)


cdef class PyObserver(object):
    '''
    Wrapper around (C++) Observer class
    '''

    # hold the C++ instance, which we're wrapping
    cdef:
        Observer *thisptr
        PyCoordGeodetic _obs_loc

    def __init__(
            self,
            double lon_deg=6.883750,
            double lat_deg=50.525,
            double alt_km=0.319
            ):
        '''
        Constructor PyObserver(double lon_deg, double lat_deg, double alt_km)
        '''

        self._obs_loc = PyCoordGeodetic(
            lon_deg=lon_deg,
            lat_deg=lat_deg,
            alt_km=alt_km
            )
        self.thisptr = new Observer(self._obs_loc._cobj)

    def __dealloc__(self):

        del self.thisptr

    def __str__(self):

        return self._obs_loc.__str__()

    def __repr__(self):

        return '<PyObserver: ' + self.__str__() + '>'

    def _get_location(self):

        return self._obs_loc

    def _set_location(self, PyCoordGeodetic loc):

        self.thisptr.SetLocation(loc._cobj)

    location = property(_get_location, _set_location, None)


cdef class PyEci(object):
    '''
    Wrapper around (C++) Eci class
    '''

    cdef:
        # hold the C++ instance, which we're wrapping
        Eci _cobj

        PyCoordGeodetic _geo_loc
        PyDateTime _dt

    def __init__(self, PyDateTime dt=None, PyCoordGeodetic geo_loc=None):
        '''
        Constructor PyEci(PyDateTime dt, PyCoordGeodetic geo_loc)
        '''

        if dt is None:

            dt = PyDateTime()

        if geo_loc is None:

            geo_loc = PyCoordGeodetic()

        self._dt = dt
        self._geo_loc = geo_loc

        self._cobj = Eci(dt._cobj, geo_loc._cobj)

    def __str__(self):

        return self._get_geo_loc().__str__() + ' ' + self._get_dt().__str__()

    def __repr__(self):

        return '<PyEci: ' + self.__str__() + '>'

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

    def _get_dt(self):

        dt = PyDateTime()
        dt._cobj = self._cobj.GetDateTime()
        return dt

    loc = property(_get_loc, None, None)
    vel = property(_get_vel, None, None)
    geo_loc = property(_get_geo_loc, None, None)
    dt = property(_get_dt, None, None)


cdef class Satellite(object):
    '''
    Calculates apparent positions of satellite for a given observer and TLE
    '''

    cdef:
        # hold C++ instances, which we're wrapping
        SGP4 *thisptr

        PyTle _tle
        PyObserver _observer

        PyDateTime _dt
        PyEci _eci
        PyCoordTopocentric _topo
        PyCoordGeodetic _geo

        double _mjd, _mjd_cache_resolution
        python_bool _pos_dirty, _tle_dirty

    def __init__(
            self,
            PyTle tle,
            PyObserver observer=None,
            object dt=None,
            double mjd_cache_resolution=MJD_RESOLUTION,
            ):
        '''
        Constructs a new Satellite object from given TLE

        if observer is None, Effelsberg location is used
        '''

        if observer is None:

            observer = PyObserver()

        self._mjd_cache_resolution = mjd_cache_resolution
        self._tle = tle  # copy reference
        self._observer = observer  # copy reference

        try:

            self.thisptr = new SGP4(deref(tle.thisptr))
            self._tle_dirty = <python_bool> False

        except:

            print('SatelliteException catched')
            self._tle_dirty = <python_bool> True


        # the following is important, otherwise self._topo and self._geo will
        # just be None after _refresh_coords()

        # initialize workspaces, otherwise we cannot assign values in
        # _refresh_coords (would produce segfault)
        # note: it is not sufficient to define these as class members (above)
        self._dt = PyDateTime(dt)  # initialize with current datetime
        self._eci = PyEci()
        self._topo = PyCoordTopocentric()
        self._geo = PyCoordGeodetic()

        self._pos_dirty = <python_bool> True

    def __dealloc__(self):

        # del self.eci_ptr
        del self.thisptr

    def _get_mjd(self):

        return self._mjd

    def _set_mjd(self, double mjd):

        assert mjd < 1000000., 'warning, make sure to use mjd'

        if fabs(self._mjd - mjd) < self._mjd_cache_resolution:
            return

        self._dt.mjd = self._mjd = mjd
        self._pos_dirty = <python_bool> True

    mjd = property(_get_mjd, _set_mjd, None, 'mjd')

    def _get_datetime(self):

        return self._dt

    def _set_datetime(self, dt):

        self._set_mjd = dt.mjd

    datetime = property(_get_datetime, _set_datetime, None, 'datetime')

    def topo_pos(self):

        if self._pos_dirty:
            self._refresh_coords()

        if self._tle_dirty:
            return None

        return self._topo

    def geo_pos(self):

        if self._pos_dirty:
            self._refresh_coords()

        if self._tle_dirty:
            return None

        return self._geo

    def eci_pos(self):

        if self._pos_dirty:
            self._refresh_coords()

        if self._tle_dirty:
            return None

        return self._eci

    def _refresh_coords(self):

        try:

            # FindPosition doesn't update ECI time, need to do manually :-/
            self._eci = PyEci(dt=self._dt)
            self._eci._cobj = self.thisptr.FindPosition(self._dt._cobj)
            self._tle_dirty = <python_bool> False

        except:

            print('SatelliteException catched')
            self._tle_dirty = <python_bool> True

            return

        self._topo._cobj = deref(self._observer.thisptr).GetLookAngle(
            self._eci._cobj
            )
        self._geo._cobj = self._eci._cobj.ToGeodetic()

        self._pos_dirty = <python_bool> False


def propagate_many(mjds, tles, observers=None):

    cdef:

        SGP4 *_sgp4_ptr
        Observer *_obs_ptr
        Eci _eci
        DateTime _dt
        CoordTopocentric _topo
        Vector _eci_pos, _eci_vel

        np.ndarray[double] mjd
        np.ndarray[object] sat
        # double[::1] mjd_v
        # Satellite[::1] sat_v
        double[::1] eci_x_v, eci_y_v, eci_z_v
        double[::1] eci_vx_v, eci_vy_v, eci_vz_v
        double[::1] az_v, el_v, dist_v, dist_rate_v
        int i, size
        bint do_topo = True

        sgp4_ptr_t* _sgp4_ptr_array
        obs_ptr_t* _obs_ptr_array

    if observers is None:
        do_topo = False

    b = np.broadcast(
        tles,
        observers if observers is not None else PyObserver()
        )
    sats = np.empty(b.shape, dtype=object)
    sats.flat = [Satellite(tle, obs) for (tle, obs) in b]

    it = np.nditer(
        [sats, mjds] + [None] * 10,
        flags=['external_loop', 'buffered', 'delay_bufalloc', 'refs_ok'],
        op_flags=[['readonly']] * 2 + [['readwrite', 'allocate']] * 10,
        op_dtypes=['object', 'float64'] + ['float64'] * 10
        )

    # it would be better to use the context manager but
    # "with it:" requires numpy >= 1.14

    it.reset()

    for itup in it:

        sat = itup[0]
        mjd = itup[1]
        eci_x_v = itup[2]
        eci_y_v = itup[3]
        eci_z_v = itup[4]
        eci_vx_v = itup[5]
        eci_vy_v = itup[6]
        eci_vz_v = itup[7]
        az_v = itup[8]
        el_v = itup[9]
        dist_v = itup[10]
        dist_rate_v = itup[11]

        size = mjd.shape[0]
        _sgp4_ptr_array = array_new[sgp4_ptr_t](size)
        _obs_ptr_array = array_new[obs_ptr_t](size)

        for i in range(size):
            # unfortunately, it is not possible in nogil loop to access
            # the cdef'ed class members; therefore, we have to maintain
            # arrays of pointers to the sgp4 and observer objects
            _sgp4_ptr_array[i] = (<Satellite> sat[i]).thisptr
            _obs_ptr_array[i] = (<Satellite> sat[i])._observer.thisptr

        for i in prange(size, nogil=True):

            _sgp4_ptr = _sgp4_ptr_array[i]
            _obs_ptr = _obs_ptr_array[i]

            # AddTicks returns a new instance...
            _dt = _dt.AddTicks(ticks_from_mjd(mjd[i]) - _dt.Ticks())

            _eci = _sgp4_ptr.FindPosition(_dt)
            if do_topo:
                _topo = _obs_ptr.GetLookAngle(_eci)

            _eci_pos = _eci.Position()
            _eci_vel = _eci.Velocity()
            eci_x_v[i] = _eci_pos.x
            eci_y_v[i] = _eci_pos.y
            eci_z_v[i] = _eci_pos.z
            eci_vx_v[i] = _eci_vel.x
            eci_vy_v[i] = _eci_vel.y
            eci_vz_v[i] = _eci_vel.z
            az_v[i] = _topo.azimuth * RAD2DEG
            el_v[i] = _topo.elevation * RAD2DEG
            dist_v[i] = _topo.distance
            dist_rate_v[i] = _topo.distance_rate

        array_delete(_sgp4_ptr_array)
        array_delete(_obs_ptr_array)

    if observers is None:
        return it.operands[2:5], it.operands[5:8]
    else:
        return it.operands[2:5], it.operands[5:8], it.operands[8:12]

