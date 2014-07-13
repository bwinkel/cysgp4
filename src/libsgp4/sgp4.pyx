#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# title                  :sgp4.pyx
# description            :Cython-powered wrapper of the sgp4lib library by Daniel Warner.
# author                 :Benjamin Winkel
#
# ####################################################################
#  Copyright (C) 2014 by Benjamin Winkel
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
# ####################################################################

#!python
#cython: cdivision=True, embedsignature=True, boundscheck=False, wraparound=False


from sgp4 cimport *
from cython.operator cimport dereference as deref, address as addr, preincrement as inc
from cpython cimport bool as python_bool
from libcpp cimport bool as cpp_bool
from libc.math cimport M_PI, floor, fabs
from datetime import datetime
cimport cython

cdef double deg2rad = M_PI / 180.
cdef double rad2deg = 180. / M_PI
cdef double MJD_RESOLUTION = 0.001 / 24. / 3600.


cdef class PyDateTime(object):
    """Wrapper around (C++) DateTime class"""

    cdef DateTime *thisptr  # holds the C++ instance being wrapped

    def __cinit__(self, object dt=None):
        """Constructor PyDateTime(datetime dt)"""
        self.thisptr = new DateTime()
        self.set_datetime(dt)

    def __dealloc__(self):
        del self.thisptr

    def set_datetime(self, object dt):
        """set from python datetime"""
        if dt is None:
            return
        self.thisptr.Initialise(
            <int> dt.year, <int> dt.month, <int> dt.day,
            <int> dt.hour, <int> dt.minute, <int> dt.second, <int> dt.microsecond
            )

    def initialize(
            self,
            int year,int month, int day,
            int hour, int minute, int second, int microsecond
            ):
        """set from python datetime"""
        self.thisptr.Initialise(
            <int> year, <int> month, <int> day,
            <int> hour, <int> minute, <int> second, <int> microsecond
            )

    #@classmethod
    #def now(cls):
        #dt = PyDateTime()
        #cdef DateTime _dt = DateTime_Now(<cpp_bool> False)
        ##dt.thisptr = deref(_dt)
        #return dt

    def asString(self):
        return str(self.thisptr.ToString())

    def __str__(self):
        return self.asString()

    def __repr__(self):
        return self.asString()



cdef class PyTle(object):
    """Wrapper around (C++) Tle class"""

    cdef Tle *thisptr  # holds the C++ instance being wrapped

    def __cinit__(self, str name, str line_one, str line_two):
        """Constructor PyTle(str line_one, str line_two)"""
        self.thisptr = new Tle(<string> name, <string> line_one, <string> line_two)

    def __dealloc__(self):
        del self.thisptr

    def asString(self):
        return str(self.thisptr.ToString())

    def __str__(self):
        return self.asString()

    def __repr__(self):
        return '\n'.join([
            str(self.thisptr.Name()),
            str(self.thisptr.Line1()),
            str(self.thisptr.Line2())
            ])


cdef class PyObserver(object):
    """Wrapper around (C++) Observer class"""

    cdef Observer *thisptr  # holds the C++ instance being wrapped

    def __cinit__(
            self,
            double latitude_deg=50.525,
            double longitude_deg=6.883750,
            double altitude_km=0.319
            ):
        """Constructor PyObserver(double lon_deg, double lat_deg, double alt_km)"""
        self.thisptr = new Observer(latitude_deg, longitude_deg, altitude_km)

    def __dealloc__(self):
        del self.thisptr

    def asString(self):
        return str(self.thisptr.GetLocation().ToString())

    def __str__(self):
        return self.asString()

    def __repr__(self):
        return self.asString()

    def _loc(self):
        g = PyCoordGeodetic()
        g.thisptr = new CoordGeodetic(self.thisptr.GetLocation())
        return g

    def _setloc(self, PyCoordGeodetic loc):
        del self.thisptr
        self.thisptr = new Observer(loc.latitude, loc.longitude, loc.altitude)

    location = property(_loc, _setloc, None, 'location property')


cdef class PyCoordGeodetic(object):
    """Wrapper around (C++) CoordGeodetic struct"""

    cdef CoordGeodetic *thisptr  # holds the C++ instance being wrapped

    def __cinit__(
            self,
            double longitude_deg=0,
            double latitude_deg=0,
            double altitude_km=0
            ):
        """Constructor PyCoordGeodetic(double lon_deg, double lat_deg, double alt_km)"""
        self.thisptr = new CoordGeodetic()
        self._setlon(longitude_deg)
        self._setlat(latitude_deg)
        self._setalt(altitude_km)

    def __dealloc__(self):
        del self.thisptr

    def __str__(self):
        return ', '.join([
            str(self.longitude),
            str(self.latitude),
            str(self.altitude),
            ])

    def __repr__(self):
        return self.__str__()

    def _lon(self):
        return rad2deg * self.thisptr.longitude

    def _setlon(self, double lon_deg):
        self.thisptr.longitude = <double> (deg2rad * lon_deg)

    def _lat(self):
        return rad2deg * self.thisptr.latitude

    def _setlat(self, double lat_deg):
        self.thisptr.latitude = <double> (deg2rad * lat_deg)

    def _alt(self):
        return self.thisptr.altitude

    def _setalt(self, double alt_km):
        self.thisptr.altitude = <double> (alt_km)

    longitude = property(_lon, _setlon, None, 'longitude property')
    latitude = property(_lat, _setlat, None, 'latitude property')
    altitude = property(_alt, _setalt, None, 'altitude property')


cdef class PyCoordTopocentric(object):
    """Wrapper around (C++) CoordTopocentric struct"""

    cdef CoordTopocentric *thisptr  # holds the C++ instance being wrapped

    def __cinit__(
            self,
            double azimuth_deg=0,
            double elevation_deg=0,
            double distance_km=0,
            double distance_rate_km_per_s=0,
            ):
        """Constructor PyCoordTopocentric(
            double az_deg, double el_deg, double dist_km, double dist_rate_km_per_s,
            )"""
        self.thisptr = new CoordTopocentric()
        self._setaz(azimuth_deg)
        self._setel(elevation_deg)
        self._setdist(distance_km)
        self._setdist_rate(distance_rate_km_per_s)

    def __str__(self):
        return ', '.join([
            str(self.azimuth),
            str(self.elevation),
            str(self.distance),
            str(self.distance_rate),
            ])

    def __repr__(self):
        return self.__str__()

    def __dealloc__(self):
        del self.thisptr

    def _az(self):
        return rad2deg * self.thisptr.azimuth

    def _setaz(self, double az_deg):
        self.thisptr.azimuth = <double> (deg2rad * az_deg)

    def _el(self):
        return rad2deg * self.thisptr.elevation

    def _setel(self, double el_deg):
        self.thisptr.elevation = <double> (deg2rad * el_deg)

    def _dist(self):
        return self.thisptr.distance

    def _setdist(self, double dist_km):
        self.thisptr.distance = <double> (dist_km)

    def _dist_rate(self):
        return self.thisptr.distance_rate

    def _setdist_rate(self, double dist_rate_km_per_s):
        self.thisptr.distance_rate = <double> (dist_rate_km_per_s)

    azimuth = property(_az, _setaz, None, 'azimuth property')
    elevation = property(_el, _setel, None, 'elevation property')
    distance = property(_dist, _setdist, None, 'distance property')
    distance_rate = property(_dist_rate, _setdist_rate, None, 'distance_rate property')



cdef class PyEci(object):
    """Wrapper around (C++) Eci class"""

    cdef Eci *thisptr  # holds the C++ instance being wrapped

    def __cinit__(self, PyDateTime dt, PyCoordGeodetic geo):
        """Constructor PyEci(PyDateTime dt, PyCoordGeodetic geo)"""
        self.thisptr = new Eci(
            <DateTime> deref((<PyDateTime> dt).thisptr),
            <CoordGeodetic> deref((<PyCoordGeodetic> geo).thisptr),
            )

    def __dealloc__(self):
        del self.thisptr

    def asString(self):
        return '\n'.join([
            str(self.thisptr.GetDateTime().ToString()),
            str(self.thisptr.ToGeodetic().ToString())
            ])

    def __str__(self):
        return self.asString()

    def __repr__(self):
        return self.asString()


cdef class Satellite(object):
    """Interface to libsgp4, calculates apparent positions of satellite for a given observer and TLE"""

    cdef Eci *eci_ptr
    cdef PyCoordTopocentric _pytopo
    cdef PyCoordGeodetic _pygeo
    cdef CoordTopocentric _topo
    cdef CoordGeodetic _geo
    cdef SGP4 *sgp4_ptr
    cdef object _tle, _observer, _dt

    cdef double _mjd
    cdef python_bool _pos_dirty, _tle_dirty

    def __cinit__(self, PyTle tle, PyObserver observer=None):
        """Constructs a new Satellite object from given TLE, if observer is None,
        Effelsberg location is used"""

        self._tle = tle  # copy reference
        self._observer= observer
        try:
            self.sgp4_ptr = new SGP4(deref((<PyTle> tle).thisptr))
            self._tle_dirty = <python_bool> False
        except:
            print 'SatelliteException catched'
            self._tle_dirty = <python_bool> True
        self._dt = PyDateTime()  # empty datetime
        self.eci_ptr = new Eci(DateTime(), CoordGeodetic())  # empty Eci, otherwise
        # we cannot assign in _refresh_coords (segfault)
        self._pytopo = PyCoordTopocentric()  # this is important, otherwise self._topo will
        # just be None after _refresh_coords()
        self._pygeo = PyCoordGeodetic()

        self._pos_dirty = <python_bool> True

    def __dealloc__(self):
        del self.eci_ptr
        del self.sgp4_ptr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def set_new_mjd(self, double _mjd):
        assert _mjd < 1000000., 'warning, make sure you use mjd instead of jd'

        if fabs(self._mjd - _mjd) < MJD_RESOLUTION:
            return
        self._mjd = _mjd
        self._refresh_datetime()
        self._pos_dirty = <python_bool> True

    def MJD(self):
        return self._mjd

    mjd = property(MJD, set_new_mjd, None, 'mjd')  # python setter decorator doesn't work directly

    @property
    def datetime(self):
        return self._dt

    def topo_pos(self, double mjd):
        self.set_new_mjd(mjd)
        if self._pos_dirty:
            self._refresh_coords()
        if self._tle_dirty:
            return None
        return self._pytopo

    def geo_pos(self, double mjd):
        self.set_new_mjd(mjd)
        if self._pos_dirty:
            self._refresh_coords()
        if self._tle_dirty:
            return None
        return self._pygeo

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def _refresh_coords(self):
        try:
            self.eci_ptr[0] = self.sgp4_ptr.FindPosition(
            deref((<PyDateTime> self._dt).thisptr)
            )
            self._tle_dirty = <python_bool> False
        except:
            print 'SatelliteException catched'
            self._tle_dirty = <python_bool> True
            return
        self._topo = deref(
            (<PyObserver> self._observer).thisptr
            ).GetLookAngle(deref(self.eci_ptr))
        self._geo = self.eci_ptr.ToGeodetic()

        del self._pytopo.thisptr
        del self._pygeo.thisptr
        self._pytopo.thisptr = new CoordTopocentric(self._topo)
        self._pygeo.thisptr = new CoordGeodetic(self._geo)
        self._pos_dirty = <python_bool> False


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def _refresh_datetime(self):
        cdef mjd = self._mjd
        cdef int mn, yr
        cdef double dy
        cdef double tmp_mjd
        cdef double d, i, f, a, b, ce, g

        mn = 0
        yr = 0
        dy = 0.0

        tmp_mjd = mjd
        tmp_mjd -= +15020L - 0.5
        d = tmp_mjd + 0.5
        i = floor(d)
        f = d - i

        if f == 1:
            f = 0
            i += 1

        if i > -115860.0:
            a = floor(i / 36524.25 + .99835726) + 14
            i += 1 + a - floor(a / 4.0)

        b = floor(i / 365.25 + .802601)
        ce = i - floor(365.25 * b + .750001) + 416
        g = floor(ce / 30.6001)
        mn = <int> g - 1
        dy = <double> ce - floor(30.6001 * g) + f
        yr = <int> b + 1899

        if g > 13.5:
            mn = <int> g - 13
        if mn < 2.5:
            yr = <int> b + 1900
        if yr < 1:
            yr -= 1

        if tmp_mjd == 0:
            mn = 12
            dy = 31.5
            yr = 1899

        cdef int i_dy, i_hh, i_mm, i_ss, i_mus

        i_dy = <int> dy
        dy -= i_dy
        dy *= 24.
        i_hh = <int> dy
        dy -= i_hh
        dy *= 60.
        i_mm = <int> dy
        dy -= i_mm
        dy *= 60.
        i_ss = <int> dy
        dy -= i_ss
        dy *= 1.e6
        i_mus = <int> dy
        self._dt.initialize(yr, mn, i_dy, i_hh, i_mm, i_ss, i_mus)
