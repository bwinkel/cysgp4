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


from sgp4 cimport *
from cython.operator cimport dereference as deref, preincrement as inc
from cpython cimport bool as python_bool
from libcpp cimport bool as cpp_bool
from libc.math cimport M_PI

cdef double deg2rad = M_PI / 180.
cdef double rad2deg = 180. / M_PI


cdef class PyDateTime(object):
    """Wrapper around (C++) DateTime class"""

    cdef DateTime *thisptr  # holds the C++ instance being wrapped

    def __cinit__(
            self,
            int year=0, int month=0, int day=0,
            int hour=0, int minute=0, int second=0
            ):
        """Constructor PyTle(str line_one, str line_two)"""
        self.thisptr = new DateTime(
            <int> year, <int> month, <int> day,
            <int> hour, <int> minute, <int> second
            )

    def __dealloc__(self):
        del self.thisptr

    @classmethod
    def now(cls):
        dt = PyDateTime()
        cdef DateTime _dt = DateTime_Now(<cpp_bool> False)
        #dt.thisptr = deref(_dt)
        return dt

    #def asString(self):
        #return str(self.thisptr.ToString())

    #def __str__(self):
        #return self.asString()

    #def __repr__(self):
        #return str(self.thisptr.Name()), str(self.thisptr.Line1()), str(self.thisptr.Line2())



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
        return str(self.thisptr.Name()), str(self.thisptr.Line1()), str(self.thisptr.Line2())


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
        self.thisptr = new Observer(longitude_deg, latitude_deg, altitude_km)

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




