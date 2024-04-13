# -*- coding: utf-8 -*-
# cython: language_level=3
# cython: cdivision=True, boundscheck=False, wraparound=False
# cython: embedsignature=True
# distutils: language = c++

# ####################################################################
#
# title                  :cysgp4.pxd
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

cimport cython
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdint cimport uint32_t, int64_t


cdef extern from 'Globals.h':

    cdef:
        const double kAE
        const double kQ0
        const double kS0
        const double kMU
        const double kXKMPER  # Earth Radius (major aka "a")
        const double kXJ2
        const double kXJ3
        const double kXJ4
        const double kXKE
        const double kCK2
        const double kCK4
        const double kQOMS2T
        const double kS
        const double kPI
        const double kTWOPI
        const double kTWOTHIRD
        const double kTHDT
        const double kF  # Earth flattening
        const double kOMEGA_E
        const double kAU
        const double kSECONDS_PER_DAY
        const double kMINUTES_PER_DAY
        const double kHOURS_PER_DAY
        # const double kEPOCH_JAN1_00H_1900
        # const double kEPOCH_JAN1_12H_1900
        # const double kEPOCH_JAN1_12H_2000
        const double kA3OVK2


cdef extern from 'Tle.h':

    cdef cppclass Tle:

        Tle(
            const string & name,
            const string & line_one,
            const string & line_two
            ) except + nogil
        Tle(const Tle& tle) except + nogil
        string ToString() nogil const
        string Name() nogil const
        string Line1() nogil const
        string Line2() nogil const
        string IntDesignator() nogil const
        unsigned int NoradNumber() nogil const
        DateTime Epoch() nogil const
        double MeanMotionDt2() nogil const
        double MeanMotionDdt6() nogil const
        double BStar() const
        double Inclination(bool in_degrees) nogil const
        double RightAscendingNode(const bool in_degrees) nogil const
        double Eccentricity() nogil const
        double ArgumentPerigee(const bool in_degrees) nogil const
        double MeanAnomaly(const bool in_degrees) nogil const
        double MeanMotion() nogil const
        unsigned int OrbitNumber() nogil const


cdef extern from 'SGP4.h':

    # Note: in order to avoid Exceptions in propagate_many, we added another
    # implementation of FindPosition, which returns NaN-filled results
    # Have to add the following to SGP4.h:

    #     Eci FindPositionNaN(const DateTime& date) const;

    # and SGP4.cpp:

    #    Eci SGP4::FindPositionNaN(const DateTime& dt) const
    #    {
    #        try {
    #            return FindPosition((dt - elements_.Epoch()).TotalMinutes());
    #        } catch(...) {
    #            const Eci nan_eci(DateTime(0), nan(""), nan(""), nan(""));
    #            return nan_eci;
    #        }
    #    }

    cdef cppclass SGP4:
        SGP4(const Tle& tle) except + nogil

        Eci FindPosition(const DateTime& date) except + nogil
        Eci FindPositionNaN(const DateTime& date) nogil


cdef extern from 'Observer.h':

    # Note, in order to create an instance of Observer on the stack, a
    # null constructor needs to be present; have to add the following
    # to Observer.h:

    #     Observer()
    #         : m_geo(0., 0., 0.),
    #         m_eci(DateTime(), m_geo)
    #     {
    #     }

    cdef cppclass Observer:

        Observer() nogil
        Observer(
            const double latitude,
            const double longitude,
            const double altitude
            ) nogil
        Observer(const CoordGeodetic &geo) nogil
        CoordTopocentric GetLookAngle(const Eci &eci) nogil
        void SetLocation(const CoordGeodetic& geo) nogil
        CoordGeodetic GetLocation() nogil const


cdef extern from 'Vector.h':

    cdef cppclass Vector:

        Vector() noexcept nogil
        Vector(
            const double arg_x, const double arg_y, const double arg_z
            ) noexcept nogil

        double x
        double y
        double z


cdef extern from 'DateTime.h':

    # Note: for Win_X86 compatibility, one needs to remove the
    #     static DateTime Now(bool microseconds = false)
    # method from DateTime; furthermore the line
    #     day = std::min(day, maxday);
    # should be replaced with
    #     if (maxday < day)
    #         day = maxday;

    cdef cppclass DateTime:

        DateTime() nogil
        DateTime(int64_t ticks) nogil
        DateTime(unsigned int year, double doy) nogil
        void Initialise(
            int year, int month, int day,
            int hour, int minute, int second, int microsecond
            ) nogil
        double ToGreenwichSiderealTime() nogil const
        double ToLocalMeanSiderealTime(const double lon) nogil const
        string ToString() nogil const
        DateTime AddTicks(int64_t ticks) nogil const
        int64_t Ticks() nogil const

        int DayOfYear(int year, int month, int day) const
        int Year() nogil const
        int Month() nogil const
        int Day() nogil const
        int Hour() nogil const
        int Minute() nogil const
        int Second() nogil const
        int Microsecond() nogil const

    # how to wrap static c++ member functions???
    # cdef DateTime_Now 'DateTime::Now' (bool microseconds)


# cdef extern from 'DateTime.h' namespace 'DateTime':
    # DateTime Now(bool microseconds)


cdef extern from 'Eci.h':

    # Note, in order to create an instance of Eci on the stack, a
    # null constructor needs to be present; have to add the following
    # to Eci.h:

    #     Eci()
    #         : m_dt(DateTime(0)),
    #         m_position(Vector())
    #     {
    #     }

    cdef cppclass Eci:

        Eci() nogil
        Eci(
            const DateTime& dt, const double latitude, const double longitude,
            const double altitude
            ) nogil
        Eci(const DateTime& dt, const CoordGeodetic& geo) nogil
        Eci(const DateTime& dt, const Vector &position) nogil
        Eci(
            const DateTime& dt,
            const Vector &position,
            const Vector &velocity
            ) nogil
        void Update(const DateTime& dt, const CoordGeodetic& geo) nogil
        DateTime GetDateTime() nogil const
        CoordGeodetic ToGeodetic() nogil const
        Vector Position() nogil const
        Vector Velocity() nogil const


cdef extern from 'CoordTopocentric.h':

    cdef cppclass CoordTopocentric:

        CoordTopocentric() nogil
        CoordTopocentric(const CoordTopocentric& topo) nogil
        CoordTopocentric(
            double az, double el, double rnge, double rnge_rate
            ) nogil
        string ToString() nogil const
        # azimuth in radians
        double azimuth
        # elevations in radians
        double elevation
        # distance in kilometers
        double distance 'range'
        # range rate in kilometers per second
        double distance_rate 'range_rate'


cdef extern from 'CoordGeodetic.h':

    cdef cppclass CoordGeodetic:

        CoordGeodetic() nogil
        CoordGeodetic(const CoordGeodetic& geo) nogil
        CoordGeodetic(
            double lat, double lon, double alt, bool is_radians=False
            ) nogil
        string ToString() nogil const
        # latitude in radians (-PI >= latitude < PI)
        double latitude
        # latitude in radians (-PI/2 >= latitude <= PI/2)
        double longitude
        # altitude in kilometers
        double altitude


# https://stackoverflow.com/questions/51006230/dynamically-sized-array-of-objects-in-cython
cdef extern from *:
    '''
    template <typename T>
    T* array_new(int n) {
        return new T[n];
    }

    template <typename T>
    void array_delete(T* x) {
        delete [] x;
    }
    '''
    T* array_new[T](int)
    void array_delete[T](T* x)


cdef int64_t ticks_from_mjd(double mjd) noexcept nogil
cdef double mjd_from_ticks(int64_t ticks) noexcept nogil
cdef DateTime datetime_from_mjd(double mjd) noexcept nogil
cdef (double, double, double) ecef_from_geo(
    double lon_rad, double lat_rad, double alt_km
    ) noexcept nogil
