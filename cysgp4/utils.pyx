# -*- coding: utf-8 -*-
# cython: language_level=3
# cython: cdivision=True, boundscheck=False, wraparound=False
# cython: embedsignature=True
# distutils: language = c++

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
cimport numpy as np
# from libc.math cimport M_PI, floor, fabs, fmod, sqrt, sin, cos
from libc.math cimport M_PI, sqrt
from libc.stdlib cimport atoi
from libc.string cimport strcpy, strcmp, memcpy
from .cysgp4 cimport kMU, kXKMPER  # aka astropy.constants.{GM_earth,R_earth}
from .cysgp4 cimport kSECONDS_PER_DAY
from .cysgp4 cimport DateTime, datetime_from_mjd

import numpy as np
from .cysgp4 import PyTle

np.import_array()


__all__ = [
    'tle_checksum', 'tle_linestrings_from_text', 'tles_from_text',
    'satellite_mean_motion', 'tle_linestrings_from_orbital_parameters',
    ]


cpdef int tle_checksum(str line):

    cdef:

        unsigned int acc = 0
        unsigned int i, u, v, w
        tmp = line.encode('UTF-8')
        const char *c_line = tmp
        char buf[2]

    strcpy(buf, b'X')  # initialize

    if len(line) != 69:
        raise ValueError('TLE line string must have 69 characters')

    acc = 0
    for i in range(68):
        memcpy(buf, &(c_line[i]), 1)
        u = atoi(buf)
        acc += atoi(buf)
        if strcmp(buf, b'-') == 0:
            acc += 1

    return acc % 10


def tle_linestrings_from_text(str tle_text):

    all_lines = tle_text.split('\r\n')
    tle_list = list(zip(*tuple(all_lines[idx::3] for idx in range(3))))

    return tle_list


def tles_from_text(str tle_text):

    return [
        PyTle(*tle)
        for tle in tle_linestrings_from_text(tle_text)
        ]


cpdef double satellite_mean_motion(
        double alt_km, double mu_km3_s2=kMU, double r_earth_km=kXKMPER
        ) nogil:
    '''
    Mean motion of satellite at altitude in Earth's gravitational field.

    See https://en.wikipedia.org/wiki/Mean_motion#Formulae
    '''

    cdef:
        double s_per_day = kSECONDS_PER_DAY
        double no = sqrt(
            4.0 * M_PI ** 2 * (alt_km + r_earth_km) ** 3 / mu_km3_s2
            ) / s_per_day

    return 1 / no


cpdef tuple tle_linestrings_from_orbital_parameters(
        str sat_name,
        int sat_nr,
        double mjd_epoch,
        double inclination_deg,
        double raan_deg,
        double mean_anomaly_deg,
        double mean_motion_per_day,
        ):
    '''
    Generate TLE strings from orbital parameters.

    Note: epoch has a very strange format: first two digits are the year, next three
    digits are the day from beginning of year, then fraction of a day is given, e.g.
    20180.25 would be 2020, day 180, 6 hours (UT?)
    '''

    # Note: RAAN = right ascention (or longitude) of ascending node

    cdef:

        DateTime dt = datetime_from_mjd(mjd_epoch)
        int year = dt.Year()
        double epoch = (
            (year % 100) * 1000. +
            dt.DayOfYear(dt.Year(), dt.Month(), dt.Day()) +
            dt.Hour() / 24. +
            dt.Minute() / 1440. +
            dt.Second() / 86400. +
            dt.Microsecond() * 1.e-6
            # mjd_epoch - <unsigned int> mjd_epoch
            )

        str tmp1 = (
        '1 {:05d}U 20001A   {:14.8f}  .00000000  00000-0  50000-4 '
        '0    0X'.format(sat_nr, epoch)
        )
        str tmp2 = (
        '2 {:05d} {:8.4f} {:8.4f} 0001000   0.0000 {:8.4f} '
        '{:11.8f}    0X'.format(
            sat_nr, inclination_deg, raan_deg,
            mean_anomaly_deg, mean_motion_per_day
        ))

        str tle1 = '{:s}{:d}'.format(tmp1[:-1], tle_checksum(tmp1))
        str tle2 = '{:s}{:d}'.format(tmp2[:-1], tle_checksum(tmp2))

    if year < 1957 or year > 2056:
        raise ValueError('Year must be between 1957 and 2056')

    return (sat_name, tle1, tle2)
