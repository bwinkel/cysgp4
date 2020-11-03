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
from libc.math cimport M_PI, sqrt, modf
from libc.stdlib cimport atoi
from libc.string cimport strcpy, strcmp, memcpy
from .cysgp4 cimport kMU, kXKMPER  # aka astropy.constants.{GM_earth,R_earth}
from .cysgp4 cimport kSECONDS_PER_DAY
from .cysgp4 cimport DateTime, datetime_from_mjd

import numpy as np
from .cysgp4 import PyTle

np.import_array()


__all__ = [
    'tle_checksum', 'tle_tuples_from_text', 'tles_from_text',
    'satellite_mean_motion', 'tle_linestrings_from_orbital_parameters',
    ]


cpdef int tle_checksum(str line):
    '''
    Compute checksum of a `TLE line string
    <https://en.wikipedia.org/wiki/Two-line_element_set>`_.

    Parameters
    ----------
    line : str
        TLE line string.

    Returns
    -------
    chksum : int
        Checksum (1-digit).

    Examples
    --------
    Assume you want to check if a TLE line string is correct::

        >>> import cysgp4

        >>> tle_line = '1 25544U 98067A   13165.59097222  .00004759  00000-0  88814-4 0    47'
        >>> chksum = cysgp4.tle_checksum(tle_line)
        >>> tle_line[-1] == str(chksum)
        True

    '''

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


def tle_tuples_from_text(str tle_text):
    '''
    Split a text containing `TLE strings
    <https://en.wikipedia.org/wiki/Two-line_element_set>`_ (including Name
    or ID, i.e., three lines per entry) into list of TLE tuples.

    Parameters
    ----------
    tle_text : str
        Text containing TLE line strings.

    Returns
    -------
    tle_list : list of (str, str, str)
        List of TLE tuples, each consisting of the line strings.

    Examples
    --------
    An example text file with TLE strings is included in `~cysgp4`
    and can be parsed like this::

        >>> import cysgp4

        >>> tle_text = cysgp4.get_example_tles()
        >>> tle_tuples = cysgp4.tle_tuples_from_text(tle_text)
        >>> tle_tuples
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

    '''

    all_lines = tle_text.replace('\r', '').split('\n')
    tle_list = list(zip(*tuple(all_lines[idx::3] for idx in range(3))))

    return tle_list


def tles_from_text(str tle_text):
    '''
    Parse a text containing `TLE strings
    <https://en.wikipedia.org/wiki/Two-line_element_set>`_ (including Name
    or ID, i.e., three lines per entry) into list of `~cysgp4.PyTle` objects.

    Parameters
    ----------
    tle_text : str
        Text containing TLE line strings.

    Returns
    -------
    tle_list : list of `~cysgp4.PyTle` objects
        List of `~cysgp4.PyTle` objects.

    Examples
    --------
    An example text file with TLE strings is included in `~cysgp4`
    and can be parsed like this::

        >>> import cysgp4

        >>> tle_text = cysgp4.get_example_tles()
        >>> tle_tuples = cysgp4.tles_from_text(tle_text)
        [<PyTle: AKEBONO (EXOS-D)        >,
         <PyTle: HST                     >,
        ...
         <PyTle: ZHANGZHENG-1 (CSES)     >]

    '''

    return [
        PyTle(*tle)
        for tle in tle_tuples_from_text(tle_text)
        ]


cpdef double satellite_mean_motion(
        double alt_km, double mu_km3_s2=kMU, double r_earth_km=kXKMPER
        ) nogil:
    '''
    Mean motion of satellite at altitude in Earth's gravitational field.

    See https://en.wikipedia.org/wiki/Mean_motion#Formulae

    Note the parameters are such that they fit to the WGS72 definition.

    Parameters
    ----------
    alt_km : float
        Satellite altitude above Earth's surface.
    mu_km3_s2 : float, optional (default: 398600.8)
        Nominal Earth mass parameter [km^3 / s^2].
    r_earth_km : float, optional (default: 6378.135)
        Earth's radius [km].

    Returns
    -------
    mean_motion : double
        Mean motion of satellite in revolutions per day [1 / day].

    Examples
    --------
    An example text file with TLE strings is included in `~cysgp4`
    and can be parsed like this::

        >>> import cysgp4

        >>> mean_motion = cysgp4.satellite_mean_motion(500.)
        >>> mean_motion
        15.219378350934464

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
        double eccentricity,
        double argument_of_perigee_deg,
        double mean_anomaly_deg,
        double mean_motion_per_day,
        ):
    '''
    Generate TLE strings from orbital parameters.

    See `TLE line strings
    <https://en.wikipedia.org/wiki/Two-line_element_set>`_ for a description
    of the two-line element format. The parameters are the `Orbital/Keplarian
    elements <https://en.wikipedia.org/wiki/Orbital_elements>`_ and Wikipedia
    has a really good description what each of them means.

    Note: The epoch (date/time) used in TLEs has a very strange format: first
    two digits are the year, next three digits are the day from beginning of
    year, then fraction of a day is given, e.g. 20180.25 would be 2020, day
    180, 6 hours (UTC).

        Parameters
    ----------
    sat_name : str
        The satellite name.
    sat_nr : str
        The satellite number.
    mjd_epoch : float
        The epoch of the orbital parameters as MJD [day].
        The function will convert this to the TLE epoch format.
    inclination_deg : float
        Inclination of the orbit [deg].
    raan_deg : float
        Right ascension of the ascending node of the orbit [deg].
    eccentricity : float
        Eccentricity of the orbit [dimensionless].
        Note that in the TLE only the decimal digits are stored.
    argument_of_perigee_deg : float
        Argument of Perigee [deg].
    mean_anomaly_deg : float
        Mean anomaly of the node [deg].
    mean_motion_per_day : float
        Mean motion of the satellite [1 / day].

    Returns
    -------
    tle_lines : (str, str, str)
        TLE tuple of the satellite, consisting of the line strings.

    Examples
    --------
    A simple use case would be like::

        >>> from datetime import datetime
        >>> import cysgp4

        >>> # Define satellite orbital parameters
        >>> sat_name, sat_nr = 'MYSAT', 1
        >>> alt_km = 3000.  # satellite altitude
        >>> mean_motion = cysgp4.satellite_mean_motion(alt_km)
        >>> inclination = 10.  # deg
        >>> raan = 35.  # deg
        >>> eccentricity = 0.0001
        >>> argument_of_perigee = 0.  # deg
        >>> mean_anomaly = 112.  # deg

        >>> # assume, the parameters are valid for the following time
        >>> dt = datetime(2019, 11, 2, 2, 5, 14)
        >>> pydt = cysgp4.PyDateTime(dt)
        >>> pydt
        <PyDateTime: 2019-11-02 02:05:14.000000 UTC>
        >>> mjd_epoch = pydt.mjd

        >>> tle_tuple = cysgp4.tle_linestrings_from_orbital_parameters(
        ...     sat_name,
        ...     sat_nr,
        ...     mjd_epoch,
        ...     inclination,
        ...     raan,
        ...     eccentricity,
        ...     argument_of_perigee,
        ...     mean_anomaly,
        ...     mean_motion,
        ...     )
        >>> tle_tuple
        ('MYSAT',
         '1 00001U 20001A   19306.08696759  .00000000  00000-0  50000-4 0    05',
         '2 00001  10.0000  35.0000 0001000   0.0000 112.0000  9.55934723    04')

    '''

    if sat_nr < 0 or sat_nr > 99999:
        raise ValueError('sat_nr must be a 5-digit number')
    if eccentricity < 0.:
        raise ValueError('Eccentricity must be >= 0.')
    if eccentricity >= 1.:
        raise ValueError('Eccentricity must be < 1.')

    cdef:

        DateTime dt = datetime_from_mjd(mjd_epoch)
        int year = dt.Year()

    if year < 1957 or year > 2056:
        raise ValueError('Year must be between 1957 and 2056')

    cdef:

        double epoch = (
            (year % 100) * 1000. +
            dt.DayOfYear(dt.Year(), dt.Month(), dt.Day()) +
            dt.Hour() / 24. +
            dt.Minute() / 1440. +
            dt.Second() / 86400. +
            dt.Microsecond() * 1.e-6
            )
        double _i
        int ecc = int(modf(eccentricity, &_i) * 1e7 + 0.5)

        str tmp1 = (
        '1 {:05d}U 20001A   {:14.8f}  .00000000  00000-0  50000-4 '
        '0    0X'.format(sat_nr, epoch)
        )
        str tmp2 = (
        '2 {:05d} {:8.4f} {:8.4f} {:07d} {:8.4f} {:8.4f} '
        '{:11.8f}    0X'.format(
            sat_nr, inclination_deg, raan_deg, ecc, argument_of_perigee_deg,
            mean_anomaly_deg, mean_motion_per_day
        ))

        str tle1 = '{:s}{:d}'.format(tmp1[:68], tle_checksum(tmp1))
        str tle2 = '{:s}{:d}'.format(tmp2[:68], tle_checksum(tmp2))

    return (sat_name, tle1, tle2)
