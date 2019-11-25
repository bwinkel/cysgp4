#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib
import os
# import requests
import pytest
import datetime
import numpy as np
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from numpy.testing import assert_equal, assert_allclose
from cysgp4 import *


TLE_ISS = (
    'ISS (ZARYA)',
    '1 25544U 98067A   13165.59097222  .00004759  00000-0  88814-4 0    47',
    '2 25544  51.6478 121.2152 0011003  68.5125 263.9959 15.50783143834295',
    )
MYSAT = (
    'MYSAT',
    '1 00001U 20001A   19306.08696759  .00000000  00000-0  50000-4 0    05',
    '2 00001  10.0000  35.0000 0001000   0.0000 112.0000  9.55934723    04',
    )


def test_tle_checksum():

    tle_line1 = TLE_ISS[1]
    tle_line2 = TLE_ISS[2]

    assert tle_line1[-1] == str(tle_checksum(tle_line1))
    assert tle_line2[-1] == str(tle_checksum(tle_line2))

    tle_text = get_example_tles()
    tle_tuples = tle_tuples_from_text(tle_text)

    for _, tle_line1, tle_line2 in tle_tuples:

        assert tle_line1[-1] == str(tle_checksum(tle_line1))
        assert tle_line2[-1] == str(tle_checksum(tle_line2))


# def test_tle_tuples_from_text():  # covered by test_tle_checksum


def test_tles_from_text():

    tle_text = get_example_tles()
    tles = tles_from_text(tle_text)

    assert repr(tles[0]) == '<PyTle: AKEBONO (EXOS-D)        >'
    assert repr(tles[-1]) == '<PyTle: ZHANGZHENG-1 (CSES)     >'


def test_satellite_mean_motion():

    assert_allclose(satellite_mean_motion(500.), 15.21937835)


def test_tle_linestrings_from_orbital_parameters():

    # Define satellite orbital parameters
    sat_name, sat_nr = 'MYSAT', 1
    alt_km = 3000.  # satellite altitude
    mean_motion = satellite_mean_motion(alt_km)
    inclination = 10.  # deg
    raan = 35.  # deg
    mean_anomaly = 112.  # deg

    # assume, the parameters are valid for the following time
    dt = datetime.datetime(2019, 11, 2, 2, 5, 14)
    pydt = PyDateTime(dt)

    tle_tuple = tle_linestrings_from_orbital_parameters(
        sat_name,
        sat_nr,
        pydt.mjd,
        inclination,
        raan,
        mean_anomaly,
        mean_motion,
        )
    assert tle_tuple == MYSAT
