#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib
import pytest
import datetime
import numpy as np
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from numpy.testing import assert_equal, assert_allclose
from cysgp4 import *


# skip over pycraf related tests, if not package present:
skip_pycraf = pytest.mark.skipif(
    importlib.util.find_spec('pycraf') is None,
    reason='"pycraf" package not installed'
    )


class TestPyDateTime:

    def setup(self):

        self.dt = datetime.datetime(2019, 1, 1, 12, 13, 14)

    def teardown(self):

        pass

    def test_constructor(self):

        t1 = PyDateTime(self.dt)
        assert str(t1) == '2019-01-01 12:13:14.000000 UTC'

    def test_times(self):

        t1 = PyDateTime(self.dt)
        assert_allclose(t1.gmst(), 4.95971515)

        assert_allclose(t1.lmst(6.88375), 5.07985925)

    def test_ticks(self):

        t1a = PyDateTime(self.dt)
        t1b = PyDateTime.from_ticks(63681941594000000)
        dt2 = datetime.datetime(1, 1, 1, 0, 0, 0)  # zero ticks
        t2a = PyDateTime(dt2)
        t2b = PyDateTime.from_ticks(0)
        dt3 = datetime.datetime(1858, 11, 17, 0, 0, 0)  # MJD-0 ticks
        t3a = PyDateTime(dt3)
        t3b = PyDateTime.from_ticks(58628880000000000)

        assert t1a.ticks == 63681941594000000
        assert t1b.ticks == 63681941594000000
        assert t2a.ticks == 0
        assert t2b.ticks == 0
        assert t3a.ticks == 58628880000000000
        assert t3b.ticks == 58628880000000000


class TestPyTle:

    def setup(self):

        self.tle_tup = (
            'ISS (ZARYA)',
            '1 25544U 98067A   13165.59097222  .00004759  00000-0  '
            '88814-4 0    47',
            '2 25544  51.6478 121.2152 0011003  68.5125 263.9959 '
            '15.50783143834295',
            )
        self.tle1_str = '''Norad Number:         25544
            Int. Designator:      98067A
            Epoch:                2013-06-14 14:10:59.999800 UTC
            Orbit Number:         83429
            Mean Motion Dt2:        0.00004759
            Mean Motion Ddt6:       0.00000000
            Eccentricity:           0.00110030
            BStar:                  0.00008881
            Inclination:           51.64780000
            Right Ascending Node: 121.21520000
            Argument Perigee:      68.51250000
            Mean Anomaly:         263.99590000
            Mean Motion:           15.50783143'''.split('\n')

    def teardown(self):

        pass

    def test_constructor(self):

        tle1 = PyTle(*self.tle_tup)
        slist = str(tle1).split('\n')
        for s1, s2 in zip(slist, self.tle1_str):
            assert s1.strip() == s2.strip()

        assert repr(tle1) == '<PyTle: ISS (ZARYA)>'


class TestPyCoordGeodetic:

    def setup(self):

        self.geo_coord = (6.88375, 50.525, 0.366)

    def teardown(self):

        pass

    def test_constructor(self):

        geo1 = PyCoordGeodetic(*self.geo_coord)
        assert str(geo1) == '6.8838d, 50.5250d, 0.3660km'


class TestPyCoordTopocentric:

    def setup(self):

        self.topo_coord = (30., 45., 900., 1.)

    def teardown(self):

        pass

    def test_constructor(self):

        topo1 = PyCoordTopocentric(*self.topo_coord)
        assert str(topo1) == '30.0000d, 45.0000d, 900.0000km, 1.0000km/s'


class TestPyEci:

    def setup(self):

        self.eci_dt = datetime.datetime(2019, 1, 1, 12, 13, 14)
        self.geo_coord = (6.88376, 50.525, 0.366)

    def teardown(self):

        pass

    def test_constructor(self):

        t1 = PyDateTime(self.eci_dt)
        geo1 = PyCoordGeodetic(*self.geo_coord)
        eci1 = PyEci(t1, geo1)
        des = '6.8838d, 50.5250d, 0.3660km 2019-01-01 12:13:14.000000 UTC'
        assert str(eci1) == des


class TestPyObserver:

    def setup(self):

        self.effbg_observer = (6.88375, 50.525, 0.366)

    def teardown(self):

        pass

    def test_constructor(self):

        geo1 = PyObserver(*self.effbg_observer)
        assert str(geo1) == '6.8838d, 50.5250d, 0.3660km'


class TestSatellite:

    def setup(self):

        self.tle_tup = (
            'ISS (ZARYA)',
            '1 25544U 98067A   13165.59097222  .00004759  00000-0  '
            '88814-4 0    47',
            '2 25544  51.6478 121.2152 0011003  68.5125 263.9959 '
            '15.50783143834295',
            )
        self.tle = PyTle(*self.tle_tup)

        self.dt_tup = (2013, 6, 15, 2, 57, 7, 200000)
        self.dt = PyDateTime(datetime.datetime(*self.dt_tup))
        self.mjd = 56458.123
        self.effbg_tup = (6.88375, 50.525, 0.366)
        self.effbg_tup_m = (6.88375, 50.525, 366.)
        self.effbg_observer = PyObserver(*self.effbg_tup)

        self.sat2 = twoline2rv(self.tle_tup[1], self.tle_tup[2], wgs72)
        self.pos2, self.vel2 = self.sat2.propagate(
            *self.dt_tup[:-2], self.dt_tup[-2] + self.dt_tup[-1] / 1e6
            )

    def teardown(self):

        pass

    def test_constructor(self):

        Satellite(self.tle, self.effbg_observer)

    def test_eci_position(self):

        sat = Satellite(self.tle, self.effbg_observer)
        # print(sat.dt)
        eci_pos = sat.eci_pos(self.mjd)

        print(eci_pos.loc)
        print(self.pos2)
        print(eci_pos.vel)
        print(self.vel2)

        assert_allclose(
            eci_pos.loc, (-4728.184444, 730.892025, 4802.515276)
            )
        assert_allclose(
            eci_pos.vel, (1.530454442, -7.065375488, 2.574384704)
            )
        assert_allclose(eci_pos.loc, self.pos2, atol=1e-2)
        assert_allclose(eci_pos.vel, self.vel2, atol=1e-5)

    def test_geo_position(self):

        sat = Satellite(self.tle, self.effbg_observer)
        # print(sat.dt)
        eci_pos = sat.eci_pos(self.mjd)
        geo_pos = sat.geo_pos(self.mjd)

        print(eci_pos.geo_loc)
        print(geo_pos.lon, geo_pos.lat, geo_pos.alt)

        assert_allclose(
            (geo_pos.lon, geo_pos.lat, geo_pos.alt),
            (-136.627536, 45.289346, 411.566712)
            )

    def test_topo_position(self):

        sat = Satellite(self.tle, self.effbg_observer)
        # print(sat.dt)
        topo_pos = sat.topo_pos(self.mjd)

        print(topo_pos.az, topo_pos.el, topo_pos.dist)

        assert_allclose(
            (topo_pos.az, topo_pos.el, topo_pos.dist),
            (334.789646, -37.384929, 8406.367773)
            )

    @skip_pycraf
    def test_topo_position_vs_pycraf(self):

        from astropy.coordinates import EarthLocation
        from astropy import time
        from pycraf import satellite

        tle_string = '\n'.join(self.tle_tup)
        location = EarthLocation.from_geodetic(*self.effbg_tup_m, 'WGS72')
        sat_obs = satellite.SatelliteObserver(location)
        dt = datetime.datetime(*self.dt_tup)
        obstime = time.Time(dt)
        gmst = obstime.sidereal_time('mean', 'greenwich').rad
        lmst = gmst + location.lon.rad
        print(obstime)
        print(gmst, lmst)
        az, el, dist = sat_obs.azel_from_sat(tle_string, obstime)
        az, el, dist = az.value, el.value, dist.value
        az %= 360.
        print(az, el, dist)

        sat = Satellite(self.tle, self.effbg_observer)
        # print(sat.dt)
        eci_pos = sat.eci_pos(self.mjd)
        topo_pos = sat.topo_pos(self.mjd)

        print(eci_pos.dt)
        print(eci_pos.dt.gmst(), eci_pos.dt.lmst(self.effbg_tup[0]))
        print(topo_pos.az, topo_pos.el, topo_pos.dist)

        assert_allclose(topo_pos.az, az, atol=1e-3)
        assert_allclose(topo_pos.el, el, atol=1e-3)
        assert_allclose(topo_pos.dist, dist, atol=2e-2)
