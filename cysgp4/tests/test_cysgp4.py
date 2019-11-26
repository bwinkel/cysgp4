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


# skip over pycraf related tests, if not package present:
skip_pycraf = pytest.mark.skipif(
    importlib.util.find_spec('pycraf') is None,
    reason='"pycraf" package not installed'
    )
# if you want to test this locally: pytest -p no:benchmark cysgp4
skip_benchmark = pytest.mark.skipif(
    importlib.util.find_spec('benchmark') is None,
    reason='"benchmark" fixture not installed'
    )

TLE_ISS = (
    'ISS (ZARYA)',
    '1 25544U 98067A   13165.59097222  .00004759  00000-0  88814-4 0    47',
    '2 25544  51.6478 121.2152 0011003  68.5125 263.9959 15.50783143834295',
    )
TLE_GPS = (
    'GPS BIIR-2  (PRN 13)',
    '1 24876U 97035A   19309.58152857 -.00000011  00000-0  00000+0 0  9996',
    '2 24876  55.4542 192.9394 0037899  66.9931 293.4794  2.00564219163301',
    )
# TLE_MMS = (  # this fails if MJD is very old
#     'MMS 1',
#     '1 40482U 15011A   19320.57120432 -.00000010  00000-0  00000+0 0  9996',
#     '2 40482  21.0630 120.0377 8869973  16.2077  55.8249  0.28561246 10609',
#     )
TLE_ECC_ERR = (  # if MJD is far in the future: Error: (e <= -0.001)
    'POLAR',
    '1 23802U 96013A   19320.59223231  .00000075  00000-0  00000+0 0  9997',
    '2 23802  78.7075 252.8203 6540600 290.3797  13.5942  1.29845581113610',
    )
TLE_ELSQ_ERR = (  # if MJD is far in the future: Error: (elsq >= 1.0)
    'CXO',
    '1 25867U 99040B   19322.60140111  .00000856  00000-0  00000-0 0  9990',
    '2 25867  65.7630 272.2977 7416434 231.1198   0.0919  0.37806068  8260',
    )
# TLE_DECAYED_ERR = (  # if MJD is far in the future: Error: Satellite decayed
#     'SWIFT',
#  '1 28485U 04047A   19321.65306734  .00000791  00000-0  23585-4 0  9996',
#  '2 28485  20.5570  74.4329 0010989 245.6346 114.2940 15.04483561821407'
#     )


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

    def test_datetimes(self):

        t1 = PyDateTime(self.dt)
        assert t1.datetime == self.dt

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

    def test_mjd(self):

        mjd = 56458.123
        t1 = PyDateTime.from_mjd(mjd)

        assert str(t1) == '2013-06-15 02:57:07.199999 UTC'
        assert_allclose(mjd, t1.mjd)

    def test_tle_epoch(self):

        # Note: this is an oddity with the TLE time format; if day is
        # zero, it is really the last day of the previous year...
        t = PyDateTime.from_tle_epoch(19000.)
        assert str(t) == '2018-12-31 00:00:00.000000 UTC'

        dt = datetime.datetime(2018, 12, 31, 0, 0, 0)
        t = PyDateTime(dt)
        print(t)
        assert_allclose(t.tle_epoch, 18365.)  # sic!

        tle_epoch = 19001.0
        t = PyDateTime.from_tle_epoch(tle_epoch)
        assert str(t) == '2019-01-01 00:00:00.000000 UTC'

        dt = datetime.datetime(2019, 1, 1, 0, 0, 0)
        t = PyDateTime(dt)
        print(t)
        assert_allclose(t.tle_epoch, tle_epoch)

        mjd = 58533.1
        tle_epoch = 19050.1
        dt = datetime.datetime(2019, 2, 19, 2, 24, 0)

        t = PyDateTime(dt)
        print(t)
        assert str(t) == '2019-02-19 02:24:00.000000 UTC'
        assert_allclose(t.tle_epoch, tle_epoch)

        t = PyDateTime.from_mjd(mjd)
        print(t)
        assert str(t) == '2019-02-19 02:23:59.999999 UTC'

        assert_allclose(t.tle_epoch, tle_epoch)

        t = PyDateTime.from_tle_epoch(tle_epoch)
        print(t)
        assert str(t) == '2019-02-19 02:24:00.000000 UTC'

        assert_allclose(t.mjd, mjd)


class TestPyTle:

    def setup(self):

        self.tle_tup = TLE_ISS
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

    def test_tle_strings(self):

        tle = PyTle(*self.tle_tup)
        for s1, s2 in zip(self.tle_tup, tle.tle_strings()):
            assert s1.strip() == s2.strip()

    def test_exception_checksum(self):
        '''
        Note, one would assume that a wrong checksum would raise an
        exception, but it doesn't in the underlying sgp4 code.
        '''

        tle_tup_wrong = (
            TLE_ISS[0],
            TLE_ISS[1][:-1] + '0',  # modify checksum
            TLE_ISS[2],
            )
        tle_wrong = PyTle(*tle_tup_wrong)
        slist = str(tle_wrong).split('\n')

        for s1, s2 in zip(slist, self.tle1_str):
            assert s1.strip() == s2.strip()

        assert repr(tle_wrong) == '<PyTle: ISS (ZARYA)>'

    # Note, there are more exceptions in the C++ code, which we could
    # test here, but it can be assumed that exception handling with
    # Cython will work similarly for other cases. We only need to test
    # that it works in general (i.e., for one example).
    def test_exception_invalid_line_length(self):

        tle_tup_wrong = (
            TLE_ISS[0],
            TLE_ISS[1][:-1],
            TLE_ISS[2],
            )
        with pytest.raises(RuntimeError) as excinfo:
            PyTle(*tle_tup_wrong)

        assert 'Invalid length for line one' in str(excinfo.value)

    def test_epoch(self):

        tle = PyTle(*self.tle_tup)

        assert_allclose(tle.epoch.mjd, 56457.590972)


class TestPyCoordGeodetic:

    def setup(self):

        self.geo_coord = (6.88375, 50.525, 0.366)

    def teardown(self):

        pass

    def test_constructor(self):

        geo1 = PyCoordGeodetic(*self.geo_coord)
        assert str(geo1) == '6.8838d, 50.5250d, 0.3660km'

    def test_ecef(self):

        geo = PyCoordGeodetic(*self.geo_coord)

        print(geo.ecef)
        assert_allclose(
            geo.ecef,
            (4033.8986696677143, 486.99458429181357, 4900.402205553537),
            atol=1.e-3
            )


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

        obs = PyObserver()
        assert str(obs) == '0.0000d, 0.0000d, 0.0000km'

        obs = PyObserver(*self.effbg_observer)
        assert str(obs) == '6.8838d, 50.5250d, 0.3660km'

    def test_location_property(self):

        obs = PyObserver(*self.effbg_observer)
        assert str(obs.loc) == '6.8838d, 50.5250d, 0.3660km'
        geo = PyCoordGeodetic(1, 2, 3)
        obs.loc = geo
        assert str(obs) == '1.0000d, 2.0000d, 3.0000km'


class TestSatellite:

    def setup(self):

        self.tle_tup = TLE_ISS
        self.tle = PyTle(*self.tle_tup)

        self.dt_tup = (2013, 6, 15, 2, 57, 7, 200000)
        self.pydt = PyDateTime(datetime.datetime(*self.dt_tup))
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

        Satellite(self.tle)
        Satellite(self.tle, self.effbg_observer)
        Satellite(self.tle, self.effbg_observer, self.pydt)
        Satellite(self.tle, self.effbg_observer, self.pydt, 1 / 86400.)

    def test_error(self):

        pydt_off = PyDateTime.from_mjd(68805.5)

        tle_ecc_err = PyTle(*TLE_ECC_ERR)
        with pytest.raises(RuntimeError) as excinfo:
            sat = Satellite(tle_ecc_err, self.effbg_observer, pydt_off)
            sat.eci_pos()

        assert 'e <= -0.001' in str(excinfo.value)

        tle_elsq_err = PyTle(*TLE_ELSQ_ERR)
        with pytest.raises(RuntimeError) as excinfo:
            sat = Satellite(tle_elsq_err, self.effbg_observer, pydt_off)
            sat.eci_pos()

        assert 'elsq >= 1.0' in str(excinfo.value)

        # tle_decayed_err = PyTle(*TLE_DECAYED_ERR)
        # with pytest.raises(RuntimeError) as excinfo:
        #     sat = Satellite(tle_decayed_err, self.effbg_observer, pydt_off)
        #     sat.eci_pos()

        # assert 'Satellite decayed' in str(excinfo.value)

    def test_error_coerce(self):

        pydt_off = PyDateTime.from_mjd(68805.5)

        tle_ecc_err = PyTle(*TLE_ECC_ERR)
        sat = Satellite(
            tle_ecc_err, self.effbg_observer, pydt_off,
            on_error='coerce_to_nan',
            )
        print(sat.eci_pos())
        assert_allclose(sat.eci_pos().loc, np.array([np.nan] * 3))

    def test_eci_position(self):

        sat = Satellite(self.tle, self.effbg_observer, self.pydt)
        # sat.mjd = self.mjd
        # print(sat.dt)
        eci_pos = sat.eci_pos()

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

        sat = Satellite(self.tle, self.effbg_observer, self.pydt)
        # sat.mjd = self.mjd
        # print(sat.dt)
        eci_pos = sat.eci_pos()
        geo_pos = sat.geo_pos()

        print(eci_pos.geo_loc)
        print(geo_pos.lon, geo_pos.lat, geo_pos.alt)

        assert_allclose(
            (geo_pos.lon, geo_pos.lat, geo_pos.alt),
            (-136.627536, 45.289346, 411.566712)
            )

    def test_topo_position(self):

        sat = Satellite(self.tle, self.effbg_observer, self.pydt)
        # sat.mjd = self.mjd
        # print(sat.dt)
        topo_pos = sat.topo_pos()

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
        location = EarthLocation.from_geodetic(
            *self.effbg_tup_m, 'WGS72'
            )
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

        sat = Satellite(self.tle, self.effbg_observer, self.pydt)
        # sat.mjd = self.mjd
        # print(sat.pydt)
        eci_pos = sat.eci_pos()
        topo_pos = sat.topo_pos()

        print(eci_pos.pydt)
        print(eci_pos.pydt.gmst(), eci_pos.pydt.lmst(self.effbg_tup[0]))
        print(topo_pos.az, topo_pos.el, topo_pos.dist)

        assert_allclose(topo_pos.az, az, atol=1e-3)
        assert_allclose(topo_pos.el, el, atol=1e-3)
        assert_allclose(topo_pos.dist, dist, atol=2e-2)

    def test_mjd_caching(self):

        # change cache resolution to 1 second
        sat = Satellite(self.tle, self.effbg_observer, self.pydt, 1. / 86400)
        topo_pos = sat.topo_pos()

        assert_allclose(
            (topo_pos.az, topo_pos.el, topo_pos.dist),
            (334.789646, -37.384929, 8406.367773)
            )

        # first half a second change must not trigger re-calculation
        sat.mjd += 0.5 / 86400.  # 0.5 s
        topo_pos = sat.topo_pos()

        print(topo_pos.az, topo_pos.el, topo_pos.dist)
        assert_allclose(
            (topo_pos.az, topo_pos.el, topo_pos.dist),
            (334.789646, -37.384929, 8406.367773)
            )

        # second half a second change must trigger re-calculation
        sat.mjd += 0.50001 / 86400.  # 0.5 s
        topo_pos = sat.topo_pos()
        print(topo_pos.az, topo_pos.el, topo_pos.dist)

        assert_allclose(
            (topo_pos.az, topo_pos.el, topo_pos.dist),
            (334.751037, -37.358432, 8402.015607)  # 4 km distance in 1 second
            )


def test_propagate_many():

    tles = PyTle(*TLE_ISS)
    observers = PyObserver(6.88375, 50.525, 0.366)
    mjds = np.linspace(56458.123, 56459.123, 4)
    result = propagate_many(mjds, tles, observers)
    eci_pos, eci_vel = result['eci_pos'], result['eci_vel']
    geo_pos = result['geo']
    topo_pos = result['topo']

    print(eci_pos.shape, topo_pos.shape)
    print(eci_pos)
    print(eci_vel)
    print(geo_pos)
    print(topo_pos)
    assert_allclose(
        eci_pos,
        np.array([
            [-4728.184444, 730.892025, 4802.515276],
            [-1147.648025, -5164.259228, 4245.417168],
            [3478.052600, -5776.901511, -858.022863],
            [4540.345522, -398.3767047, -5052.505474],
            ]),
        atol=1.e-5
        )

    assert_allclose(
        eci_vel,
        np.array([
            [1.530454, -7.065375, 2.574385],
            [5.324882, -4.172626, -3.615363],
            [3.719869, 3.104738, -5.930719],
            [-1.505247, 7.243366, -1.925367],
            ]),
        atol=1.e-5
        )

    assert_allclose(
        geo_pos,
        np.array([
            [-136.627536, 45.289346, 411.566712],
            [-170.697920, 38.923546, 413.350097],
            [112.553254, -7.296947, 419.680645],
            [46.159782, -48.125882, 438.181917],
            ]),
        atol=1.e-5
        )

    assert_allclose(
        topo_pos,
        np.array([
            [334.789646, -37.384929, 8406.367770, -4.351749],
            [358.119800, -43.467272, 9374.294560, 3.058513],
            [82.269620, -51.326038, 10485.152100, 4.306016],
            [153.995031, -50.639241, 10376.500700, 3.230162],
            ]),
        atol=1.e-5
        )


def test_propagate_many_broadcast():

    tles = np.array([PyTle(*TLE_ISS), PyTle(*TLE_GPS)])[:, np.newaxis]
    observers = PyObserver(6.88375, 50.525, 0.366)
    mjds = np.linspace(56458.123, 56459.123, 4)[np.newaxis, :]
    result = propagate_many(mjds, tles, observers)

    assert result['eci_pos'][..., 0].shape == (2, 4)


def test_propagate_many_pos_switches():

    tles = PyTle(*TLE_ISS)
    mjds = np.linspace(56458.123, 56459.123, 4)

    result = propagate_many(
        mjds, tles,
        do_eci_pos=True, do_eci_vel=True, do_geo=True, do_topo=True,
        )
    assert all(k in result for k in ['eci_pos', 'eci_vel', 'geo', 'topo'])

    result = propagate_many(
        mjds, tles,
        do_eci_pos=False, do_eci_vel=True, do_geo=True, do_topo=True,
        )
    assert ('eci_pos' not in result) and ('eci_vel' in result)


def test_propagate_many_raises_error():

    mjd_off = 68805.5
    obs = PyObserver(6.88375, 50.525, 0.366)

    tle_ecc_err = PyTle(*TLE_ECC_ERR)
    with pytest.raises(RuntimeError) as excinfo:
        propagate_many(mjd_off, tle_ecc_err, obs)

    assert 'e <= -0.001' in str(excinfo.value)

    tle_elsq_err = PyTle(*TLE_ELSQ_ERR)
    with pytest.raises(RuntimeError) as excinfo:
        propagate_many(mjd_off, tle_elsq_err, obs)

    assert 'elsq >= 1.0' in str(excinfo.value)

    tle_ecc_err = PyTle(*TLE_ECC_ERR)
    with pytest.raises(RuntimeError) as excinfo:
        propagate_many_slow(mjd_off, tle_ecc_err, obs)

    assert 'e <= -0.001' in str(excinfo.value)

    tle_elsq_err = PyTle(*TLE_ELSQ_ERR)
    with pytest.raises(RuntimeError) as excinfo:
        propagate_many_slow(mjd_off, tle_elsq_err, obs)

    assert 'elsq >= 1.0' in str(excinfo.value)


def test_propagate_many_suppress_error():

    mjd_off = 68805.5
    obs = PyObserver(6.88375, 50.525, 0.366)

    tle_ecc_err = PyTle(*TLE_ECC_ERR)
    res = propagate_many(mjd_off, tle_ecc_err, obs, on_error='coerce_to_nan')

    assert_allclose(res['eci_pos'], np.array([np.nan] * 3))

    res = propagate_many_slow(
        mjd_off, tle_ecc_err, obs, on_error='coerce_to_nan'
        )

    assert_allclose(res['eci_pos'], np.array([np.nan] * 3))


def propagate_many_sgp4(
        mjds, tles, observers=None,
        do_eci_pos=True, do_eci_vel=True,
        # bint do_geo=True, bint do_topo=True,
        ):
    '''
    This is an sgp4-based variant (non-parallelized, Python-looping) of
    `~cysgp4.propagate_many` that is meant for testing and benchmarking
    only. It has the same interface.
    '''

    pnum = 0
    out_dts = []
    if do_eci_pos:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    if do_eci_vel:
        out_dts.append(np.dtype(('float64', 3)))
        pnum += 1
    # if do_geo:
    #     out_dts.append(np.dtype(('float64', 3)))
    #     pnum += 1
    # if do_topo:
    #     out_dts.append(np.dtype(('float64', 4)))
    #     pnum += 1

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

            dt_tup = PyDateTime.from_mjd(mjd[i]).get_datetime_tuple()
            line1, line2 = tle[i].tle_strings()[1:]
            sat = twoline2rv(line1, line2, wgs72)
            pos, vel = sat.propagate(
                *dt_tup[:-2], dt_tup[-2] + dt_tup[-1] / 1e6
                )

            eci_pos_x, eci_pos_y, eci_pos_z = pos
            eci_vel_x, eci_vel_y, eci_vel_z = vel

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

            # if do_geo:
            #     itup[n][i][0] = geo_pos.lon
            #     itup[n][i][1] = geo_pos.lat
            #     itup[n][i][2] = geo_pos.alt
            #     n += 1

            # if do_topo:
            #     itup[n][i][0] = topo_pos.az
            #     itup[n][i][1] = topo_pos.el
            #     itup[n][i][2] = topo_pos.dist
            #     itup[n][i][3] = topo_pos.dist_rate
            #     n += 1

    result = {}

    n = 3
    if do_eci_pos:
        result['eci_pos'] = it.operands[n]
        n += 1

    if do_eci_vel:
        result['eci_vel'] = it.operands[n]
        n += 1

    # if do_geo:
    #     result['geo'] = it.operands[n]
    #     n += 1

    # if do_topo:
    #     result['topo'] = it.operands[n]
    #     n += 1

    return result


def _propagate_prepare():

    # url = 'http://celestrak.com/NORAD/elements/science.txt'
    # ctrak_science = requests.get(url)
    # all_lines = ctrak_science.text.split('\r\n')

    this_dir, this_filename = os.path.split(__file__)
    fname = os.path.join(this_dir, 'data', 'science.txt')
    with open(fname, 'r') as f:
        all_lines = [s.strip() for s in f.readlines()]

    tle_list = list(zip(*tuple(
        all_lines[idx::3] for idx in range(3)
        )))
    tles = np.array([
        PyTle(*tle) for tle in tle_list
        ])[np.newaxis, np.newaxis, :20]
    observers = np.array([
        PyObserver(6.88375, 50.525, 0.366),
        PyObserver(16.88375, 50.525, 0.366),
        ])[np.newaxis, :, np.newaxis]
    mjds = np.linspace(58805.5, 58806.5, 1000)[:, np.newaxis, np.newaxis]

    return mjds, tles, observers


def _propagate_many_cysgp4():

    return propagate_many(*_propagate_prepare())


def _propagate_many_cysgp4_slow():

    return propagate_many_slow(*_propagate_prepare())


def _propagate_many_sgp4():

    return propagate_many_sgp4(*_propagate_prepare())


def test_propagate_many_cysgp4_vs_many_cysgp4_slow():

    res_many = _propagate_many_cysgp4()
    res_many_slow = _propagate_many_cysgp4_slow()

    for k in ['eci_pos', 'eci_vel', 'geo', 'topo']:

        assert_allclose(res_many[k], res_many_slow[k], atol=1.e-5)


@pytest.mark.xfail
def test_propagate_many_cysgp4_vs_many_sgp4():
    '''
    Currently, very distant satellites can have a relatively large deviation
    of almost 1%. The problem with "rtol" is, that it doesn't work with
    positive vs. negative numbers (which happens when the satellite crosses
    zero in x, y or z)
    '''

    res_many = _propagate_many_cysgp4()
    res_many_sgp4 = _propagate_many_sgp4()

    # idxs = np.where(np.abs(res_many['eci_pos'] - res_many_sgp4['eci_pos']) > 1.e-2)

    # sat_idxs = np.unique(idxs[2])
    # for sat_idx in sat_idxs:
    #     mask = sat_idx == idxs[2]
    #     tidxs = tuple(idx[mask] for idx in idxs)
    #     print(sat_idx)
    #     print(res_many['eci_pos'][tidxs])
    #     print(res_many_sgp4['eci_pos'][tidxs])

    assert_allclose(res_many['eci_pos'], res_many_sgp4['eci_pos'], rtol=1.e-2)
    assert_allclose(res_many['eci_vel'], res_many_sgp4['eci_vel'], rtol=1.e-5)


@skip_benchmark
def test_propagate_many_cysgp4_benchmark(benchmark):

    benchmark(_propagate_many_cysgp4)


@skip_benchmark
def test_propagate_many_cysgp4_slow_benchmark(benchmark):

    benchmark(_propagate_many_cysgp4_slow)


@skip_benchmark
def test_propagate_many_sgp4_benchmark(benchmark):

    benchmark(_propagate_many_sgp4)
