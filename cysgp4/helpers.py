#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pkg_resources import resource_stream
import numpy as np
from .cysgp4 import PyDateTime, Satellite, PyObserver


__all__ = ['get_example_tles', 'propagate_many_slow', 'propagate_many_sgp4']


def get_example_tles():
    '''
    Return a text with example TLE entries for testing purposes.

    The example is the content of the file "science.txt", which was downloaded
    from `Celestrak <https://celestrak.com/>`_ on Nov 25 2019.

    Returns
    -------
    tle_text : str
        Text containing TLE line strings.


    Examples
    --------
    Use like this::

        >>> import cysgp4

        >>> tle_text = cysgp4.get_example_tles()
        >>> tle_text[:200]
        'AKEBONO (EXOS-D)        \\r\\n1 19822U 89016A   19321.49921565  .00016421  94291-6  28704-3 0  9992\\r\\n2 19822  75.0304 327.7460 1766579 276.4058  63.9734 12.00957719 13956\\r\\nHST                     \\r\\n1 2058'

        >>> tles = cysgp4.tle_tuples_from_text(tle_text)
        >>> tles  # doctest: +SKIP
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

    with resource_stream('cysgp4', 'tests/data/science.txt') as f:
        text = f.read()

    return text.decode('utf-8')


def propagate_many_slow(
        mjds, tles, observers=None,
        do_eci_pos=True, do_eci_vel=True,
        do_geo=True, do_topo=True,
        on_error='raise',
        ):
    '''
    This is a slow (non-parallelized, Python-looping) version of
    `~cysgp4.propagate_many` that is meant for testing and benchmarking
    only. It has the same interface.
    '''

    assert on_error in ['raise', 'coerce_to_nan']

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

            sat = Satellite(
                tle[i], obs[i], PyDateTime.from_mjd(mjd[i]),
                on_error=on_error
                )
            eci = sat.eci_pos()

            eci_pos_x, eci_pos_y, eci_pos_z = eci.loc
            eci_vel_x, eci_vel_y, eci_vel_z = eci.vel

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

            if do_geo:
                geo_pos = sat.geo_pos()
                itup[n][i][0] = geo_pos.lon
                itup[n][i][1] = geo_pos.lat
                itup[n][i][2] = geo_pos.alt
                n += 1

            if do_topo:
                topo_pos = sat.topo_pos()
                itup[n][i][0] = topo_pos.az
                itup[n][i][1] = topo_pos.el
                itup[n][i][2] = topo_pos.dist
                itup[n][i][3] = topo_pos.dist_rate
                n += 1

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

    return result


def propagate_many_sgp4(
        mjds, tles, observers=None,
        do_eci_pos=True, do_eci_vel=True,
        ):
    '''
    This is an `sgp4`_-based variant (non-parallelized, Python-looping) of
    `~cysgp4.propagate_many` that is meant for testing and benchmarking
    only. It has (almost) the same interface, except that `do_geo`, `do_topo`,
    and `on_error` are not handled.

    Note: `sgp4 <https://pypi.org/project/sgp4/>`_ Package by Brandon Rhodes
    '''

    try:
        from sgp4.earth_gravity import wgs72
        from sgp4.io import twoline2rv
    except ImportError as e:
        import warnings
        warnings.warn('''
            Package sgp4 needed to use this function. Can install with
            "pip install sgp4" or "conda install -c conda-forge sgp4".
            ''')
        raise e

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
