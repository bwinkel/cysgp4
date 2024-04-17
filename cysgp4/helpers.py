#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pkg_resources import resource_stream
import numpy as np
from .cysgp4 import PyDateTime, Satellite, PyObserver
from .cysgp4 import eci_to_geo, geo_to_eci, lookangles, _propagate_many_cysgp4


__all__ = ['get_example_tles', 'propagate_many_slow', 'propagate_many']


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


# def propagate_many_sgp4(
#         mjds, tles, observers=None,
#         do_eci_pos=True, do_eci_vel=True,
#         do_geo=True, do_topo=True,
#         do_obs_pos=False, do_sat_azel=False,
#         sat_frame='zxy',  # on_error='raise',
#         ):
#     '''
#     This is an `sgp4`_-based variant (non-parallelized, Python-looping) of
#     `~cysgp4.propagate_many` that is meant for testing and benchmarking
#     only. It has the same interface.

#     Note: `sgp4 <https://pypi.org/project/sgp4/>`_ Package by Brandon Rhodes
#     '''

#     try:
#         import sgp4
#     except ImportError as e:
#         import warnings
#         warnings.warn('''
#             Package sgp4 needed to use this function. Can install with
#             "pip install sgp4" or "conda install -c conda-forge sgp4".
#             ''')
#         raise e

#     # check whether sgp4 v1 or v2 are installed

#     try:
#         from sgp4 import api
#         version = 2
#     except ImportError:
#         version = 1

#     if version == 1:
#         from sgp4.earth_gravity import wgs72
#         from sgp4.io import twoline2rv
#     elif version == 2:
#         from sgp4.api import Satrec, SGP4_ERRORS

#     if observers is None:
#         observers = PyObserver()

#     out_dts = [np.dtype(('float64', 3))]

#     it = np.nditer(
#         [tles, observers, mjds] + [None] * 2,
#         flags=['external_loop', 'buffered', 'delay_bufalloc', 'refs_ok'],
#         op_flags=[['readonly']] * 3 + [['readwrite', 'allocate']] * 2,
#         op_dtypes=['object', 'object', 'float64'] + out_dts + out_dts
#         )

#     it.reset()
#     for itup in it:

#         tle = itup[0]
#         mjd = itup[2]

#         size = mjd.shape[0]

#         # TODO code path for api.accelerated == True

#         for i in range(size):

#             line1, line2 = tle[i].tle_strings()[1:]

#             if version == 1:
#                 dt_tup = PyDateTime.from_mjd(mjd[i]).get_datetime_tuple()
#                 sat = twoline2rv(line1, line2, wgs72)
#                 pos, vel = sat.propagate(
#                     *dt_tup[:-2], dt_tup[-2] + dt_tup[-1] / 1e6
#                     )
#             elif version == 2:
#                 sat = Satrec.twoline2rv(line1, line2)
#                 _f, _i = np.modf(mjd[i])
#                 err_code, pos, vel = sat.sgp4(_i + 2400000.5, _f)

#                 if err_code:
#                     raise ValueError(
#                         'Satellite propagation error', err_code,
#                         '({})'.format(SGP4_ERRORS[err_code])
#                         )

#             eci_pos_x, eci_pos_y, eci_pos_z = pos
#             eci_vel_x, eci_vel_y, eci_vel_z = vel

#             itup[3][i][0] = eci_pos_x
#             itup[3][i][1] = eci_pos_y
#             itup[3][i][2] = eci_pos_z

#             itup[4][i][0] = eci_vel_x
#             itup[4][i][1] = eci_vel_y
#             itup[4][i][2] = eci_vel_z

#     eci_pos = it.operands[3]
#     eci_vel = it.operands[4]
#     eci_pos_x, eci_pos_y, eci_pos_z = (eci_pos[..., i] for i in range(3))
#     eci_vel_x, eci_vel_y, eci_vel_z = (eci_vel[..., i] for i in range(3))

#     # why is this not working???
#     # obs_geo_lon, obs_geo_lat, obs_geo_alt = np.vectorize(
#     #     lambda o: np.array([o.loc.lon, o.loc.lat, o.loc.alt]),
#     #     otypes=[np.float64, np.float64, np.float64]
#     #     )(observers)

#     obs_geo_lon = np.empty(observers.shape, dtype=np.float64)
#     obs_geo_lon.flat[:] = [o.loc.lon for o in observers.flat]
#     obs_geo_lat = np.empty(observers.shape, dtype=np.float64)
#     obs_geo_lat.flat[:] = [o.loc.lat for o in observers.flat]
#     obs_geo_alt = np.empty(observers.shape, dtype=np.float64)
#     obs_geo_alt.flat[:] = [o.loc.alt for o in observers.flat]

#     (
#         tles,
#         obs_geo_lon, obs_geo_lat, obs_geo_alt,
#         mjds,
#         eci_pos_x, eci_pos_y, eci_pos_z,
#         eci_vel_x, eci_vel_y, eci_vel_z,
#         ) = np.broadcast_arrays(
#             tles,
#             obs_geo_lon, obs_geo_lat, obs_geo_alt,
#             mjds,
#             eci_pos_x, eci_pos_y, eci_pos_z,
#             eci_vel_x, eci_vel_y, eci_vel_z,
#             )

#     result = {}

#     if do_eci_pos:
#         result['eci_pos'] = np.stack(
#             [eci_pos_x, eci_pos_y, eci_pos_z], axis=-1
#             )

#     if do_eci_vel:
#         result['eci_vel'] = np.stack(
#             [eci_vel_x, eci_vel_y, eci_vel_z], axis=-1
#             )

#     if do_geo:
#         geo_lon, geo_lat, geo_alt = eci_to_geo(
#             eci_pos_x, eci_pos_y, eci_pos_z, mjds
#             )
#         result['geo'] = np.stack([geo_lon, geo_lat, geo_alt], axis=-1)

#     if do_topo or do_sat_azel:
#         obs_az, obs_el, sat_az, sat_el, dist, dist_rate = lookangles(
#             eci_pos_x, eci_pos_y, eci_pos_z,
#             eci_vel_x, eci_vel_y, eci_vel_z,
#             mjds, observers, sat_frame
#             )
#         topo = np.stack([obs_az, obs_el, dist, dist_rate], axis=-1)
#         sat_azel = np.stack([sat_az, sat_el, dist], axis=-1)

#     if do_topo:
#         result['topo'] = topo

#     if do_sat_azel:
#         result['sat_azel'] = sat_azel

#     if do_obs_pos:
#         obs_x, obs_y, obs_z = geo_to_eci(
#             obs_geo_lon, obs_geo_lat, obs_geo_alt, mjds
#             )
#         result['obs_pos'] = np.stack([obs_x, obs_y, obs_z], axis=-1)

#     return result


def _sat_eci_sgp4_nonaccelerated(tles, observers, mjds):

    from sgp4.api import Satrec, SGP4_ERRORS

    out_dts = [np.dtype(('float64', 3))]

    it = np.nditer(
        [tles, observers, mjds] + [None] * 2,
        flags=['external_loop', 'buffered', 'delay_bufalloc', 'refs_ok'],
        op_flags=[['readonly']] * 3 + [['readwrite', 'allocate']] * 2,
        op_dtypes=['object', 'object', 'float64'] + out_dts + out_dts
        )

    it.reset()
    for itup in it:

        tle = itup[0]
        mjd = itup[2]

        size = mjd.shape[0]

        # TODO code path for api.accelerated == True

        for i in range(size):

            line1, line2 = tle[i].tle_strings()[1:]
            sat = Satrec.twoline2rv(line1, line2)
            _f, _i = np.modf(mjd[i])
            err_code, pos, vel = sat.sgp4(_i + 2400000.5, _f)

            if err_code:
                raise ValueError(
                    'Satellite propagation error', err_code,
                    '({})'.format(SGP4_ERRORS[err_code])
                    )

            eci_pos_x, eci_pos_y, eci_pos_z = pos
            eci_vel_x, eci_vel_y, eci_vel_z = vel

            itup[3][i][0] = eci_pos_x
            itup[3][i][1] = eci_pos_y
            itup[3][i][2] = eci_pos_z

            itup[4][i][0] = eci_vel_x
            itup[4][i][1] = eci_vel_y
            itup[4][i][2] = eci_vel_z

    eci_pos = it.operands[3]
    eci_vel = it.operands[4]

    return eci_pos, eci_vel


def _sat_eci_sgp4_accelerated(tles, observers, mjds):
    '''
    Note, on a single core, this version is not really faster than
    the non-accelerated version. Probably the overhead of the python
    iteration can mostly be neglected.
    '''

    from sgp4.api import Satrec, SGP4_ERRORS

    unique_tles = np.unique(tles)
    tles, observers, mjds = np.broadcast_arrays(tles, observers, mjds)

    eci_pos = np.empty(tles.shape + (3,), dtype=np.float64)
    eci_vel = np.empty(tles.shape + (3,), dtype=np.float64)

    for u_tle in unique_tles:

        # doesn't work, as array doesn't know tle_strings method
        # mask = u_tle == tles
        mask = np.array(
            [u_tle == t for t in tles.flat], dtype=bool
            ).reshape(tles.shape)

        mjd = mjds[mask]

        line1, line2 = u_tle.tle_strings()[1:]

        sat = Satrec.twoline2rv(line1, line2)
        _f, _i = np.modf(mjd)
        err_code, pos, vel = sat.sgp4_array(_i + 2400000.5, _f)

        if np.any(err_code):
            raise ValueError(
                'Satellite propagation error', err_code,
                '({})'.format(SGP4_ERRORS[err_code])
                )

        eci_pos[mask] = pos
        eci_vel[mask] = vel

    return eci_pos, eci_vel


def _propagate_many_sgp4(
        mjds, tles, observers=None,
        do_eci_pos=True, do_eci_vel=True,
        do_geo=True, do_topo=True,
        do_obs_pos=False, do_sat_azel=False,
        do_sat_rotmat=False,
        sat_frame='zxy',  # on_error='raise',
        ):
    '''
    This is an `sgp4`_-based variant (Python-looping) of
    `~cysgp4.propagate_many` that is meant for testing and benchmarking
    only. It has the same interface.

    Note: `sgp4 <https://pypi.org/project/sgp4/>`_ Package by Brandon Rhodes
    '''

    try:
        from sgp4 import api
    except ImportError as e:
        import warnings
        warnings.warn('''
            Package sgp4 2.x needed to use this function. Can install with
            "pip install sgp4" or "conda install -c conda-forge sgp4".
            ''')
        raise e

    if observers is None:
        observers = PyObserver()

    if api.accelerated:
        eci_pos, eci_vel = _sat_eci_sgp4_accelerated(tles, observers, mjds)
    else:
        eci_pos, eci_vel = _sat_eci_sgp4_nonaccelerated(tles, observers, mjds)

    eci_pos_x, eci_pos_y, eci_pos_z = (eci_pos[..., i] for i in range(3))
    eci_vel_x, eci_vel_y, eci_vel_z = (eci_vel[..., i] for i in range(3))

    # why is this not working???
    # obs_geo_lon, obs_geo_lat, obs_geo_alt = np.vectorize(
    #     lambda o: np.array([o.loc.lon, o.loc.lat, o.loc.alt]),
    #     otypes=[np.float64, np.float64, np.float64]
    #     )(observers)

    obs_geo_lon = np.empty(observers.shape, dtype=np.float64)
    obs_geo_lon.flat[:] = [o.loc.lon for o in observers.flat]
    obs_geo_lat = np.empty(observers.shape, dtype=np.float64)
    obs_geo_lat.flat[:] = [o.loc.lat for o in observers.flat]
    obs_geo_alt = np.empty(observers.shape, dtype=np.float64)
    obs_geo_alt.flat[:] = [o.loc.alt for o in observers.flat]

    (
        tles,
        obs_geo_lon, obs_geo_lat, obs_geo_alt,
        mjds,
        eci_pos_x, eci_pos_y, eci_pos_z,
        eci_vel_x, eci_vel_y, eci_vel_z,
        ) = np.broadcast_arrays(
            tles,
            obs_geo_lon, obs_geo_lat, obs_geo_alt,
            mjds,
            eci_pos_x, eci_pos_y, eci_pos_z,
            eci_vel_x, eci_vel_y, eci_vel_z,
            )

    result = {}

    if do_eci_pos:
        result['eci_pos'] = np.stack(
            [eci_pos_x, eci_pos_y, eci_pos_z], axis=-1
            )

    if do_eci_vel:
        result['eci_vel'] = np.stack(
            [eci_vel_x, eci_vel_y, eci_vel_z], axis=-1
            )

    if do_geo:
        geo_lon, geo_lat, geo_alt = eci_to_geo(
            eci_pos_x, eci_pos_y, eci_pos_z, mjds
            )
        result['geo'] = np.stack([geo_lon, geo_lat, geo_alt], axis=-1)

    if do_topo or do_sat_azel or do_sat_rotmat:
        res = lookangles(
            eci_pos_x, eci_pos_y, eci_pos_z,
            eci_vel_x, eci_vel_y, eci_vel_z,
            mjds, observers, sat_frame,
            do_sat_rotmat=do_sat_rotmat,
            )

        if do_sat_rotmat:
            obs_az, obs_el, sat_az, sat_el, dist, dist_rate, *sat_bas = res
        else:
            obs_az, obs_el, sat_az, sat_el, dist, dist_rate = res

    if do_topo:
        result['topo'] = np.stack([obs_az, obs_el, dist, dist_rate], axis=-1)

    if do_sat_azel:
        result['sat_azel'] = np.stack([sat_az, sat_el, dist], axis=-1)

    if do_sat_rotmat:
        new_sh = sat_bas[0].shape + (3, 3)
        result['sat_rotmat'] = np.swapaxes(
            np.stack(sat_bas, axis=-1).reshape(new_sh),
            -2, -1,
            )

    if do_obs_pos:
        obs_x, obs_y, obs_z = geo_to_eci(
            obs_geo_lon, obs_geo_lat, obs_geo_alt, mjds
            )
        result['obs_pos'] = np.stack([obs_x, obs_y, obs_z], axis=-1)

    return result


def propagate_many(
        mjds, tles, observers=None,
        do_eci_pos=True, do_eci_vel=True,
        do_geo=True, do_topo=True,
        do_obs_pos=False, do_sat_azel=False,
        do_sat_rotmat=False,
        sat_frame='zxy', on_error='raise',
        method='dwarner'
        ):
    '''
    Calculate positions of many satellites at a various times at once.

    This is an array interface to the sgp4 calculations, which allows to
    perform calculations for different satellite TLEs, observers and times
    in a parallelized manner. `~numpy` broadcasting rules apply.

    With the Boolean parameters, `do_eci_pos`, `do_eci_vel`, `do_geo`,
    `do_topo`, `do_obs_pos`, and `do_sat_azel` the user can decide, which
    position frames are returned in the output dictionary.

    Satellite are defined via TLEs (see `~cysgp4.PyTle`). The
    `~.cysgp4.propagate_many` function works with a single (scalar) PyTle
    object or a list/array of them. The same is true for the `mjds` and
    optional `observers` parameters, which must be a scalar (a double and
    `~cysgp4.PyObserver` instance, respectively) or a list/array of them.

    As in most use cases, a large number of positions is probably queried,
    the returned values do not use the `~cysgp4.PyEci` or
    `~cysgp4.PyCoordTopocentric` classes, but pack everything in arrays.
    This is not only faster, but makes it easier to further process the
    results.

    For parallelization, `OpenMP <https://www.openmp.org/>`_ is utilized.
    To change the number of CPU cores that are used, one can either set
    the environment variable `OMP_NUM_THREADS` or use
    `~cysgp4.set_num_threads`.

    Parameters
    ----------
    mjds : `~numpy.ndarray`, `~list`, or scalar of float
        Modified Julian Date.
    tles : `~numpy.ndarray`, `~list`, or scalar of `~cysgp4.PyTle`
        TLE instance of the satellite of interest.
    observers : `~numpy.ndarray`, `~list`, or scalar of `~cysgp4.PyObserver` or None (default: None)
        Observer instance. If `None` then the observer location is set to
        (0 deg, 0 deg, 0 km).
    do_eci_pos : Boolean, optional (default: True)
        Whether to include ECI cartesian positions in the result.
    do_eci_vel : Boolean, optional (default: True)
        Whether to include ECI cartesian velocity in the result.
    do_geo : Boolean, optional (default: True)
        Whether to include geographic/geodetic positions in the result.
    do_topo : Boolean, optional (default: True)
        Whether to include topocentric positions in the result.
    do_obs_pos : Boolean, optional (default: True)
        Whether to include the observer ECI position in the results.
    do_sat_azel : Boolean, optional (default: True)
        Whether to include the observer position as seen by the satellite
        (distance/azimuth/elevation) in the results.
    do_sat_rotmat : Boolean, optional (default: False)
        Whether to include the rotation matrix that converts the
        (moving and rotated) satellite frame (in cartesian) into
        cartesian ECI-aligned coordinates in the results. This can be useful
        for cases where the user needs to transform additional vectors
        between both frames (and is not only interested the observer
        position in the satellite frame as returned by `do_sat_azel`).
    sat_frame : 'zxy' or 'xyz', optional (default: 'zxy')
        How the moving satellite frame is defined. Two options are
        implemented, 'zxy' and 'xyz'. If 'zxy' is chosen, the moving
        satellite frame is constructed such that the `z` axis is
        aligned with the satellite motion vector. The `y` axis is lies
        perpendicularly to the plane defined by the motion vector and
        the ECI zero point (aka the Earth centre). The resulting `x`
        axis, which is orthogonal to the `y` and `z` axes, is then
        approximately pointing towards nadir. Alternatively, if the
        frame is set as `xyz`, the `x` axis is the motion vector, `y`
        has the same meaning (but points into the opposite direction)
        and `z` is approximately pointing towards the nadir. The
        definition of the output polar angles is different for the two
        reference frames, see Returns.
    on_error : str, optional (either 'raise' or 'coerce_to_nan', default: 'raise')
        If the underlying SGP C++ library throws an error (which often
        happens if one works with times that are strongly deviating from
        the TLE epoch), a Python RuntimeError is usually thrown. For batch
        processing, this is not always desirable. If
        `on_error = 'coerce_to_nan'` then C++ errors will be suppressed and
        the resulting position vectors will be converted to NaN-values
        instead.

        Note: this option is only relevant for `method='dwarner'`!
    method : str, optional (either 'dwarner' or 'vallado')
        This option can be used to change from the default SGP C++ library
        (by Daniel Warner) to the Vallado model. For the latter, the
        external `python-sgp4` library is used. 'vallado' may be slower
        due to the additional boiler-plate code and differences in the
        model itself and implementation details.

    Returns
    -------
    result : dict
        Resulting positions for each requested frame:

        - `eci_pos` : `~numpy.ndarray` of float

          Satellites ECI cartesian positions. Last dimension has length 3,
          one for each of `x`, `y`, and `z`.

        - `eci_vel` : `~numpy.ndarray` of float

          Satellites ECI cartesian velicities. Last dimension has length 3,
          one for each of `v_x`, `v_y`, and `v_z`.

        - `geo` : `~numpy.ndarray` of float

          Satellites Geodetic positions. Last dimension has length 3, one
          for each of `lon`, `lat`, and `alt`.

        - `topo` : `~numpy.ndarray` of float

          Satellites Topocentric positions. Last dimension has length 4,
          one for each of `az`, `el`, `dist`, and `dist_rate`.

        - `obs_pos` : `~numpy.ndarray` of float

          Observer positions in ECI frame (Cartesian). Last dimension has
          length 3, one for each of `x`, `y`, and `z`.

        - `sat_azel` : `~numpy.ndarray` of float

          If `sat_frame` is 'zxy', `z` lies in the direction of motion,
          `y` perpendicular to the z-axis and the Earth center, `x` is
          pointing approximately towards nadir, see also `sat_frame`
          parameter description. The Observer positions in the
          (co-moving) satellite frame are given as azimuth, elevation
          in the specified reference frame, and distance (`az`, `el`,
          `dist`). `az` is the angle between the projection of the vector
          towards the Observer onto the xy-plane and the x-axis. -180
          deg < `az` < 180 deg. `el` is the angle between the normal
          vector and the xy-plane. -90 deg < `el` < 90 deg.

          If `sat_frame` is 'xyz', `x` lies in the direction of motion,
          `y` is perpendicular to `z` and the Earth center, `z` is pointing
          approximately towards nadir, see also `sat_frame` parameter
          description. The Observer positions in the (moving)
          satellite frame are given as azimuth and polar angle in the
          specified reference frame, and distance (`az`, `theta`, `dist`). `az`
          is the angle between the projection of the vector towards
          the observer onto the xy-plane and the x-axis. -180 deg < `az`
          < 180 deg. `theta` is the angle between the normal vector and
          the z-axis. -90 deg < `theta` < 90 deg.
          
        - `sat_rotmat` : `~numpy.ndarray` of float

          Rotation matrices which would transform a vector defined in the
          (moving and rotated) satellite frames (in cartesian) to the
          cartesian ECI-aligned basis frame. It is noted that the origin of
          this ECI-aligned frame is still at the satellite center.

          This can be useful for cases where the user needs to transform
          additional vectors between both frames (and is not only interested
          the observer position in the satellite frame as returned by
          `do_sat_azel`).

          Likewise, the inverse of these rotation matrices (aka the
          transposed) can be used to rotate any vector from ECI-aligned
          satellite basis frame to the satellite frame.

        In all cases the first dimensions are determined by the
        (broadcasted) shape of the inputs `mjd`, `tles`, and `observers`.

    Examples
    --------
    The following demonstrates how to use the `~cysgp4.propagate_many`
    function::

        >>> import numpy as np
        >>> from cysgp4 import PyTle, PyObserver, propagate_many
        >>> from cysgp4 import get_example_tles, tles_from_text

        >>> tle_text = get_example_tles()
        >>> tles = np.array(
        ...     tles_from_text(tle_text)
        ...     )[np.newaxis, np.newaxis, :20]  # use first 20 TLEs
        >>> observers = np.array([
        ...     PyObserver(6.88375, 50.525, 0.366),
        ...     PyObserver(16.88375, 50.525, 0.366),
        ...     ])[np.newaxis, :, np.newaxis]
        >>> mjds = np.linspace(
        ...     58805.5, 58806.5, 1000  # 1000 time steps
        ...     )[:, np.newaxis, np.newaxis]

        >>> result = propagate_many(mjds, tles, observers)
        >>> print(sorted(result.keys()))
        ['eci_pos', 'eci_vel', 'geo', 'topo']

        >>> # shapes are as follows
        >>> print(np.broadcast(mjds, tles, observers).shape)
        (1000, 2, 20)
        >>> print(result['eci_pos'].shape, result['topo'].shape)
        (1000, 2, 20, 3) (1000, 2, 20, 4)

        >>> result = propagate_many(
        ...     mjds, tles, observers,
        ...     do_eci_pos=False, do_eci_vel=False, do_geo=False, do_topo=True
        ...     )
        >>> print(sorted(result.keys()))
        ['topo']

    '''

    if method not in ['dwarner', 'vallado']:
        raise ValueError('"method" most be one of ["dwarner", "vallado"]')

    if on_error not in ['raise', 'coerce_to_nan']:
        raise ValueError(
            '"on_error" most be one of ["raise", "coerce_to_nan"]'
            )

    if sat_frame not in ['zxy', 'xyz']:
        raise ValueError('"sat_frame" most be one of ["zxy", "xyz"]')

    if method == 'dwarner':

        return _propagate_many_cysgp4(
            mjds, tles, observers=observers,
            do_eci_pos=do_eci_pos, do_eci_vel=do_eci_vel,
            do_geo=do_geo, do_topo=do_topo,
            do_obs_pos=do_obs_pos, do_sat_azel=do_sat_azel,
            do_sat_rotmat=do_sat_rotmat,
            sat_frame=sat_frame, on_error=on_error,
            )

    elif method == 'vallado':

        return _propagate_many_sgp4(
            mjds, tles, observers=observers,
            do_eci_pos=do_eci_pos, do_eci_vel=do_eci_vel,
            do_geo=do_geo, do_topo=do_topo,
            do_obs_pos=do_obs_pos, do_sat_azel=do_sat_azel,
            do_sat_rotmat=do_sat_rotmat,
            sat_frame=sat_frame,
            )
