#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pkg_resources import resource_stream


__all__ = ['get_example_tles']


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
