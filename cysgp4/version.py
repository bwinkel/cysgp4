#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution(__name__.split('.')[0]).version
except DistributionNotFound:
    # package is not installed
    __version__ = None

__all__ = ['__version__']
