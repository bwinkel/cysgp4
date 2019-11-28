#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest


def pytest_collection_modifyitems(items, config):
    '''
    Skip tests marked with 'benchmark' if pytest-benchmark is not installed.
    '''

    if config.pluginmanager.hasplugin('benchmark'):
        return

    for item in items:
        if item.get_closest_marker('benchmark'):
            item.add_marker(
                pytest.mark.skip(reason='pytest-benchmark is not installed')
            )
