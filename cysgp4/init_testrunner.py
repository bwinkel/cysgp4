#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
__all__ = []

try:

    # Create the test function for self test
    from astropy.tests.runner import TestRunner
    test = TestRunner.make_test_runner_in(os.path.dirname(__file__))
    test.__test__ = False
    __all__ += ['test']

except ImportError:
    pass
