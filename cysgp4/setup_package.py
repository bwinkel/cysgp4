#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Allow cythonizing of our pyx files and provide custom compiler options.
'''

import os
import glob
from setuptools.extension import Extension
from extension_helpers import get_compiler, add_openmp_flags_if_available
import platform
import numpy as np
# Note: importing numpy from here won't work, see:
# http://docs.astropy.org/en/stable/development/ccython.html#using-numpy-c-headers
# import numpy as np
# 'include_dirs': [np.get_include()], --> 'include_dirs': 'numpy'

PYXDIR = os.path.relpath(os.path.dirname(__file__))
PYXFILES = [
    'cysgp4.pyx', 'utils.pyx'
    ]
CPPSOURCES = glob.glob('cextern/sgp4-f5cb54b3/libsgp4/*.cc')


def get_extensions():

    print('Using compiler', get_compiler())

    comp_args = {
        'extra_compile_args': ['-O3', '-std=c++11'],
        'language': 'c++',
        'include_dirs': [
            np.get_include(),
            'cysgp4/',
            'cextern/sgp4-f5cb54b3/libsgp4',
            ],
        }

    if platform.system().lower() == 'windows':

        comp_args = {
            'include_dirs': [
                np.get_include(),
                'cysgp4/',
                'cextern/sgp4-f5cb54b3/libsgp4',
                ],
            }

    elif 'darwin' in platform.system().lower():

        # from subprocess import getoutput

        mac_version = '10.7'

        extra_compile_args = [
            '-fopenmp', '-O3', '-std=c++11',
            '-mmacosx-version-min={:s}'.format(mac_version),
            ]

        # if 'clang' in getoutput('gcc -v'):
        if 'clang' in get_compiler():
            extra_compile_args += ['-stdlib=libc++', ]

        comp_args['extra_compile_args'] = extra_compile_args

    ext_list = []
    for pyx in PYXFILES:
        ext = Extension(
            name='cysgp4.{}'.format(pyx.replace('.pyx', '')),
            sources=[os.path.join(PYXDIR, pyx)] + CPPSOURCES,
            **comp_args
            )
        add_openmp_flags_if_available(ext)
        ext_list.append(ext)

    return ext_list

