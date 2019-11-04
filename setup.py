#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy
import platform
import glob

EX_COMP_ARGS = []
if 'mac' in platform.system().lower():
    EX_COMP_ARGS += ['-stdlib=libc++', '-mmacosx-version-min=10.7']


def get_compile_args():

    comp_args = {
        'extra_compile_args': ['-O3', '-std=c++11'],
        'language': 'c++',
        'include_dirs': [
            numpy.get_include(),
            'cysgp4/',
            'cextern/sgp4-05d8cc2fc596/libsgp4/',
            ],
        }

    if 'darwin' in platform.system().lower():

        import sys
        py3 = sys.version_info[0] >= 3
        mac_version = '10.7' if py3 else '10.6'

        if py3:
            from subprocess import getoutput
        else:
            from subprocess import check_output

            def getoutput(s):
                return check_output(s.split())

        extra_compile_args = [
            '-O3', '-std=c++11',
            '-mmacosx-version-min={:s}'.format(mac_version),
            ]
        if ('clang' in getoutput('gcc -v')) and all(
                'command not found' in getoutput('gcc-{:d} -v'.format(d))
                for d in [6, 7, 8]
                ):
            extra_compile_args += ['-stdlib=libc++', ]
        comp_args['extra_compile_args'] = extra_compile_args

    return comp_args


CPPSOURCES = glob.glob('cextern/sgp4-05d8cc2fc596/libsgp4/*.cpp')
# print(cppsources)

SGP_EXT = Extension(
    'cysgp4.cysgp4',
    ['cysgp4/cysgp4.pyx'] + CPPSOURCES,
    **get_compile_args()
    )


setup(
    name='cysgp4',
    version='0.1.0',
    author='Benjamin Winkel',
    author_email='bwinkel@mpifr.de',
    description=(
        'cysgp4: a wrapper around the SGP4 package, for sat TLE calculations'
        ),
    long_description='''cysgp4 ... Cython-powered wrapper of the
    sgp4lib (Daniel Warner) library to compute satellite positions
    from two-line elements (TLE).''',
    install_requires=[
        'setuptools',
        'cython>=0.27.1',
        'numpy>=1.11',
        # 'astropy>=1.0',
        ],
    packages=['cysgp4'],
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        SGP_EXT,
        ],
    zip_safe=False,
    # url='https://github.com/bwinkel/cyaatm/',
    # download_url='https://github.com/bwinkel/cyaatm/tarball/0.1.0',
    # keywords=['astronomy']
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    )
