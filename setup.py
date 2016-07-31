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

CPPSOURCES = glob.glob('cextern/sgp4-05d8cc2fc596/libsgp4/*.cpp')
# print(cppsources)

SGP_EXT = Extension(
    'cysgp4.cysgp4',
    ['cysgp4/cysgp4.pyx'] + CPPSOURCES,
    extra_compile_args=['-O3'] + EX_COMP_ARGS,
    # extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
        'cysgp4/',
        'cextern/sgp4-05d8cc2fc596/libsgp4/'
        # 'cextern/deprecated/'
    ]
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
        'cython>=0.20.2',
        'numpy>=1.8',
        # 'astropy>=1.0',
    ],
    packages=['cysgp4'],
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        SGP_EXT,
    ],
    # url='https://github.com/bwinkel/cyaatm/',
    # download_url='https://github.com/bwinkel/cyaatm/tarball/0.1.0',
    # keywords=['astronomy']
)
