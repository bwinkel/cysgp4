#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import os
import numpy
import platform
import glob
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser

# Produce annotated html files
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True


# Get some values from the setup.cfg
conf = ConfigParser()
conf.read(['setup.cfg'])
metadata = dict(conf.items('metadata'))

PACKAGENAME = metadata.get('package_name', 'cysgp4')
DESCRIPTION = metadata.get('description', 'Cysgp4: a wrapper around the SGP4 package, for sat TLE calculations')
LONG_DESCRIPTION = metadata.get('long_description', '')
AUTHOR = metadata.get('author', 'Benjamin Winkel')
AUTHOR_EMAIL = metadata.get('author_email', '')
LICENSE = metadata.get('license', 'unknown')
URL = metadata.get('url', 'https://github.com/bwinkel/cysgp4')
__minimum_python_version__ = metadata.get("minimum_python_version", "3.5")

if os.path.exists('README.rst'):
    with open('README.rst') as f:
        LONG_DESCRIPTION = f.read()

EX_COMP_ARGS = []
if 'mac' in platform.system().lower():
    EX_COMP_ARGS += ['-stdlib=libc++', '-mmacosx-version-min=10.7']


def get_compile_args():

    comp_args = {
        'extra_compile_args': ['-fopenmp', '-O3', '-std=c++11'],
        'extra_link_args': ['-fopenmp'],
        'language': 'c++',
        'include_dirs': [
            numpy.get_include(),
            'cysgp4/',
            'cextern/sgp4-05d8cc2fc596/libsgp4/',
            ],
        }

    if platform.system().lower() == 'windows':
        comp_args['extra_compile_args'] = ['/openmp']
        del comp_args['extra_link_args']

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
            '-fopenmp', '-O3', '-std=c++11',
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

UTILS_EXT = Extension(
    'cysgp4.utils',
    ['cysgp4/utils.pyx'] + CPPSOURCES,
    **get_compile_args()
    )

for e in [SGP_EXT, UTILS_EXT]:
    e.cython_directives = {'language_level': "3"}  # all are Python-3


print('Cython.Compiler.Options.annotate', Cython.Compiler.Options.annotate)

# NOTE: for github pages, put an empty .nojekyll into the root dir of
# the web directory (gh-pages branch root)

setup(
    name=PACKAGENAME,
    # version='0.1.0',
    use_scm_version=True,  # provides version
    # use  the following to read version in package files
    # from pkg_resources import get_distribution, DistributionNotFound
    # try:
    #     __version__ = get_distribution(__name__).version
    # except DistributionNotFound:
    #     # package is not installed
    #     pass
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    license=LICENSE,
    install_requires=[
        'setuptools',
        'numpy>=1.13.1',
        ],
    tests_require=['pytest', 'numpy>=1.13.1', 'sgp4'],
    packages=find_packages(),
    package_data={
        PACKAGENAME: ['tests/data/science.txt']
        },
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        UTILS_EXT, SGP_EXT,
        ],
    zip_safe=False,
    # url='https://github.com/bwinkel/cyaatm/',
    # download_url='https://github.com/bwinkel/cyaatm/tarball/0.1.0',
    # keywords=['astronomy']
    python_requires='>={}'.format(__minimum_python_version__),
    )
