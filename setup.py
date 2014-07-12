from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


ext_module_cysgp4 = Extension(
    'cysgp4',
    [
        'src/libsgp4/sgp4.pyx',
        'src/libsgp4/CoordGeodetic.cpp',
        'src/libsgp4/CoordTopocentric.cpp',
        'src/libsgp4/DateTime.cpp',
        'src/libsgp4/Eci.cpp',
        'src/libsgp4/Globals.cpp',
        'src/libsgp4/Observer.cpp',
        'src/libsgp4/OrbitalElements.cpp',
        'src/libsgp4/SGP4.cpp',
        'src/libsgp4/SolarPosition.cpp',
        'src/libsgp4/TimeSpan.cpp',
        'src/libsgp4/Tle.cpp',
        'src/libsgp4/Util.cpp',
        'src/libsgp4/Vector.cpp',
        ],
    language='c++',
    extra_compile_args=['-O2']
    #extra_compile_args=['-fopenmp','-O3'],
    #extra_link_args=['-fopenmp'],
    #libraries=['m'],
)


setup(
    name = 'cysgp4',
    version = '0.1',
    description = 'cysgp4',
    author = 'Benjamin Winkel',
    author_email = 'bwinkel@mpifr.de',
    url = 'http://www.astro.uni-bonn.de/~bwinkel',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext_module_cysgp4],
    #package_data = {'ralib' : ['rttools/data/receivers.csv', 'rttools/data/horizons.txt'] },
    long_description = """cysgp4 ... Cython-powered wrapper of the
    sgp4lib (Daniel Warner) library to compute satellite positions
    from two-line elements (TLE)."""
) 
