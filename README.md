# Introduction #

- *Version*: 0.1
- *Authors*: Benjamin Winkel

# Purpose#

`cysgp4` is a thin wrapper around the sgp4lib (C++) library, which calculates
positions of satellites from two-line elements (TLE)

# Features

* Exposes a (minimal) set of sgp4lib classes:
    PyDateTime, PyTle, PyObserver,
    PyCoordGeodetic, PyCoordTopocentric, PyEci
* Adds some helper classes/functions to conveniently compute satellite position

# Usage #

### Installation ###

TODO

### Dependencies ###

Dependencies are kept as minimal as possible. The following packages are
required:
* `numpy 1.10` or later
* `cython 0.23.4` or later
Furthermore, a C++ compiler is needed.

### Examples ###

TODO

### Who do I talk to? ###

If you encounter any problems or have questions, do not hesitate to raise an
issue or make a pull request. Moreover, you can contact the dev directly:

* <bwinkel@mpifr.de>

### Contributions ###

Contributions are more than welcome.
