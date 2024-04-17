0.3.6 (2024-04-17)
=======================

New Features
------------
- Add an option to allow querying the rotation matrices that would convert any
  vector from the (moving) satellite frame to the ECI-aligned frame to the
  `propagate_many` function. It is noted that the origin of this ECI-aligned
  frame is still at the satellite center (See function description.)

Trivia
~~~~~~

- Migrated to a new package structure (`pyproject.toml`-based). [#29]

0.3.5 (2022-11-19)
=======================
This is a pure maintenance release. No features have been added.

0.3.4 (2022-09-07)
=======================

New Features
------------
- Add an option to choose the satellite coordinate frame in the
  `propagate_many` function. The two available choices differ in the order
  of the Cartesian axes (how x/y/z are mapped to motion and nadir directions)
  and the choice of angles. (See function description.)

0.3.3 (2020-05-13)
=======================

Bugfixes
~~~~~~~~~~
- `PyDateTime` was using wrong time zone if current time was requested (e.g.,
  by calling it with empty init arguments, or setting the `datetime` property
  to None). The wrong behavior had its origin in a call to `datetime.now()`
  while the correct call must be `datetime.utcnow` (as `PyDateTime` is
  defined in the UTC frame. [#13]

0.3.2 (2020-03-26)
=======================

New Features
------------
- Add observer coordinates (in eci frame and moving satellite frame) to
  propagate_many function. This can be useful if one wants to calculate
  satellite communication antenna gains w.r.t. observer on Earth. [#12]

Bugfixes
~~~~~~~~~~
- In the `tle_linestrings_from_orbital_parameters` two parameters were missing
  (eccentricity and argument of perigee), which are necessary for a complete
  orbit description. [#11]

0.3.1 (2020-01-07)
=======================

New Features
------------
- Migrate to the latest `libsgp4` version (master, hash: f5cb54b3), which can
  be downloaded from Dan Warners repository: https://github.com/dnwrnr/sgp4.
  Thanks to Cees Bassa for pointing out that the ancient version that was
  previously in use was having accuracy issues.

0.3.0 (2019-11-28)
=======================
Add a benchmark section to the User Manual.

Bugfixes
~~~~~~~~~~
- The pytest benchmarks were not working properly (so we had them skipped),
  because it was not clear how to check if the pytest-benchmark plugin was
  installed in the first place. This is now clear: can add a config routine
  to `conftest.py`.

0.2.20 (2019-11-27)
=======================
Finish the User Manual.

0.2.19 (2019-11-27)
=======================

New Features
------------
- Add two example tutorial notebooks.
- Introduce several helper/utitlity functions, e.g. to create your own TLEs
  from orbital parameters.

Bugfixes
~~~~~~~~~~
- Include tests in (installed) package, such that test runner can find them.
- Fix github pages rendering, where css/js wasn't found in static directories
  with a leading underscore (an empty .nojekyll is necessary!)

0.2.18 (2019-11-24)
=======================

Bugfixes
~~~~~~~~~~
- Really show README on PyPI.

0.2.17 (2019-11-24)
=======================

Bugfixes
~~~~~~~~~~
- Show README on PyPI.
- Twine upload didn't work for Python3.5.

0.2.16 (2019-11-24)
=======================
Add wheels for MacOS and Windows.


0.2.0 (2019-11-21)
=======================
Initial release
