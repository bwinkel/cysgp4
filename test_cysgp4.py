#!/usr/bin/python


import cysgp4
from datetime import datetime
from ralib.cyastrometry import astrometry



tle_name = 'ASTRA 2F'
tle_line1 = '1 38778U 12051A   12288.95265372  .00000136  00000-0  00000+0 0   217'
tle_line2 = '2 38778 000.0698 254.6769 0000479 231.1384 284.5280 01.00269150   226'
#tle_name = 'UK-DMC 2'
#tle_line1 = '1 35683U 09041C   12289.23158813  .00000484  00000-0  89219-4 0  5863'
#tle_line2 = '2 35683  98.0221 185.3682 0001499 100.5295 259.6088 14.69819587172294'
TLEstring = '\n'.join([tle_name, tle_line1, tle_line2])

TLEstring = """ISS
1 25544U 98067A   14191.54936479  .00016717  00000-0  10270-3 0  9011
2 25544  51.6452 343.7596 0003472 238.6843 121.3969 15.50506927 14952"""
tle_name, tle_line1, tle_line2 = TLEstring.split('\n')


tle = cysgp4.PyTle(tle_name, tle_line1, tle_line2)
print tle

eff_observer = cysgp4.PyObserver()

mysat = cysgp4.Satellite(tle, eff_observer)
#mysat.mjd = 54444.123
#print mysat.datetime

mjdnow = astrometry.astrometry().current_mjd()

mysat.topo_pos(mjdnow)
mysat.geo_pos(mjdnow)

from ralib.astrometry import satellite_from_tle
sO = satellite_from_tle.satelliteObserver()
sO.obtainAzElFromTLE(TLEstring, mjdnow)





t = cysgp4.PyCoordTopocentric()

now = datetime.now()
dt = cysgp4.PyDateTime(now)
dt.set_datetime(datetime.now())


eci = cysgp4.PyEci(cysgp4.PyDateTime(datetime.now()), cysgp4.PyObserver().location)

# -----------------------------------
import cysgp4
from datetime import datetime
from ralib.cyastrometry import astrometry
from ralib.astrometry import satellite_from_tle

sO = satellite_from_tle.satelliteObserver()

def get_azel_from_tle(tle_string, mjd):
    tle_name, tle_line1, tle_line2 = tle_string.split('\n')
    eff_observer = cysgp4.PyObserver()
    tle = cysgp4.PyTle(tle_name, tle_line1, tle_line2)
    mysat = cysgp4.Satellite(tle, eff_observer)
    res = mysat.topo_pos(mjd)
    return tle_name, res.azimuth, res.elevation


TLEstring = """ISS
1 25544U 98067A   14191.54936479  .00016717  00000-0  10270-3 0  9011
2 25544  51.6452 343.7596 0003472 238.6843 121.3969 15.50506927 14952"""
mjdnow = astrometry.astrometry().current_mjd()
sO.obtainAzElFromTLE(TLEstring, mjdnow)
get_azel_from_tle(TLEstring, mjdnow)

#% timeit -r 3 -n 3 sO.obtainAzElFromTLE(TLEstring, mjdnow)
#% timeit -r 3 -n 3 get_azel_from_tle(TLEstring, mjdnow)
