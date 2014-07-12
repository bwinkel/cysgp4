import cysgp4

tle_name = 'UK-DMC 2'
tle_line1 = '1 35683U 09041C   12289.23158813  .00000484  00000-0  89219-4 0  5863'
tle_line2 = '2 35683  98.0221 185.3682 0001499 100.5295 259.6088 14.69819587172294'

tle = cysgp4.PyTle(tle_name, tle_line1, tle_line2)
print tle


eff_observer = cysgp4.PyObserver()

t = cysgp4.PyCoordTopocentric()
