import ephem

obs = ephem.Observer()
obs.lon = "-105.2630"
obs.lat = "40.00371"
obs.elev = 1653
obs.date = "2018/08/03 00:00:00"

entry = "2002QF15, e, inclination, ascending node, argument of perihelion, semimajor axis,, eccentricity, mean anomaly, epoch for mean anomaly, 2000,,"

asteroid = ephem.readdb(entry)

asteroid.compute(obs)

print(asteroid.a_ra, asteroid.a_dec)

