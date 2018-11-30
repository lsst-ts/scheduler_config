import lsst.ts.schedulerConfig.instrument.park
assert type(config)==lsst.ts.schedulerConfig.instrument.park.Park, 'config is of type %s.%s instead of import lsst.ts.schedulerConfig.instrument.park.Park' % (type(config).__module__, type(config).__name__)
# Telescope azimuth (units=degrees) in the park position.
config.telescope_azimuth=45.0

# Dome azimuth (units=degrees) in the park position.
config.dome_azimuth=45.0

# Dome altitude (units=degrees) in the park position.
config.dome_altitude=90.0

# Camera filter for the park position.
config.filter_position='z'

# Telescope altitude (units=degrees) in the park position.
config.telescope_altitude=86.5

# Telescope rotator angle (units=degrees) in the park position.
config.telescope_rotator=0.0

