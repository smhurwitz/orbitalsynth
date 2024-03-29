from orbitalsynth.orbits.onebody import OneBody
import numpy as np
from orbitalsynth.config import *

# Make a circular orbit using each constructor and plot using both methods...
#==============================================================================

# circular = OneBody.circular(2, 3)
# circular = OneBody.from_m_ε_a_r(1, 0, 1, [1, 0])
# circular = OneBody.from_m_r_v(1, [0,1], [np.sqrt(6.67430 * 1e-11),0])
# circular = OneBody.from_T_ε_a_r(769089.7201971824, 0, 1, [1, 0])
# circular.plot_orbit()
# circular.track_orbit()

# Test some other orbits...
#==============================================================================

# OneBody.from_m_ε_a_r(2, 0, 1, [1, 0]).track_orbit()

# Clockwise vs. counterclockwise...
#==============================================================================
orbit = OneBody.from_m_r_v(masses['sun'], radii['earth'] * e['right'], velocities['earth'] * e['down'])
orbit.track_orbit()
