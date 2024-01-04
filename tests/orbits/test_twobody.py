from orbitalsynth.orbits.twobody import TwoBody
import numpy as np

# Make a circular orbit using each constructor and plot using both methods...
#==============================================================================

# TwoBody.circular(2, 1, 1).track_orbit()
# TwoBody.from_m_ε_a1_R_r1(1, 1, 0, 1, [0,0], [1,0]).track_orbit()
# TwoBody.from_m_R_r1_v1(1, 1, [0,0], [0, 1], [np.sqrt(6.67430 * 1e-11 / 4), 0]).track_orbit()

# Test some other orbits...
#==============================================================================

# TwoBody.from_m_R_r1_v1(2, 1, [0,0], [0, 1], [np.sqrt(6.67430 * 1e-11 / 2), 0]).plot_orbit()
# TwoBody.from_m_R_r1_v1(2, 1, [0,0], [0, 1], [np.sqrt(6.67430 * 1e-11 / 2), 0]).track_orbit()