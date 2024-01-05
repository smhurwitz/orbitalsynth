from orbitalsynth.orbits.twobody import TwoBody
from orbitalsynth.config import *
import numpy as np

# Test that realistic orbital parameters give expected results...
#==============================================================================
sun_earth = TwoBody.from_m1_T_ε_a1_R_r1(masses['earth'], periods['earth'], 0, radii['earth'], [0,0], [radii['earth'], 0])
print(f'The expected mass of the sun is {masses["sun"]:.3E}, and the calculated mass is {sun_earth.m0:.3E}.')

# Make a circular orbit using each constructor and plot using both methods...
#==============================================================================

# TwoBody.circular(2, 1, 1).track_orbit()
# TwoBody.from_m_ε_a1_R_r1(1, 1, 0, 1, [0,0], [1,0]).track_orbit()
# TwoBody.from_m_R_r1_v1(1, 1, [0,0], [0, 1], [np.sqrt(6.67430 * 1e-11 / 4), 0]).track_orbit()

# Test some other orbits...
#==============================================================================

# TwoBody.from_m_R_r1_v1(2, 1, [0,0], [0, 1], [np.sqrt(6.67430 * 1e-11 / 2), 0]).plot_orbit()
# TwoBody.from_m_R_r1_v1(2, 1, [0,0], [0, 1], [np.sqrt(6.67430 * 1e-11 / 2), 0]).track_orbit()

# Test that errors are raised when non-physical parameters are calculated.
#==============================================================================
# TwoBody.from_m0_T_ε_a1_R_r1(1, 1, 0, 1, [0, 0], [1, 0]) #should throw error
# TwoBody.from_m0_T_ε_a1_R_r1(5.92e11, 1, 0, 1, [0, 0], [1, 0]) #should not throw error
# TwoBody.from_m0_T_ε_a1_R_r1(1, 7.7e5, 0, 1, [0, 0], [1, 0]) #should not throw error
# TwoBody.from_m0_T_ε_a1_R_r1(1, 1, 0, 1.18e-4, [0, 0], [1, 0]) #should not throw error
# TwoBody.from_m0_T_ε_a1_R_r1(5.91e11, 1, 0, 1, [0, 0], [1, 0]) #should throw error
# TwoBody.from_m0_T_ε_a1_R_r1(1, 7.68e5, 0, 1, [0, 0], [1, 0]) #should throw error
# TwoBody.from_m0_T_ε_a1_R_r1(1, 1, 0, 1.2e-4, [0, 0], [1, 0]) #should throw error

