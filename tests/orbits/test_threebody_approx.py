from orbitalsynth.orbits.threebody_approx import ThreeBodyApprox
from orbitalsynth.orbits.twobody import TwoBody
import numpy as np
from orbitalsynth.config import *

# Test that realistic orbital parameters give expected results...
#==============================================================================

# CASE 1: sun-earth-moon
# m2 = masses['moon']; T1 = periods['moon']; a2 = radii['moon']
# T0 = periods['earth']; a1 = radii['earth']
# twobody = TwoBody.from_m1_T_ε_a1_R_r1(m2, T1, 0, a2, [0,0], [a2, 0])
# threebody = ThreeBodyApprox.from_twobody1_T0_ε0_A_R(twobody, T0, 0, a1, [a1,0])
# threebody()
# print(f'We expect that m0 = {masses["sun"]:.3E}.')

# Test that we get expected orbital plots...
#==============================================================================

# CASE 1: sun-earth-moon
# m2 = masses['moon']; T1 = periods['moon']; a2 = radii['moon']
# T0 = periods['earth']; a1 = radii['earth']
# twobody = TwoBody.from_m1_T_ε_a1_R_r1(m2, T1, 0, a2, [0,0], [a2, 0])
# threebody = ThreeBodyApprox.from_twobody1_T0_ε0_A_R(twobody, T0, 0, a1, [a1,0])
# threebody.plot_orbit()

# CASE 2: integer period ratio, check that it rotates properly
# m2 = masses['moon']; T1 = periods['earth'] / 10; a2 = radii['earth'] / 10
# T0 = periods['earth']; a1 = radii['earth']
# twobody = TwoBody.from_m1_T_ε_a1_R_r1(m2, T1, 0, a2, [0,0], [a2, 0])
# threebody = ThreeBodyApprox.from_twobody1_T0_ε0_A_R(twobody, T0, 0, a1, [a1,0])
# threebody.plot_orbit()
# threebody.track_orbit(save = True)
# twobody = TwoBody.from_m1_T_ε_a1_R_r1(m2, T1, 0, a2, [0,0], [-a2, 0])
# threebody = ThreeBodyApprox.from_twobody1_T0_ε0_A_R(twobody, T0, 0, a1, [a1,0])
# threebody.plot_orbit()
# threebody.track_orbit()

# CASE 3: two equal-mass planets
# m0 = masses['sun']; m1 = masses['earth']*10000; m2 = m1
# a0 = radii['earth']; a2 = radii['earth'] / 20
# twobody = TwoBody.from_m_ε_a1_R_r1(m1, m2, 0, a2, [0,0], [a2,0])
# threebody = ThreeBodyApprox.from_twobody1_m0_ε0_A_R(twobody, m0, 0, a0, [a0,0])
# threebody.plot_orbit()
# threebody.track_orbit(save=True)

# CASE 4: square
# threebody = ThreeBodyApprox.square(masses['moon'], radii['earth'])
# threebody.track_orbit()