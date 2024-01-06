from orbitalsynth.synthesizer.waveform import Waveform
from orbitalsynth.orbits.onebody import OneBody
from orbitalsynth.orbits.twobody import TwoBody
from orbitalsynth.orbits.threebody_approx import ThreeBodyApprox
from orbitalsynth.config import *


#tests: plotting, how rotating works plotting, plotting a single wave, check fundamental frequnecy

# A circular orbit should make a sin wave, no matter how rotated
#==============================================================================
# circle = OneBody.circular(1, 1)
# sin = Waveform(circle, 0)
# sin.save(440, filename='circular')

# ellipse = OneBody.from_m_ε_a_r(1, 0.9, 1, [1,0])
# Waveform(ellipse, 0).save(440, filename='ellipse_0deg')

# square = ThreeBodyApprox.square()
# Waveform(square, 2).save(440, filename='square')

m0 = masses['sun']; m1 = masses['earth']*10000; m2 = m1
a0 = radii['earth']; a2 = radii['earth'] / 20
twobody = TwoBody.from_m_ε_a1_R_r1(m1, m2, 0, a2, [0,0], [a2,0])
spiral = ThreeBodyApprox.from_twobody1_m0_ε0_A_R(twobody, m0, 0, a0, [a0,0])
Waveform(spiral, 1).save(440, filename='spiral')
# spiral.track_orbit()

# m0 = masses['sun']; m1 = masses['earth']*10000; m2 = m1/3
# a0 = radii['earth']; a2 = radii['earth'] / 10
# twobody = TwoBody.from_m_ε_a1_R_r1(m1, m2, 0.7, a2, [0,0], [a2,0])
# spiral = ThreeBodyApprox.from_twobody1_m0_ε0_A_R(twobody, m0, 0, a0, [a0/1.4,a0/1.4])
# spiral.track_orbit()
# Waveform(spiral, 1).save(440, filename='funny1')

# m0 = masses['sun']; m1 = masses['earth']*10000; m2 = m1/1.5
# a0 = radii['earth']; a2 = radii['earth'] / 10
# twobody = TwoBody.from_m_ε_a1_R_r1(m1, m2, 0.7, a2, [0,0], [a2,0])
# spiral = ThreeBodyApprox.from_twobody1_m0_ε0_A_R(twobody, m0, 0.5, a0, [a0/1.4,a0/1.4])
# Waveform(spiral, 1).save(440, filename='funny2')
# spiral.track_orbit()

# m0 = masses['sun']; m1 = masses['earth']*30000; m2 = m1
# a0 = radii['earth']; a2 = radii['earth'] / 10
# twobody = TwoBody.from_m_ε_a1_R_r1(m1, m2, 0.7, a2, [0,0], [a2,0])
# spiral = ThreeBodyApprox.from_twobody1_m0_ε0_A_R(twobody, m0, 0.5, a0, [a0/1.4,a0/1.4])
# spiral.track_orbit(save = True)
# Waveform(spiral, 1).save(440, filename='funny3')

