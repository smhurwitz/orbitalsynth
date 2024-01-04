from orbitalsynth.orbits.threebody_approx import ThreeBodyApprox
from orbitalsynth.orbits.twobody import TwoBody
import numpy as np

outer = TwoBody.from_m0_T_ε_a1_R_r1(1, 1e6, 0, 1, [0,0], [1,0])
m1 = outer.m1
inner = TwoBody.from_m0_T_ε_a1_R_r1(m1/2, 0.5 * 1e6, 0, 0.25, [0,0], [0.25,0])

ThreeBodyApprox(outer, inner).plot_orbit()


# t, r, epsilon, (a)

# a2,T1->m1+m2
# a1,T0->m0+m1+m2