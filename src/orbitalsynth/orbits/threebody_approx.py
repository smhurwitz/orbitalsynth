import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import CirclePolygon
from matplotlib.animation import FuncAnimation
from .onebody import OneBody
from .twobody import TwoBody
from .kerplerian import Keplerian

class ThreeBodyApprox(Keplerian):
    r"""
    Class to calculate the elliptical motion of three massive bodies under the
    assumptions that m0≫m1≫m2 and |r0-r1|≫|r1-r2|. The constructors for this 
    method are overloaded and take 
    the following parameters in various combinations:

    Parameters
    ==========
    m0 : The mass of the first body.
    m1 : The mass of the second body.
    R : The 2D Cartesian center-of-mass of both bodies.
    orbit_1B : The OneBody describing the reduction of the two-body problem.
    ε : The eccentricity of the ellipses.
    a1 : The length of the semi-major axis of the second body's ellipse.
    r1 : The 2D Cartesian position vector for the second body at some time t. 
    v1 : The 2D Cartesian velocity vector for the second body at some time t.
    """

    #==========================================================================
    #  CONSTRUCTORS
    #==========================================================================

    def __init__(self, orbit_2B_out, orbit_2B_in):
        self.orbit_2B_out = orbit_2B_out #2B for m0, m1+m_2
        self.orbit_2B_in = orbit_2B_in #2Bfor m1, m2

        #TODO: raise exception if masses on both don't align correctly
        super().__init__()

    #==========================================================================
    #  ORBITAL CALCULATIONS
    #==========================================================================

    def r_from_ν(self, ν_out): 
        ν_in = self.orbit_2B_in.T * ν_out / self.orbit_2B_out.T
        r0, R = self.orbit_2B_out.r_from_ν(ν_out)
        r1_minus_R, r2_minus_R = self.orbit_2B_in.r_from_ν(ν_in)
        r1 = r1_minus_R + R; r2 = r2_minus_R + R
        return r0, r1, r2

    # def r_from_t(self, t, rtol=1e-3):
    #     r = self.orbit_1B.r_from_t(t, rtol)
    #     r0 = self.R[:, None] - self.m1 / (self.m0 + self.m1) * r
    #     r1 = self.R[:, None] + self.m0 / (self.m0 + self.m1) * r
    #     return r0, r1
    
    #==========================================================================
    #  VISUALIZATION
    #==========================================================================

    def plot_orbit(self, N = 100):
        ν = np.linspace(0, 2 * np.pi, N)
        r0, r1, r2 = self.r_from_ν(ν)
        all_x = np.concatenate((r0[0], r1[0], r2[0]))
        all_y = np.concatenate((r0[1], r1[1], r2[1]))
        
        scale = 1.1; scale_dif = 0.5 * (scale - 1)
        dx = scale * (np.max(all_x) - np.min(all_x))
        dy = scale * (np.max(all_y) - np.min(all_y))
        x_min =  np.min(all_x) - scale_dif * dx; x_max = x_min + dx
        y_min =  np.min(all_y) - scale_dif * dx; y_max = y_min + dy
        axis_len = scale * max(dx, dy) 
        # COM = CirclePolygon((self.R[0], self.R[1]), 0.025 * axis_len, 
        #                     resolution=3, color='b')

        plt.figure()
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        # axis.add_patch(COM)
        plt.plot(r0[0], r0[1], color='b')
        plt.plot(r1[0], r1[1], color='b')
        plt.plot(r2[0], r2[1], color='b')
        plt.show()

    # def track_orbit(self, duration=10, fps=30, save=False):
    #     N = math.floor(fps * duration) #number of frames
    #     evals_per_frame = math.floor(max([1000, N]) / N)
    #     num_evals = evals_per_frame * N
    #     t = np.linspace(0, self.T, num_evals)
    #     r0, r1 = self.r_from_t(t)
    #     all_x = np.concatenate((r0[0], r1[0]))
    #     all_y = np.concatenate((r0[1], r1[1]))

    #     scale = 1.1; scale_dif = 0.5 * (scale - 1)
    #     dx = scale * (np.max(all_x) - np.min(all_x))
    #     dy = scale * (np.max(all_y) - np.min(all_y))
    #     x_min =  np.min(all_x) - scale_dif * dx; x_max = x_min + dx
    #     y_min =  np.min(all_y) - scale_dif * dx; y_max = y_min + dy 
    #     axis_len = scale * max(dx, dy)
    #     COM = CirclePolygon((self.R[0], self.R[1]), 0.025 * axis_len, 
    #                         resolution=3, color='b')
        
    #     fig = plt.figure()  
    #     axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
    #     axis.set_aspect('equal')
    #     axis.set_axis_off()
    #     axis.add_patch(COM)
    #     ball0 = Circle((0, 0), 0.01 * axis_len); ball1 = Circle((0, 0), 0.01 * axis_len);
    #     axis.add_patch(ball0); axis.add_patch(ball1)
    #     line0, = axis.plot([], [], lw = 2); line1, = axis.plot([], [], lw=2)
    #     def animate(i):  
    #         eval_no = i * evals_per_frame
    #         line0.set_data(r0[0, :eval_no], r0[1, :eval_no])  
    #         line1.set_data(r1[0, :eval_no], r1[1, :eval_no])  
    #         ball0.set_center((r0[0, eval_no], r0[1, eval_no]))
    #         ball1.set_center((r1[0, eval_no], r1[1, eval_no]))
    #         return line0, line1, ball0, ball1
    
    #     anim = FuncAnimation(fig=fig, func=animate, frames=N, interval=fps)
    #     if not save: plt.show()
    #     else: anim.save('2B_orbit.gif')

    #==========================================================================
    #  PRE-DEFINED ORBITS
    #==========================================================================    

    # @classmethod
    # def circular(cls, m0, m1, a1): 
    #     return cls.from_m_ε_a1_R_r1(m0, m1, 0, a1, [0, 0], [1, 0])