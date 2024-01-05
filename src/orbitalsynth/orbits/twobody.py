import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import CirclePolygon
from matplotlib.animation import FuncAnimation
from scipy.optimize import newton
from .onebody import OneBody
from .kerplerian import Keplerian
from ..config import *

class TwoBody(Keplerian):
    r"""
    Class to calculate the elliptical motion of two massive bodies about their
    center-of-mass. The constructors for this method are overloaded and take 
    the following parameters in various combinations:

    Parameters
    ==========
    m0 : The mass of the first body.
    m1 : The mass of the second body.
    R : The 2D Cartesian center-of-mass of both bodies.
    onebody : The OneBody describing the reduction of the two-body problem.
    ε : The eccentricity of the ellipses.
    a1 : The length of the semi-major axis of the second body's ellipse.
    r1 : The 2D Cartesian position vector for the second body at some time t. 
    v1 : The 2D Cartesian velocity vector for the second body at some time t.
    """

    #==========================================================================
    #  CONSTRUCTORS
    #==========================================================================

    def __init__(self, m0, m1, R, onebody):
        self.m0 = m0
        self.m1 = m1
        self.R = np.array(R)
        self.T = onebody.T
        self.ε = onebody.ε
        self.θi = onebody.θi
        self.a0 = m1 * onebody.a / (m0 + m1)
        self.a1 = m0 * onebody.a / (m0 + m1)
        self.ν_init = onebody.ν_from_E(onebody.E_from_t(0))
        self.onebody = onebody
        self.clockwise = self.onebody.clockwise

        super().__init__()

    def closed(self): True
    
    @classmethod
    def from_m_ε_a1_R_r1(cls, m0, m1, ε, a1, R, r1, clkwise=True):
        R = np.array(R); r1 = np.array(r1)
        r0 = ((m0 + m1) * R - m1 * r1) / m0
        r = r1 - r0
        m = m0 + m1
        a = (m0 + m1) * a1 / m0
        onebody = OneBody.from_m_ε_a_r(m, ε, a, r, clkwise)
        return cls(m0, m1, R, onebody)

    @classmethod
    def from_m_R_r1_v1(cls, m0, m1, R, r1, v1):
        R = np.array(R); r1 = np.array(r1); v1 = np.array(v1)
        r0 = ((m0 + m1) * R - m1 * r1) / m0
        v0 = - m1 * v1 / m0
        r = r1 - r0
        v = v1 - v0
        m = m0 + m1
        orbit = OneBody.from_m_r_v(m, r, v)
        return cls(m0, m1, R, orbit)
    
    @classmethod
    def from_m1_T_ε_a1_R_r1(cls, m1, T, ε, a1, R, r1):
        f = lambda m0: 4 * π**2 * (m0 + m1)**2 * a1**3 - G * T**2 * m0**3
        guess = 4 * π**2 * a1**3 / (G * T**2)
        m0 = newton(f, guess, rtol=1e-3)
        clkwise = True if T > 0 else False
        if m0 <= 0: 
            raise ValueError('The parameters chosen yield a non-physical value of m1. Try increasing m0, increasing T, or decreasing a1.')

        return cls.from_m_ε_a1_R_r1(m0, m1, ε, a1, R, r1, clkwise)

    #==========================================================================
    #  ORBITAL CALCULATIONS
    #==========================================================================

    def r_from_ν(self, ν): 
        r = self.onebody.r_from_ν(ν)
        r0 = self.R[:, None] - self.m1 / (self.m0 + self.m1) * r
        r1 = self.R[:, None] + self.m0 / (self.m0 + self.m1) * r
        return r0, r1

    def r_from_t(self, t, rtol=1e-3):
        r = self.onebody.r_from_t(t, rtol)
        r0 = self.R[:, None] - self.m1 / (self.m0 + self.m1) * r
        r1 = self.R[:, None] + self.m0 / (self.m0 + self.m1) * r
        return r0, r1
    
    #==========================================================================
    #  VISUALIZATION
    #==========================================================================

    def plot_orbit(self, N = 100):
        ν = np.linspace(0, 2 * π, N)
        r0, r1 = self.r_from_ν(ν)
        all_x = np.concatenate((r0[0], r1[0]))
        all_y = np.concatenate((r0[1], r1[1]))
        
        scale = 1.1; scale_dif = 0.5 * (scale - 1)
        dx = scale * (np.max(all_x) - np.min(all_x))
        dy = scale * (np.max(all_y) - np.min(all_y))
        x_min =  np.min(all_x) - scale_dif * dx; x_max = x_min + dx
        y_min =  np.min(all_y) - scale_dif * dx; y_max = y_min + dy
        axis_len = scale * max(dx, dy) 
        COM = CirclePolygon((self.R[0], self.R[1]), 0.025 * axis_len, 
                            resolution=3, color='b')

        plt.figure()
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        axis.add_patch(COM)
        plt.plot(r0[0], r0[1], color='b')
        plt.plot(r1[0], r1[1], color='b')
        plt.show()

    def track_orbit(self, duration=10, fps=30, save=False):
        N = math.floor(fps * duration) #number of frames
        evals_per_frame = math.floor(max([1000, N]) / N)
        num_evals = evals_per_frame * N
        t = np.linspace(0, self.T, num_evals)
        r0, r1 = self.r_from_t(t)
        all_x = np.concatenate((r0[0], r1[0]))
        all_y = np.concatenate((r0[1], r1[1]))

        scale = 1.1; scale_dif = 0.5 * (scale - 1)
        dx = scale * (np.max(all_x) - np.min(all_x))
        dy = scale * (np.max(all_y) - np.min(all_y))
        x_min =  np.min(all_x) - scale_dif * dx; x_max = x_min + dx
        y_min =  np.min(all_y) - scale_dif * dx; y_max = y_min + dy 
        axis_len = scale * max(dx, dy)
        COM = CirclePolygon((self.R[0], self.R[1]), 0.025 * axis_len, 
                            resolution=3, color='b')
        
        fig = plt.figure()  
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        axis.set_axis_off()
        axis.add_patch(COM)
        ball0 = Circle((0, 0), 0.01 * axis_len); ball1 = Circle((0, 0), 0.01 * axis_len);
        axis.add_patch(ball0); axis.add_patch(ball1)
        line0, = axis.plot([], [], lw = 2); line1, = axis.plot([], [], lw=2)
        def animate(i):  
            eval_no = i * evals_per_frame
            line0.set_data(r0[0, :eval_no], r0[1, :eval_no])  
            line1.set_data(r1[0, :eval_no], r1[1, :eval_no])  
            ball0.set_center((r0[0, eval_no], r0[1, eval_no]))
            ball1.set_center((r1[0, eval_no], r1[1, eval_no]))
            return line0, line1, ball0, ball1
    
        anim = FuncAnimation(fig=fig, func=animate, frames=N, interval=fps)
        if not save: plt.show()
        else: anim.save('outputs/2B_orbit.gif')

    #==========================================================================
    #  PRE-DEFINED ORBITS
    #==========================================================================    

    @classmethod
    def circular(cls, m0, m1, a1): 
        return cls.from_m_ε_a1_R_r1(m0, m1, 0, a1, [0, 0], [1, 0])