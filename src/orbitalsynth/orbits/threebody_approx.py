import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import CirclePolygon
from matplotlib.animation import FuncAnimation
from .onebody import OneBody
from .twobody import TwoBody
from .kerplerian import Keplerian
from ..config import *

class ThreeBodyApprox(Keplerian):
    r"""
    Class to calculate the elliptical motion of three massive bodies under the
    assumptions that m0≫m1≫m2 and |r0-r1|≫|r1-r2|. The constructors for this 
    method are overloaded and take 
    the following parameters in various combinations:

    Parameters
    ==========
    mi : The mass of the ith body.
    Ti : The period of the ith system
    a0 : The semi-major axis in the 0 system of the 0th body
    A : The semi-major axis in the 0 system of the 1st & 2nd bodies
    a1 : The semi-major axis in the 1 system of the 1st body
    a2 : The semi-major axis in the 1 system of the 2nd body
    εi : The eccentricity of the ith system
    θii : The inclination angle of the ellipse in the ith system
    νi_init : The initial ν in the ith system
    twobodyi : The ith system
    """

    #==========================================================================
    #  CONSTRUCTORS
    #==========================================================================

    def __init__(self, twobody0, twobody1):
        self.m0 = twobody0.m0; self.m1 = twobody1.m0; self.m2 = twobody1.m1
        self.T0 = twobody0.T; self.T1 = twobody1.T
        self.a0 = twobody0.a0; self.A = twobody0.a1
        self.a1 = twobody1.a0; self.a2 = twobody1.a1
        self.ε0 = twobody0.ε; self.ε1 = twobody1.ε
        self.θi0 = twobody0.θi; self.θi1 = twobody1.θi
        self.ν0_init = twobody0.ν_init; self.ν1_init = twobody1.ν_init
        self.twobody0 = twobody0; self.twobody1 = twobody1
        self.clkwise0 = twobody0.clockwise; self.clkwise1 = twobody1.clockwise
        if not math.isclose(twobody0.m1, self.m1 + self.m2, rel_tol=1e-5):
            raise ValueError('''The parameters must be set such that the masses 
                             of the planet-moon system are identical across 
                             both two-body orbits.''')
        super().__init__()

    def __call__(self):
        print(f'[m0, m1, m2] = [{self.m0:.3E}, {self.m1:.3E}, {self.m2:.3E}]')
        print(f'[a0, A]  = [{self.a0:.3E}, {self.A:.3E}]')
        print(f'[a1, a2] = [{self.a1:.3E}, {self.a2:.3E}]')
        print(f'[T0, T1] = [{self.T0:.3E}, {self.T1:.3E}]')
        print(f'[ε0, ε1] = [{self.ε0:.3E}, {self.ε1:.3E}]')

    @classmethod
    def from_twobody1_T0_ε0_A_R(cls, twobody1, T0, ε0, A, R):
        M = twobody1.m0 + twobody1.m1
        twobody0 = TwoBody.from_m1_T_ε_a1_R_r1(M, T0, ε0, A, [0, 0], R)
        return cls(twobody0, twobody1)
    
    @classmethod
    def from_twobody1_m0_ε0_A_R(cls, twobody1, m0, ε0, A, R, clkwise0=True):
        M = twobody1.m0 + twobody1.m1
        twobody0 = TwoBody.from_m_ε_a1_R_r1(m0, M, ε0, A, [0,0], R, clkwise0)
        return cls(twobody0, twobody1)


    #==========================================================================
    #  ORBITAL CALCULATIONS
    #==========================================================================

    def r_from_ν(self, ν0): 
        ν1 = self.ν1_init + self.T0 / self.T1 * (ν0 - self.ν0_init)
        r0, R = self.twobody0.r_from_ν(ν0)
        r1, r2 = self.twobody1.r_from_ν(ν1) + R
        return r0, r1, r2

    def r_from_t(self, t, rtol=1e-3):
        r0, R = self.twobody0.r_from_t(t)
        # print(self.twobody0.clockwise); exit()
        r1, r2 = self.twobody1.r_from_t(t) + R
        return r0, r1, r2
    
    #==========================================================================
    #  VISUALIZATION
    #==========================================================================

    def plot_orbit(self, N = 1000):
        ν = np.linspace(0, 2 * π, N)
        r0, r1, r2 = self.r_from_ν(ν)
        all_x = np.concatenate((r0[0], r1[0], r2[0]))
        all_y = np.concatenate((r0[1], r1[1], r2[1]))
        
        scale = 1.1; scale_dif = 0.5 * (scale - 1)
        dx = scale * (np.max(all_x) - np.min(all_x))
        dy = scale * (np.max(all_y) - np.min(all_y))
        x_min =  np.min(all_x) - scale_dif * dx; x_max = x_min + dx
        y_min =  np.min(all_y) - scale_dif * dx; y_max = y_min + dy

        plt.figure()
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        plt.plot(r0[0], r0[1], color='b')
        plt.plot(r1[0], r1[1], color='b')
        plt.plot(r2[0], r2[1], color='b')
        plt.show()

    def track_orbit(self, duration=10, fps=30, save=False):
        N = math.floor(fps * duration) #number of frames
        evals_per_frame = math.floor(max([1000, N]) / N)
        num_evals = evals_per_frame * N
        t = np.linspace(0, self.T0, num_evals)
        r0, r1, r2 = self.r_from_t(t)
        all_x = np.concatenate((r0[0], r1[0], r2[0]))
        all_y = np.concatenate((r0[1], r1[1], r2[1]))

        scale = 1.1; scale_dif = 0.5 * (scale - 1)
        dx = scale * (np.max(all_x) - np.min(all_x))
        dy = scale * (np.max(all_y) - np.min(all_y))
        x_min =  np.min(all_x) - scale_dif * dx; x_max = x_min + dx
        y_min =  np.min(all_y) - scale_dif * dx; y_max = y_min + dy 
        axis_len = scale * max(dx, dy)
        
        fig = plt.figure()  
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        axis.set_axis_off()
        ball0 = Circle((0, 0), 0.01 * axis_len); axis.add_patch(ball0); 
        ball1 = Circle((0, 0), 0.01 * axis_len); axis.add_patch(ball1)
        ball2 = Circle((0, 0), 0.01 * axis_len); axis.add_patch(ball2)
        line0, = axis.plot([], [], lw=2)
        line1, = axis.plot([], [], lw=2)
        line2, = axis.plot([], [], lw=2)

        def animate(i):  
            eval_no = i * evals_per_frame
            line0.set_data(r0[0, :eval_no], r0[1, :eval_no])  
            line1.set_data(r1[0, :eval_no], r1[1, :eval_no])  
            line2.set_data(r2[0, :eval_no], r2[1, :eval_no])  
            ball0.set_center((r0[0, eval_no], r0[1, eval_no]))
            ball1.set_center((r1[0, eval_no], r1[1, eval_no]))
            ball2.set_center((r2[0, eval_no], r2[1, eval_no]))
            return line0, line1, line2, ball0, ball1, ball2
    
        anim = FuncAnimation(fig=fig, func=animate, frames=N, interval=fps)
        if not save: plt.show()
        else: anim.save('outputs/3B_orbit.gif')

    #==========================================================================
    #  PRE-DEFINED ORBITS
    #==========================================================================    

    @classmethod 
    def fourier(cls, T0, T1, m2, A, a2, θA, θ2):
        R = A * np.array([np.cos(θA), np.sin(θA)])
        r2 = a2 * np.array([np.cos(θ2), np.sin(θ2)])
        twobody1 = TwoBody.from_m1_T_ε_a1_R_r1(m2, T1, 0, a2, [0,0], [a2,0])
        return cls.from_twobody1_T0_ε0_A_R(twobody1, T0, 0, A, [A,0])

    @classmethod
    def square(cls, m2=masses['moon'], A=radii['earth']):
        T0 = 1
        T1 = -1.0 / 3
        a2 = A * 1.0 / 9
        θA = 0; θ2 = 0
        return cls.fourier(T0, T1, m2, A, a2, θA, θ2)
