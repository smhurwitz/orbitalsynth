import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
from scipy.optimize import brentq
from .kerplerian import Keplerian

class OneBody(Keplerian):
    r"""
    Class to calculate the elliptical motion of one massive body about another 
    stationary massive body. The constructors for this method are overloaded 
    and take the following parameters in various combinations:

    Parameters
    ==========
    m : The mass of the stationary body.
    ε : The eccentricity of the ellipse.
    a : The length of the semi-major axis of the ellipse.
    tp : The time at perihelion.
    θi : The rotation of the semi-major axis relative to the Cartesian x-axis.
    r : The 2D Cartesian heliocentric displacement at some time t.
    v : The 2D Cartesian heliocentric velocity at some time t.
    """

    #==========================================================================
    #  CONSTRUCTORS
    #==========================================================================

    def __init__(self, m, ε, a, θi, tp):
        self.m = m
        self.ε = ε
        self.a = a
        self.θi = θi
        self.tp = tp

        self.β = ε / (1 + np.sqrt(1 - ε**2))
        self.p = a * (1 - ε**2)
        self.T = np.sqrt(4 * np.pi**2 * a**3 / (super().G * m))

        if ε >= 1: raise ValueError("The parameters chosen result in ε≧1.")
        super().__init__()

    @classmethod
    def from_m_ε_a_r(cls, m, ε, a, r):
        if ε == 0: return cls(m, ε, a, 0, 0) 
        T = np.sqrt(4 * np.pi**2 * a**3 / (super().G * m))
        E = np.arccos((1 - np.linalg.norm(r) / a) / ε)
        β = ε / (1 + np.sqrt(1 - ε**2))
        ν = E + 2 * np.arctan2(β * np.sin(E), 1 - β * np.cos(E))
        M = E - ε * np.sin(E) % (2 * np.pi)
        tp = -T * M / (2 * np.pi)
        θi = -ν + np.arctan2(r[1], r[0])
        return cls(m, ε, a, θi, tp) 
    
    @classmethod
    def from_m_r_v(cls, m, r, v):
        r_3D = np.concatenate([r, [0]])
        v_3D = np.concatenate([v, [0]])
        print(np.cross(v_3D, np.cross(r_3D, v_3D)) / (super().G * m))
        ε = np.linalg.norm(np.cross(v_3D, np.cross(r_3D, v_3D)) / (super().G * m) 
                           - r_3D / np.linalg.norm(r_3D)) 
        a = (2 / np.linalg.norm(r) - np.linalg.norm(v)**2 / (super().G * m))**(-1)
        return cls.from_m_ε_a_r(m, ε, a, r)
    
    @classmethod
    def from_T_ε_a_r(cls, T, ε, a, r):
        m = 4 * np.pi**2 * a**3 / (T**2 * super().G)
        return cls.from_m_ε_a_r(m, ε, a, r)

    #==========================================================================
    #  ORBITAL CALCULATIONS
    #==========================================================================

    def r_from_ν(self, ν): 
        direction = np.array([np.cos(ν + self.θi), np.sin(ν + self.θi)])
        return self.p / (1 + self.ε * np.cos(ν)) * direction

    def r_from_t(self, t, rtol=1e-3):
        t = np.atleast_1d(t)
        r = []
        for t in t:
            M = 2 * np.pi * (t - self.tp) / self.T
            divs, Mdiv = divmod(M, 2 * np.pi)
            kepler = lambda E: E - self.ε * np.sin(E) - Mdiv
            E = brentq(kepler, 0, 2 * np.pi, rtol=rtol) + 2 * np.pi * divs
            ν = E + 2 * np.arctan2(self.β * np.sin(E), 1 - self.β * np.cos(E))
            direction = np.array([np.cos(ν + self.θi), np.sin(ν + self.θi)])
            r.append(self.a * (1 - self.ε * np.cos(E)) * direction)
        return np.array(r).transpose()

    #==========================================================================
    #  VISUALIZATION
    #==========================================================================
    
    def plot_orbit(self, N=100):
        ν = np.linspace(0, 2 * np.pi, N)
        r = self.r_from_ν(ν)

        scale = 1.1; scale_dif = 0.5 * (scale - 1)
        dx = scale * (np.max(r[0]) - np.min(r[0]))
        dy = scale * (np.max(r[1]) - np.min(r[1]))
        x_min = np.min(r[0]) - scale_dif * dx; x_max = x_min + dx
        y_min = np.min(r[1]) - scale_dif * dy; y_max = y_min + dy
        axis_len = max(dx, dy)

        plt.figure()  
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        axis.add_patch(Circle((0,0), 0.025 * axis_len))
        plt.plot(r[0], r[1])
        plt.show()

    def track_orbit(self, duration=10, fps=30, save=False):
        N = math.floor(fps * duration) #number of frames
        evals_per_frame = math.floor(max([1000, N]) / N)
        num_evals = evals_per_frame * N
        t = np.linspace(0, self.T, num_evals)
        r = self.r_from_t(t)

        scale = 1.1; scale_dif = 0.5 * (scale - 1)
        dx = scale * (np.max(r[0]) - np.min(r[0]))
        dy = scale * (np.max(r[1]) - np.min(r[1]))
        x_min = np.min(r[0]) - scale_dif * dx; x_max = x_min + dx
        y_min = np.min(r[1]) - scale_dif * dy; y_max = y_min + dy
        axis_len = max(dx, dy)
        planet = Circle((0, 0), 0.01 * axis_len)

        fig = plt.figure()
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        axis.set_axis_off()
        axis.add_patch(Circle((0,0), 0.025 * axis_len))
        axis.add_patch(planet)
        line, = axis.plot([], [], lw = 2)
        def animate(i):  
            eval_no = i * evals_per_frame
            line.set_data(r[0, :eval_no], r[1, :eval_no])  
            planet.set_center((r[0, eval_no], r[1, eval_no]))
            return line, planet
    
        anim = FuncAnimation(fig=fig, func=animate, frames=N, interval=fps)
        if not save: plt.show()
        else: anim.save('1B_orbit.gif')

    #==========================================================================
    #  PRE-DEFINED ORBITS
    #==========================================================================    

    @classmethod
    def circular(cls, m, a): return cls.from_m_ε_a_r(m, 0, a, [0, a])