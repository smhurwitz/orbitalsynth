import numpy as np
import math
import scipy
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import CirclePolygon
import matplotlib.animation as animation
from .onebody import OneBody

class TwoBody():
    r"""
    Class to calculate the elliptical motion of two massive bodies about their
    center-of-mass. The constructors for this method are overloaded and take 
    the following parameters in various combinations:

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

    def __init__(self, m0, m1, R, orbit_1B):
        self.m0 = m0
        self.m1 = m1
        self.R = R
        self.T = orbit_1B.T
        self.orbit_1B = orbit_1B

        super().__init__()

    @classmethod
    def from_m_ε_a1_R_r1(cls, m0, m1, ε, a1, R, r1):
        r0 = ((m0 + m1) * R - m1 * r1) / m0
        r = r1 - r0
        m = m0 + m1
        a = (m0 + m1) * a1 / m0
        orbit = OneBody.from_m_ε_a_r(m, ε, a, r)
        cls(m0, m1, R, orbit)

    @classmethod
    def from_m_R_r1_v1(cls, m0, m1, R, r1, v1):
        r0 = ((m0 + m1) * R - m1 * r1) / m0
        v0 = - m1 * r1 / m0
        r = r1 - r0
        v = v1 - v0
        m = m0 + m1
        orbit = OneBody.from_m_r_v(m, r, v)
        cls(m0, m1, R, orbit)

    #==========================================================================
    #  ORBITAL CALCULATIONS
    #==========================================================================

    def r_from_ν(self, ν): 
        r = self.orbit.r_from_ν(ν)
        r0 = self.R[:, None] - self.m1 / (self.m0 + self.m1) * r
        r1 = self.R[:, None] + self.m0 / (self.m0 + self.m1) * r
        return r0, r1

    def r_from_t(self, t, rtol=1e-3):
        r = self.orbit.r_from_t(t, rtol)
        r0 = self.R[:, None] - self.m1 / (self.m0 + self.m1) * r
        r1 = self.R[:, None] + self.m0 / (self.m0 + self.m1) * r
        return r0, r1
    
    #==========================================================================
    #  VISUALIZATION
    #==========================================================================

    def plot_orbit(self, N = 100):
        ν = np.linspace(0, 2 * np.pi, N)
        r0, r1 = self.r_ν(ν)
        all_x = np.concatenate((r0[0], r1[0]))
        all_y = np.concatenate((r0[1], r1[1]))
        
        scale = 1.1; scale_dif = 0.5 * (scale - 1)
        dx = scale * (np.max(all_x) - np.min(all_x))
        dy = scale * (np.max(all_y) - np.min(all_y))
        axis_len = scale * max(dx, dy)
        x_min =  np.min(all_x) - scale_dif * dx; x_max = x_min + dx
        y_min =  np.min(all_y) - scale_dif * dx; y_max = y_min + dy 
        
        COM = CirclePolygon((self.R[0], self.R[1]), 0.025 * axis_len, resolution=3)

        plt.figure()
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        axis.add_patch(COM)
        plt.plot(r0[0], r0[1])
        plt.plot(r1[0], r1[1])
        plt.show()

    def track_orbit(self, vid_len = 10, save = False):
        fps = 30
        N = math.floor(fps * vid_len)
        vals_per_frame = math.floor(max([1000, N]) / N)
        t = np.linspace(0, self.T, vals_per_frame * N)
        r0, r1 = self.r_t(t)

        axis_len = 1.1 * max(np.max(np.concatenate((r0[0], r1[0]))) - np.min(np.concatenate((r0[0], r1[0]))), 
                             np.max(np.concatenate((r0[1], r1[1]))) - np.min(np.concatenate((r0[1], r1[1]))))
        dx = 1.1 * (np.max(np.concatenate((r0[0], r1[0]))) - np.min(np.concatenate((r0[0], r1[0]))))
        dy = 1.1 * (np.max(np.concatenate((r0[1], r1[1]))) - np.min(np.concatenate((r0[1], r1[1]))))
        x_min =  np.min(np.concatenate((r0[0], r1[0]))) - 0.05 * dx; x_max = x_min + dx
        y_min =  np.min(np.concatenate((r0[1], r1[1]))) - 0.05 * dx; y_max = y_min + dy   
        fig = plt.figure()  
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        axis.set_axis_off()
        axis.add_patch(CirclePolygon((self.R[0],self.R[1]), 0.025 * axis_len, resolution = 3))
        ball0 = Circle((0, 0), 0.01 * axis_len); ball1 = Circle((0, 0), 0.01 * axis_len);
        axis.add_patch(ball0); axis.add_patch(ball1)
        line0, = axis.plot([], [], lw = 2); line1, = axis.plot([], [], lw=2)
        def animate(i):  
            line0.set_data(r0[0,:i * vals_per_frame], r0[1,:i * vals_per_frame])  
            line1.set_data(r1[0,:i * vals_per_frame], r1[1,:i * vals_per_frame])  
            ball0.set_center((r0[0,i * vals_per_frame], r0[1,i * vals_per_frame]))
            ball1.set_center((r1[0,i * vals_per_frame], r1[1,i * vals_per_frame]))
            return line0, line1, ball0, ball1
    
        anim = animation.FuncAnimation(fig=fig, func=animate, frames=N, interval=30)
        if save == False: plt.show()
        else: anim.save('2B_orbit.gif')

    # # ================================ STANDARD ORBITS ================================
    
    @classmethod
    def circular(self, m0, m1, a1): return TwoBody(m0, m1, [0, 0], [0, a1], 0, a1)
    
    # # ================================= VISUALIZATION =================================

    def plot_orbit(self, N = 100):
        ν = np.linspace(0, 2 * np.pi, N)
        r0, r1 = self.r_ν(ν)
        
        axis_len = 1.1 * max(np.max(np.concatenate((r0[0], r1[0]))) - np.min(np.concatenate((r0[0], r1[0]))), 
                             np.max(np.concatenate((r0[1], r1[1]))) - np.min(np.concatenate((r0[1], r1[1]))))
        dx = 1.1 * (np.max(np.concatenate((r0[0], r1[0]))) - np.min(np.concatenate((r0[0], r1[0]))))
        dy = 1.1 * (np.max(np.concatenate((r0[1], r1[1]))) - np.min(np.concatenate((r0[1], r1[1]))))
        x_min =  np.min(np.concatenate((r0[0], r1[0]))) - 0.05 * dx; x_max = x_min + dx
        y_min =  np.min(np.concatenate((r0[1], r1[1]))) - 0.05 * dx; y_max = y_min + dy        
        plt.figure()
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        axis.add_patch(CirclePolygon((self.R[0],self.R[1]), 0.025 * axis_len, resolution = 3))
        plt.plot(r0[0], r0[1])
        plt.plot(r1[0], r1[1])
        plt.show()

    def track_orbit(self, vid_len = 10, save = False):
        fps = 30
        N = math.floor(fps * vid_len)
        vals_per_frame = math.floor(max([1000, N]) / N)
        t = np.linspace(0, self.T, vals_per_frame * N)
        r0, r1 = self.r_t(t)

        axis_len = 1.1 * max(np.max(np.concatenate((r0[0], r1[0]))) - np.min(np.concatenate((r0[0], r1[0]))), 
                             np.max(np.concatenate((r0[1], r1[1]))) - np.min(np.concatenate((r0[1], r1[1]))))
        dx = 1.1 * (np.max(np.concatenate((r0[0], r1[0]))) - np.min(np.concatenate((r0[0], r1[0]))))
        dy = 1.1 * (np.max(np.concatenate((r0[1], r1[1]))) - np.min(np.concatenate((r0[1], r1[1]))))
        x_min =  np.min(np.concatenate((r0[0], r1[0]))) - 0.05 * dx; x_max = x_min + dx
        y_min =  np.min(np.concatenate((r0[1], r1[1]))) - 0.05 * dx; y_max = y_min + dy   
        fig = plt.figure()  
        axis = plt.axes(xlim = (x_min, x_max), ylim = (y_min, y_max))
        axis.set_aspect('equal')
        axis.set_axis_off()
        axis.add_patch(CirclePolygon((self.R[0],self.R[1]), 0.025 * axis_len, resolution = 3))
        ball0 = Circle((0, 0), 0.01 * axis_len); ball1 = Circle((0, 0), 0.01 * axis_len);
        axis.add_patch(ball0); axis.add_patch(ball1)
        line0, = axis.plot([], [], lw = 2); line1, = axis.plot([], [], lw=2)
        def animate(i):  
            line0.set_data(r0[0,:i * vals_per_frame], r0[1,:i * vals_per_frame])  
            line1.set_data(r1[0,:i * vals_per_frame], r1[1,:i * vals_per_frame])  
            ball0.set_center((r0[0,i * vals_per_frame], r0[1,i * vals_per_frame]))
            ball1.set_center((r1[0,i * vals_per_frame], r1[1,i * vals_per_frame]))
            return line0, line1, ball0, ball1
    
        anim = animation.FuncAnimation(fig=fig, func=animate, frames=N, interval=30)
        if save == False: plt.show()
        else: anim.save('2B_orbit.gif')