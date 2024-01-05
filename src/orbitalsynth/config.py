import numpy as np

#==============================================================================
#  CONSTANTS
#==============================================================================

π = np.pi
G = 6.67430 * 1e-11

masses = {
    "sun": 1.988e30,
    "earth": 5.972e24,
    "moon": 7.3459e22
}

radii = {
    "earth": 1.495978707e11,
    "moon": 3.84e8
}

periods = {
    "earth": 3.154e7,
    "moon": 2.36e6
}

velocities = {
    "earth": 29.78e3,
    "moon": 1.022e3
}

e = {
    "up": np.array([0, 1]),
    "down": np.array([0, -1]),
    "left": np.array([-1, 0]),
    "right": np.array([1, 0])
}