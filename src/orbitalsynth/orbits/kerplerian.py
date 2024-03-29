from abc import ABC, abstractmethod, abstractproperty
from .orbit import Orbit
class Keplerian(Orbit):
    r"""
    Abstract class to calculate the Keplerian motion of various bodies under gravity.
    """

    #==========================================================================
    #  ORBITAL CALCULATIONS
    #==========================================================================

    @abstractmethod
    def r_from_ν(self, ν): 
        r"""
        Returns the 2D orbital position vector for some true anomaly ν.
        """
        pass
    
    #==========================================================================
    #  VISUALIZATION
    #==========================================================================

    @abstractmethod
    def plot_orbit(self, N):
        r"""
        Plots one full period of orbit.
        """
        pass