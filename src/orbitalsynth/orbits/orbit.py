from abc import ABC, abstractmethod, abstractproperty

class Orbit(ABC): 
    r"""
    Abstract class to calculate the motion of various bodies under gravity.
    """

    #==========================================================================
    #  PARAMETERS
    #==========================================================================

    @abstractproperty
    def closed(self): pass

    #==========================================================================
    #  ORBITAL CALCULATIONS
    #==========================================================================

    @abstractmethod
    def r_from_t(self, t): 
        r"""
        Returns the 2D orbital position vector for some time t.
        """
        pass
    
    #==========================================================================
    #  VISUALIZATION
    #==========================================================================

    @abstractmethod
    def track_orbit(self, t, vid_len, save):
        r"""
        Creates a movie tracking the motion of the orbiting bodies.
        """
        pass