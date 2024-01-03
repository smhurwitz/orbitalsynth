class Keplerian():
    r"""
    Abstract class to calculate the motion of various bodies under gravity.
    """

    G = 6.67430 * 1e-11

    #==========================================================================
    #  ORBITAL CALCULATIONS
    #==========================================================================

    def r_from_ν(self, ν): 
        r"""
        Returns the 2D orbital position vector for some true anomaly ν.
        """
        raise NotImplementedError('I need to be implemented!')
    
    def r_from_t(self, t): 
        r"""
        Returns the 2D orbital position vector for some time t.
        """
        raise NotImplementedError('I need to be implemented!')
    
    #==========================================================================
    #  VISUALIZATION
    #==========================================================================

    def plot_orbit(self, N):
        r"""
        Plots one full period of orbit.
        """
        raise NotImplementedError('I need to be implemented!')


    def track_orbit(self, vid_len, save):
        r"""
        Creates a movie tracking the motion of the orbiting bodies.
        """
        raise NotImplementedError('I need to be implemented!')