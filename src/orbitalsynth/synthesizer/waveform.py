import numpy as np
# import playsound
from ..orbits.onebody import OneBody
from ..orbits.twobody import TwoBody
from ..orbits.threebody_approx import ThreeBodyApprox
import scipy.io.wavfile as wav

class Waveform():
    def __init__(self, orbit, level, sample_rate=44100):
        self.orbit = orbit
        self.level = level #object no to track
        self.sample_rate = sample_rate

    def __call__(self, f, τ, φ): self.play(f, τ, φ)

    def output(self, func, τ, gain):
        samples = self.sample_rate * τ
        output = np.zeros((samples))
        for n in range(samples):
            t = n / self.sample_rate
            output[n] = func(t)
            if n %10000 == 0: print(self.quasiperiod_ε(output[:n]))
        return output * 10**(0.05 * gain) / np.max(np.absolute(output))

    def play(self, f, τ, φ):
        raise NotImplementedError("TODO")
    
    def plotwave(self):
        raise NotImplementedError("TODO")
    
    def quasiperiod_ε(self, wave):
        return np.var(wave) / np.mean(wave)
    
    def save(self, f, τ=3, φ=0, gain=0, filename='sound'):
        filename = "outputs/sounds" + filename + ".wav"
        func = lambda t: self.wave_func(t, f, φ)
        output = self.output(func, τ, gain)
        wav.write(filename, self.sample_rate, output.astype(np.float32))

    def wave_func(self, t, f, φ):
        def r_proj(t):
            if isinstance(self.orbit, OneBody):
                r = self.orbit.r_from_t(t)[:,0]
            else:
                r = self.orbit.r_from_t(t)[self.level][:,0]
            line_of_projection = [np.cos(φ), np.sin(φ)]
            r_proj = np.dot(r, line_of_projection)
            return r_proj
        
        #obtain fundamental freq, use FFT for full 3body.
        if isinstance(self.orbit, OneBody) or isinstance(self.orbit, TwoBody): 
            fundamental = 1 / self.orbit.T
        elif isinstance(self.orbit, ThreeBodyApprox): 
            fundamental = 1 / self.orbit.T0
        else: raise NotImplementedError("This needs to be implemented!")
        t_scaled = (f / fundamental) * t
        
        return r_proj(t_scaled)