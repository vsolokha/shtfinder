import freegs
import numpy as np
from scipy import interpolate
from scipy.signal import butter, lfilter

class MirnovCoil():
    area    = 1.0
    n       = 1.0 
    lowcut  = 10.0
    highcut = 50.0

    def __init__(self, RZ, at=0.0, an=0.0):
        self.R  = RZ[0]
        self.Z  = RZ[1]
        self.at = at
        self.an = an 
    
    def get_psi(self, eq):
        return eq.psiRZ(self.R, self.Z)

    def set_exp_data(self, shot_data, lowcut, highcut):
        return self.butter_bandpass_filter(shot_data, lowcut, highcut)

    def butter_bandpass(self, lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype='band')
        return b, a

    def butter_bandpass_filter(self, data, lowcut, highcut, fs, order=5):
        b, a = self.butter_bandpass(lowcut, highcut, fs, order)
        y = lfilter(b, a, data)
        return y

    def get_voltage(self, times, eqs, n_spl=1e3):
        psis = []
        for eq in eqs:
            psi = self.get_psi(eq)
            psis.append(psi)
        psi_spl   = interpolate.UnivariateSpline(times, psis, k=1, s=0)
        volt_spl  = psi_spl.derivative()
        t         = np.linspace(times[0], times[-1], int(n_spl))
        fs        = 1./np.diff(times).mean()
        filt_volt = self.butter_bandpass_filter(volt_spl(t), self.lowcut, self.highcut, fs)
        return filt_volt, t, volt_spl, psi_spl
    
