from scipy.fftpack import fftshift, fftfreq
from variables import speed_of_light, np

# Defining a class for the simulation parameters
# # Class for holding info about the simulation params
class SIM_config:
    def __init__(self,N,dt,frequency0,correction_factor,wavelength0):
        self.wavelength0=wavelength0
        self.correction_factor=correction_factor
        self.number_of_points=N
        self.time_step=dt
        t = np.linspace(0,N*dt,N)                                                                                   # Time step array
        t = t-np.mean(t)
        self.t=t
        self.tmin=self.t[0]                                                                                         # Minimum time
        self.tmax=self.t[-1]                                                                                        # Maximum time
        f=fftshift(fftfreq(N,d=dt))                                                                                 # Frequency step array
        self.f=f
        self.frequency0=frequency0
        self.f_abs=self.f + self.frequency0                                                                          # Generates array of frequencies centered at the carrier frequency.
        self.frequency0_without_correction=self.frequency0/self.correction_factor
        self.f_abs_plot=self.f + self.frequency0_without_correction
        self.fmin=self.f[0]                                                                                          # Minimum frequency
        self.fmax=self.f[-1]                                                                                         # Maximum frequency
        self.freq_step=self.f[1]-self.f[0]
        assert np.min(self.f_abs) >= 0, f"ERROR! Lowest frequency of {np.min(self.f_abs):.3f} is below 0. Consider increasing the center frequency!"                  
        self.wavelength=speed_of_light/self.f                                                                        # Wavelength step array
        self.wavelength_abs=speed_of_light/self.f_abs + self.wavelength0 - self.wavelength0/self.correction_factor   # absolute Wavelength step array

# Class for holding info about the fiber
class Fiber_config:
    def __init__(self,nsteps,L,gamma,beta2,beta3,alpha_dB_per_m,self_steepening):
        self.nsteps=nsteps
        self.ntraces=self.nsteps+1                                           # NOTE: If we want to do 100 steps, we will get 101 calculated pulses
        self.Length=L
        self.dz=L/nsteps
        self.zlocs_array=np.linspace(0,L,self.ntraces)                       # Locations of each calculated pulse
        self.gamma=gamma
        self.beta2=beta2
        self.beta3=beta3
        self.alpha_dB_per_m=alpha_dB_per_m
        self.alpha_Np_per_m=alpha_dB_per_m*np.log(10)/10.0                   # The loss coefficient is usually specified in dB/km, but Nepers/km is more useful for calculations
        self.self_steepening=self_steepening