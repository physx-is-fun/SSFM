from libraries import *

# constants
pi = np.pi
speed_of_light=3e8 # Speed of light [m/s]

# Defining parameters for the simulation
# Initialize Gaussian pulse parameters (OCTAVIUS-85M-HP from THORLABS)
wavelength0=800*1e-9                                        # Pulse central wavelengt [m]
frequency0_without_correction=speed_of_light/wavelength0    # Pulse central frequency (without correction) [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm
correction_factor=16                                        # Correction factor
frequency0=correction_factor*frequency0_without_correction  # Pulse central frequency (with correction) [Hz]
duration=20*1e-15                                           # Pulse duration in FWHM [s]
assert duration < 20, f"ERROR! 20 fs is the lowest duration which can be simulated without any problem!"
repetition_frequency=85*1e6                                 # Pulse repetition frequency [Hz]
average_power=600*1e-3                                      # Pulse average power [W]
pulse_energy=average_power/repetition_frequency             # Pulse energy [J]
peak_power=pulse_energy/duration                            # Pulse peak power [W]
amplitude = np.sqrt(peak_power)                             # Electrical field strength amplitude in units of sqrt(W)
# Initialize frequency, time and wavelength domain
N=2**15                                                     # Number of points                                                    
dt = 1/frequency0                                           # Time resolution [s]

# Defining the parameters of the fiber
# Define fiberulation parameters
Length=5*1e-2                                                               # Fiber length in m
nsteps=2**10                                                                # Number of steps we divide the fiber into
effective_mode_diameter=5e-6                                                # Effective mode diameter [m] for 780HP single mode fiber @ 850 nm from THORLABS
effective_mode_area=(pi/4)*effective_mode_diameter**2                       # Effective mode area [m^2]
nonlinear_refractive_index=2.7e-20                                          # Nonlinear refractive index [m^2/W] of fused silica @ 800 nm from https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-26-37940&id=424534
gamma=(2*pi*nonlinear_refractive_index)/(wavelength0*effective_mode_area)   # Nonlinear parameter [1/(W*m)]
beta2=36.16                                                                 # Dispersion in fs^2/mm (units typically used when referring to beta2) of fused silica @ 800nm from https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
beta2*=(1e-27)                                                              # Convert fs^2 to s^2 so everything is in SI units of fused silica @ 800nm
beta3=27.47                                                                 # Dispersion in fs^3/mm (units typically used when referring to beta3) of fused silica @ 800nm from https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
beta3*=(1e-42)                                                              # Convert f3^2 to s^3 and mm to m so everything is in SI units of fused silica @ 800nm
alpha_dB_per_m=0.2e-3                                                       # Power attenuation coeff in decibel per m. Usual value at 1550 nm is 0.2 dB/km
# NOTE: beta2>0 is normal dispersion with red light pulling ahead, causing a negative leading chirp
# NOTE: beta2<0 is anomalous dispersion with blue light pulling ahead, causing a positive leading chirp

self_steepening = False

# Define parameters for fused silica
B1  =   0.696166300
B2  =   0.407942600
B3  =   0.897479400	
C1  =   4.67914826e-3 * 1e-12 # in m^2
C2  =   1.35120631e-3 * 1e-12 # in m^2
C3  =   97.9340025    * 1e-12 # in um^2

