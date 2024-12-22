from variables import *
from libraries import *
from classes import *

# Defining fuctions to calculating the phase, chirp and wavelet
def getPhase(pulse):
    phi=np.unwrap(np.angle(pulse)) #Get phase starting from 1st entry
    phi=phi-phi[int(len(phi)/2)]   #Center phase on middle entry
    return phi    

def getChirp(time,pulse):
    phi=getPhase(pulse)
    dphi=np.diff(phi ,prepend = phi[0] - (phi[1]  - phi[0]  ),axis=0) #Change in phase. Prepend to ensure consistent array size 
    dt  =np.diff(time,prepend = time[0]- (time[1] - time[0] ),axis=0) #Change in time.  Prepend to ensure consistent array size
    return -1/(2*pi)*dphi/dt

def wavelet(t,duration_s,frequency_Hz):
    wl = np.exp(-1j*2*pi*frequency_Hz*t)*np.sqrt(np.exp(-0.5*(t/duration_s)**2 )/np.sqrt(2*pi)/duration_s)
    return wl

# Defining Functions to simulate a Gaussian pulse
# # Function returns pulse power or spectrum PSD
def getPower(amplitude):
    return np.abs(amplitude)**2

# Function gets the energy of a pulse or spectrum by integrating the power
def getEnergy(time_or_frequency,amplitude):
    return np.trapz(getPower(amplitude),time_or_frequency)

def GaussianPulseTime(time,amplitude,duration):
    return amplitude*np.exp(-2*np.log(2)*((time)/(duration))**2)*(1+0j)

def GaussianPulseFrequency(frequency,amplitude,duration):
    return 2*amplitude*duration*np.sqrt(pi/(8*np.log(2)))*np.exp(-((duration**2)/(8*np.log(2)))*(2*pi*frequency)**2)*(1+0j)

# Getting the spectrum based on a given pulse
def getSpectrumFromPulse(time,frequency,pulse_amplitude):
    pulseEenergy=getEnergy(time,pulse_amplitude) # Get pulse energy
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    spectrumEnergy=getEnergy(frequency,spectrum_amplitude) # Get spectrum energy
    err=np.abs((pulseEenergy/spectrumEnergy-1))
    #assert( err<1e-7 ), f'ERROR = {err}: Energy changed when going from Pulse to Spectrum!!!'
    return spectrum_amplitude

def getPulseFromSpectrum(time,frequency,spectrum_aplitude):
    spectrumEnergy=getEnergy(frequency,spectrum_aplitude)
    dt=time[1]-time[0]
    pulse=ifft(ifftshift(spectrum_aplitude))/dt
    pulseEnergy=getEnergy(time,pulse)
    err=np.abs((pulseEnergy/spectrumEnergy-1))
    #assert( err<1e-7 ), f'ERROR = {err}: Energy changed when going from Spectrum to Pulse!!!'
    return pulse

# Equivalent function for generating a Gaussian spectrum
def GaussianSpectrum(time,frequency,amplitude,bandwidth):
    return getSpectrumFromPulse(time,frequency,GaussianPulseTime(time,amplitude,1/bandwidth))

# Getting FWHM based on a given pulse
# Find the FWHM of the frequency/time domain of the signal
def FWHM(X, Y):
    deltax = X[1] - X[0]
    half_max = max(Y) / 2.
    l = np.where(Y > half_max, 1, 0)
    return np.sum(l) * deltax

# The first linear term
def GVD_term(fiber,sim):
    return fiber.beta2/2*(2*pi*sim.f)**2

# The second linear term
def TOD_term(fiber,sim):
    return fiber.beta3/6*(2*pi*sim.f)**3

# The third linear term
def loss_term(fiber,sim):
    return -fiber.alpha_dB_per_m/2

# The SPM nonlinear term
def SPM_term(pulse):
    pulse_power = getPower(pulse)
    return pulse_power

# The self-steepening nonlinear term
def self_steepening_term(sim,fiber,pulse,):
    pulse_power = getPower(pulse)
    output = 1j/2/pi/(sim.frequency0_without_correction) / (pulse+np.sqrt(np.max(pulse_power))/1e6*(1+0j))*np.gradient(pulse_power*pulse, sim.t)
    return output if fiber.self_steepening else 0

# Defining the Split-Step Fourier Method's function
def SSFM(fiber:Fiber_config,sim:SIM_config,pulse):
    # Initialize arrays to store pulse and spectrum throughout fiber
    pulseMatrix=np.zeros((fiber.nsteps+1,len(sim.t)),dtype=np.complex_)
    spectrumMatrix=np.copy(pulseMatrix)
    pulseMatrix[0,:]=pulse
    spectrumMatrix[0,:]=getSpectrumFromPulse(sim.t,sim.f,pulse)

    # Pre-calculate effect of dispersion and loss as it's the same everywhere
    disp_and_loss=np.exp((1j*(GVD_term(fiber,sim)+TOD_term(fiber,sim))+loss_term(fiber,sim))*fiber.dz)

    # Precalculate constants for nonlinearity
    nonlinearity=1j*fiber.gamma*fiber.dz

    for n in range(fiber.nsteps):
        pulse*=np.exp(nonlinearity*(SPM_term(pulse)+self_steepening_term(sim,fiber,pulse))) # Apply nonlinearity
        spectrum=getSpectrumFromPulse(sim.t,sim.f,pulse)*disp_and_loss # Go to spectral domain and apply disp and loss
        pulse=getPulseFromSpectrum(sim.t,sim.f,spectrum) # Return to the time domain
        
        # Store results and repeat
        pulseMatrix[n+1,:]=pulse
        spectrumMatrix[n+1,:]=spectrum

    # return results
    return pulseMatrix, spectrumMatrix

def savePlot(fileName):
    if not os.path.isdir('results/'):
        os.makedirs('results/')
    plt.savefig('results/%s.png'%(fileName))

def plotFirstAndLastPulse(matrix, sim:SIM_config):
  t=sim.t*1e15
  plt.figure()
  plt.title("Initial pulse and final pulse")
  power = getPower(matrix[0,:])
  maximum_power=np.max(power)
  plt.plot(t,getPower(matrix[0,:])/maximum_power,label="Initial Pulse")
  plt.plot(t,getPower(matrix[-1,:])/maximum_power,label="Final Pulse")
  plt.axis([-5*duration*1e15,5*duration*1e15,0,1])
  plt.xlabel("Time [fs]")
  plt.ylabel("Power [a.u.]")
  plt.legend()
  #savePlot('initial and final pulse')
  plt.show()

def plotPulseMatrix2D(matrix,fiber:Fiber_config,sim:SIM_config):
  #Plot pulse evolution throughout fiber in normalized lin/log scale
  fig, ax = plt.subplots()
  ax.set_title('Distance-time pulse evolution (a.u.)')
  t=sim.t*1e15
  z = fiber.zlocs_array 
  T, Z = np.meshgrid(t, z)
  P=getPower(matrix[:,:])/np.max(getPower(matrix[:,:]))
  P[P<1e-100]=1e-100
  surf=ax.contourf(T, Z, P,levels=40)
  ax.set_xlabel('Time [fs]')
  ax.set_ylabel('Distance [m]')
  cbar=fig.colorbar(surf, ax=ax)
  ax.set_xlim(left=-5*duration*1e15)
  ax.set_xlim(right=5*duration*1e15)
  savePlot('distance-time pulse evolution')
  plt.show()

def plotFirstAndLastSpectrum(matrix,sim:SIM_config,FWHM_frequency_final):
    f=sim.f_abs_plot/1e15
    frequency0_without_correction=sim.frequency0_without_correction
    plt.figure()
    plt.title("Initial spectrum and final spectrum")
    power = getPower(matrix[0,:])
    maximum_power=np.max(power)
    plt.plot(f,getPower(matrix[0,:])/maximum_power,label="Initial Spectrum")
    plt.plot(f,getPower(matrix[-1,:])/maximum_power,label="Final Spectrum")
    plt.axis([frequency0_without_correction/1e15-FWHM_frequency_final*1e-15,frequency0_without_correction/1e15+FWHM_frequency_final*1e-15,0,1])
    plt.xlabel("Frequency [PHz]")
    plt.ylabel("Power spectral density [a.u.]")
    plt.legend()
    savePlot('initial and final spectrum')
    plt.show()

def PSD_wavelength(matrix,sim:SIM_config):
    wavelength=sim.wavelength_abs*1e9
    wavelength0=sim.wavelength0*1e9
    power=getPower(matrix[0,:])*2*pi*speed_of_light/wavelength**2
    maximum_power=np.max(power)
    plt.plot(wavelength,(getPower(matrix[0,:])*2*pi*speed_of_light/wavelength**2)/maximum_power,label="Initial Spectrum")
    plt.plot(wavelength,(getPower(matrix[-1,:])*2*pi*speed_of_light/wavelength**2)/maximum_power,label="Final Spectrum")
    plt.title('Power spectral density as function of the wavelength')
    plt.axis([wavelength0-0.001*wavelength0,wavelength0+0.001*wavelength0,0,1])
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Power spectral density [a.u.]")
    plt.legend()
    savePlot('power spectral density in function of wavelength')
    plt.show()

def plotSpectrumMatrix2D_1(matrix,fiber:Fiber_config,sim:SIM_config,FWHM_frequency_final):
  #Plot pulse evolution throughout fiber in normalized lin/log scale
  frequency0=sim.frequency0_without_correction/1e15
  fig, ax = plt.subplots()
  ax.set_title('Distance-spectrum evolution (a.u.)')
  f=sim.f_abs_plot/1e15
  z = fiber.zlocs_array 
  F, Z = np.meshgrid(f, z)
  Pf=getPower(matrix[:,:])/np.max(getPower(matrix[:,:]))
  Pf[Pf<1e-100]=1e-100
  surf=ax.contourf(F, Z, Pf,levels=40)
  ax.set_xlabel('Frequency [PHz]')
  ax.set_ylabel('Distance [m]')
  ax.set_xlim(left=frequency0 - FWHM_frequency_final*1e-15)
  ax.set_xlim(right=frequency0 + FWHM_frequency_final*1e-15)
  cbar=fig.colorbar(surf, ax=ax)
  savePlot('distance-spectrum evolution') 
  plt.show()

def plot_spectrogram(sim_config:SIM_config, pulse, nrange_pulse, nrange_spectrum, time_resolution_s:float, label=None):
    t=sim_config.t
    fc=sim_config.frequency0_without_correction
    f = sim_config.f_abs_plot
    Nmin_pulse = np.max([int(sim_config.number_of_points / 2 - nrange_pulse), 0])
    Nmax_pulse = np.min([int(sim_config.number_of_points / 2 + nrange_pulse),sim_config.number_of_points - 1,])
    Nmin_spectrum = np.max([int(sim_config.number_of_points / 2 - nrange_spectrum), 0])
    Nmax_spectrum = np.min([int(sim_config.number_of_points / 2 + nrange_spectrum),sim_config.number_of_points - 1,])
    t=t=sim_config.t[Nmin_pulse:Nmax_pulse]
    pulse = pulse[Nmin_pulse:Nmax_pulse]
    f_rel = sim_config.f[Nmin_spectrum:Nmax_spectrum]
    result_matrix = np.zeros((len(f_rel),len(t)))*1j
    for idx, f_Hz in enumerate(f_rel):
        current_wavelet = lambda time: wavelet(time,time_resolution_s,f_Hz)
        result_matrix[idx,:] = signal.fftconvolve(current_wavelet(t), pulse, mode='same')
    Z = np.abs(result_matrix) ** 2
    Z /= np.max(Z)
    fig, ax = plt.subplots(dpi=300)
    ax.set_title('Wavelet transform (spectrogram) of the %s pulse'%(label))
    T, F = np.meshgrid(t, f[Nmin_spectrum:Nmax_spectrum])
    surf = ax.contourf(T*1e15, F/1e15, Z , levels=40)
    ax.set_xlabel(f"Time [fs]")
    ax.set_ylabel(f"Frequency [PHz]")
    tkw = dict(size=4, width=1.5)
    #ax.yaxis.label.set_color('b')
    n_ticks = len(ax.get_yticklabels())-2
    norm=plt.Normalize(0,1)
    cbar = fig.colorbar(surf, ax=ax)
    text='spectrogram of the ' + label + 'pulse'
    savePlot(text) 
    plt.show()

    
# Function to calculate the refractive index using Sellmeier's formula in the frequency domain.
def getRefractiveIndex(frequency,B1,B2,B3,C1,C2,C3):
    first_term = 1
    second_term = (B1*speed_of_light**2)/(speed_of_light**2-C1*frequency**2)
    third_term = (B2*speed_of_light**2)/(speed_of_light**2-C2*frequency**2)
    fourth_term = (B3*speed_of_light**2)/(speed_of_light**2-C3*frequency**2)
    return np.sqrt(first_term + second_term + third_term + fourth_term)

# Function to calculate the extintion coefficient.
def getExtintionCoefficient(refractiveIndex):
    return 1.0 + hilbert(refractiveIndex) / pi

# Function to calculate the absorption coefficient.
def getAbsorptionCoefficient(ExtintionCoefficient,frequency):
    return 4*pi*frequency*ExtintionCoefficient/speed_of_light

def plotRefractiveIndex(x,y):
    plt.figure()
    plt.title(f"Refractive index in the frequency domain")
    plt.plot(x/1e15,y)
    plt.xlabel("Frequency [PHz]")
    plt.ylabel("Refractive index")
    plt.legend()
    plt.show()

def plotExtintionCoefficient(x,y):
    plt.figure()
    plt.title(f"Extintion coefficient in the frequency domain")
    plt.plot(x/1e15,y)
    plt.xlabel("Frequency [PHz]")
    plt.ylabel("Extintion coefficient")
    plt.legend()
    plt.show()

def plotAbsorptionCoefficient(x,y):
    plt.figure()
    plt.title(f"Absorption coefficient in the frequency domain")
    plt.plot(x/1e15,y*8.69)
    plt.xlabel("Frequency [PHz]")
    plt.ylabel("Absorption coefficient (dB/m)")
    plt.legend()
    plt.show()

# Calculating the GVD from Sellmeier formula
# Function to calculate the first derivative of refractive index using Sellmeier's formula in the frequency domain.
def getFirstDerivateveOfRefractiveIndex(frequency,refractiveIndex,B1,B2,B3,C1,C2,C3):
    first_term = B1*C1/((speed_of_light/frequency)**2-C1)**2
    second_term = B2*C2/((speed_of_light/frequency)**2-C2)**2
    third_term = B3*C3/((speed_of_light/frequency)**2-C3)**2
    return -(speed_of_light/(refractiveIndex*frequency))*(first_term + second_term + third_term)

# Function to calculate the second derivative of refractive index using Sellmeier's formula in the frequency domain.
def getSecondDerivateveOfRefractiveIndex(frequency,refractiveIndex,firstDerivateveOfRefractiveIndex,B1,B2,B3,C1,C2,C3):
    first_term = B1*C1*(3*(speed_of_light/frequency)**2+C1)/((speed_of_light/frequency)**2-C1)**3
    second_term = B2*C2*(3*(speed_of_light/frequency)**2+C2)/((speed_of_light/frequency)**2-C2)**3
    third_term = B3*C3*(3*(speed_of_light/frequency)**2+C3)/((speed_of_light/frequency)**2-C3)**3
    return (1/refractiveIndex)*(first_term + second_term + third_term - (firstDerivateveOfRefractiveIndex)**2)

# Function to calculate the group velocity dispersion using Sellmeier's formula in the frequency domain.
def getGroupVelocityDispersion(frequency,secondDerivateveOfRefractiveIndex):
    return (1/(2*pi*frequency))*(speed_of_light**2/frequency**2)*(1/speed_of_light)*secondDerivateveOfRefractiveIndex

def plotFirstDerivateveOfRefractiveIndex(x,y):
    plt.figure()
    plt.title(f"First derivative of refractive index in the frequency domain")
    plt.plot(x/1e15,y)
    plt.xlabel("Frequency [PHz]")
    plt.ylabel("First derivative of refractive index [1/m]")
    plt.legend()
    plt.show()

def plotSecondDerivateveOfRefractiveIndex(x,y):
    plt.figure()
    plt.title(f"Second derivative of refractive index in the frequency domain")
    plt.plot(x/1e15,y)
    plt.xlabel("Frequency [PHz]")
    plt.ylabel("Second derivative of refractive index [1/m^2]")
    plt.legend()
    plt.show()

def plotGroupVelocityDispersion(x,y):
    plt.figure()
    plt.title(f"Group velocity dispersion in the frequency domain")
    plt.plot(x/1e15,y*10e24)
    plt.xlabel("Frequency [PHz]")
    plt.ylabel("Group velocity dispersion [fs^2/m]")
    plt.legend()
    plt.show()