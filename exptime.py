#! /usr/bin/python
import numpy as np
from astropy import units as u
from astropy import constants as c
#from astropy.constants import c
#from astropy.constants import h

'''
Exposure time calculator for a list of target stars, desired measurement quality, and optical system specifications. Currently a proof of concept derived from HARPS, but capable to some degree of considering stellar properties.

Makes simplifications, and is currently only valid in optical (-), but can be expanded to red (-) and near IR (-).

Notably, this uses an analytic approximation that SNR is constant over the entire wavelength range. This is because time from SNR is easy, but SNR from V_RMS is hard.
'''

# Beatty spectra (simulated): 100 A chunks; Solar metallicity;
# R = 100k (4000-6500 A), 165000 (6500-10000 A), 88000 (10000A - 25000 A)
# In general need to specify wavelength range, as corrections vary with the one chosen (optical, red, near IR).
# Table from http://www.personal.psu.edu/tgb15/beattygaudi/table1.dat
beatty = np.genfromtxt("table1.dat", dtype=None)
beatty.dtype.names = ('Angstrom', 'Teff', 'Uncertaintykms')

# HARPS
telescope = 3.566 * u.m
#area = np.pi * (telescope/2)**2 # telescope area, central obstruction ignored
area = 8.8564 * u.m * u.m
efficiency = 0.02 #could be ~1-3%, depending on what is measured
λ_min = 3800
λ_max = 6800
Δλ = 100
BeattyWaves = np.arange(λ_min+Δλ/2, λ_max+Δλ/2, Δλ)
R = 110000 #110k or 120k, depending on source
gain = 1/1.42 # ADU/e-, assume 1:1 photon to electron generation
read_noise = 7.07 #RMS of ± spurious electrons
dark_current = 4 / 3600 # (4 e-/hour converted to per second)
n_pix = 4

# Alpha Cen B
Teff = 5214 # 5214 K
logg = 4.37
FeH = 0.23
rstar = 0.865 * u.solRad
dstar = 1.3 * u.pc
vsini = 2*np.pi*rstar/(38.7*u.day) * u.s/u.km
theta_rot = 1.13*vsini.si

BTSettl = np.genfromtxt(str(round(Teff,-2)), dtype=float)
#print("Teff:", Teff, " BT-Settl:", round(Teff,-2), " Beatty RV:", round(Teff/200)*200)

# BT Settl spectra are labled by Teff, and available every 100 K.
# These spectra have a wavelength (Angstroms), and a flux (1e8 erg/s/cm^2/Angstrom) column
# sum of flux(λ)*Δλ(λ)/1e8 == total flux (in power/area) emitted. 
# Multiply bt stellar_radius^2/distance^2 for recieved.
# Spectra downloaded from other sources will use different units!

#De-duping, this is not necessary if spectra are downloaded from another source.
model = [np.array([0,0])]
for x in BTSettl:
	if np.array_equal(x,model[-1]) == False:
		model.append(x)
photons = 0
#Counting up photons for SNR -> exposure time calc
for i in np.arange(0, len(model)-1):
	if((model[i][0] >= λ_min) and (model[i][0] <= λ_max)):
		power = (model[i][1]*1e-8*u.erg/u.cm**2/u.s/u.angstrom) * ((model[i+1][0]-model[i][0]) * u.angstrom)
		power *= rstar**2/dstar**2 #rescaling emitted spectrum based on stellar surface area and distance from us
		#Opacity is treated as universal, but is not. Need to calc per wavelength later.
		#power *= opacity #rescaling because of atmosphere.
		power = power.si
		photons += (power * model[i][0] * u.angstrom / (c.h * c.c)) * area * efficiency *u.s
print("photons/s", photons)
Σλ = 0 # number of wavelength bands RV info is summed over

t = round(Teff/200)*200
signal_v2 = 0
for b in beatty:
	if ((b[1] == t) and (b[0] >= λ_min) and (b[0] <= λ_max)):
		Σλ += 1
		signal_v2 += b[2]**(2) # (km/s)**2
print("RMS RV info", np.sqrt(signal_v2))
print("Wavelength chunks", Σλ)
sigma_v = 5e-5 # target single measurement photon noise precision in km/s

# General Stellar error/precision effects:
# Metallicity effects on number of lines and their depth.
v_FeH = 10**(-0.27*FeH) # within 15% for near-solar (-2 to 0.5), biggest diff at [Fe/H] = -2

dTeff = Teff/5800 - 1

#log(g), or pressure broadening on linewidth (and velocity precision)
#m_opt = -0.27505*(1 - 1.22211*dTeff - 4.17622*dTeff**2) #4500-6500 A
#m_red = -0.33507*(1 - 1.41362*dTeff - 4.63727*dTeff**2) #6500-10000 A
#m_nir = -0.43926*(1 - 1.12505*dTeff - 4.53938*dTeff**2) #10000-25000 A
m = -0.27505*(1 - 1.22211*dTeff - 4.17622*dTeff**2) #optical, 4500-6500 A
v_logg = m*(logg-4.5)+1 # #f(log g), m varies with Teff and wavelength. Eqn only good for 4.0 to 5.0

#Effective temperature effects on number of lines and their depth.
#v_teff_opt = 1 + 2.04515*dTeff + 3.13362*dTeff**2 + 4.23845*dTeff**3
#v_teff_red = 1 + 2.18311*dTeff + 4.00361*dTeff**2 + 5.62077*dTeff**3
#v_teff_nir = 1 + 1.62418*dTeff + 2.62018*dTeff**2 + 5.01776*dTeff**3
v_Teff = 1 + 2.04515*dTeff + 3.13362*dTeff**2 + 4.23845*dTeff**3 #Optical

#theta0 also varies with wavelength choice
#theta0 is the inherent width of the Voigt profile
#Θ0_opt = 5.10521*(1-0.6395*dTeff) #4000 to 6500 A
#Θ0_red = 3.73956*(1-0.1449*dTeff) #6500 to 10000 A
#Θ0_nir = 6.42622*(1-0.2737*dTeff ) #10000 to 25000 A
theta_0 = 5.10521*(1-0.6395*dTeff) #optical, 4000 to 6500 A

theta_R = 299792.458/R #c/R in km/s

# Macroturbulence error
if (Teff > 5000) and (Teff < 7600):
	v_mac = 1.976 + 16.14*dTeff + 19.713*dTeff**2
elif (Teff < 5000):
	v_mac = 0.51
theta_mac = np.sqrt(2*np.log(2))*v_mac #Need to get macroturbulence velocity empirically.

# Plugging it all in to get SNR
SNR = (((0.5346*theta_0 + np.sqrt(0.2166*theta_0**2+theta_R**2+0.518*theta_rot**2+theta_mac**2))/theta_0)**1.5 * v_Teff * v_logg * v_FeH / sigma_v) * np.sqrt(299792.458 / signal_v2 / (n_pix * R))

# getting time from SNR
#time1 = (photons*gain + np.sqrt((photons*gain)**2 + 4*(photons/SNR)**2 * ((2*2.2*gain*read_noise)**2 + dark_current**2)))/(2*(photons/SNR)**2)
# 2.2 pixels per res element, factor of 2 vertical spread
#time2 = (photons*gain - np.sqrt((photons*gain)**2 + 4*(photons/SNR)**2 * ((2*2.2*gain*read_noise)**2 + dark_current**2)))/(2*(photons/SNR)**2)
A = (gain*photons/SNR)**2
B = gain*photons + n_pix*gain*dark_current
C = 2*2.2*n_pix*read_noise**2
time = (B + np.sqrt(B**2 + 4*A*C))/(2*A)
print("V_RMS (km/s), SNR, exptime")
print(sigma_v, SNR, time)
