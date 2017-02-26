#! /usr/bin/python
import numpy as np
from astropy import units as u
from astropy import constants as c

'''
Radial Velocity uncertainty based off of signal in a somewhat idealized spectrograph, and stellar type. Equations and data from Beatty and Gaudi (2015). DOI: 10.1086/684264

Also uses http://www.aanda.org/articles/aa/full/2001/29/aa1316/aa1316.right.html for reference.

Current limitations: Does not consider actual detector SNR, or vsini (mostly).
'''

# Beatty spectra (simulated): 100 A chunks; Solar metallicity;
# R = 100k (4000-6500 A), 165000 (6500-10000 A), 88000 (10000A - 25000 A)
# In general need to specify wavelength range, as corrections vary with the one chosen (optical, red, near IR).
# All velocities are in km/s
# Table from http://www.personal.psu.edu/tgb15/beattygaudi/table1.dat but reformatted as a CSV.
beatty = np.genfromtxt("table1.csv", dtype=None, delimiter=",", names=True)

# Currently assuming a spectrometer entirely in the 'optical' range.
λ_min = 4000
λ_max = 6500
Δλ = 100
BeattyWaves = np.arange(4050, 6550, Δλ) #Optical range wavelength bits

Teff = 2600 #may need to interpolate between 2 later, as Beatty's table goes in 200 K chunks.
FeH = 0.0 #solar metallicity, can load from elsewhere
logg = 4.5 #solar surface gravity, can load from elsewhere
R = 100000 #actually varies with wavelength band, but...
vsini = 2.0 #solarish placeholder, will use actual values if possible
theta_rot = 1.13*vsini # Rough approximation, but rotational effects are near-linear, no matter limb-darkening.
# need stellar diameter/distance (or apparent magnitude) correction
mag = 12 #Are BT-Settl spectra for mag 0 stars?)

Q = 0 #Quality factor from summing up weights, ignoring SNR

exptime = 300 * u.s #5 minutes
efficiency = 0.1 # total optical system efficiency, currently includes CCD quantum efficiency
telescope = 3.5 * u.m #diameter
area = np.pi * (telescope/2)**2

BTSettl = np.genfromtxt(str(Teff), dtype=float)
# BT Settl spectra are labled by Teff
# These spectra have a wavelength (Angstroms), and a flux (1e8 erg/s/cm^2/Angstrom) column
# sum of flux(λ)*Δλ(λ)/1e8 == total flux (in power/area) emitted. Need radius^2/distance^2 for recieved.

#De-duping
model = [np.array([0,0])]
for x in BTSettl:
	if np.array_equal(x,model[-1]) == False:
		model.append(x)

I_0 = np.zeros((len(BeattyWaves), 2)) #Wavelength bin, and intensity at that velocity/wavelength element in terms of photons.
for λ in np.arange(0, len(BeattyWaves)): #need some C-style array traversal.
	for i in np.arange(0, len(model)-1):
		if ((model[i][0] >= BeattyWaves[λ]) and (model[i][0] < BeattyWaves[λ]+Δλ) and (model[i][0] >= λ_min) and (model[i][0] <= λ_max)):
			# fix: multiply by radius^2/dist^2 for brightness mod
			# alternative rescaling, sum up power over eg: V-band -> PV
			# PV/(k*10^(-0.4mv) = flux rescaling needed
			# k is relevant conversion factor because Vega mags are terrible
			power = (model[i][1]*1e-8*u.erg/u.cm**2/u.s/u.angstrom) * ((model[i+1][0]-model[i][0]) * u.angstrom)
			photons = power * model[i][0] * u.angstrom / (c.h * c.c)
			I_0[λ][1] += photons * area * efficiency * exptime
			I_0[λ][0] = BeattyWaves[λ]
print(I_0)

#need resolution element to pixel conversion
'''
# per pixel
gain = 4 # electrons generated per photon
read_noise = 5 #RMS of +- spurious electrons
dark_noise = 0
efficiency = 0.1 # total optical system efficiency, currently includes CCD quantum efficiency
n_pix = 5 # width of spectrum in pixels. Need length per 100 A chunk.
SNR1 = photons*time1*gain / math.sqrt(photons*time1*gain + (gain*read_noise*2.2*2)**2 + dark_noise**2)
# For a given element i, I_0 == SNR_pix^2 * n_pix / dV_chunk, need per-pixel SNR, number of pixels, and velocity span of the 100 Angstrom velocity chunk
# dV_chunk == dlambda*c/R? (ie: ~3e-8 km/s)
# SNR_pix should be SNR1
# n_pix - need width (5?) and length per 100 A chunk. Presumably 1e-3 m.
#dV_chunk = 
#I_1 = SNR1**2 * n_pix/dV_chunk
'''

for b in beatty: #each line of b is a tuple with wavelength, teff, uncertainty
#if b[0] is in wavelength range, and b[1] is in Teff range, use the uncertainty
	if ((b[0] >= λ_min) and (b[0] <= λ_max) and b[1] == Teff):
		for i in I_0:
			if i[0] == b[0]:
				Q += i[1]/b[2]**2 # Summation to get RMS velocity error over the wavelength range, given that it is known in each bin.
		# need to consider SNR in each 100 A bin, especially per pixel. This is currently 1 photon per velocity element.
		# 
Q = 1/np.sqrt(Q)
print("QualityV_RMS (km/s):", Q)

# Metallicity effects on number of lines and their depth.
v_FeH = 10**(-0.27*FeH) # within 15% for near-solar (-2 to 0.5), biggest diff at [Fe/H] = -2

dTeff = Teff/5800 - 1

#v_logg, or f(log g)
#m_opt = -0.27505*(1 - 1.22211*dTeff - 4.17622*dTeff**2) #4500-6500 A
#m_red = -0.33507*(1 - 1.41362*dTeff-4.63727*dTeff**2) #6500-10000 A
#m_nir = -0.43926*(1 - 1.12505*dTeff - 4.53938*dTeff**2) #10000-25000 A
m = -0.27505*(1 - 1.22211*dTeff - 4.17622*dTeff**2) #optical, 4500-6500 A
v_logg = m*(logg-4.5)+1 # #f(log g), m varies with Teff and wavelength. Eqn only good for 4.0 to 5.0

#v_teff, or f(Teff), effective temperature effects on number of lines and their depth.
#v_teff_opt = 1 + 2.04515*dTeff+ 3.13362*dTeff**2+4.23845*dTeff**3
#v_teff_red = 1 + 2.18311*dTeff+ 4.00361*dTeff**2+5.62077*dTeff**3
#v_teff_nir = 1 + 1.62418*dTeff+ 2.62018*dTeff**2+5.01776*dTeff**3
v_Teff = 1 + 2.04515*dTeff+ 3.13362*dTeff**2+4.23845*dTeff**3 #Optical

#theta0 also varies with wavelength choice
#Θ0_opt = 5.10521*(1-0.6395*dTeff) #4000 to 6500 A
#Θ0_red = 3.73956*(1-0.1449*dTeff) #6500 to 10000 A
#Θ0_nir = 6.42622*(1-0.2737d*Teff ) #10000 to 25000 A
theta_0 = 5.10521*(1-0.6395*dTeff) #optical, 4000 to 6500 A

theta_R = 299792.458/R #c/R in km/s

# Macroturbulence error
if (Teff > 5000) and (Teff < 7600):
	v_mac = 1.976 + 16.14*dTeff + 19.713*dTeff**2
elif (Teff < 5000):
	v_mac = 0.51
theta_mac = np.sqrt(2*np.log(2))*v_mac #Need to get macroturbulence velocity emperically.

sigma_v = Q * ((0.5346*theta_0 + np.sqrt(0.2166*theta_0**2+theta_R**2+0.518*theta_rot**2+theta_mac**2))/theta_0)**1.5 * v_Teff * v_logg * v_FeH

#debugging
#print("theta_0, theta_R, theta_rot, theta_mac")
#print(theta_0, theta_R, theta_rot, theta_mac)
#print("v_Teff, v_logg, v_FeH")
#print(v_Teff, v_logg, v_FeH)
print("sigma_v", sigma_v)
