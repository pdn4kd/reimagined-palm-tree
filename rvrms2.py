#! /usr/bin/python
import numpy as np
from astropy import units as u
from astropy import constants as c

'''
Python 2 compatable variant (no greek in variable names)

Calculates SNR and RV precision, given telescope and stellar parameters. Radial Velocity uncertainty based off of signal in a somewhat idealized spectrograph, and stellar type. Stellar spectra are BT-Settl, if in a somewhat non-standard format. Other equations and data (http://www.personal.psu.edu/tgb15/beattygaudi/table1.dat) from Beatty and Gaudi (2015). DOI: 10.1086/684264

Also uses http://www.aanda.org/articles/aa/full/2001/29/aa1316/aa1316.right.html for reference.

Current limitations: Detector properties not closely based on actual hardware, microturbulence fixed at ~1 km/s; granulation/starspots/other jitter not considered; fixed assumptions made of Teff (200 K steps) and spectrograph R values made based on wavelength range.
'''

# Beatty spectra (simulated): 100 A chunks; Solar metallicity;
# R = 100k (4000-6500 A), 165000 (6500-10000 A), 88000 (10000A - 25000 A)
# In general need to specify wavelength range, as corrections vary with the one chosen (optical, red, near IR).
# Table from http://www.personal.psu.edu/tgb15/beattygaudi/table1.dat
beatty = np.genfromtxt("table1.dat", dtype=None)
beatty.dtype.names = ('Angstrom', 'Teff', 'Uncertaintykms')

# Currently assuming a spectrometer entirely in the 'optical' range.
# Optical is basically Johnson V-band, so m = 0 means 3460 Jy
lam_min = 4000
lam_max = 6500
dellam = 100
BeattyWaves = np.arange(4050, 6550, dellam) #Optical range wavelength bits

Teff = 5800 #may need to interpolate between 2 later, as Beatty's table goes in 200 K chunks.
FeH = 0.0 #solar metallicity, can load from elsewhere
logg = 4.5 #solar surface gravity, can load from elsewhere
R = 100000 #actually varies with wavelength band, but...
vsini = 2.0 #solarish placeholder, will use actual values if possible
theta_rot = 1.13*vsini # Rough approximation, but rotational effects are near-linear, no matter limb-darkening.
# need stellar diameter/distance (or apparent magnitude in wavelength range)
# for correct brightness
rstar = 1.0 * u.solRad
dstar = 10.0 * u.pc

exptime = 300 * u.s
efficiency = 0.1 # total optical system efficiency, currently includes CCD quantum efficiency
telescope = 3.5 * u.m #diameter
area = np.pi * (telescope/2)**2

# Alpha Cen B test with HARPS
# How much do R and wavelength range matter?
exptime = 10 * u.minute
telescope = 3.566 * u.m
area = 8.8564 * u.m * u.m
efficiency = 0.02 #could be ~1-3%, depending on what is measured
lam_min = 3800
lam_max = 6800
BeattyWaves = np.arange(3850, 6850, dellam)
R = 110000 #110k or 120k, depending on source
Teff = 5200 # 5214 K
logg = 4.37
FeH = 0.23
rstar = 0.865 * u.solRad
dstar = 1.3 * u.pc
vsini = 2*np.pi*rstar/(38.7*u.day) * u.s/u.km
theta_rot = 1.13*vsini.si
print("vsini", vsini.si)
# Detectors: fast and slow modes


BTSettl = np.genfromtxt(str(Teff), dtype=float)
# BT Settl spectra are labled by Teff
# These spectra have a wavelength (Angstroms), and a flux (1e8 erg/s/cm^2/Angstrom) column
# sum of flux(lam)*dellam(lam)/1e8 == total flux (in power/area) emitted. 
# Multiply bt stellar_radius^2/distance^2 for recieved.
# Spectra downloaded from other sources will use different units!

#De-duping, this is not necessary if spectra are downloaded from another source.
model = [np.array([0,0])]
for x in BTSettl:
	if np.array_equal(x,model[-1]) == False:
		model.append(x)

#I_0 = np.zeros((len(BeattyWaves), 4)) #Wavelength bin, intensity at that velocity/wavelength element in terms of power and photons, overall SNR.
I_0 = np.zeros(len(BeattyWaves), dtype=[('wavelength','f4'),('intensity','f8'),('photons','f8'),('SNR','f8'),('real photons','f8')])
for lam in np.arange(0, len(BeattyWaves)): #need some C-style array traversal.
	for i in np.arange(0, len(model)-1):
		if ((model[i][0] >= BeattyWaves[lam]) and (model[i][0] < BeattyWaves[lam]+dellam) and (model[i][0] >= lam_min) and (model[i][0] <= lam_max)):
			# alternative rescaling, sum up power over eg: V-band -> PV
			# PV/(k*10^(-0.4mv)) = flux rescaling needed
			# k is relevant conversion factor because Vega mags are terrible
			power = (model[i][1]*1e-8*u.erg/u.cm**2/u.s/u.angstrom) * ((model[i+1][0]-model[i][0]) * u.angstrom)
			power *= rstar**2/dstar**2 #rescaling emitted spectrum based on stellar surface area and distance from us
			I_0[lam][1] += power.si * u.m**2/u.watt
			photons = power * model[i][0] * u.angstrom / (c.h * c.c)
			I_0[lam][2] += photons * area * efficiency * exptime
			I_0[lam][0] = BeattyWaves[lam]

gain = 1.0 # electrons generated per photon
read_noise = 1.0 #RMS of +- spurious electrons
dark_current = 4 / u.hour #electrons per time
#Resolution elements: approximately 5 pix wide by 5 pix long for first guess -- 25 in total
n_pix = 25

# HARPS
#Jasmin/slow: 0.63 e/ADU, 2.87e- RON
#Linda/slow: N/A?
#Jasmin/fast: 1.42 e-/ADU, 7.07e- RON
#Linda/fast: 1.4 e-/ADU, 6.08 e- RON
#Normal/slow: 2.00 ADU/e-, 2.9e- RON
#Normal/fast: 0.75 ADU/e-, 4.9e- RON
#HARPS
gain = 1/1.42 # ADU/e-, assume 1:1 photon to electron generation
read_noise = 7.07 #RMS of +- spurious electrons
dark_current = 4 / u.hour
n_pix = 4

# A resolution element is not 100 A long!
# R = l/dl, so 1 resolution element at wavelength lambda is lambda/R wide
# Or, a bin contains 100 A/(lambda/R) = 100A*R/lambda resolution elements
for lam in np.arange(0, len(BeattyWaves)):
	pix = n_pix*dellam*R/I_0[lam][0]
	# For now, assume equal signal per pixel. A gaussian would be better.
	I_0[lam][3] = I_0[lam][2]/pix*gain / np.sqrt(I_0[lam][2]/pix*gain + (gain*read_noise*2.2*2)**2 + (gain*dark_current*exptime)**2)
	I_0[lam][4] = I_0[lam][3]**2 * pix / ((c.c/(u.km/u.s)).si * dellam/I_0[lam][0])

print(I_0)
#print(sum(I_0[:,1]), "W/m^2")
print(sum(I_0['intensity']), "W/m^2")
I_SNR = sum(I_0['photons'])
print("Total SNR", I_SNR/pix*gain / np.sqrt(I_SNR/pix*gain + (gain*read_noise*2.2*2)**2 + (dark_current*exptime)**2), "Average SNR", np.mean(I_0['SNR']))


Q = 0 # "Quality factor" -- summation of intensity/signal over default velocity uncertainty in each bin. Ignores error sources that are considered later.
for b in beatty: #each line of b is a tuple with wavelength, teff, uncertainty
#if b[0] is in wavelength range, and b[1] is in Teff range, use the uncertainty
	if ((b[0] >= lam_min) and (b[0] <= lam_max) and b[1] == Teff):
		for i in I_0:
			if i[0] == b[0]:
				#print(b[0], b[2], 1/np.sqrt(i[4]/b[2]**2))
				Q += i[4]/b[2]**2 # Summation to get RMS velocity error over the wavelength range, given that it is known in each bin.
Q = 1/np.sqrt(Q)#Quality factor from summing up weights, ignoring other noise sources
print("Noise/Feh/logg-Free V_RMS (km/s):", Q)

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
#theta0_opt = 5.10521*(1-0.6395*dTeff) #4000 to 6500 A
#theta0_red = 3.73956*(1-0.1449*dTeff) #6500 to 10000 A
#theta0_nir = 6.42622*(1-0.2737*dTeff ) #10000 to 25000 A
theta_0 = 5.10521*(1-0.6395*dTeff) #optical, 4000 to 6500 A

theta_R = 299792.458/R #c/R in km/s

# Macroturbulence error
if (Teff > 5000) and (Teff < 7600):
	v_mac = 1.976 + 16.14*dTeff + 19.713*dTeff**2
elif (Teff < 5000):
	v_mac = 0.51
theta_mac = np.sqrt(2*np.log(2))*v_mac #Need to get macroturbulence velocity empirically.

# "Final" value.
sigma_v = Q * ((0.5346*theta_0 + np.sqrt(0.2166*theta_0**2+theta_R**2+0.518*theta_rot**2+theta_mac**2))/theta_0)**1.5 * v_Teff * v_logg * v_FeH


#debugging
print("theta_0, theta_R, theta_rot, theta_mac")
print(theta_0, theta_R, theta_rot, theta_mac)
print("v_Teff, v_logg, v_FeH")
print(v_Teff, v_logg, v_FeH)
print("Stellar Noise V_RMS", sigma_v)

#print("No thetas:", Q*v_Teff*v_logg*v_FeH)
#print("theta_0", Q*v_Teff*v_logg*v_FeH*((0.5346*theta_0+np.sqrt(0.2166*theta_0**2))/theta_0)**1.5)
v_logg = 1
v_FeH = 1
v_Teff = 1
print(Q * ((0.5346*theta_0 + np.sqrt(0.2166*theta_0**2+theta_R**2+0.518*theta_rot**2+theta_mac**2))/theta_0)**1.5 * v_Teff * v_logg * v_FeH)