#! /usr/bin/python
import numpy as np
from astropy import units as u
from astropy import constants as c
import rvrms

'''
Exposure time calculator for a list of target stars, desired measurement quality, and optical system specifications. Currently a proof of concept derived from HARPS, but capable to some degree of considering stellar properties.

Makes simplifications, and is currently only valid in optical (-), but can be expanded to red (-) and near IR (-).

Guesses at exposure time based on saturating the detector, then scales with t^-0.5 to get the precision. Also considers readout time.
'''

# Beatty spectra (simulated): 100 A chunks; Solar metallicity;
# R = 100k (4000-6500 A), 165000 (6500-10000 A), 88000 (10000A - 25000 A)
# In general need to specify wavelength range, as corrections vary with the one chosen (optical, red, near IR).
# Table from http://www.personal.psu.edu/tgb15/beattygaudi/table1.dat
beatty = np.genfromtxt("table1.dat", dtype=None)
beatty.dtype.names = ('Angstrom', 'Teff', 'Uncertaintykms')


def λ_peak(Teff, λ_min=0 * u.angstrom, λ_max=1e13 * u.angstrom):
	λ = (c.b_wien / (Teff * u.K)).si # Star is assumed to be a perfect blackbody
	# Correcting for if peak wavelength is outside of detector, angstroms assumed
	if (λ > λ_max):
		λ = λ_max
	if (λ < λ_min):
		λ = λ_min
	return λ

def extinction(λ):
	# best guess, based on measurements at Texas A&M Univeristy, and CFHT in Mauka Kea
	# http://www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/extinction
	# http://adsabs.harvard.edu/abs/1994IAPPP..57...12S
	# Values should not be considered trustworthy outside of ~3000-10000 Angstroms, and really 3500-9500 at that.
	return (0.09 + (3080.0/λ)**4)

def airmass(zenith_angle, site_elevation=0, scale_height=8400):
	# number is scaled to sea-level, not site
	return np.exp(-site_elevation/scale_height)/np.cos(zenith_angle)

def time_guess(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_peak):
	'''Generating an initial guess on exposure time. We assume a saturated detector is appropriate, thus putting one in the photon noise limited regime.'''
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
	
	#adding up photons in 100 A chunk centered on peak wavelength:
	for i in np.arange(0, len(model)-1):
		if((model[i][0]>= λ_peak/u.angstrom-50) and (model[i][0] <= λ_peak/u.angstrom+15)):
			power = (model[i][1]*1e-8*u.erg/u.cm**2/u.s/u.angstrom) * ((model[i+1][0]-model[i][0]) * u.angstrom)
			power *= rstar**2/dstar**2 #rescaling emitted spectrum based on stellar surface area and distance from us
			opacity = np.exp(-extinction(model[i][0]) * atmo)
			power *= opacity #rescaling because of atmospheric scattering/absorption.
			power = power.si
			photons += (power * model[i][0] * u.angstrom / (c.h * c.c)) * area * efficiency
	photons /= (100*u.angstrom / (λ_peak/R)) #100 Angstrom range -> per resolution element
	photons /= n_pix #per resolution element -> per pixel
	#print("photons/s/pixel", photons)
	time_guess = well_depth / photons
	#print("Estimated exposure time", time_guess, time_guess.si)
	return time_guess.si

def time_actual(sigma_v, Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, exptime, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_min, λ_max, Δλ):
	rv_actual = rvrms.rvcalc(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, exptime, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_min, λ_max, Δλ)
	time_actual = exptime * (sigma_v/rv_actual)**2
	time_actual += read_time * np.ceil(time_actual/exptime)
	return time_actual

if __name__ == '__main__':
	# HARPS
	telescope = 3.566 * u.m
	#area = np.pi * (telescope/2)**2 # telescope area, central obstruction ignored
	area = 8.8564 * u.m * u.m
	efficiency = 0.02 #could be ~1-3%, depending on what is measured
	λ_min = 3800 * u.angstrom
	λ_max = 6800 * u.angstrom
	Δλ = 100 * u.angstrom
	#BeattyWaves = np.arange((λ_min+Δλ/2)/u.angstrom, (λ_max+Δλ/2)/u.angstrom, Δλ/u.angstrom)
	R = 110000 #110k or 120k, depending on source
	gain = 1/1.42 # ADU/e-, assume 1:1 photon to electron generation
	read_noise = 7.07 #RMS of ± spurious electrons
	dark_current = 4 / u.hour
	n_pix = 4
	well_depth = 30000 #electrons (or photons)
	read_time = 10 * u.s
	
	# Alpha Cen B
	Teff = 5214 # 5214 K
	logg = 4.37
	FeH = 0.23
	rstar = 0.865 * u.solRad
	dstar = 1.3 * u.pc
	vsini = 2*np.pi*rstar/(38.7*u.day) * u.s/u.km
	theta_rot = 1.13*vsini.si
	
	sigma_v = 5e-5 # target single measurement photon noise precision in km/s
	
	# Atmo
	zenith_angle = 0 # Can be read in from stellar observations!
	scale_height = 8400 * u.m #Can make arguments for `8500 m to ~7500 m, but the higher temperatures typical of the lower atmosphere are more relevant here.
	site_elevation = 2016 * u.m #Kitt Peak/WIYN
	atmo = airmass(zenith_angle, site_elevation, scale_height)

	λ_peak = λ_peak(Teff, λ_min, λ_max)
	#BeattyWaves = np.arange(λ_min/u.angstrom, λ_max/u.angstrom, Δλ/u.angstrom)
	guess = time_guess(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_peak)
	actual = time_actual(sigma_v, Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, guess, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_min, λ_max, Δλ)
	print("guess, actual times")
	print(guess, actual)
