#! /usr/bin/python
import numpy as np
from astropy import units as u
from astropy import constants as c

'''
Calculates SNR and RV precision, given telescope and stellar parameters. Radial Velocity uncertainty based off of signal in a somewhat idealized spectrograph, and stellar type. Stellar spectra are BT-Settl, if in a somewhat non-standard format. Other equations and data (http://www.personal.psu.edu/tgb15/beattygaudi/table1.dat) from Beatty and Gaudi (2015). DOI: 10.1086/684264

Also uses http://www.aanda.org/articles/aa/full/2001/29/aa1316/aa1316.right.html for reference.

Current limitations: Detector properties not closely based on actual hardware, microturbulence fixed at ~1 km/s; granulation/starspots/other jitter not considered; fixed assumptions made of Teff (200 K steps) and spectrograph R values made based on wavelength range.
'''
'''
starlist = np.genfromtxt("stars.csv", delimiter=",", names=True)
	τ = 0.11 # Should vary per wavelength bin.
	for star in starlist:
	starobs = np.genfromtxt(star['pl_hostname'], delimiter=",", names=True)
		for obs in starobs:
			if (obs['altitude'] > 0):
				# get v_rms for every observation
				zenith_angle = (90.0-altitude)*np.pi/180
				opacity = np.exp(-τ/np.cos(zenith_angle))
				rvcalc(star['Teff'], star['FeH'], star['logg'], star['vsini'], star['theta_rot'], star['rstar'], star['dstar'], opacity, obs['duration'], star['efficiency'], star['area'], star['R'], star['gain'], star['read_noise'], star['dark_current'], star['n_pix']):
'''


def extinction(λ):
	# best guess, based on measurements at Texas A&M Univeristy, and CFHT in Mauka Kea
	# http://www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/extinction
	# http://adsabs.harvard.edu/abs/1994IAPPP..57...12S
	# Values should not be considered trustworthy outside of ~3000-10000 Angstroms, and really 3500-9500 at that.
	return (0.09 + (3080.0/λ)**4)

def rvcalc(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, airmass, exptime, efficiency, area, R, gain, read_noise, dark_current, n_pix, λmin, λmax, dλ):
	# Need to fix to make less unit dependent
	λ_max = λmax.to(u.angstrom).value
	λ_min = λmin.to(u.angstrom).value
	Δλ = dλ.to(u.angstrom).value
	
	BeattyWaves = np.arange((λ_min+Δλ/2), (λ_max+Δλ/2), Δλ)
	
	BTSettl = np.genfromtxt(str(int(round(Teff,-2))), dtype=float)
	# BT Settl spectra are labled by Teff, and available every 100 K.
	# These spectra have a wavelength (Angstroms), and a flux (1e8 erg/s/cm^2/Angstrom) column
	# sum of flux(λ)*Δλ(λ)/1e8 == total flux (in power/area) emitted. 
	# Multiply bt stellar_radius^2/distance^2 for recieved.
	# Spectra downloaded from other sources will use different units!
	#
	# Note that Beatty RV information is only available every 200 K.
	
	#De-duping, this is not necessary if spectra are downloaded from another source.
	model = [np.array([0,0])]
	for x in BTSettl:
		if np.array_equal(x,model[-1]) == False:
			model.append(x)
	
	#I_0 = np.zeros((len(BeattyWaves), 4)) #Wavelength bin, intensity at that velocity/wavelength element in terms of power and photons, overall SNR.
	I_0 = np.zeros(len(BeattyWaves), dtype=[('wavelength','f4'),('intensity','f8'),('photons','f8'),('SNR','f8'),('real photons','f8')])
	for λ in np.arange(0, len(BeattyWaves)): #need some C-style array traversal.
		for i in np.arange(0, len(model)-1):
			if ((model[i][0] >= BeattyWaves[λ]) and (model[i][0] < BeattyWaves[λ]+Δλ) and (model[i][0] >= λ_min) and (model[i][0] <= λ_max)):
				power = (model[i][1]*1e-8*u.erg/u.cm**2/u.s/u.angstrom) * ((model[i+1][0]-model[i][0]) * u.angstrom)
				power *= rstar**2/dstar**2 #rescaling emitted spectrum based on stellar surface area and distance from us
				τ = extinction(model[i][0])
				opacity = np.exp(-τ*airmass)
				power *= opacity #rescaling because of atmosphere.
				I_0[λ][1] += power.si * u.m**2/u.watt
				photons = power * model[i][0] * u.angstrom / (c.h * c.c)
				I_0[λ][2] += photons * area * efficiency * exptime
				I_0[λ][0] = BeattyWaves[λ]
	
	
	# A resolution element is not 100 A long!
	# R = l/dl, so 1 resolution element at wavelength lambda is lambda/R wide
	# Or, a bin contains 100 A/(lambda/R) = 100A*R/lambda resolution elements
	for λ in np.arange(0, len(BeattyWaves)):
		pix = n_pix*Δλ*R/I_0[λ][0]
		# For now, assume equal signal per pixel. A gaussian would be better.
		I_0[λ][3] = I_0[λ][2]/pix*gain / np.sqrt(I_0[λ][2]/pix*gain + (gain*read_noise*2.2*2)**2 + (gain*dark_current*exptime)**2)
		I_0[λ][4] = I_0[λ][3]**2 * pix / ((c.c/(u.km/u.s)).si * Δλ/I_0[λ][0])
	
	#print(I_0)
	#print(sum(I_0[:,1]), "W/m^2")
	#print(sum(I_0['intensity']), "W/m^2")
	I_SNR = sum(I_0['photons'])
	#print("Total SNR", I_SNR/pix*gain / np.sqrt(I_SNR/pix*gain + (gain*read_noise*2.2*2)**2 + (dark_current*exptime)**2), "Average SNR", np.mean(I_0['SNR']))
	
	
	Q = 0 # "Quality factor" -- summation of intensity/signal over default velocity uncertainty in each bin. Ignores error sources that are considered later.
	# Should this table really be opened in here and not passed from elsewhere/loaded at startup?
	beatty = np.genfromtxt("table1.dat", dtype=None)
	beatty.dtype.names = ('Angstrom', 'Teff', 'Uncertaintykms')
	for b in beatty: #each line of b is a tuple with wavelength, teff, uncertainty
	#Teff is rounded to nearest 200 because of limitations in Beatty RVs
	#if b[0] is in wavelength range, and b[1] is in Teff range, use the uncertainty
		if ((b[0] >= λ_min) and (b[0] <= λ_max) and b[1] == round(Teff/200)*200):
			for i in I_0:
				if i[0] == b[0]:
					#print(b[0], b[2], 1/np.sqrt(i[4]/b[2]**2))
					Q += i[4]/b[2]**2 # Summation to get RMS velocity error over the wavelength range, given that it is known in each bin.
	Q = 1/np.sqrt(Q)#Quality factor from summing up weights, ignoring other noise sources
	#print("Noise/Feh/logg-Free V_RMS (km/s):", Q)
	
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
	
	# Macroturbulence error estimate, only really applies for log(g) > 4.0
	# Teff between 5000 K and 6500 K most accurate, but up to 7600 K okay.
	# Empirical numbers are better if available.
	if (Teff > 5000) and (Teff < 7600):
		v_mac = 1.976 + 16.14*dTeff + 19.713*dTeff**2
	elif (Teff <= 5000): # exactly 5000 K not specified by Beatty, but this is more conservative.
		v_mac = 0.51
	theta_mac = np.sqrt(2*np.log(2))*v_mac
	
	# "Final" value.
	sigma_v = Q * ((0.5346*theta_0 + np.sqrt(0.2166*theta_0**2+theta_R**2+0.518*theta_rot**2+theta_mac**2))/theta_0)**1.5 * v_Teff * v_logg * v_FeH
	#print("V_RMS", sigma_v)
	return sigma_v

def expest():
	'''quick exposure time estimator, given target SNR and detector properties'''

if __name__ == '__main__':
	'''
	starlist = np.genfromtxt("stars.csv", delimiter=",", names=True)
		for star in starlist:
		starobs = np.genfromtxt(star['pl_hostname'], delimiter=",", names=True)
			for obs in starobs:
		# get v_rms for every observation
	'''
	site_elevation = 2016.0 # m

	# Beatty spectra (simulated): 100 A chunks; Solar metallicity;
	# R = 100k (4000-6500 A), 165000 (6500-10000 A), 88000 (10000A - 25000 A)
	# In general need to specify wavelength range, as corrections vary with the one chosen (optical, red, near IR).
	# Table from http://www.personal.psu.edu/tgb15/beattygaudi/table1.dat
	beatty = np.genfromtxt("table1.dat", dtype=None)
	beatty.dtype.names = ('Angstrom', 'Teff', 'Uncertaintykms')
	
	# Currently assuming a spectrometer entirely in the 'optical' range.
	# Optical is basically Johnson V-band, so m = 0 means 3460 Jy
	λ_min = 3800 * u.angstrom
	λ_max = 6800 * u.angstrom
	Δλ = 100 * u.angstrom
	BeattyWaves = np.arange(3850, 6850, 100) #Optical range wavelength bits
	
	#atmospheric opacity is, uh, problematic. An innacurate, but better than nothing assumption of optical depth being ~0.11/airmass is used.
	# A better assuption uses a baseline (0.09/airmass), along with rayleigh scattering factor, calibrated to be 1 at 6080 Angstroms. This has some basis in data.
	# Airmass amount uses sec(zenith_angle) because this is accurate over the ranges that telescopes actually point. (up to 60 degrees off zenith)
	zenith_angle = 0 # Can be read in from stellar observations!
	#τ = 0.11 # Should vary per wavelength bin.
	#opacity = np.exp(-τ/np.cos(zenith_angle))
	airmass = np.exp(-site_elevation/8400)/np.cos(zenith_angle)
	
	# Alpha Cen B
	Teff = 5214 # 5214 K. Note that there is rounding to the nearest 100 or 200 K in places.
	# Sun is 5800 K
	logg = 4.37 # Surface gravity, sun is 4.5
	FeH = 0.23 # Metallicity, sun is 0
	# need stellar diameter/distance (or apparent magnitude in wavelength range)
	# for correct brightness
	rstar = 0.865 * u.solRad
	dstar = 1.3 * u.pc
	vsini = 2*np.pi*rstar/(38.7*u.day) * u.s/u.km # Sun is ~2 km/s
	theta_rot = 1.13*vsini.si # Rough approximation, but rotational effects are near-linear, no matter limb-darkening.
	
	# HARPS instrument on ??? telescope (nominally in La Silla)
	telescope = 3.566 * u.m #diameter
	#area = np.pi * (telescope/2)**2 # telescope area, central obstruction ignored
	area = 8.8564 * u.m * u.m
	efficiency = 0.02 #total optical system efficiency. Could be ~1-3% here, depending on what is measured
	λ_min = 3800 * u.angstrom
	λ_max = 6800 * u.angstrom
	Δλ = 100 * u.angstrom
	#BeattyWaves = np.arange((λ_min+Δλ/2)/u.angstrom, (λ_max+Δλ/2)/u.angstrom, Δλ/u.angstrom)
	R = 110000 #110k or 120k, depending on source
	gain = 1/1.42 # ADU/e-, assume 1:1 photon to electron generation
	read_noise = 7.07 #RMS of ± spurious electrons
	dark_current = 4 / u.hour
	n_pix = 4 # resolution element total (length x width) pixel count
	well_depth = 30000 #electrons (or photons)
	read_time = 10 * u.s
	exptime = 76 * u.s

	print("Radial Velocity test with HARPS on Alpha Cen B. Expected precision 5 cm/s (5e-5 km/s):")
	print(rvcalc(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, airmass, exptime, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_min, λ_max, Δλ))
