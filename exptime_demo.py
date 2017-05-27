#! /usr/bin/python
import math
from astropy import units as u

'''
Exposure time calculator for a list of target stars, desired measurement quality, and optical system specifications. Currently a proof of concept derived from ESPRESSO performance on fixed stellar types.
'''

mags = [0,1,2,3,4,5,6,7,8,9,10,11,12]

telescope = 3.5 * u.meter # telescope diameter
area = math.pi * (telescope/2)**2 # telescope area, central obstruction ignored
λ = 5450 * u.angstrom #central wavelength of spectrograph region
R = 100000 # spectrograph resolution
gain = 4 # electrons generated per photon
read_noise = 5 #RMS of +- spurious electrons
dark_noise = 0
efficiency = 0.1 # total optical system efficiency, currently includes CCD quantum efficiency
f0 = 2.518021002e-8 * u.watt / u.meter / u.meter # Approximate flux for an apparent bolometric magnitude 0 star, defined in IAU 2015 meeting
from astropy.constants import c
from astropy.constants import h
SNR = 600 #arbitrary target number from ESPRESSO, real ones will vary with stars/instruments
time_min = 300.0 #seconds
time_max = 3600.0 #seconds

print("Mag\t|\tTime to SNR "+str(SNR)+" (s)")
print("--------+---------------------------")
for mag in mags:
	photons = f0 * 10**(-0.4*mag) * λ.si/R * area * efficiency / h / c * u.s #seconds here to fix units later
	# tends to ~1267 photons/s/cm^2/angstrom because bolometric, need to assume more like 995-1005 for Johnson V-band
	time1 = (photons*gain + math.sqrt((photons*gain)**2 + 4*(photons/SNR)**2 * ((2*2.2*gain*read_noise)**2 + dark_noise**2)))/(2*(photons/SNR)**2)
# 2.2 pixels per res element, factor of 2 vertical spread
	#time2 = (photons*gain - math.sqrt((photons*gain)**2 + 4*(photons/SNR)**2 * ((2*gain*read_noise)**2 + dark_noise**2)))/(2*(photons/SNR)**2)
	print(mag, '\t|\t', time1)
#Need to do something if t < 300 s, or t > 3600 s, capping obs times and noting SNR
	if (time1 < time_min):
		time1 = time_min
		SNR1 = photons*time1*gain / math.sqrt(photons*time1*gain + (gain*read_noise*2.2*2)**2 + dark_noise**2)
		print("Below min time("+str(time_min)+"), actual SNR:", SNR1.value)
	elif (time1 > time_max):
		time1 = time_max
		SNR1 = photons*time1*gain / math.sqrt(photons*time1*gain + (gain*read_noise*2.2*2)**2 + dark_noise**2)
		print("Above max time("+str(time_max)+"), actual SNR:", SNR1.value)
#ESPRESSO says precision varies with 1/SNR, and SNR 600 yields 0.1 m/s
#ie: 600*0.1/SNR = RV precision in m/s
#realistic figures will be done elsewhere.
#Additional systematic errors (jitter?) are added in quadrature to these RV precisions
#multi-epoch observations beat down the noise with precision = single_epoch/sqrt(N_obs)
#Unsure about systematics for multi-epoch. Add afterwards?
