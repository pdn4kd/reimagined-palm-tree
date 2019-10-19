'''
Using functionality of rvrms.py to get a bunch of RV precisions from existing observation logs. Also throws in some Gaussian noise.
'''
import numpy as np
import rvrms
from astropy import units as u

#load up full simulation bits
import simulation
sim = simulation.simulation('simulation.ini') #metadata is assumed not to change
import datetime
now = str((datetime.datetime.now()).isoformat())

#opening a folder with results from a previous run
import glob
import os
simnumber = input('Enter 5 digit sim number (including any leading zeros): ')
simpath = glob.glob('./results/*.'+simnumber+'/')[0]
os.stat(simpath)
print(simpath)

Δλ = 100 * u.angstrom
λ_max = sim.instruments[0].λmax * u.angstrom
λ_min = sim.instruments[0].λmin * u.angstrom
area = sim.telescopes[0].area * u.m * u.m
dark_current = sim.instruments[0].dark_current / u.hour
instrms = sim.instruments[0].general_noise

#target_list = np.genfromtxt("targetstars.csv", delimiter=",", dtype=("U11", "U10", "U9", float, float, float, float, int, float, float, float), names=True) #bad, should have it in observation run metadata
target_list = np.genfromtxt("targetstars.csv", delimiter=",", dtype=(float, float, "U15", float, float, float, float, float, float, float, float, float, "U30"), names=True)
#target_list = np.genfromtxt("targetstars.csv", delimiter=",", names=True)
for target in target_list:
	if (os.access(simpath+target['Name']+".txt", os.R_OK)):#with current bad setup, target listiing might not include all stars
		star = np.genfromtxt(simpath+target['Name']+".txt", delimiter=",", names=True)
		if (star.shape != ()):
			# if we have actual observations, not just the test/setup one, we can calculate RVs
			print(target['Name'])
			star_rms = open(simpath+target['Name']+"_rv.txt", 'w')
			star_rms.write("obs_start,obs_end,duration,altitude,azimuth,exposures,photonprec,instprec,rvprec,rvmeas\n")
			Teff = target['K']
			FeH = target['Sun']
			logg = target['cms']
			vsini = target['kms'] * u.km / u.s
			theta_rot = 1.13 * target['kms']
			rstar = target['solRad'] * u.solRad
			dstar = target['pc'] * u.pc
			vmac = target['Vmac']
			for i in star:
				exptime = i['duration'] * u.minute
				zenith_angle = (90-i['altitude'])*np.pi/180.0
				airmass = np.exp(-sim.elevation/8400)/np.cos(zenith_angle)
				n_expose = i['exposures']
				photonrms = 1000*rvrms.rvcalc(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, vmac, airmass, exptime, sim.instruments[0].efficiency, area, sim.instruments[0].R, sim.instruments[0].gain, sim.instruments[0].read_noise, dark_current, sim.instruments[0].n_pix, λ_min, λ_max, Δλ, n_expose)
				vrms = np.sqrt(photonrms**2+instrms**2)
				vmeas = np.random.normal(0.0, vrms)
				line = str(i['obs_start'])+","+str(i['obs_end'])+","+str(i['duration'])+","+str(i['altitude'])+","+str(i['azimuth'])+","+str(i['exposures'])+","+str(photonrms)+","+str(instrms)+","+str(vrms)+","+str(vmeas)+"\n"
				star_rms.write(line)
				print(vrms, vmeas)
			star_rms.close()
