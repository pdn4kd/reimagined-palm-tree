'''
Using functionality of rvrms.py to get a bunch of RV precisions from existing observation logs. Also throws in some Gaussian noise.
'''
import numpy as np
import rvrms
from astropy import units as u

#load up full simulation bits, will drop later
import simulation
sim = simulation.simulation('simulation.ini') #not necessary if metadata is included
import datetime
now = str((datetime.datetime.now()).isoformat())

#opening a folder with results from a previous run
#borrowed from vis.py
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
			star_rms.write("obs_start,obs_end,duration,altitude,azimuth,photonprec,instprec,rvprec,rvmeas\n")
			Teff = target['K']
			FeH = target['Sun']
			logg = target['cms']
			vsini = target['kms'] * u.km / u.s
			theta_rot = 1.13 * target['kms']
			rstar = target['solRad'] * u.solRad
			dstar = target['pc'] * u.pc
			vmac = target['Vmac']
			for i in star:
				exptime = i['duration'] * u.s
				zenith_angle = (90-i['altitude'])*np.pi/180.0
				airmass = np.exp(-sim.elevation/8400)/np.cos(zenith_angle)
				photonrms = 1000*rvrms.rvcalc(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, vmac, airmass, exptime, sim.instruments[0].efficiency, area, sim.instruments[0].R, sim.instruments[0].gain, sim.instruments[0].read_noise, dark_current, sim.instruments[0].n_pix, λ_min, λ_max, Δλ)
				vrms = np.sqrt(photonrms**2+instrms**2)
				vmeas = np.random.normal(0.0, vrms)
				line = str(i['obs_start'])+","+str(i['obs_end'])+","+str(i['duration'])+","+str(i['altitude'])+","+str(i['azimuth'])+","+str(photonrms)+","+str(instrms)+","+str(vrms)+","+str(vmeas)+"\n"
				star_rms.write(line)
				print(vrms, vmeas)
		star_rms.close()

'''
# How the other batch file was going to grab star info
	starlist = np.genfromtxt("stars.csv", delimiter=",", dtype=None, names=True)
	for star in starlist:
		starname = str(star['pl_hostname'])[2:-1]
		starobs = np.genfromtxt(starname + ".txt", delimiter=",", names=True)
		starmeas = open(starname + "rv.txt", 'w')
		starmeas.write("obs_time,V_RMS\n")
		for obs in starobs:
			if (obs['altitude'] > 0):
				# get v_rms for every observation
				zenith_angle = (90.0-obs['altitude'])*np.pi/180
				airmass = np.exp(-site_elevation/8400)/np.cos(zenith_angle)
				RMS = rvcalc(star['Teff'], star['FeH'], star['logg'], 1.13*star['vsini'], star['rstar']*u.solRad, star['dstar']*u.pc, airmass, obs['duration']*u.s, efficiency, area, R, gain, read_noise, dark_current, n_pix)
				starmeas.write(str(.5*(obs['obs_start']+obs['obs_end'])) + "," + str(RMS) + "\n")
		starmeas.close()
'''
if __name__ == '__main__':
	#"manually" load up modules for stars. Also, how does the visualization pick run #?
	'''
	starlist = np.genfromtxt("stars.csv", delimiter=",", names=True)
		for star in starlist:
		starobs = np.genfromtxt(star['pl_hostname'], delimiter=",", names=True)
			for obs in starobs:
		# get v_rms for every observation
	'''
'''
#from vis.py -- a way to grab an output
if __name__ == '__main__':
    #ipdb.set_trace()
    # get the full target list
    target_list = simbad_reader.read_simbad('./secret/eta_list.txt') #need to switch to output list
    simnumber = input('Enter 5 digit sim number (including any leading zeros): ')
    simpath = glob.glob('./results/*.'+simnumber+'/')[0]
    os.stat(simpath)
    print(simpath)
    
    all_ras = []
    all_decs = []
    max_num_obs = 0
    min_num_obs = 3000
    num_stars_obs=0
    for target in target_list:

        all_ras.append(target['ra'])
        all_decs.append(target['dec'])
        try:
            tdays,ttimes,alts = get_target(simpath,target['name'])
            target['num_obs'] = len(ttimes)
            num_stars_obs+=1
            if target['num_obs']>max_num_obs:
                max_num_obs = target['num_obs']
            if target['num_obs']<min_num_obs:
                min_num_obs = target['num_obs']
        except:
            target['num_obs'] = 0
    print(num_stars_obs)
'''
'''
#Instrument will fake being HARPS, albeit installed at Kitt Peak where the WIYN should be (2096 m)
site_elevation = 2096.0 # m

#HD 114783
Teff = 5135 #SPOCS
FeH = 0.12 #SPOCS
logg = 4.53 #SPOCS
vsini = 0.9 #SPOCS
theta_rot = 1.13*vsini # Rough approximation, but rotational effects are near-linear, no matter limb-darkening.
rstar = 0.783 * u.solRad #SPOCS
dstar = 20.83 * u.pc #SPOCS

# HARPS For the WIYN!
telescope = 3.566 * u.m
#area = np.pi * (telescope/2)**2 # telescope area, central obstruction ignored
area = 8.8564 * u.m * u.m
efficiency = 0.02 #could be ~1-3%, depending on what is measured
λ_min = 3800 * u.angstrom
λ_max = 6800 * u.angstrom
Δλ = 100 * u.angstrom
BeattyWaves = np.arange(3850, 6850, 100) #Optical range wavelength bits
#BeattyWaves = np.arange((λ_min+Δλ/2)/u.angstrom, (λ_max+Δλ/2)/u.angstrom, Δλ/u.angstrom)
R = 110000 #110k or 120k, depending on source
gain = 1/1.42 # ADU/e-, assume 1:1 photon to electron generation
read_noise = 7.07 #RMS of ± spurious electrons
dark_current = 4 / u.hour
n_pix = 4
exptime = 3600 * u.s
'''
