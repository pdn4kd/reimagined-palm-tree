'''
Calculating all exposure times given a stellar input list. Requires a "targetstars.csv" in the same folder, which it uses for Teff, Fe/H, log(g), v*sin(i), distance, and stellar radius. Instrumental parameters are taken from config/instrument.ini. Output is eta_list.txt, which is formatted to be useable by the dispatch scheduler (move to secret/).

This is a somewhat ad-hoc script, and you absolutely should edit it if you find other choices of headings clearer. Likewise sigma_v for differing precisions.

This particular variant finds a minimum p-mode timescale based off of stellar mass and radius.
'''

import numpy as np
import exptime
from astropy import units as u
from astropy import coordinates as coord
import simulation
import datetime
import re
now = str((datetime.datetime.now()).isoformat())
sim = simulation.simulation('simulation.ini')

Δλ = 100 * u.angstrom
λ_max = sim.instruments[0].λmax * u.angstrom
λ_min = sim.instruments[0].λmin * u.angstrom
#t_min = 300 * u.s # 5 minutes in seconds
zenith_angle = 10*np.pi/180 #obviously should be set based on typical object altitude
atmo = exptime.airmass(zenith_angle,sim.elevation,8400)
area = sim.telescopes[0].area * u.m * u.m

sigma_v = 27e-5 # km/s
starlist = np.genfromtxt("targetstars.csv", delimiter=",", dtype=None, names=True)
eta_list = open("eta_list.txt", 'w')
eta_csv = open("eta_list.csv", 'w')

eta_csv.write("name,Exp time,Sky Time,Single Exposure,Number of Exposures,SNR_approx,RV_approx\n")
eta_list.write("Generated on "+now+" with a target precision of "+str(sigma_v)+" m/s or SNR "+str(sim.instruments[0].SNR)+"\n\n")
eta_list.write(" IdentList\n----------\n\n")
eta_list.write(" # \ttyped ident\t  coord1 (ICRS,J2000/2000)     \tMag V\t      spec type  \tExp Time\tSky Time\tSingle Exposure\tNumber of Exposures\n")
eta_list.write("--- ----------- -----------------------------   ------  --------------- ------------- ------- ---- --\n")
for star in starlist:
    Teff = star['K']
    λ_peak = exptime.λ_peak(Teff, λ_min, λ_max)
    FeH = star['Sun']
    logg = star['cms']
    try:
        vsini = float(star['kms']) * u.km / u.s
        theta_rot = 1.13 * float(star['kms'])
    except:
        print("Warning: v * sin(i) not found. Assuming 2 km/s")
        vsini = 2.0 * u.km / u.s
        theta_rot = 2.26
    if np.isnan(theta_rot):
        print("Warning: v * sin(i) not found. Assuming 2 km/s")
        vsini = 2.0 * u.km / u.s
        theta_rot = 2.26
    rstar = star['solRad'] * u.solRad
    dstar = star['pc'] * u.pc
    try:
        vmac = star['Vmac']
        vmac + 1.0
    except:
        print("Warning: Macroturbulence not found. Estimating from other properties.")
        vmac = float('nan')
    t_min = 300 * u.s * np.sqrt(star['solRad']**3/star['Mstar']) # p-mode oscillation timescale, scaled by the dynamical timescale of the star. This defaults to 5 minutes for the Sun.
    dark_current = sim.instruments[0].dark_current / u.hour
    read_time = sim.instruments[0].read_time * u.s
    guess, SNR_sat = exptime.time_guess(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, vmac, atmo, sim.instruments[0].efficiency, area, sim.instruments[0].R, sim.instruments[0].gain, sim.instruments[0].read_noise, dark_current, sim.instruments[0].n_pix, λ_peak, sim.instruments[0].well_depth)
    actual, readout, exposure, number, SNR_actual, RV_actual, SNR_final, RV_final = exptime.time_actual(sigma_v, Teff, FeH, logg, vsini, theta_rot, rstar, dstar, vmac, atmo, guess, sim.instruments[0].efficiency, area, sim.instruments[0].R, sim.instruments[0].gain, sim.instruments[0].read_noise, dark_current, sim.instruments[0].n_pix, λ_min, λ_max, Δλ, read_time, t_min, sim.instruments[0].SNR, SNR_sat)
    # We don't know if we're getting decimal degrees or sexigesimal formatted coordinates
    try:
        RADEC = str(star['hms'])[2:-1]+' '+str(star['dms'])[2:-1]
    except ValueError:
        coords = coord.SkyCoord(ra=star['_RAJ2000']*u.degree, dec=star['_DEJ2000']*u.degree)
        RADEC = coords.to_string('hmsdms')
        # Unfortunately, Astropy as of the writing of this (2019-09-09) does not convert to
        # Simbad style format (spaces between hours, minutes, etc) but with letters in between.
        # So we fix this with a few regexes.
        RADEC = re.sub('[dhm]', ' ', RADEC)
        RADEC = re.sub('s', '', RADEC)
    try:
        MK = str(star['MK'])[2:-1] #Spectra type is nice to have
    except:
        MK = "XXX"
    if MK == "":
        MK = "XXX"
    line = '0'+'\t'+str(star['Name'])[2:-1]+'\t'+RADEC+'\t'+str(star['mag'])+"\t"+MK+"\t"+str((actual+readout).to(u.minute).value)+'\t'+str(actual.to(u.minute).value)+'\t'+str(exposure.to(u.minute).value)+'\t'+str(int(number))+'\n'
    print(line)
    eta_list.write(line)
    line = str(star['Name'])[2:-1]+","+str((actual+readout).to(u.minute).value)+','+str(actual.to(u.minute).value)+','+str(exposure.to(u.minute).value)+','+str(int(number))+','+str(SNR_final)+','+str(RV_final)+'\n'
    eta_csv.write(line)
    #print("guess, actual exposure, readout(s), total")
    #print(guess, actual, readout.to(u.minute), actual+readout)
eta_list.close()
eta_csv.close()
'''
def time_guess(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_peak):
actual, readout = time_actual(sigma_v, Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, guess, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_min, λ_max, Δλ)
'''
