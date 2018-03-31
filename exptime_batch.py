'''
Calculating all exposure times given a stellar input list. Requires a "targetstars.csv" in the same folder, which it uses for Teff, Fe/H, log(g), v*sin(i), distance, and stellar radius. Instrumental parameters are taken from config/instrument.ini. Output is eta_list.txt, which is formatted to be useable by the dispatch scheduler (move to secret/).
'''

import numpy as np
import exptime
from astropy import units as u
import simulation
import datetime
now = str((datetime.datetime.now()).isoformat())
sim = simulation.simulation('simulation.ini')

Δλ = 100 * u.angstrom
λ_max = sim.instruments[0].λmax * u.angstrom
λ_min = sim.instruments[0].λmin * u.angstrom
sigma_v = 1e-5 # km/s
zenith_angle = 10*np.pi/180 #obviously should be set based on typical object altitude
atmo = exptime.airmass(zenith_angle,sim.elevation,8400)
area = sim.telescopes[0].area * u.m * u.m

starlist = np.genfromtxt("targetstars.csv", delimiter=",", dtype=None, names=True)
eta_list = open("eta_list.txt", 'w')
eta_list.write("Automatically generated list by exptime_batch.py on "+now+"\n\n")
eta_list.write(" IdentList\n----------\n\n")
eta_list.write(" # \ttyped ident\t  coord1 (ICRS,J2000/2000)     \tMag V\t      spec type  \tExp Time\n")
eta_list.write("--- ----------- -----------------------------   ------  --------------- -------------\n")
for star in starlist:
    Teff = star['K']
    λ_peak = exptime.λ_peak(Teff, λ_min, λ_max)
    FeH = star['Sun']
    logg = star['cms']
    vsini = star['kms'] * u.km / u.s
    theta_rot = 1.13 * star['kms']
    rstar = star['solRad'] * u.solRad
    dstar = star['pc'] * u.pc
    dark_current = sim.instruments[0].dark_current / u.hour
    read_time = sim.instruments[0].read_time * u.s
    guess = exptime.time_guess(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, sim.instruments[0].efficiency, area, sim.instruments[0].R, sim.instruments[0].gain, sim.instruments[0].read_noise, dark_current, sim.instruments[0].n_pix, λ_peak, sim.instruments[0].well_depth)
    actual, readout = exptime.time_actual(sigma_v, Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, guess, sim.instruments[0].efficiency, area, sim.instruments[0].R, sim.instruments[0].gain, sim.instruments[0].read_noise, dark_current, sim.instruments[0].n_pix, λ_min, λ_max, Δλ, read_time)
    line = '0'+'\t'+str(star['Name'])[2:-1]+'\t'+str(star['hms'])[2:-1]+' '+str(star['dms'])[2:-1]+'\t'+str(star['mag'])+"\tDQZ\t"+str((actual+readout).to(u.minute).value)+'\t'+str(actual.to(u.minute).value)+'\n'
    print(line)
    eta_list.write(line)
    #print("guess, actual exposure, readout(s), total")
    #print(guess, actual, readout.to(u.minute), actual+readout)
eta_list.close()
'''
def time_guess(Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_peak):
actual, readout = time_actual(sigma_v, Teff, FeH, logg, vsini, theta_rot, rstar, dstar, atmo, guess, efficiency, area, R, gain, read_noise, dark_current, n_pix, λ_min, λ_max, Δλ)
'''