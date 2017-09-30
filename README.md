# reimagined-palm-tree
Spectrograph performance, both in terms of an exposure time calculation for a target SNR, and RMS velocity error for input signal and stellar parameters.

All wavelengths are in ångströms, and all velocities are in km/s.
Also, distances are in parsecs, and star sizes are in solar radii. Usage of astropy.units should limit the need to keep careful track of what is used where.

Required libraries: numpy, astropy

##### exptime.py
Calculates exposure time to reach a target SNR, given telescope parameters and target apparent magnitude. Also acts as a library for exptime_batch.py

#### exptime_batch.py
Generates exposure times, and creates the target list of stars (eta_list.txt) that dispatch_scheduler (https://github.com/pdn4kd/dispatch_scheduler) needs, given an existing CSV with the relevant properties (targetstars.csv). Relies on the dispatch_scheduler for telescope, site, etc. parameters. To use, put exptime_batch, exptime, and rvrms in the dispatch_scheduler base directory.

As this uses rvrms.py, it suffers from that file's limitations. Also, macroturbulence really should be using observational values where possible.

#### exptime_demo.py
Simple example of calculating exposure times across a range of magnitudes. Demonstrates that you can analytically get an exposure time if you have a target SNR. (Target radial velocities are rather more complex)

##### rvrms.py
Calculates SNR and RV precision, given telescope and stellar parameters. Radial Velocity uncertainty based off of signal in a somewhat idealized spectrograph, and stellar type. Stellar spectra are BT-Settl, if in a somewhat non-standard format. Other equations and data (http://www.personal.psu.edu/tgb15/beattygaudi/table1.dat) from Beatty and Gaudi (2015). DOI: 10.1086/684264

Also uses http://www.aanda.org/articles/aa/full/2001/29/aa1316/aa1316.right.html for reference.

Current limitations: Detector properties not that closely based on actual hardware (especially sensitivity vs wavelength), macroturbulence assumed based on effective temperature instead of observations, microturbulence fixed at ~1 km/s; granulation/starspots/other jitter not considered; fixed assumptions made of Teff (200 K steps) and spectrograph R values made based on wavelength range.

The spectrograph issues bother pdn4kd greatly, and a way to not have to use the table (with its assumptions on R and wavelength ranges) would be greatly appreciated.

##### rvrms2.py
A Python 2 compatable version of rvrms.py (all unicode replaced with ASCII). This is probably silly, given pdn4kd's tendency to use unicode in other functions, including the modified version of dispatch_scheduler.

#### rvrms_batch.py
Convert observations found in dispatch_scheduler to radial velocities. Currently more proof of concept than anything else. Will require the dispatch_scheduler, as it will be primarily used to generate radial velocities from observations done in it.
