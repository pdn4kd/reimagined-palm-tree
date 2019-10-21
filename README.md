# reimagined-palm-tree
Spectrograph performance, both in terms of an exposure time calculation for a target SNR, and RMS velocity error for input signal and stellar parameters.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

All wavelengths are in ångströms, and all velocities are in m/s if unlabled (If lableled, they are likely to be in km/s).
Also, distances are in parsecs, and star sizes are in solar radii. Usage of astropy.units should limit the need to keep careful track of what is used where.

Required libraries: numpy, astropy, re

Requires Python 3, notably due to extensive unicode variable names.

##### exptime.py
Calculates exposure time to reach a target SNR, given telescope parameters and target apparent magnitude. Also acts as a library for exptime_batch.py. Depends on rvrms.py

#### exptime_batch.py
Generates exposure times, and creates the target list of stars (eta_list.txt) that dispatch_scheduler (https://github.com/pdn4kd/dispatch_scheduler) needs, given an existing CSV with the relevant properties (targetstars.csv). Secondarily, it generates a CSV (eta_list.csv) with the object names and outputs for data analysis purposes. Relies on the dispatch_scheduler for telescope, site, etc. parameters. To use, put exptime_batch, exptime, and rvrms in the dispatch_scheduler base directory.

As this uses rvrms.py, it suffers from that file's limitations. Also, macroturbulence really should be using observational values where possible.

''input file information''
The input file (targetstars.csv by default) is a CSV that requires the following parameters (with the default column names in parentheses): Effective Temperature (K), Metallicity (FeH), Distance (pc), Radius (solRadius), log(g) (cms). Optional: Macroturbulence (Vmac, estimated from temperature if unavailable), V*sin(i) (kms, assumed to be 2 km/s if unavailable)
While not needed for the RV precision and exposure time calculations, several additional columns are needed for the output: an identifying name (Name); right ascension (_RAJ2000); declination (_DEJ2000); a magnitude (mag); and a spectral type (MK).
These are somewhat baroque/nonstandard names, which hopefully means that the user has closely checked their data.
An example targetstars.csv is included for clarity.

''output file information''
The file is in a quasi-Simbad format, with tabs between the major values (except for RA/Dec, which are seperated by a space). The "columns" are: # (an arbitrary number, set to 0 by default. This can be changed if desired); typed ident (a name for the target object); coord1 (ICRS,J2000/2000) (RA/DEC in HMS/DMS format); Mag V (V-band magnitude, usually apparent); spec type (Spectral Type)
ExpTime is total time spent in terms of imaging and readouts that the telescope spends pointed at a given target. ("wall-clock" time)

SkyTime is the time with the shutter open staring at a given target (SingleExposure * N).

SingleExposure is how long an individual exposure (excluding readout) as part of single epoch measurement takes.

#### exptime_demo.py
Simple example of calculating exposure times across a range of magnitudes. Demonstrates that you can analytically get an exposure time if you have a target SNR. (Target radial velocities are rather more complex)

##### rvrms.py
Calculates SNR and RV precision, given telescope and stellar parameters. Radial Velocity uncertainty based off of signal in a somewhat idealized spectrograph, and stellar type. Stellar spectra are BT-Settl, if in a somewhat non-standard format. Other equations and data (http://www.personal.psu.edu/tgb15/beattygaudi/table1.dat) from Beatty and Gaudi (2015). DOI: 10.1086/684264

Example spectra are included (named 1700-9800, every 100 K below 7000 K and every 200 K above). Also uses http://www.aanda.org/articles/aa/full/2001/29/aa1316/aa1316.right.html for reference.

Current limitations: Detector properties not that closely based on actual hardware (especially sensitivity vs wavelength), microturbulence fixed at ~1 km/s; granulation/starspots/other jitter not considered; fixed assumptions made of Teff (200 K steps for RV information, must be between 2600 K and 7600 K inclusive, Stellar spectra in 100 K steps); wavelength ranges (400-2500 nm); spectrograph R values made based on wavelength range.

The spectrograph issues bother pdn4kd greatly, and a way to not have to use the table (with its assumptions on R and wavelength ranges) would be greatly appreciated.

#### rvrms_batch.py
Convert observations found in dispatch_scheduler to radial velocities. Needs additional work due to fragility/assumptions. Requires the dispatch_scheduler repo, and it is primarily used to generate radial velocities from observations done there.
