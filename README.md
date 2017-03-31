# reimagined-palm-tree
Spectrograph performance, both in terms of an exposure time calculation for a target SNR, and RMS velocity error for input signal and stellar parameters.

All wavelengths are in ångströms, and all velocities are in km/s.
Also, distances are in parsecs, and star sizes are in solar radii.

##### exptime.py
Calculates exposure time to reach a target SNR, given telescope parameters and target apparent magnitude.


##### rvrms.py
Calculates SNR and RV precision, given telescope and stellar parameters. Radial Velocity uncertainty based off of signal in a somewhat idealized spectrograph, and stellar type. Stellar spectra are BT-Settl, if in a somewhat non-standard format. Other equations and data from Beatty and Gaudi (2015). DOI: 10.1086/684264

Also uses http://www.aanda.org/articles/aa/full/2001/29/aa1316/aa1316.right.html for reference.

Current limitations: Detector properties not closely based on actual hardware, microturbulence fixed at ~1 km/s, granulation/starspots/other jitter not considered.
