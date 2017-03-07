from __future__ import print_function
from builtins import range
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass

import numpy as n

# let's treat this galaxy as if it were many galaxies, with many redshifts, dust extinctions,
# and fluxnorms, and calculate the magnitudes in a single bandpass from those many galaxies

# pretend we got a list (or array of strings) of galaxy SED names
galaxykeys = []
for i in range(0, 10):
    galaxykeys.append('exampleSED.dat')

# pretend we got fluxnorm from somewhere
fluxnorm_min = 3.81e-12
fluxnorm_max = 1.2e-13
fluxnormstep = (fluxnorm_max - fluxnorm_min) / float(len(galaxykeys))
fluxnorm = n.arange(fluxnorm_min, fluxnorm_max, fluxnormstep)

# and similarly, that we got an array of EBV values for the galaxy internal extinction
ebvmin = 0
ebvmax = 2
ebvstep = (ebvmax - ebvmin)/float(len(galaxykeys))
ebv_gal = n.arange(ebvmin, ebvmax, ebvstep)

# and that we got an array of EBV values for the milky way line of sight extinction
ebvmin = 0
ebvmax = 2
ebvstep = (ebvmax - ebvmin)/float(len(galaxykeys))
ebv_mw = n.arange(ebvmin, ebvmax, ebvstep)

# and that we have a range of values for redshift
redshiftmin = 0.5
redshiftmax = 1.5
redshiftstep = (redshiftmax - redshiftmin)/float(len(galaxykeys))
redshifts = n.arange(redshiftmin, redshiftmax, redshiftstep)

# okay, instantiate the bandpass used to calculate the magnitude 
rband = Bandpass()
rband.readThroughput("exampleBandpass.dat")

# instantiate gals into dictionary
gals = {}
rootdir = "./"
for key in galaxykeys:
    gals[key] = Sed()
    gals[key].readSED_flambda(rootdir + key)

# we know the gals have the same wavelength array (gals[i].wavelen)
# (because they were the same SED!)  BUT, in general we may not know this is true.
# so check .. (this is not strictly necessary, but does speed up some of the calculations
#  below. so, for many gals, do it.). For GALAXIES we want to be sure to base this on
# a wide enough range, so use a galaxy that is good for the wavelen_base.

wavelen_base = gals[key].wavelen
for key in galaxykeys:
    if gals[key].needResample(wavelen_match = wavelen_base):
        gals[key].resampleSED(wavelen_match=wavelen_base)

# now for galaxies, we have to do (a) multiply the fluxnorm, (b) add the dust internal to
#  the galaxy, (c) redshift the galaxy (d) add the dust from the milky way.
# Note that step (a) can be done at any time, but b-c-d must be done in order. 

# create the a,b arrays for all the gals (because we resampled the gals onto the
#  same wavelength range we can just calculate a/b once, and this is slow)

a_gal, b_gal = gals[galaxykeys[0]].setupCCMab()

# pretend we want to read mags into an array .. you could just as easily put it into a
# dictionary or list, with small variations in the code
mags = n.empty(len(galaxykeys), dtype='float')

for i in range(len(galaxykeys)):
    # make a copy of the original SED if you want to 'reuse' the SED for multiple magnitude
    # calculations with various fluxnorms, redshifts and dusts
    tmpgal = Sed(wavelen=gals[galaxykeys[i]].wavelen, flambda=gals[galaxykeys[i]].flambda)
    # add the dust internal to the distant galaxy
    tmpgal.addCCMDust(a_gal, b_gal, ebv=ebv_gal[i])
    # redshift the galaxy
    tmpgal.redshiftSED(redshifts[i], dimming=False)
    # add the dust from our milky way - have to recalculate a/b because now wavelenghts
    # for each galaxy are *different*
    a_mw, b_mw = tmpgal.setupCCMab()
    tmpgal.addCCMDust(a_mw, b_mw, ebv=ebv_mw[i])
    tmpgal.multiplyFluxNorm(fluxnorm[i])
    mags[i] = tmpgal.calcMag(rband)


# show results
print("#sedname      fluxnorm     redshift  ebv_gal   ebv_mw  magnitude ")
for i in range(len(galaxykeys)):
    print("%s %.5g  %.3f %.5f %.5f %.5f" %(galaxykeys[i], fluxnorm[i], redshifts[i], ebv_gal[i], ebv_mw[i], mags[i]))
    
