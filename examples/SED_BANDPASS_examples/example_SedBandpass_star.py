from __future__ import print_function
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass

import numpy as n

# let's treat this star as if it were many stars, with many fluxnorms and calculate the magnitude
# in a single bandpass from those many stars

# pretend we got a list (or array of strings) of star SED names
starskeys = []
for i in range(0, 10):
    starskeys.append('exampleSED.dat')

# pretend we got fluxnorm from somewhere
fluxnorm_min = 3.81e-12
fluxnorm_max = 1.2e-13
fluxnormstep = (fluxnorm_max - fluxnorm_min) / float(len(starskeys))
fluxnorm = n.arange(fluxnorm_min, fluxnorm_max, fluxnormstep)

# and similarly, that we got an array of EBV values
ebvmin = 0
ebvmax = 2
ebvstep = (ebvmax - ebvmin)/float(len(starskeys))
ebv = n.arange(ebvmin, ebvmax, ebvstep)

# okay, instantiate the bandpass used to calculate the magnitude 
rband = Bandpass()
rband.readThroughput("exampleBandpass.dat")

# instantiate stars into dictionary
stars = {}
rootdir = "./"
for key in starskeys:
    stars[key] = Sed()
    stars[key].readSED_flambda(rootdir + key)

# we know the stars have the same wavelength array (stars[i].wavelen)
# (because they were the same SED!)  BUT, in general we may not know this is true.
# so check .. (this is not strictly necessary, but does speed up some of the calculations
#  below. so, for many stars, do it.)

wavelen_base = rband.wavelen
for key in starskeys:
    if stars[key].needResample(wavelen_match = wavelen_base):
        stars[key].resampleSED(wavelen_match=wavelen_base)

# create the a,b arrays for all the stars (because we resampled the stars onto the
#  same wavelength range we can just calculate a/b once, and this is slow)

a, b = stars[starskeys[0]].setupCCMab()

# pretend we want to read mags into an array .. you could just as easily put it into a
# dictionary or list, with small variations in the code
mags = n.empty(len(starskeys), dtype='float')
mags2 = n.empty(len(starskeys), dtype='float')
sedlist = []
for i in range(len(starskeys)):
    # make a copy of the original SED *if* you want to 'reuse' the SED for multiple magnitude
    # calculations with various fluxnorms and dust applications (otherwise just use the object
    # you instantiated above)
    tmpstar = Sed(wavelen=stars[starskeys[i]].wavelen, flambda=stars[starskeys[i]].flambda)
    tmpstar.addCCMDust(a, b, ebv=ebv[i])
    tmpstar.multiplyFluxNorm(fluxnorm[i])
    mags[i] = tmpstar.calcMag(rband)
    # This is for showing an example of the manyMagCalc function on bandpass
    sedlist.append(tmpstar)


# Now, pretend we're actually wanting to calculate the magnitude in multiple bandpasses.
# (ignore the part that uses exampleBandpass.dat as the source .. you would replace that with
#  rootdir + "total_" + filter where filter is a member of lsstfilterlist
lsstfilterlist = ['u', 'g', 'r', 'i', 'z', 'y']
lsst = {}
rootdir = "./"
for filter in lsstfilterlist:
    lsst[filter] = Bandpass()
    lsst[filter].readThroughput(rootdir+"exampleBandpass.dat")
    # you have to do this now - sbToPhi to use the multi-mag calc
    lsst[filter].sbTophi()

# I *do* know that all bandpasses are using the same wavlength array, so let's take it
# for granted that bandpass.wavelen is the same for each.

# Now, we can calculate the magnitudes in multiple bandpasses, for each SED.

# make the bandpass list
bplist = []
for filter in lsstfilterlist:
    bplist.append(lsst[filter])
phiArray, dlambda = tmpstar.setupPhiArray(bplist)
# store the values in a 2-d array (could do in a dictionary of arrays too)
mags = n.empty((len(starskeys), len(bplist)), dtype='float')
for i in range(len(starskeys)):
    # make a copy of the original SED *if* you want to 'reuse' the SED for multiple magnitude
    # calculations with various fluxnorms and dust applications (otherwise just use the object
    # you instantiated above)
    tmpstar = Sed(wavelen=stars[starskeys[i]].wavelen, flambda=stars[starskeys[i]].flambda)
    tmpstar.addCCMDust(a, b, ebv=ebv[i])
    tmpstar.multiplyFluxNorm(fluxnorm[i])
    mags[i] = tmpstar.manyMagCalc(phiArray, dlambda)

print("#sedname mag_u  mag_g   mag_r   mag_i   mag_z   mag_y")
for i in range(len(starskeys)):
    print("%s  %.4f %.4f %.4f %.4f %.4f %.4f" %(starskeys[i], mags[i][0], mags[i][1],
                                                mags[i][2], mags[i][3], mags[i][4],
                                                mags[i][5]))

