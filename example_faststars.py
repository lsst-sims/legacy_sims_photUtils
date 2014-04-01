#

# CONCLUSION: don't  bother using bandpass's manyMagCalc function. Use Sed's instead.
# the memory management required for bandpass's function seems to make it *slower* than doing the magnitude calculation
# in the straight-forward manner. (actually, Bandpass's manyMagCalc function has been removed).  

# Otherwise, the manyMagCalc function speeds up magnitude calculation by about 0.7 times
#  (if the filter-by-filter magnitude calculation was done in a smart way, using the same wavelength grid
#  for all stars when they were read in and pre-calculating the MW dust extinction).
#  Could calculate 6 magnitudes with dust extinction for 10,000 stars in 24 s without manyMagCalc,
#    and in 17 s using manyMagCalc.

# Note: if you did not pre-calculate the dust extinction and did not place the stars seds onto
# the same wavelength grid as the bandpass before starting any calculations, this would increase
# the time requirements. 

import os

# Set up timing information and simple routine. 
import time
def dtime(time_prev):
   return (time.time() - time_prev, time.time())

# Import other modules and classes. 
import numpy
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass

# Initialize starting time.
t = time.time()

# Set the wavelength step for bandpasses. 
wavelen_step = 0.1
#wavelen_step = 0.25

# Read in LSST bandpasses.
bpdir = os.getenv("LSST_THROUGHPUTS_BASELINE")
filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
lsstbp = {}
for f in filterlist:
    lsstbp[f] = Bandpass()
    lsstbp[f].readThroughput(os.path.join(bpdir, 'total_'+f+'.dat'), wavelen_step=wavelen_step)
dt, t = dtime(t)
print "Reading %d filters took %f s" %(len(filterlist), dt)

# Read in a set of star seds. 
# Replace stardir with your root galaxy sed directory.
stardir = os.path.join(os.getenv('SED_DATA'), 'starSED/kurucz')
stars = {}
starlist = os.listdir(stardir)
for star in starlist:
    stars[star] = Sed()
    stars[star].readSED_flambda(os.path.join(stardir,star))

dt, t = dtime(t)
print "Reading %d star seds took %f s" %(len(starlist), dt)

# Check on resampling - want all seds to have the same wavelength range and it can match the
# bandpass wavelength grid because we're not redshifting the stars.
resampled = False
for star in starlist:
    if stars[star].needResample(wavelen_match = lsstbp['g'].wavelen):
        stars[star].resampleSED(wavelen_match = lsstbp['g'].wavelen)
        resampled = True
dt, t = dtime(t)
if resampled:
    print "Checking for (and doing) resampling took %f s" %(dt)
else:
    print "Checking for (but not doing any) resampling took %f s" %(dt)

# Okay, now decide how many stars to include in this 'chunk'.
# We will assign the star seds above randomly to each star within this chunk (to simulate getting
# an sed name from the database), then add dust, and calculate the magnitude of the star in all
# ugrizy bandpasses.  This is similar to the example_fastgals.py test, but will use
# Bandpass.manyMagCalc as well as Sed.manyMagCalc (in two different magnitude calculations, as a test).
# This assumes that you need to access the result as an array of u magnitudes, then an array of g magnitudes, etc.

num_star = 10000
ebv_mw = numpy.random.rand(num_star)  # random 0-1
star_name = numpy.random.rand(num_star) * len(starlist)  # 'assign' sed names
star_name = numpy.array(star_name, dtype='int')
fluxnorm = (numpy.random.rand(num_star) + 2) * 1e-14

dt, t = dtime(t)
print "Picking random numbers for dust/fluxnorm/sed id  took %f s" %(dt)

# First - start 'regular' (but slightly slower) method of calculating magnitudes for stars, for comparison.
# If you're only calculating magnitudes for a few stars, this might actually be just as fast and
# would likely be easier to code/read, so is useful as a comparison.
# Actual difference in timing between this method and the next (more optimized) method can be determined
# by simply running 'python example_faststars.py', but on my mac it was about 2.5 times faster optimized, with
# wavelen_step = 0.1 nm. (wavelen_step will have an impact on the speed difference).
# There will be a difference in the speedup reported here vs. in example_fastgals.py primarily because the
# galaxies require more calculations (which aren't part of manyMagCalc) and partly because we're using a different
# optimization (bandpass.manyMagCalc instead of sed.manyMagCalc). 

# Start 'regular' magnitude calculation.

# Calculate dust extinction a_x and b_x vectors.
a_mw, b_mw = stars[starlist[0]].setupCCMab()
# Set up dictionary + arrays to hold calculated magnitude information. 
mags1 = {}
for f in filterlist:
    mags1[f] = numpy.zeros(num_star, dtype='float')
# For each star (in num_star's), apply apply MW dust, fluxnorm & calculate mags. 
for i in range(num_star):
    starname = starlist[star_name[i]]
    tmpstar = Sed(wavelen=stars[starname].wavelen, flambda=stars[starname].flambda)
    tmpstar.addCCMDust(a_mw, b_mw, ebv=ebv_mw[i])
    tmpstar.multiplyFluxNorm(fluxnorm[i])
    # Note that the stars have already been matched to the bandpass wavelength grid.
    # Just want to be sure that fnu is already calculated.
    tmpstar.flambdaTofnu()
    for f in filterlist:
        mags1[f][i] = tmpstar.calcMag(lsstbp[f])
dt, t = dtime(t)
print "Calculating dust/fluxnorm/%d magnitudes with some smart usage for %d stars took %f s" \
      %(len(filterlist), num_star, dt)


# Test Sed.manyMagCalc : 

# First: (re) calculate internal a/b on wavelength range required for dust extinction.
a_mw, b_mw = stars[starlist[0]].setupCCMab()  
# Also: set up phi for each bandpass - ahead of time. And set up a list of bandpasses, then create phiarray 
# and dlambda to set up for manyMagCalc method.
bplist = []
for f in filterlist:
    lsstbp[f].sbTophi()
    bplist.append(lsstbp[f])
phiarray, dlambda = stars[starlist[0]].setupPhiArray(bplist)
# Set up dictionary + arrays to hold calculated magnitude information. 
mags2 = {}
for f in filterlist:
    mags2[f] = numpy.zeros(num_star, dtype='float')
# For each star (in num_star's), apply MW dust, fluxnorm, and then calculate mags using manyMagCalc. 
for i in range(num_star):
    starname = starlist[star_name[i]]
    tmpstar = Sed(wavelen=stars[starname].wavelen, flambda=stars[starname].flambda)
    tmpstar.addCCMDust(a_mw, b_mw, ebv=ebv_mw[i])
    tmpstar.multiplyFluxNorm(fluxnorm[i])
    tmpmags = tmpstar.manyMagCalc(phiarray, dlambda)
    j = 0
    for f in filterlist:
        mags2[f][i] = tmpmags[j]
        j = j+1
dt, t = dtime(t)
print "Calculating dust/fluxnorm/%d magnitudes for %d stars optimized with Sed.manyMagCalc took %f s" \
      %(len(filterlist), num_star, dt)


"""
This test turned out so badly that I removed these methods from Bandpass.
Basically, it was slower than just the normal calculation above, and took up too much memory.

Code kept here as example. Replace the next two def's into Bandpass if you need to use this for
some bizarre reason.

## Bonus, many-magnitude calculation for many SEDs with a single bandpass
##  Note: testing of these functions indicate that this is NOT a good optimization.
##  Use these if you have to, but it is SLOWER than just calculating a magnitude a normal way, and
##  definitely slower than using Sed's manyMagCalc function. 

    def setupFnuArray(self, seddict, sedlist):
        ""Set up 2-d array of sed FNU values and zeropoint suitable for input to Bandpass's manyMagCalc. 

        Input a *dictionary* of seds.
        This is intended to be used once on a list of Sed's to prepare them for several calls
        to Bandpass's manyMagCalc, such as before calculating many magnitudes in different bands.
        It will check if the seds are resampled to this bandpass's wavelength grid and if not, will
        carry out that resampling on the fnu array values (only here - not affecting the sed itself). 
        Returns 2-d array of sed fnu values.
        ATTENTION - testing has proved that this is slower than any normal mag calculation.
        Memory management most likely the culprit. Use Sed.manyMagCalc function instead. ""
        fnuArray = n.empty((len(sedlist), len(self.wavelen)), dtype='float')
        i = 0
        zp = seddict[sedlist[0]].zp
        for s in sedlist:
            # Check that zeropoint is the same for this Sed as all previous.
            # This will always be the case for UW Seds as they are all in the same units. 
            if seddict[s].zp != zp:
                raise Exception(""Sed %s has a zeropoint of %f (others have %f). Please check your sed units.""\
                                % (s, seddict[s].zp, zp))
            # Check that fnu has already been calculated.
            if seddict[s].fnu == None:
                seddict[s].flambdaTofnu()
            # Check that fnu is on the necessary wavelength grid. (remember self=bandpass).
            if seddict[s].needResample(wavelen=seddict[s].wavelen, wavelen_match=self.wavelen):
                wavelen, fnu = self.resampleSED(seddict[s].wavelen, seddict[s].fnu, wavelen_match=self.wavelen)
                fnuArray[i] = fnu
            else:
                fnuArray[i] = seddict[s].fnu
            i = i + 1
        return fnuArray, zp
    
    def manyMagCalc(self, fnuArray, zp):
        "" Calculate many magnitudes for many seds using a single bandpass.

        So, testing has shown that this method is slow, slower than just calculating a normal
        magnitude (see examples/example_faststars.py for proof).
        Thus I am now going to delete it from Bandpass. Keeping a copy in SVN though, just in
        case it's ever useful. 

        This method assumes that there will be flux within a particular bandpass
        (could return '-Inf' for a magnitude if there is none), and use setupFnuArray first, on this
        bandpass itself or on another bandpass with the same wavelength grid (i.e. here it is assumed that
        the wavelength grid is correct).
        Also assumes that phi has already been calculated for this bandpass (using self.sbTophi()). 
        These assumptions are to avoid error checking within this function (for speed), but could lead
        to errors if method is used incorrectly.
        ATTENTION - testing has proved that this is slower than any normal mag calculation.
        Memory management most likely the culprit. Use Sed.manyMagCalc function instead. ""
        mags = n.empty(len(fnuArray), dtype='float')
        dlambda = self.wavelen[1] - self.wavelen[0]
        mags = -2.5*n.log10(n.sum(self.phi*fnuArray, axis=1)*dlambda) - zp
        return mags

###

# Test Bandpass.manyMagCalc

# First: (re) calculate internal a/b on wavelength range required for dust extinction.
a_mw, b_mw = stars[starlist[0]].setupCCMab()
# Set up all the stars in memory (so set num_stars / chunksize carefully)
allstars = {}
allstarlist = []
for i in range(num_star):
    starname = starlist[star_name[i]]
    this_starname = starname + "%d" %(i)
    allstarlist.append(this_starname)
    allstars[this_starname] = Sed(wavelen=stars[starname].wavelen, flambda=stars[starname].flambda)
    allstars[this_starname].addCCMDust(a_mw, b_mw, ebv=ebv_mw[i])
    allstars[this_starname].multiplyFluxNorm(fluxnorm[i])
# and then set up fnu array
fnuArray, zp = lsstbp['g'].setupFnuArray(allstars, allstarlist)
dt1, t = dtime(t)
print ""Setting up fnu array took %f s"" \
      %( dt1)
# Now release the other memory ... (hopefully)
stars = None
allstars = None
mags3 = {}
for f in filterlist:
    mags3[f] = lsstbp[f].manyMagCalc(fnuArray, zp)
dt, t = dtime(t)
print ""Calculating dust/fluxnorm/%d magnitudes for %d stars optimized with Bandpass.manyMagCalc took %f s"" \
      %(len(filterlist), num_star, (dt+dt1))
"""

# Check for differences in magnitudes.
import pylab
pylab.figure()
colors = ['m', 'g', 'r', 'b', 'k', 'y']
i = 0
diff1 = {}
#diff2 = {}
for f in filterlist:
    flags = numpy.isfinite(mags1[f])
    if (flags.any() == 'False'):
        print "Found %d finite magnitudes out of %d in non-optimized mags %s" %(len(flags[0]), len(mags1[f]), f)
    flags = numpy.isfinite(mags2[f])
    if (flags.any() == 'False'):
        print "Found %d finite magnitudes out of %d in optimized mags2 %s" %(len(flags[0]), len(mags2[f]), f)
    #flags = numpy.isfinite(mags3[f])
    #if (flags.any() == 'False'):
    #    print "Found %d finite magnitudes out of %d in optimized mags3 %s" %(len(flags[0]), len(mags3[f]), f)
    diff1[f] = numpy.zeros(num_star, dtype='float')
    #diff2[f] = numpy.zeros(num_star, dtype='float')
    diff1[f] = numpy.abs(mags1[f] - mags2[f])
    #diff2[f] = numpy.abs(mags1[f] - mags3[f])
    print "In ", f, " mags2 max :", diff1[f].max()#, " mags3 max: ", diff2[f].max()
    #condition  = (diff1[f]>0.01)
    #print f, diff1[f][condition], mags1[f][condition], mags2[f][condition]
    #condition  = (diff2[f]>0.01)
    #print f, diff2[f][condition], mags1[f][condition], mags3[f][condition]
    pylab.plot(mags1[f], mags2[f], colors[i]+'+')
    #pylab.plot(mags1[f], mags3[f], colors[i]+"x")
    i = i + 1
x = numpy.arange(10, 35, 1)
pylab.plot(x, x, 'k-')
pylab.show()
