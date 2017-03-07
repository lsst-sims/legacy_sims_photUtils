from __future__ import print_function
# ljones@astro.washington.edu

# This is an example case of using the Sed.py and Bandpass.py classes. 
# Only some of the (probably most often used) functionality is illustrated here,
#  so please see the comments in Sed and Bandpass for more information.

import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass

# Note that Bandpass and Sed are set up for NANOMETERS units on input data

# Instantiate a Bandpass object to hold the throughput information
rband = Bandpass()
# Read the bandpass data file into that object -
#  (note exampleBandpass.dat contains THROUGHPUT information, not PHI)
bpfile = 'exampleBandpass.dat'
rband.readThroughput(bpfile)

# Instantiate a Sed object to hold the spectral energy distribution information
star = Sed()
# Read the Sed data file into that object -
#  (note that exampleSED.dat contains wavelength/F_lambda .. but is possible to also
#    read in wavelength/F_nu data using readSED_fnu)
sedfile = 'exampleSED.dat'
star.readSED_flambda(sedfile)

print("")
print("Read %s into Bandpass and %s into Sed." %(bpfile, sedfile))

# Simply calculate the magnitude of this source in the bandpass.
mag = star.calcMag(rband)

print("")
print("Without any scaling of the SED, the magnitude is %.4f" %(mag))

# That was probably pretty small, right? Maybe we actually know what
# magnitude we expect this source to have in this bandpass, and then want to scale
# the SED to that appropriate magnitude (and then calculate the magnitudes once it's
# scaled properly, in other bandpasses).

mag_desired = 24.5
print("Now going to apply a scaling factor to the SED to set magnitude to %.4f" %(mag_desired))

# Calculate the scaling factor.
fluxnorm = star.calcFluxNorm(mag_desired, rband)
# Apply the scaling factor. 
star.multiplyFluxNorm(fluxnorm)

# Try the magnitude calculation again. 
mag = star.calcMag(rband)
print("After scaling, magnitude of SED is now %.4f (desired magnitude was %.4f)" %(mag, mag_desired))

# And let's calculate what the expected photon counts for LSST would be.
counts = star.calcADU(rband, expTime=30)
print("This would correspond to roughly %f counts in the LSST focal plane, in a 30s exposure." %(counts))

# For fun, let's see what else can happen.

ebv = 0.5
print("")
print("Let's try adding %.2f E(B-V) dust extinction to this star." %(ebv))
a, b = star.setupCCMab()
# You can use addCCMDust on the 'star' object itself, but I'm illustrating here how you could also
#  do otherwise - preserve the original 'star' object as is, and create a new Sed object that does
#  include the effects of dust ('dustystar').  Star's data will be unchanged by the dust. 
dustywavelen, dustyflambda = star.addCCMDust(a, b, ebv=ebv, wavelen=star.wavelen, flambda=star.flambda)
dustystar = Sed(wavelen=dustywavelen, flambda=dustyflambda)
magdust = dustystar.calcMag(rband)
print("With this dust, the magnitude of the star in this bandpass is now %.4f." %(magdust))

redshift = 0.2
print("What if this star was at a redshift of %f?" %(redshift))
# Here (unlike above with the dust), I'm applying the redshift to the 'star' object itself.
#  Star's data will be changed by the redshifting. 
star.redshiftSED(redshift, dimming=True)
magredshift = star.calcMag(rband)
print("")
print("Redshifted to %.2f, and adding cosmological dimming, the magnitude is now %.4f" \
      %(redshift, magredshift))
print(" (this was a pretty hot star, so redshifting brings more flux into this particular bandpass.)")
print("")

# There is more functionality in the class. Please see the code. 
