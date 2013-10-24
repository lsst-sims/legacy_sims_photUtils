# 11/15/2010
# lynne jones
# simple example of calculating the colors of a set of varying temperature
#  black bodies.  (Zeljko was looking for a table with this info for SDSS and LSST,
#  and I thought it might serve as an example of using Sed/Bandpass in a method 
#  other than previously shown in example_SedBandpass, etc.)


import os
import numpy

import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass


lsstdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
sdssdir = os.getenv("SDSS_THROUGHPUTS")

sdssflist = ('u', 'g', 'r', 'i', 'z')
lsstflist = ('u', 'g', 'r', 'i', 'z', 'y')

sdss = {}
for f in sdssflist:
    sdss[f] = Bandpass()
    sdss[f].readThroughput(os.path.join(sdssdir, "sdss_"+f+".dat"))

lsst = {}
for f in lsstflist:
    lsst[f] = Bandpass()
    lsst[f].readThroughput(os.path.join(lsstdir, "total_"+f+".dat"))



LIGHTSPEED = 299792458       # speed of light, = 2.9979e8 m/s
PLANCK = 6.626068e-27        # planck's constant, = 6.626068e-27 ergs*seconds
NM2M = 1.00e-9               # nanometers to meters conversion = 1e-9 m/nm
BOLTZMAN = 1.3806504e-16     # boltzman constant, = 1.38e-16 ergs/K

wavelen = numpy.arange(300, 1200, 0.1)
nu = LIGHTSPEED/(wavelen*NM2M)

# set up header
writestring = "## Temperature  SDSS:"
for i in range(len(sdssflist)-1):
    col = "%s-%s " %(sdssflist[i], sdssflist[i+1])
    writestring = writestring + col
writestring = writestring + " LSST:"
for i in range(len(lsstflist)-1):
    col = "%s-%s " %(lsstflist[i], lsstflist[i+1])
    writestring = writestring + col
print writestring

# calculate colors
temperatures = numpy.arange(3000, 30000, 1000, dtype='float')
for temperature in temperatures:
    # blackbody: B_lambda dl = (2*h*c^2)/(wavelen^5)*(e^(h*c/(wavelen*k*T))-1)^-1 dl
    # blackbody: B_nu dnu = (2*h/c^2)*nu^3 * (e^(h*nu/k*T)-1)^-1 dnu
    # Set up SED object with blackbody spectrum.
    fnu = ((2.0*PLANCK/LIGHTSPEED**2)*(nu**3) /
           (numpy.e**(PLANCK*nu/BOLTZMAN/temperature)-1))
    blackbody = Sed(wavelen, fnu=fnu)
    # wien's law: wavelen_max(m) = 0.0029/T
    #max = 0.0029/temperature/NM2M
    writestring = "%d " %(temperature)
    # Calculate magnitudes.
    sdssmag = {}
    for f in sdssflist:
        sdssmag[f] = blackbody.calcMag(sdss[f])
    lsstmag = {}
    for f in lsstflist:
        lsstmag[f] = blackbody.calcMag(lsst[f])
    # Calculate colors. 
    for i in range(len(sdssflist)-1):
        col = sdssmag[sdssflist[i]] - sdssmag[sdssflist[i+1]]
        writestring = writestring + "%.3f " %(col)
    for i in range(len(lsstflist)-1):
        col = lsstmag[lsstflist[i]] - lsstmag[lsstflist[i+1]]
        writestring = writestring + "%.3f " %(col)    
    print writestring


