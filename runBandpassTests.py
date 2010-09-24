"""
  Questions or comments, email : ljones.uw@gmail.com

$Id$

 Hopefully this is helpful to the engineering types.

"""

import os
import numpy
import pylab
import Sed
import Bandpass
import BandpassSet

# Step one: 
# Set rootdir = to the directory where the throughput curves being tested are located.
rootdir = os.getenv("LSST_THROUGHPUTS_BASELINE")
# or you could set rootdir to a particular directory on your system ..
# rootdir = "/Users/rhiannonjones/throughputs/"

# Step two:
# What kind of bandpass tests are you doing?
# Are you only looking to evaluate a single throughput curve but for multiple filters?
# If so, leave the next line uncommented and fill in the information with what you want to use. 
single = True
if single == True:
    # this is the list of filters you wish to test
    filterlist = ('u', 'g', 'r', 'i', 'z','y')
    # this is the root of the file name (at the start of the filename, before the 'filter' part)
    rootname = "total_"
    # this is the suffix of the file name (at the end of the filename, after the 'filter' part)
    rootsuffix = ".dat"
    # will thus result in evaluating a series of files called total_u.dat, total_g.dat, etc.

# Or do you need to build a final throughput curve from various components, so you can swap
# those components for different tests?  (uncomment the next line and place the relevant info afterwards).
#single = False
if single == False:
    filterlist = ('u', 'g', 'r','i', 'z', 'y')
    # these files will be used for each separate filter component.
    separate_filter_complist = ('filter_u.dat', 'filter_g.dat', 'filter_r.dat',
                                'filter_i.dat', 'filter_z.dat', 'filter_y.dat')
    # these files will be used for the components which apply to all filters. 
    all_filter_complist = ('detector.dat', 'lens1.dat', 'lens2.dat',
                           'lens3.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                           'atmos.dat')



# Okay, let's read this information and then do the tests. 
testSet = BandpassSet.BandpassSet()
if single:
    testSet.setThroughputs_SingleFiles(filterlist=filterlist, rootdir=rootdir, rootname=rootname,
                                       rootsuffix=rootsuffix)
else:
    testSet.setThroughputs_ComponentList(filterlist=filterlist,
                                         separate_filter_complist = separate_filter_complist,
                                         all_filter_complist = all_filter_complist,
                                         rootdir = rootdir)

# And let's do some tests.

# First we'll calculate the easy things: effective wavelengths and telescope zeropoints.
testSet.calcFilterEffWave()
testSet.calcZeroPoints(gain=1.0)  # note - the zeropoint calculation assumes a 6.5m effective area and 15s exp.

# Now let's calculate the 'edges' of the filters.
# First let's do this for 50% of max/50% absolute throughput values.
testSet.calcFilterEdges(drop_peak=0.5, drop_percent=0.5)
# And then also for 10% of max/10% absolute throughput values.
testSet.calcFilterEdges(drop_peak=0.1, drop_percent=0.1)

# Now let's make a plot of the filter throughputs.
# You can change xlim to other wavelength limits for extra plot clarify if desired.
xlim = [300,1100]
# There are lots of options to play with in plotting the filters.
if (xlim[0] > 500) | (xlim[1]<900):
    filter_tags='side'
else:
    filter_tags='normal'
testSet.plotFilters(xlim=xlim, throughput=True, ploteffsb=True, savefig=True, 
                    figroot='bandpass', atmos=False, leg_tag=None, filter_tags=filter_tags)


# Okay, and now let's calculate the filter leaks.
testSet.calcFilterLeaks(makeplot=True, savefig=True, figroot = "bandpass")


pylab.show()
