from __future__ import print_function
import sys
import os
import copy
import numpy
import scipy.optimize
import pylab
from .Sed import Sed
from .SedSet import SedSet
from .Bandpass import Bandpass
from .BandpassSet import BandpassSet


# Set some values for overall test. 
filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
# Set all filter 'edges' (95% throughput boundaries). 
filterEdges = [311, 398, 550, 691, 813, 918, 1050]
# Set values that we will use to generate different 'versions' of the filters.
# taper = distance between the top of the filter and the bottom (0).
tapers = [0,2,4,8,16,32,64]
# offset = distance between original filter 'edge' and the top of the filter
# 2* offset = gap between filters 
offsets = [0,2,4,8,16,32,64] 

def generateFilter(filterEdges, taper, offset):
    # generate 'filter only'
    lsst_filter = {}
    # Generate the actual 'filter' for this set of offset/taper values.
    i = 0
    for f in filterlist:
        wavelen = [300.0, filterEdges[i]-taper+offset, filterEdges[i]+offset, 
                   filterEdges[i+1]-offset, filterEdges[i+1]-offset+taper, 1200.0]
        sb = [0, 0, 0.99, 0.99, 0, 0]
        wavelen = numpy.array(wavelen)
        sb = numpy.array(sb)
        lsst_filter[f] = Bandpass.Bandpass(wavelen=wavelen, sb=sb)
        i = i + 1
    filterBpSet = BandpassSet.BandpassSet()
    filterBpSet.setBandpassSet(lsst_filter, filterlist)
    return filterBpSet

def generateFullBandpass(filterBpSet):
    # include atmosphere, detector, etc.
    lsst_stat = BandpassSet.BandpassSet()
    rootdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
    complist = ('detector.dat', 'lens1.dat', 'lens2.dat', 'lens3.dat', 'm1.dat', 'm2.dat', 'm3.dat', 'atmos.dat')
    lsst_stat.setThroughputs_ComponentFiles(filterlist=filterlist, all_filter_complist = complist,
                                            rootdir = rootdir, verbose=False)
    bpSet = lsst_stat.multiplyBandpassSets(filterBpSet)
    return bpSet

def calcColorDistanceStats(sedSet1, sedSet2):
    colors_for_distance = ('ug', 'gr')
    dmin = 100
    dave = 0
    dmax = -100
    dnum = 0
    for i in range(len(sedSet1.sedlist)):
        for j in range(len(sedSet2.sedlist)):
            d2 = 0
            for c in colors_for_distance:
                d2 = d2 + (sedSet1.colors[c][i] - sedSet2.colors[c][j])**2
            d2 = numpy.sqrt(d2)
            dave = dave + d2
            dnum = dnum + 1
            if d2 < dmin:
                dmin = d2
            if d2 > dmax:
                dmax = d2            
    dave = dave / float(dnum)
    if dmin==100:
        dmin=0
    if dmax==-100:
        dmax=0
    return dmin, dave, dmax


if __name__ == "__main__":

    # Read in data for SEDS.
    kurucz = SedSet.SedSet(rootdir= "/Users/rhiannonjones/seds/kurucz_r", sedtype='kurucz', verbose=False)
    whitedwarf = SedSet.SedSet(rootdir= "/Users/rhiannonjones/seds/white_dwarfs_r/H",
                               sedtype='white_dwarf', verbose=False)
    whitedwarfHe = SedSet.SedSet(rootdir= "/Users/rhiannonjones/seds/white_dwarfs_r/He",
                               sedtype='white_dwarf', verbose=False)
    quasars = SedSet.SedSet(rootdir="/Users/rhiannonjones/seds/quasar", sedtype='quasar')
    quasars.redshiftSEDS(redshiftlim=[0, 2.51], redshiftstep=0.1)
    print("Read %d kurucz models, %d whitedwarf models, and %d quasars" %(len(kurucz.sedlist), len(whitedwarf.sedlist),
                                                                          len(quasars.sedlist)))
    
    savefigs=True
    figformat='png'
    # loop through versions of the filters.
    #offsets=[0, 32]
    #tapers = [0, 16, 32]
    dcolormin = numpy.zeros((len(offsets), len(tapers)), dtype='float')
    dcolormax = numpy.zeros((len(offsets), len(tapers)), dtype='float')
    dcolorave = numpy.zeros((len(offsets), len(tapers)), dtype='float')
    metgrid = numpy.zeros((len(offsets), len(tapers)), dtype='float')
    o = 0
    for offset in offsets:
        t = 0
        for taper in tapers:
            figroot = "filters_%do_%dt" %(offset, taper)
            # Set up filter only part.
            bpSet = generateFilter(filterEdges, taper, offset)
            # Add in other components.
            bpSet = generateFullBandpass(bpSet)
            # Plot this filter set.
            bpSet.plotFilters(atmos=False, rootdir=os.getenv("LSST_THROUGHPUTS_DEFAULT"),
                              title = "Offset : %f  Taper : %f" %(offset, taper),
                              filter_tags='normal', ploteffsb=False, savefig=savefigs, figroot=figroot)
            # Calculate colors in this bandpassSet for the stars, white dwarfs and quasars.
            kurucz.calcAllColors(bpSet.bandpass, bpSet.filterlist)
            whitedwarf.calcAllColors(bpSet.bandpass, bpSet.filterlist)
            whitedwarfHe.calcAllColors(bpSet.bandpass, bpSet.filterlist)
            quasars.calcAllColors(bpSet.bandpass, bpSet.filterlist)
            # QUASAR STAR separation
            # Plot separation in color-color space, u-g/g-r plot.
            xcolor = 'ug'
            ycolor = 'gr'
            kurucz.plotColorColor(xcolor, ycolor, ptcolor='r', linestyle='.', withlines=False,
                                  newfig=True, savefig=False)
            whitedwarf.plotColorColor(xcolor, ycolor, ptcolor='k', linestyle='.', withlines=False,
                                      newfig=False, savefig=False)
            whitedwarfHe.plotColorColor(xcolor, ycolor, ptcolor='k', linestyle='.', withlines=False,
                                        newfig=False, savefig=False)            
            quasars.plotColorColor(xcolor, ycolor, ptcolor='g', linestyle='.', withlines=True,
                                   newfig=False, savefig=False)
            pylab.xlim(-0.5, 2.5)
            pylab.ylim(-0.5, 1.5)
            pylab.title("Offset : %f  Taper : %f" %(offset, taper))
            pylab.legend(numpoints=1, loc='lower right')
            if savefigs:
                pylab.savefig(figroot+"_colors." + figformat, format=figformat)                
            # Calculate minimum distance between each of these populations (for separation)
            dmin1, dave1, dmax1 = calcColorDistanceStats(kurucz, quasars)
            dmin2, dave2, dmax2 = calcColorDistanceStats(whitedwarf, quasars)
            dmin3, dave3, dmax3 = calcColorDistanceStats(whitedwarfHe, quasars)
            dcolormin[o][t] = min(dmin1, dmin2, dmin3)
            dcolormax[o][t] = max(dmax1, dmax2, dmax3)
            dcolorave[o][t] = min(dave1, dave2, dmax3)
            # METALLICITY determination
            # select stars from kurucz models with same g-r color
            pylab.figure()
            grbins = numpy.arange(-0.3, 1.5, 0.5)
            colors = ['b', 'g', 'r', 'c', 'm']
            ci = 0
            metslope = {}
            for grcol in grbins:
                condition = (abs(kurucz.colors['gr']-grcol)<0.2)
                # pull out ug and metallicity values for the stars with this color
                ug = kurucz.colors['ug'][condition]
                mets = kurucz.data['FeH'][condition]
                if len(ug) < 2:
                    metslope[grcol] = 0
                    continue
                pylab.plot(mets, ug, colors[ci]+'x', label='%.2f' %(grcol))
                # Fit least squares line to metallicity/ug values
                fitfunc = lambda p, x: p[0] + p[1]*x # target function
                errfunc = lambda p, x, y: fitfunc(p, x) - y  # distance to target function
                p0 = [0, 0] # initial guess for parameters
                p1, success = scipy.optimize.leastsq(errfunc, p0, args=(mets, ug))
                x = numpy.arange(mets.min(), mets.max(), 0.1)
                y = fitfunc(p1, x)
                metslope[grcol] = p1[1]
                pylab.plot(x, y, color=colors[ci])
                ci = ci+1
            # find grbin closest to 0.3, which is what we want to use to save the color info
            dmin = 5
            for grcol in grbins:
                if (abs(grcol - 0.3)<dmin):
                    dmin = abs(grcol-0.3)
                    grbin = grcol
            metgrid[o][t] = metslope[grcol]
            pylab.legend(numpoints=1, loc=(0.8, 0.3))
            pylab.xlabel("Metallicity")
            pylab.ylabel("ug")
            pylab.ylim(-1, 4)
            pylab.title("Offset : %f  Taper : %f" %(offset, taper))
            if savefigs:
                pylab.savefig(figroot+"_met." + figformat, format=figformat)                
            t = t +1
        o = o + 1

    # Create attempts at summary plots for quasar-star separation
    pylab.figure()
    pylab.imshow(dcolormin, interpolation='nearest',
                 extent=([min(tapers), max(tapers), max(offsets), min(offsets)]))
    pylab.colorbar()
    pylab.xlabel("Taper")
    pylab.ylabel("Offset")
    pylab.title("Color distance, minimum")
    if savefigs:
        pylab.savefig("grid_colordist_min." + figformat, format=figformat)                
    pylab.figure()
    pylab.imshow(dcolormax, interpolation='nearest', aspect='equal',
                 extent=([min(tapers), max(tapers), max(offsets), min(offsets)]))
    pylab.colorbar()
    pylab.xlabel("Taper")
    pylab.ylabel("Offset")
    pylab.title("Color distance, maximum")
    if savefigs:
        pylab.savefig("grid_colordist_max." + figformat, format=figformat)                
    pylab.figure()
    pylab.imshow(dcolorave, interpolation='nearest',
                 extent=([min(tapers), max(tapers), max(offsets), min(offsets)]))
    pylab.colorbar()
    pylab.xlabel("Taper")
    pylab.ylabel("Offset")
    pylab.title("Color distance, average")
    if savefigs:
        pylab.savefig("grid_colordist_ave." + figformat, format=figformat)
    pylab.figure()
    pylab.imshow(metgrid, interpolation='nearest',
                 extent=([min(tapers), max(tapers), max(offsets), min(offsets)]))
    pylab.colorbar()
    pylab.xlabel("Taper")
    pylab.ylabel("Offset")
    pylab.title("Delta(ug)/Delta(FeH) at gr=0.3")
    if savefigs:
        pylab.savefig("grid_metgrid." + figformat, format=figformat)                


    
    #pylab.show()
