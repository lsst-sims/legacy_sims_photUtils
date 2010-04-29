"""
  Written by: Lynne Jones - UW - 4/28/10
  Questions or comments, email : ljones.uw@gmail.com

Similarly to the explanation given in BandpassSet, the main 
point of this class is to provide a way to group many seds 
together (into one object which has some useful methods: like 
calculating magnitudes for all seds or colors, or plotting the 
color-color diagram, or plotting all of the seds, or just 
reading all of the seds in a particular directory and also handing you
a list of the dictionary keys). 

There are several functions to read a bunch of seds from a particular 
directory: 
 (the init is the primary method to get to these)
 __init__, readSeds, readKurucz, readWhiteDwarfs, readRedDwarfs, 
 readAGBStars, readAsteroids, readLymanGals, readGalaxies, readQuasar. 
Then there is an interesting method to redshift all of a set of base
seds you have (say if you have a bunch of galaxy seds and then need
to replicate these at many redshifts .. use redshiftSEDS. 

Then things for plotting - plotSeds, normalizeSEDS is useful for helping
put things on the same scale to be plotted. 

calcEffObjLambdasDict - calculate effective wavelengths for each of your
seds in a BandpassSet (ie effective wavelength for each sed, for each filter).

calcMagArray - calculate magnitudes for each of your seds in one bandpass,
but returns a numpy array of the magnitudes. 
calcMagsDict - similar to above, but returns a magnitude dictionary (keyed 
to sed name). 

The calcColor/calcDeltaMag array and dict methods are similar to the mag
Array/Dict methods above, but for colors. 

Then some useful shortcuts to plotting color-color diagrams: 
plotColorColor / plotLSSTAllColor (subpanels appropriate for LSST filters)
plotSDSSAllColor (subpanels appropriate for SDSS filters). 

"""


import os
import numpy as n
import pylab as pyl
import Bandpass
import Sed
import BandpassSet

figformat = 'png'

class SedSet:
    """ Set up a dictionary of a bunch of seds."""
    def __init__(self, sedlist=None, rootdir='seds', type=None, verbose=True):
        """Initialize Seds with in list sedlist, in directory dir with root name root"""
        if sedlist == None:
            if rootdir!='':
                rootdir = rootdir + "/"
            if verbose:
                print "No files specified."
                print "Will try to read all seds known to be of type %s from %s%s" \
                    %(type, rootdir, type)
            if type == 'kurucz':
                rootdir = rootdir
                self.readKurucz(rootdir=rootdir, verbose=verbose)
            elif (type == 'white_dwarf_H') | (type=='wdH'):
                rootdir = rootdir + 'white_dwarfs/H'
                self.readWhiteDwarfs(rootdir=rootdir, verbose=verbose)
            elif (type == 'white_dwarf_He') | (type=='wdHe'):
                rootdir = rootdir + 'white_dwarfs/He'
                self.readWhiteDwarfs(rootdir=rootdir,  verbose=verbose)
            elif type == 'reddwarf':
                rootdir = rootdir + 'mlt'
                self.readRedDwarfs(rootdir=rootdir, verbose=verbose)
            elif type == 'AGB':
                rootdir = rootdir + 'AGB'
                self.readAGBStars(rootdir=rootdir, verbose=verbose)
            elif (type=='asteroid') | (type=='asteroids'):
                rootdir = rootdir + 'asteroids'
                self.readAsteroids(rootdir=rootdir, verbose=verbose)                
            elif ((type == 'lyman_gal') | (type=='lyman')) :
                rootdir = rootdir + 'lyman_gal'
                self.readLymanGals(rootdir=rootdir, verbose=verbose)
            elif (type=='galaxy') | (type=='galaxies'):
                rootdir = rootdir + 'galaxies'
                self.readGalaxies(rootdir=rootdir,  verbose=verbose)
            elif (type=='quasar'):
                rootdir= rootdir + 'quasar'
                self.readQuasars(rootdir=rootdir, alpha=0.0, verbose=verbose)
            else :
                print "Will just try to read everything in directory %s" %(rootdir)
                sedlist = os.listdir(rootdir)
                self.sedlist, self.seds = self.readSeds(sedlist, rootdir=rootdir, verbose=verbose)
        else:
            self.type = type
            self.sedlist, self.seds = self.readSeds(sedlist, rootdir=rootdir, verbose=verbose)
        return

    def readSeds(self, sedlist, rootdir='seds', verbose=True):
        """Read sim object sed files specified in sedlist, returns sedlist and seds."""
        # Read the files, set up the sed dictionary.
        seds = {}
        if rootdir != '':
            rootdir = rootdir + "/"
        for sedname in sedlist:
            filename = rootdir  + sedname
            # read seds from sedlist in rootdir
            if verbose:
                print "Reading sed from %s" %(filename)
            seds[sedname] = Sed.Sed()
            seds[sedname].readSED_flambda(filename)
        return sedlist, seds

    def readKurucz(self, rootdir='', verbose=True):
        """ Initialize all kurucz models seds"""
        # not all metallicities and temperatures have corresponding kurucz models 
        tmplist = os.listdir(rootdir)
        kuruczlist = []
        # check for files that match the 'kurucz' filename template k(m10)_Temp_g*.fits[_Temp]
        # not guaranteed 100% but should be okay
        for tmp in tmplist:
            if (tmp[:1] == 'k') & (tmp[:5][-1:] == '_'):
                kuruczlist.append(tmp)
        self.sedlist, self.seds = self.readSeds(kuruczlist, rootdir=rootdir, verbose=verbose)
        self.type = 'kurucz'
        return
    
    def readWhiteDwarfs(self, rootdir='', verbose=True):
        """Initialize all white dwarf seds"""
        tmplist = os.listdir(rootdir)
        wdlist = []
        # Look for files which match 'bergeron_*'
        for tmp in tmplist:
            if (tmp[:8]=='bergeron'):
                wdlist.append(tmp)
        self.sedlist, self.seds = self.readSeds(wdlist, rootdir=rootdir, verbose=verbose)
        self.type = 'white_dwarf'
        return
        
    def readRedDwarfs(self, rootdir='', verbose=True):
        """Initialize all cool red MLT dwarf seds"""
        sedlist=[]
        mdwarflist = []
        ldwarflist = []
        bdwarflist = []
        tmplist = os.listdir(rootdir)
        for tmp in tmplist:
            if tmp[:1]=='m':
                mdwarflist.append(tmp)
            if (tmp[:1]=='l') | (tmp[:1]=="L"):
                ldwarflist.append(tmp)
            if tmp[:1]=='t':
                tdwarflist.append(tmp)
            if tmp[:4]=='burr':
                bdwarflist.append(tmp)
        sedlist = mdwarflist + ldwarflist  + bdwarflist
        self.sedlist, self.seds = self.readSeds(sedlist, rootdir=rootdir,verbose=verbose)
        self.type = 'reddwarf'
        self.mdwarflist = mdwarflist
        self.ldwarflist = ldwarflist
        self.bdwarflist = bdwarflist
        return 

    def readAGBStars(self, rootdir='', verbose=True):
        """ Initialize all AGB stars """
        sedlist = []
        Cstarlist = []
        Ostarlist = []
        tmplist = os.listdir(rootdir)
        for tmp in tmplist:
            if tmp[-4:]=='.dat':
                sedlist.append(tmp)
                if tmp[:2]=='C_':
                    Cstarlist.append(tmp)
                if tmp[:2]=='O_':
                    Ostarlist.append(tmp)
        self.sedlist, self.seds = self.readSeds(sedlist, rootdir=rootdir, verbose=verbose)
        self.type='AGB'
        self.Cstarlist = Cstarlist
        self.Ostarlist = Ostarlist
        return

    def readAsteroids(self, rootdir='', verbose=True):
        """ Initialize all asteroid seds"""
        sedlist = []
        tmplist = os.listdir(rootdir)
        for tmp in tmplist:
            if tmp[-4:]=='.dat':
                sedlist.append(tmp)
        self.sedlist, self.seds = self.readSeds(sedlist, rootdir=rootdir, verbose=verbose)
        self.type='asteroid'
        return

    def readLymanGals(self, rootdir='',  verbose=True):
        """Initialize lyman break galaxy seds"""
        # there are only a few of these, so just specify here
        lymanlist = ('spec1.dat', 'spec2.dat', 'spec3.dat', 'spec4.dat')
        self.sedlist, self.seds = self.readSeds(lymanlist, rootdir=rootdir, verbose=verbose)
        self.type = 'lyman_gal'
        return

    def readGalaxies(self, rootdir='', verbose=True):
        """Initialize galaxy seds."""
        galaxylist = []
        tmplist = os.listdir(rootdir)
        for tmp in tmplist:
            if tmp[-4:] == "spec":
                galaxylist.append(tmp)
        self.sedlist, self.seds = self.readSeds(galaxylist, rootdir=rootdir, verbose=verbose)
        self.type = 'galaxy'
        return

    def readQuasars(self, alpha=0.0, rootdir='', verbose=True):
        """ Initialize quasar seds from QSO composite spectra """
        seds = {}
        if rootdir != '':
            rootdir = rootdir + "/"
        filename = rootdir + 'quasar.sed'
        quasar = Sed.Sed()
        # NOTE NEED TO FIX Quasar sed (not nm, not flambda)
        quasar.readSED_fnu(filename) 
        quasar.wavelen = quasar.wavelen/10.0
        # add 'color' variation
        quasar.flambda = quasar.flambda* n.power((quasar.wavelen/4000), alpha)
        quasar.flambdatofnu()
        seds = {}
        sedlist = ('quasar%.2f' %(alpha),)
        self.sedlist = sedlist
        seds[sedlist[0]] = quasar
        self.seds = seds
        self.type = 'quasar'
        return

    def redshiftSEDS(self, redshiftlim=(0, 8), redshiftstep=0.1, dimming=False):
        """Add redshifted versions of all base seds.

        Key seds by sedname, but sedname has redshift added (origname_%.2f) %(redshift)"""
        # now add redshifted galaxy seds
        if verbose:
            print "Now adding redshifts to these galaxies"
        redshifts = n.arange(redshiftlim[0], redshiftlim[1]+redshiftstep, redshiftstep)
        newsedlist = []
        basesedlist = self.sedlist
        newseds = {}
        for basesedname in basesedlist:
            basesed = self.seds[basesedname]
            for redshift in redshifts:
                sedname = basesed + "_%.2f" %(redshift)
                galaxylist.append(sedname)
                wavelen, flambda = basesed.redshiftSED(redshift, wavelen=basesed.wavelen, 
                                                       flambda=basesed.flambda)
                newseds[sedname] = Sed.Sed(wavelen, flambda)
        self.redshifts = redshifts
        self.sedlist = newsedlist
        self.seds = newseds
        return

    def plotSeds(self, sedlist, fnu=False, flambda=True, lambdaflambda=False, 
                 xlim=(300, 1100), ylimfnu=None, ylimflambda=None, ylabel=None, 
                 sed_xtags=0.68, sed_ytags=0.85, linestyle='-' , 
                 newfig=True, savefig=False, figroot='seds',):
        """ Plot the sed fnu and/or flambdas, with limits xlim/ylimfnu/ylimflambda """        
        try:
            seds = self.seds
        except AttributeError:
            print "There don't seem to be any seds defined."
            return
        # note that sedlist doesn't have to be the full set of seds in self (b/c passed in)
        # set up colors for plot output
        colors = ('c', 'k', 'r', 'y', 'g', 'b', 'm')
        if (fnu):
            if newfig:
                pyl.figure()
            # plot fnus
            colorindex = 0
            for sedname in sedlist:
                color = colors[colorindex]
                colorindex = colorindex+1
                if colorindex == len(colors):
                    colorindex = 0 
                pyl.plot(seds[sedname].wavelen, seds[sedname].fnu, color+linestyle)
            # add names to filter throughputs
            # x location specified by sed_xtags (can be turned off by setting sed_xtags to None)
            colorindex = 0
            ytag = sed_ytags
            if sed_xtags != None:
                for sedname in sedlist:
                    color = colors[colorindex]
                    colorindex = colorindex + 1
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.figtext(sed_xtags, ytag, sedname, color=color, va='top')
                    ytag = ytag - 0.04
            # set x/y limits
            pyl.xlim(xmin=xlim[0], xmax=xlim[1])
            if ylimfnu != None:
                pyl.ylim(ymin=ylimfnu[0], ymax=ylimfnu[1])
            pyl.xlabel("Wavelength (A)")
            if ylabel!= None:
                if ylabel=='auto':
                    ylabel = "F_nu (Jansky)  "
                pyl.ylabel(ylabel)
            if savefig:
                figname = figroot + "_fnu.png"
                pyl.savefig(figname, format=figformat)   
        if (flambda):
            if newfig:
                pyl.figure()
            # plot flambda
            colorindex = 0
            for sedname in sedlist:
                color = colors[colorindex]
                colorindex = colorindex+1
                if colorindex == len(colors):
                    colorindex = 0 
                if lambdaflambda:
                    pyl.plot(seds[sedname].wavelen, seds[sedname].flambda*seds[sedname].wavelen, 
                             color+linestyle)
                else:
                    pyl.plot(seds[sedname].wavelen, seds[sedname].flambda, color+linestyle)
            # add names to filter throughputs
            # x location specified by sed_xtags (can be turned off by setting sed_xtags to None)
            colorindex = 0
            if sed_xtags != None:
                ytag = sed_ytags
                for sedname in sedlist:
                    color = colors[colorindex]
                    colorindex = colorindex + 1
                    if colorindex == len(colors):
                        colorindex = 0
                    pyl.figtext(sed_xtags, ytag, sedname, color=color, va='top')
                    ytag = ytag - 0.04
            # set x/y limits
            pyl.xlim(xmin=xlim[0], xmax=xlim[1])
            if ylimflambda != None:
                pyl.ylim(ymin=ylimflambda[0], ymax=ylimflambda[1])
            pyl.xlabel("Wavelength (A)")
            if ylabel!= None:
                if (ylabel=='F_nu (Jansky)  ') | (ylabel=='auto'):
                    ylabel="F_lambda (ergs/s/cm^2/A)"
                pyl.ylabel(ylabel)
            if savefig:
                figname = figroot + "_flambda.png"
                pyl.savefig(figname, format=figformat)      
        return
   
    def normalizeSEDS(self, mags, bandpass):
        """ Given a set of seds, calculate and apply normalization constants.

        so that they have the magnitudes given in mags, in filter 'bandpass'
        mags can be one number or a numpy array of len equal to number of seds"""
        # check if mags is an array or a single number
        try:
            magnum = len(mags)
        except TypeError: # single number
            magnum = 0
        # check length of magnitude list, return if not single number or sedlist len
        if ((magnum != 0 ) & (magnum!=len(self.sedlist))):
            print "'mags' must be a single number or have length equal to sedlist"
            return
        seds = self.seds
        i = 0
        for sedname in self.sedlist:
            if magnum ==0:
                magmatch = mags
            else:
                magmatch = mags[i]
                i = i + 1
            norm = self.seds[sedname].calcFluxNorm(magmatch, bandpass)
            self.seds[sedname].multiplyFluxNorm(norm)
        return

    def calcEffObjLambdasDict(self, bandpassSet):
        """Calculate the effective wavelengths for these seds, in each of the filters.
        
        Returns dictionary of effective wavelengths for each SED (keyed by sedname) and keys"""
        effobjlambda = {}
        for sedname in self.sedlist:
            effobjlambda[sedname] = {}
            if self.seds[sedname].fnu.sum() == 0:
                # then there was nothing in this data file
                for filter in bandpassSet.filterlist:
                    effobjlambda[sedname][filter] = 0
                continue
            for filter in bandpassSet.filterlist:
                tmp1 = (bandpassSet.filters[filter].wavelen * bandpassSet.filters[filter].phi 
                        * self.seds[sedname].fnu )
                tmp2 = bandpassSet.filters[filter].phi * self.seds[sedname].fnu
                effobjlambda[sedname][filter] = tmp1.sum() / tmp2.sum()
        return effobjlambda, self.sedlist

    def calcMagsDict(self, bandpassSet, verbose=False):
        """ Given a set of seds, calculate the magnitudes using a set of filters"""
        """ Returns magnitude dictionary (keyed sed / bandpass)  """
        mags = {}
        if verbose:
            # print out magnitudes header
            printstring = "#      "
            for filter in bandpassSet.filterlist:
                printstring = printstring + "  "  + filter + " "
            print printstring
        for sedname in self.sedlist:
            mags[sedname] = {}
            for filter in bandpassSet.filterlist:
                mags[sedname][filter] = self.seds[sedname].calcMag(bandpassSet.bandpass[filter])
            if verbose:
                # print out magnitudes
                printstring = sedname + " "
                for filter in bandpassSet.filterlist:
                    printstring = printstring + " %0.3f" %(mags[sedname][filter])
                print printstring
        self.mags = mags
        return mags, bandpassSet.filterlist

    def calcColorsDict(self, bandpassSet, colorfilter='r', verbose=False):
        """ Calculates colors, relative to colfilter (i.e. g-r, u-r...) of seds.
        Returns color dictionary (keyed sed/color) """
        try:
            self.mags
        except AttributeError:
            self.calcMagsDict(bandpassSet, verbose=False)
        mags = self.mags
        # set up list of colors (keys)
        colorlist = []
        for filter in bandpassSet.filterlist:
            color = filter + colorfilter
            colorlist.append(color)
        # calculate colors for all seds
        colormags = {}
        # if verbose, print colorheader
        if verbose:
            printstring = "#     "
            for color in colorlist:
                printstring = printstring + " " + color + "  "
            print printstring
        for sedname in self.sedlist:
            colormags[sedname] = {}
            for filter in bandpassSet.filterlist:
                color = filter+colorfilter
                colormags[sedname][color] = mags[sedname][filter] - mags[sedname][colorfilter]
            if verbose:
                # print out colors 
                printstring = sedname + " "
                for color in colorlist:
                    printstring = printstring+  " %0.5f" %(colormags[sedname][color])
                print printstring
        self.colormags = colormags
        return colormags, colorlist

    def calcMagArray(self, bandpass):
        """ Calculate magnitudes in bandpass and returns numpy array"""
        mags = n.empty(len(self.sedlist), float)
        i = 0
        for sedname in self.sedlist:
            mags[i] = self.seds[sedname].calcMag(bandpass)
            i = i +1
        return mags

    def calcDeltaMagArray(self, sedlist, bandpass1, bandpass2):
        """ Calculates the difference in magnitudes between filter (bandpass) 1 and 2 """
        """  Returns numpy array """
        deltamags = n.empty(len(sedlist), float)
        i = 0
        for sedname in sedlist:
            deltamags[i] = (self.seds[sedname].calcMag(bandpass1) - 
                            self.seds[sedname].calcMag(bandpass2))
            i = i +1
        return deltamags
                         

    def plotColorColor(self, sedlist, bandpass1, bandpass2, bandpass3, bandpass4, 
                       color='r', linestyle='.', withlines=True, xlim=None, ylim=None, 
                       xlabel=None, ylabel=None, title=None, grid=False, 
                       leg_xtag=0.6, leg_ytag=0.8, leg_text="", 
                       newfig=True, savefig=False, figroot='color-color'):
        """ Plot color-color diagram for self """
        if newfig:
            pyl.figure()
        # calculate the values for the color color plot 
        # bandpass1 - bandpass2 = color1 = X axis
        # bandpass3 - bandpass4 = color2 = Y axis
        color1 = self.calcDeltaMagArray(sedlist, bandpass1, bandpass2)
        color2 = self.calcDeltaMagArray(sedlist, bandpass3, bandpass4)
        color = color
        if withlines:
            pyl.plot(color1, color2, color)
        pyl.plot(color1, color2, color+linestyle)
        pyl.figtext(leg_xtag, leg_ytag, leg_text, color=color)
        if xlim!=None:
            pyl.xlim(xlim[0], xlim[1])
        if ylim!=None:
            pyl.ylim(ylim[0], ylim[1])
        if grid:
            pyl.grid()
        if xlabel!=None:
            pyl.xlabel(xlabel)
        if ylabel!=None:
            pyl.ylabel(ylabel)
        if title!=None:
            pyl.title(title)
        if savefig:
            figname = figroot + '.' + figformat
            pyl.savefig(figname, format=figformat)
        return

    def plotLSSTAllColor(self, sedlist, lsst, color='r', linestyle='.', 
                         newfig=True, savefig=False, figroot='allcolor'):
        """ Plot all colors, expecting LSST bands"""
        if newfig: 
            pyl.figure()
        # assume the filters are lsst filters and want g-r/u-g, r-i/g-r, i-z/r-i, z-y/i-z        
        pyl.subplot(221)
        self.plotColorColor(sedlist, lsst.bandpass['u'], lsst.bandpass['g'], 
                            lsst.bandpass['g'], lsst.bandpass['r'], 
                            xlabel='u-g', ylabel='g-r', color=color, linestyle=linestyle, 
                            withlines=False, grid=True, newfig=False, savefig=False)
        pyl.subplot(222)
        self.plotColorColor(sedlist, lsst.bandpass['g'], lsst.bandpass['r'], 
                            lsst.bandpass['r'], lsst.bandpass['i'], 
                            xlabel='g-r', ylabel='r-i', color=color, linestyle=linestyle, 
                            withlines=False, grid=True, newfig=False, savefig=False)
        pyl.subplot(223)
        self.plotColorColor(sedlist, lsst.bandpass['r'], lsst.bandpass['i'], 
                            lsst.bandpass['i'], lsst.bandpass['z'], 
                            xlabel='r-i', ylabel='i-z', color=color, linestyle=linestyle, 
                            withlines=False, grid=True, newfig=False, savefig=False)
        pyl.subplot(224)
        self.plotColorColor(sedlist, lsst.bandpass['i'], lsst.bandpass['z'], 
                            lsst.bandpass['z'], lsst.bandpass['y'], 
                            xlabel='i-z', ylabel='z-y', color=color, linestyle=linestyle, 
                            withlines=False, grid=True, newfig=False, savefig=False)
        if savefig:
            figname = figroot + '.' + figformat
            pyl.savefig(figname, format=figformat)
        return

    def plotSDSSAllColor(self, sedlist, sdss, color='r', linestyle='.', 
                         newfig=True, savefig=False, figroot='allcolor'):
        """ Plot all colors, expecting SDSS bands"""
        if newfig: 
            pyl.figure()
        # assume the filters are sdss filters and want g-r/u-g, r-i/g-r, i-z/r-i
        pyl.subplot(221)
        self.plotColorColor(sedlist, sdss.bandpass['u'], sdss.bandpass['g'], 
                            sdss.bandpass['g'], sdss.bandpass['r'], 
                            xlabel='u-g', ylabel='g-r', color=color, 
                            linestyle=linestyle, withlines=False, grid=True, 
                            newfig=False, savefig=False)
        pyl.subplot(222)
        self.plotColorColor(sedlist, sdss.bandpass['g'], sdss.bandpass['r'],
                            sdss.bandpass['r'], sdss.bandpass['i'], 
                            xlabel='g-r', ylabel='r-i', color=color, linestyle=linestyle, 
                            withlines=False, grid=True, newfig=False, savefig=False)
        pyl.subplot(223)
        self.plotColorColor(sedlist, sdss.bandpass['r'], sdss.bandpass['i'],
                            sdss.bandpass['i'], sdss.bandpass['z'], 
                            xlabel='r-i', ylabel='i-z', color=color, linestyle=linestyle, 
                            withlines=False, grid=True, newfig=False, savefig=False)
        if savefig:
            figname = figroot + '.' + figformat
            pyl.savefig(figname, format=figformat)
        return
