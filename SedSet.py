"""
  Questions or comments, email : ljones.uw@gmail.com

$Id$

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
 __init__, readSeds,

Then there is an interesting method to redshift all of a set of base
seds you have (say if you have a bunch of galaxy seds and then need
to replicate these at many redshifts .. use redshiftSEDS. 

Then things for plotting - plotSeds, normalizeSEDS is useful for helping
put things on the same scale to be plotted. 

calcEffLambdasSed - calculate effective wavelengths for each of your
seds in a BandpassSet (ie effective wavelength for each sed, for each filter).

calcMag - calculate magnitudes for each of your seds in one bandpass,
but returns a numpy array of the magnitudes. 

calcDeltaMag - similar to above, but give it two bandpasses so it can calculate colors.

calcAllMags - calculate magnitudes for each of the SEDS, in all bandpasses in a bandpassDictionary.
calcAllColors - similar to above, but calculates adjacent colors for all bandpasses in the dictionary.

"""


import os
import copy
import numpy as n
import pylab as pyl
import Bandpass
import Sed
import photUtils
import BandpassSet

figformat = 'png'

class SedSet:
    """ Set up a dictionary of a bunch of seds."""
    def __init__(self, sedlist=None, rootdir='.', sedtype=None, verbose=True):
        """Initialize Seds.

        Give root directory for SED location (rootdir). Will read all files from that directory unless you - 
        Optional: give list of SEDS from that directory (sedlist).
        Optional: give 'type' of SEDS, which will attempt to pull more identifing information (T/FeH/g, etc.), and
                  also constrain the list of SEDS from the rootdir."""
        # Check and set sedtype.
        sedtypes = ('white_dwarf', 'kurucz', 'red_dwarf', 'asteroid', 'AGB', 'galaxy', 'quasar')
        if sedtype not in sedtypes:
            print "sedtype is not recognized: defaulting to 'None'. Please use one of %s" %(sedtypes)
            self.sedtype = None
        else:
            self.sedtype = sedtype
        # Comment on source directory.
        if verbose:
            print "Reading files from %s" %(rootdir)
        # Set up sedlist. 
        if sedlist == None:
            if verbose:
                print "No files specified. Will scan directory."
            self.sedlist = os.listdir(rootdir)
        else:
            self.sedlist = sedlist
        if self.sedtype != None:                
            if verbose:
                print "Using sedtype of %s."
            self.setupSedType(rootdir)
        # Read seds.
        if self.sedtype == 'quasar':
            if verbose:
                print "Quasar is a special case and has already read in SED in setupSedType."
        else:
            self.readSeds(rootdir=rootdir, verbose=verbose)
        # Convert data into numpy arrays so masking/conditions can work. 
        for key in self.data.keys():
            self.data[key] = n.array(self.data[key])
        return

    def readSeds(self, rootdir, verbose=True):
        """Read sim object sed files specified in self.sedlist, removing invalid sed files from list."""
        # Set up the dictionary and good sed list.
        self.seds = {}
        tmpsed = Sed.Sed()
        validsedlist = []
        i = 0 
        # Attempt to read seds from sedlist.
        for sedname in self.sedlist:
            filename = os.path.join(rootdir, sedname)
            if verbose:
                print "Reading sed from %s" %(filename)
            try:
                tmpsed.readSED_flambda(filename)
                validsedlist.append(sedname)
                # Store the SED in the self.seds dictionary. 
                self.seds[sedname] = Sed.Sed(wavelen=tmpsed.wavelen, flambda=tmpsed.flambda)
                # Set up fnu while we're here.
                self.seds[sedname].flambdaTofnu()
                i = i + 1
            except ValueError:
                print "Oops - file %s was not a valid sedfile. Removing from sedlist and other data info."
                del self.sedlist[i]
                for key in self.data.keys():
                    del self.data[key][i]
        self.sedlist = validsedlist                
        return 

    def setupSedType(self, rootdir):
        """Use self.sedlist to trim sedlist to only seds known to be of this 'type' and gather
        additional data appropriate to type (such as Teff/FeH/logg/asteroid type)."""
        self.data = {}
        # Kurucz
        if self.sedtype == 'kurucz':
            kuruczlist = []
            self.data['Teff'] = []
            self.data['FeH'] = []
            self.data['logg'] = []
            # Check for files that match the 'kurucz' filename template k(m10)_Temp.fits_g(40)[_Temp]
            # And gather relevant information.
            for sedname in self.sedlist:
                tmp = sedname.split('_')
                if (tmp[0][0] == 'k'):
                    kuruczlist.append(sedname)
                    # Teff will either be tmp[1] or tmp[3]
                    if len(tmp) > 3:
                        self.data['Teff'].append(int(tmp[3]))
                    else:
                        self.data['Teff'].append(int(tmp[1].split('.')[0]))
                    # logg should be tmp[2], minus the 'g'
                    self.data['logg'].append(float(tmp[2][1:])/10.0)
                    # FeH should be tmp[0], after accounting for 'm' or 'p'
                    metsign = 1
                    if tmp[0][1:][:1] == 'm':
                        metsign = -1
                    self.data['FeH'].append(metsign*float(tmp[0][2:])/10.0)
            self.sedlist = kuruczlist
        # White Dwarf
        elif self.sedtype == 'white_dwarf':
            wdlist = []
            self.data['Teff'] = []
            self.data['logg'] = []
            # Look for files which match 'bergeron_Temp_logg.dat[_Temp]'
            for sedname in self.sedlist:
                tmp = sedname.split('_')
                if (tmp[0]=='bergeron'):
                    wdlist.append(sedname)
                    if tmp[1] == "He":
                        if len(tmp)>4:
                            self.data['Teff'].append(int(tmp[4]))
                        else:
                            self.data['Teff'].append(int(tmp[2]))
                        self.data['logg'].append(tmp[3].split('.')[0])
                    else:                        
                        if len(tmp)>3:
                            self.data['Teff'].append(int(tmp[3]))
                        else:
                            self.data['Teff'].append(int(tmp[1]))
                        self.data['logg'].append(float(tmp[2].split('.')[0])/10.0)
            self.sedlist = wdlist
        # Red dwarf (mlty stars)
        elif self.sedtype == 'red_dwarf':
            self.data['mlist'] = []
            self.data['llist'] = []
            self.data['tlist'] = []
            self.data['blist'] = []
            for sedname in self.sedlist:
                if sedname[:1]=='m':
                    self.data['mlist'].append(sedname)
                if (sedname[:1]=='l') | (sedname[:1]=="L"):
                    self.data['llist'].append(sedname)
                if sedname[:1]=='t':
                    self.data['tlist'].append(sedname)
                if sedname[:4]=='burr':
                    self.data['blist'].append(sedname)
            self.sedlist = self.data['mlist'] + self.data['llist'] + self.data['tlist'] + self.data['blist']
        # AGB stars
        elif self.sedtype == 'AGB':
            self.data['Clist'] = []
            self.data['Olist']= []
            for sedname in self.sedlist:
                if sedname[:2]=='C_':
                    self.data['Clist'].append(tmp)
                if sedname[:2]=='O_':
                    self.data['Olist'].append(tmp)
            self.sedlist = self.data['Clist'] + self.data['Olist']
        # Asteroids
        elif self.sedtype == 'asteroid':
            sedlist = []
            self.data['taxtype'] = []
            for sedname in self.sedlist:
                if sedname[-4:]=='.dat':
                    sedlist.append(sedname)
                self.data['taxtype'].append(sedname.split('.')[0])
            self.sedlist = sedlist
        # Galaxies
        elif self.sedtype == 'galaxy':
            galaxylist = []
            self.data['galtype'] = []
            self.data['age'] = []
            for sedname in self.sedlist:
                tmp = sedname.split('.')
                if tmp[-1:] == "spec":
                    galaxylist.append(sedname)
                    self.data['galtype'] = tmp[0]
                    self.data['age'] = tmp[3]
            self.sedlist = galaxylist
        # Quasars
        elif self.sedtype == 'quasar':
            # This is a special case .. there is only one quasar spectrum, but we can 'shift' the slope.
            quasar = Sed.Sed()
            quasar.readSED_flambda(os.path.join(rootdir, 'quasar.dat'))
            self.data['alpha'] = n.arange(-0.5, 0.6, 0.5, dtype='float')
            self.seds = {}
            self.sedlist = []
            for alpha in self.data['alpha']:
                dictkey = "q_" + "%.1f" %(alpha)
                self.sedlist.append(dictkey)
                self.seds[dictkey] = Sed.Sed(wavelen=quasar.wavelen, flambda=quasar.flambda)
                # Add 'color' variation.
                self.seds[dictkey].flambda = self.seds[dictkey].flambda * \
                                             n.power(self.seds[dictkey].wavelen/400.0, alpha)
        return

    def redshiftSEDS(self, redshiftlim=(0, 8), redshiftstep=0.2, dimming=False, verbose=True):
        """Add redshifted versions of all base seds.

        Key seds by sedname, but sedname has redshift added (origname_%.2f) %(redshift)"""
        # This is really only applicable to galaxies or SN or quasars. 
        if verbose:
            print "Now adding redshifts to these seds"
        redshifts = n.arange(redshiftlim[0], redshiftlim[1]+redshiftstep, redshiftstep)
        newsedlist = []
        basesedlist = self.sedlist
        newseds = {}
        newalpha = []
        newredshift = []
        i = 0
        for basesedname in basesedlist:
            basesed = self.seds[basesedname]
            alpha = self.data['alpha'][i]
            i = i + 1
            for redshift in redshifts:
                sedname = basesedname + "_%.2f" %(redshift)
                newsedlist.append(sedname)
                newalpha.append(alpha)
                newredshift.append(redshift)
                wavelen, flambda = basesed.redshiftSED(redshift, wavelen=basesed.wavelen, 
                                                       flambda=basesed.flambda)
                newseds[sedname] = Sed.Sed(wavelen=wavelen, flambda=flambda)
                newseds[sedname].flambdaTofnu()
        self.data['redshift'] = newredshift
        self.data['alpha'] = newalpha
        self.basesedlist = basesedlist
        self.sedlist = newsedlist
        self.seds = newseds
        return

    def plotSeds(self, sedlist=None, fnu=False, flambda=True, lambdaflambda=False, 
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
        if sedlist == None:
            sedlist = self.sedlist
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
            # add names to fnu values
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
        """ Given a set of seds, calculate and apply offsets so all seds have the same mag in this bandpass.

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
            fluxnorm = self.seds[sedname].calcFluxNorm(magmatch, bandpass)
            self.seds[sedname].multiplyFluxNorm(fluxnorm)
            self.seds[sedname].flambdaTofnu()
        return

    def calcEffLambdasSed(self, bandpassDict, filterlist, verbose=True):
        """Calculate the effective wavelengths for these seds, in each of the filters.

        Give this a dictionary of bandpasses and the list of filter names. 
        Sets efflambda_sed dictionary (key: sed, then filter) for self."""
        efflambda_sed = {}
        for sedname in self.sedlist:
            efflambda_sed[sedname] = {}
            if self.seds[sedname].fnu.sum() == 0:
                # then there was nothing in this data file
                for f in filterlist:
                    efflambda_sed[sedname][f] = 0
                continue
            for f in filterlist:
                tmp1 = (bandpassDict[f].wavelen * bandpassDict[f].phi 
                        * self.seds[sedname].fnu )
                tmp2 = bandpassDict[f].phi * self.seds[sedname].fnu
                efflambda_sed[sedname][f] = tmp1.sum() / tmp2.sum()
        self.efflambda_sed = efflambda_sed
        if verbose:
            writestring = "#SEDNAME efflambda:  "
            for f in filterlist:
                writestring = writestring + " %s " %(f)
            print writestring
            for sedname in self.sedlist:
                writestring = sedname
                for f in filterlist:
                    writestring = writestring + " %.2f " %(self.efflambda_sed[sedname][f])
            print writestring
        return 

    def calcMags(self, bandpass):
        """ Calculate magnitudes in bandpass and returns numpy array"""
        # Please note - this is not the *fastest* way to do this calculation, but it is
        #  easy to read and quick to write like this. But look at manyMagsCalc versions in Sed and Bandpass.
        mags = n.empty(len(self.sedlist), float)
        i = 0
        for sedname in self.sedlist:
            mags[i] = self.seds[sedname].calcMag(bandpass)
            i = i +1
        return mags    

    def calcDeltaMag(self, bandpass1, bandpass2):
        """ Calculates the difference in magnitudes between filter (bandpass) 1 and 2 """
        """  Returns numpy array """
        deltamags = n.empty(len(self.sedlist), float)
        i = 0
        for sedname in self.sedlist:
            deltamags[i] = (self.seds[sedname].calcMag(bandpass1) - 
                            self.seds[sedname].calcMag(bandpass2))
            i = i +1
        return deltamags
                         
    def calcAllMags(self, bandpassDict, filterlist, verbose=False):
        """Calculate magnitudes for all seds in all bandpasses. Stores result in self as dictionary:
                 mags[filter]=array."""
        # Set up storage dictionary for resulting magnitudes.
        self.mags = {}
        for f in filterlist:
            self.mags[f] = n.zeros(len(self.sedlist), dtype='float')
        # First set up bandpass list (in filterlist order).
        phiarray, dlambda = photUtils.setupPhiArray_dict(bandpassDict, filterlist)
        # Now get set to resample SEDs onto the bandpass wavelength grid.
        wavelen_min = bandpassDict[filterlist[0]].wavelen.min()
        wavelen_max = bandpassDict[filterlist[0]].wavelen.max()
        wavelen_step = bandpassDict[filterlist[0]].wavelen[1] - bandpassDict[filterlist[0]].wavelen[0]
        i = 0
        for sedname in self.sedlist:
            # Set up for fast mag calculation. Use copy of sed because we'll be resampling. 
            tmpsed = copy.deepcopy(self.seds[sedname])
            # Resample and set up fnu.
            tmpsed.synchronizeSED(wavelen_min=wavelen_min, wavelen_max=wavelen_max,
                                  wavelen_step=wavelen_step)
            # Calculate the magnitudes. 
            tmpmags = tmpsed.manyMagCalc(phiarray, dlambda)
            # Assign list of magnitudes returned to individual elements. 
            j = 0
            for f in filterlist:
                self.mags[f][i] = tmpmags[j]
                j = j + 1
            i = i + 1        
        if verbose:
            writestring = "# SEDname "
            for f in filterlist:
                writestring = writestring + " %s " %(f)
            print writestring
            for i in range(len(self.sedlist)):
                writestring = self.sedlist[i]
                for  f in filterlist:
                    writestring = writestring + " %.3f" %(self.mags[f][i])
                print writestring
        return
    
    def calcAllColors(self, bandpassDict, filterlist, verbose=False):
        # First calculate all the magnitudes.
        self.calcAllMags(bandpassDict, filterlist, verbose=False)
        # Calculate all the colors.
        self.colors = {}
        i = 0
        colornames = []
        for i in range(len(filterlist)-1):
            colorname = filterlist[i]+filterlist[i+1]
            colornames.append(colorname)
            self.colors[colorname] = self.mags[filterlist[i]] - self.mags[filterlist[i+1]]
        self.colornames = colornames
        if verbose:
            writestring = "# SEDname "
            for f in filterlist:
                writestring = writestring + " %s " %(f)
            print writestring
            for i in range(len(self.sedlist)):
                writestring = self.sedlist[i]
                for c in colornames:
                    writestring = writestring + " %.3f" %(self.colors[c][i])
                print writestring
        return

    def plotColorColor(self, xcolor, ycolor, ptcolor='r', linestyle='.', withlines=True,
                       xlim=None, ylim=None, newfig=True, savefig=False, figroot='colorcolor'):
        if newfig:
            pyl.figure()
        try:
            self.colors[xcolor]
        except AttributeError:
            raise Exception("You haven't calculated this %s color yet - calculate and add to colors dictionary." \
                            %(xcolor))
        try:
            self.colors[ycolor]
        except AttributeError:
            raise Exception("You haven't calculated this %s color yet - calculate and add to colors dictionary." \
                            %(ycolor))
        if withlines:
            # figure out what to connect.
            if self.sedtype == 'kurucz':
                # connect stars with the same metallicity, same logg but different temperatures.
                mets = n.unique(self.data['FeH'])
                loggs = n.unique(self.data['logg'])
                for met in mets:
                    for logg in loggs:
                        condition = ((self.data['FeH'] == met) & (self.data['logg']==logg))
                        pyl.plot(self.colors[xcolor][condition], self.colors[ycolor][condition], ptcolor+"-")
            elif self.sedtype == 'quasar':
                # connect stars with the same alpha, but different redshifts
                alphas = n.unique(self.data['alpha'])
                for alpha in alphas:
                    condition = (self.data['alpha'] == alpha)
                    pyl.plot(self.colors[xcolor][condition], self.colors[ycolor][condition], ptcolor+"-")
            elif self.sedtype == 'white_dwarf':
                # connect stars with the same logg but different temperatures
                loggs = n.unique(self.data['logg'])
                for logg in loggs:
                    condition = (self.data['logg'] == logg)
                    pyl.plot(self.colors[xcolor][condition], self.colors[ycolor][condition],  ptcolor+"-")
            else:
                pyl.plot(self.colors[xcolor], self.colors[ycolor], color=ptcolor, linestyle='-')
        pyl.plot(self.colors[xcolor], self.colors[ycolor], ptcolor+linestyle, label=self.sedtype)
        pyl.xlabel(xcolor)
        pyl.ylabel(ycolor)
        if savefig:
            figname = figroot + '.' + figformat
            pyl.savefig(figname, format=figformat)
        return

