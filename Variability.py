import numpy
import linecache
import math
import os
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

class Variability(object):
    """Variability class for adding temporal variation to the magnitudes of
    objects in the base catalog.  All methods should return a dictionary of
    magnitude offsets.
    """

    def __init__(self, baseDir=None, cache=False):
        self.lcCache = {}
        self.cache = cache
        if baseDir is None:
            try:
                self.datadir = os.path.join(os.environ.get("CAT_SHARE_DATA"),"data","LightCurves")
            except:
                raise("No directory specified and $CAT_SHARE_DATA is undefined")
            

    def applyStdPeriodic(self, params, keymap, expmjd, inPeriod=None,
            inDays=True, interpFactory=None):
        expmjd = numpy.asarray(expmjd)
        filename = params[keymap['filename']]
        toff = float(params[keymap['t0']])
        epoch = expmjd - toff
        if self.lcCache.has_key(filename):
            splines = self.lcCache[filename]['splines']
            period = self.lcCache[filename]['period']
        else:
            print filename
            lc = numpy.loadtxt(self.datadir+"/"+filename, unpack=True, comments='#')
            if inPeriod is None:
                dt = lc[0][1] - lc[0][0]
                period = lc[0][-1] + dt
            else:
                period = inPeriod
                
            if inDays:    
                lc[0] /= period

            splines  = {}
            if interpFactory is not None:
                splines['u'] = interpFactory(lc[0], lc[1])
                splines['g'] = interpFactory(lc[0], lc[2])
                splines['r'] = interpFactory(lc[0], lc[3])
                splines['i'] = interpFactory(lc[0], lc[4])
                splines['z'] = interpFactory(lc[0], lc[5])
                splines['y'] = interpFactory(lc[0], lc[6])
                if self.cache:
                    self.lcCache[filename] = {'splines':splines, 'period':period}
            else:
                splines['u'] = interp1d(lc[0], lc[1])
                splines['g'] = interp1d(lc[0], lc[2])
                splines['r'] = interp1d(lc[0], lc[3])
                splines['i'] = interp1d(lc[0], lc[4])
                splines['z'] = interp1d(lc[0], lc[5])
                splines['y'] = interp1d(lc[0], lc[6])
                if self.cache:
                    self.lcCache[filename] = {'splines':splines, 'period':period}

        phase = epoch/period - epoch//period
        magoff = {}
        for k in splines:
            magoff[k] = splines[k](phase)
        return magoff

    def applyMflare(self, params, expmjd):
        params['lcfilename'] = "mflare/"+params['lcfilename'][:-5]+"1.dat"
        keymap = {'filename':'lcfilename', 't0':'t0'}
        magoff = self.applyStdPeriodic(params, keymap, expmjd, inPeriod=params['length'])
        for k in magoff.keys():
            magoff[k] = -magoff[k]
        return magoff

    def applyRRly(self, params, expmjd):
        keymap = {'filename':'filename', 't0':'tStartMjd'}
        return self.applyStdPeriodic(params, keymap, expmjd,
                interpFactory=InterpolatedUnivariateSpline)

    def applyEb(self, params, expmjd):
        keymap = {'filename':'lcfile', 't0':'t0'}
        dMags = self.applyStdPeriodic(params, keymap, expmjd, 
                interpFactory=InterpolatedUnivariateSpline)
        for k in dMags.keys():
            dMags[k] = -2.5*numpy.log10(dMags[k])
        return dMags

    def applyMicrolens(self, params, expmjd):
        expmjd = numpy.asarray(expmjd)
        epochs = expmjd - params['t0']
        dMags = {}
        u = numpy.sqrt(params['umin']**2 + ((2.0*epochs/params['that'])**2))
        magnification = (u**2+2.0)/(u*numpy.sqrt(u**2+4.0))
        dmag = -2.5*numpy.log10(magnification)
        dMags['u'] = dmag
        dMags['g'] = dmag
        dMags['r'] = dmag
        dMags['i'] = dmag
        dMags['z'] = dmag
        dMags['y'] = dmag
        return dMags


    def applyAgn(self, params, expmjd):
        dMags = {}
        expmjd = numpy.asarray(expmjd)
        toff = params['t0_mjd']
        seed = int(params['seed'])
        sfint = {}
        sfint['u'] = params['agn_sfu']
        sfint['g'] = params['agn_sfg']
        sfint['r'] = params['agn_sfr']
        sfint['i'] = params['agn_sfi']
        sfint['z'] = params['agn_sfz']
        sfint['y'] = params['agn_sfy']
        tau = params['agn_tau']
        epochs = expmjd - toff
        if epochs.min() < 0:
            raise("WARNING: Time offset greater than minimum epoch.  Not applying variability")
        endepoch = epochs.max()

        dt = tau/100.
        nbins = int(math.ceil(endepoch/dt))
        dt = (endepoch/nbins)/tau
        sdt = math.sqrt(dt)
        numpy.random.seed(seed=seed)
        es = numpy.random.normal(0., 1., nbins)
        for k in sfint.keys():
            dx = numpy.zeros(nbins+1)
            dx[0] = 0.
            for i in range(nbins):
                #The second term differs from Zeljko's equation by sqrt(2.)
                #because he assumes stdev = sfint/sqrt(2)
                dx[i+1] = -dx[i]*dt + sfint[k]*es[i]*sdt + dx[i]
            x = numpy.linspace(0, endepoch, nbins+1)
            intdx = interp1d(x, dx)
            magoff = intdx(epochs)
            dMags[k] = magoff
        return dMags

    def applyMicrolensing(self, params, expmjd):
        dMags = {}
        epochs = numpy.asarray(expmjd) - params['t0']
        u = numpy.sqrt(params['umin']**2 + (2.0 * epochs /\
            params['that'])**2)
        magnification = (u + 2.0) / (u * numpy.sqrt(u**2 + 4.0))
        dmag = -2.5 * numpy.log10(magnification)
        dMags['u'] = dmag
        dMags['g'] = dmag
        dMags['r'] = dmag
        dMags['i'] = dmag
        dMags['z'] = dmag
        dMags['y'] = dmag 
        return dMags

    def applyAmcvn(self, params, expmjd):
        dMag = {}
        epochs = numpy.asarray(expmjd)
        # get the light curve of the typical variability
        uLc   = params['amplitude']*numpy.cos((epochs - params['t0'])/params['period'])
        gLc   = uLc
        rLc   = uLc
        iLc   = uLc
        zLc   = uLc
        yLc   = uLc

        # add in the flux from any bursting
        if params['does_burst']:
            adds = np.zeros(len(epochs))
            for o in np.linspace(params['t0'] + params['burst_freq'],\
                                 params['t0'] + maxyears*365.25, \
                                 np.ceil(maxyears*365.25/params['burst_freq'])):
                tmp = np.exp( -1*(epochs - o)/params['burst_scale'])/np.exp(-1.)
                adds -= params['amp_burst']*tmp*(tmp < 1.0)  ## kill the contribution 

            ## add some blue excess during the outburst
            uLc += adds +  2.0*params['color_excess_during_burst']*adds/min(adds)
            gLc += adds + params['color_excess_during_burst']*adds/min(adds)
            rLc += adds + 0.5*params['color_excess_during_burst']*adds/min(adds)
            iLc += adds
            zLc += adds
            yLc += adds
                      
        self.dMag['u'] = uLc
        self.dMag['g'] = gLc
        self.dMag['r'] = rLc
        self.dMag['i'] = iLc
        self.dMag['z'] = zLc
        self.dMag['y'] = yLc
        return dMag

    def applyBHMicrolens(self, params, expmjd):
        expmjd = numpy.asarray(expmjd)
        filename = params['filename']
        toff = float(params['t0'])
        epoch = expmjd - toff
        lc = numpy.loadtxt(self.datadir+"/"+filename, unpack=True, comments='#')
        dt = lc[0][1] - lc[0][0]
        period = lc[0][-1]
        #BH lightcurves are in years
        lc[0] *= 365.
        minage = lc[0][0]
        maxage = lc[0][-1]
        #I'm assuming that these are all single point sources lensed by a
        #black hole.  These also can be used to simulate binary systems.
        #Should be 8kpc away at least.
        splines  = {}
        magnification = InterpolatedUnivariateSpline(lc[0], lc[1])

        magoff = {}
        moff = []
        for ep in epoch:
            if ep < minage or ep > maxage:
                moff.append(1.)
            else:
                moff.append(magnification(ep))
        moff = numpy.asarray(moff)
        for k in ['u','g','r','i','z','y']:
            magoff[k] = -2.5*numpy.log(moff)
        return magoff
