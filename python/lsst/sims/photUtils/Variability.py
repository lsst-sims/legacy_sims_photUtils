import numpy
import linecache
import math
import os
import inspect
import json as json
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from lsst.sims.catalogs.measures.instance import compound

def variabilityRegistration(ff):
    """
    A decorator to indicate which methods need to be added to the
    register of variability methods
    """

    def decoratedFunction(*args,**kwargs):
        args[0].logIt = True
        return ff(*args,**kwargs)
    
    return decoratedFunction


class Variability(object):
    """Variability class for adding temporal variation to the magnitudes of
    objects in the base catalog.  All methods should return a dictionary of
    magnitude offsets.
    """
    
    variabilityInitialized = False
    variabilityMethods = {}
    
    def initializeVariability(self,doCache=False):
        self.variabilityInitialized=True        
        listOfMembers=inspect.getmembers(self)
        for methodName, actualMethod in listOfMembers:
            self.logIt = False
            if methodName != "initializeVariability": 
            #so that we don't recursively call this function
                try:
                    actualMethod()
                except:
                    pass
                
                if self.logIt == True:
                    self.variabilityMethods[methodName] = actualMethod
            
        #below are variables to cache the light curves of variability models
        self.variabilityLcCache = {}
        self.variabilityCache = doCache
        try:
            self.variabilityDataDir = os.path.join(os.environ.get("CAT_SHARE_DATA"),"data","LightCurves")
        except:
            print "No directory specified and $CAT_SHARE_DATA is undefined"
            raise
    
    
    @compound('lsst_u_var','lsst_g_var','lsst_r_var','lsst_i_var',
    'lsst_z_var','lsst_y_var')
    def get_stellar_variability(self):
        uu = self.column_by_name('lsst_u')
        gg = self.column_by_name('lsst_g')
        rr = self.column_by_name('lsst_r')
        ii = self.column_by_name('lsst_i')
        zz = self.column_by_name('lsst_z')
        yy = self.column_by_name('lsst_y')
        
        varParams = self.column_by_name('varParamStr')
        
        uuout = []
        ggout = []
        rrout = []
        iiout = []
        zzout = []
        yyout = []
        
        i=0
        for vv in varParams:
            if vv != numpy.unicode_("None"):
                deltaMag = self.applyVariability(vv)
                uuout.append(uu[i]+deltaMag['u'])
                ggout.append(gg[i]+deltaMag['g'])
                rrout.append(rr[i]+deltaMag['r'])
                iiout.append(ii[i]+deltaMag['i'])
                zzout.append(zz[i]+deltaMag['z'])
                yyout.append(yy[i]+deltaMag['y'])
            else:
                uuout.append(uu[i])
                ggout.append(gg[i])
                rrout.append(rr[i])
                iiout.append(ii[i])
                zzout.append(zz[i])
                yyout.append(yy[i])
            i+=1
            
        return numpy.array([uuout,ggout,rrout,iiout,zzout,yyout])        
    
    @compound('uRecalc_var', 'gRecalc_var', 'rRecalc_var', 'iRecalc_var',
          'zRecalc_var', 'yRecalc_var',
          'uAgn_var', 'gAgn_var', 'rAgn_var', 'iAgn_var', 'zAgn_var', 'yAgn_var')
    def get_galaxy_variability(self):
        
        uTotal = self.column_by_name("uRecalc")
        gTotal = self.column_by_name("gRecalc")
        rTotal = self.column_by_name("rRecalc")
        iTotal = self.column_by_name("iRecalc")
        zTotal = self.column_by_name("zRecalc")
        yTotal = self.column_by_name("yRecalc")
        
        uBulge = self.column_by_name("uBulge")
        gBulge = self.column_by_name("gBulge")
        rBulge = self.column_by_name("rBulge")
        iBulge = self.column_by_name("iBulge")
        zBulge = self.column_by_name("zBulge")
        yBulge = self.column_by_name("yBulge")
        
        uDisk = self.column_by_name("uDisk")
        gDisk = self.column_by_name("gDisk")
        rDisk = self.column_by_name("rDisk")
        iDisk = self.column_by_name("iDisk")
        zDisk = self.column_by_name("zDisk")
        yDisk = self.column_by_name("yDisk")
        
        uAgn = self.column_by_name("uAgn")
        gAgn = self.column_by_name("gAgn")
        rAgn = self.column_by_name("rAgn")
        iAgn = self.column_by_name("iAgn")
        zAgn = self.column_by_name("zAgn")
        yAgn = self.column_by_name("yAgn")
        
        varParams = self.column_by_name("varParamStr")
        
        uTotalOut = []
        gTotalOut = []
        rTotalOut = []
        iTotalOut = []
        zTotalOut = []
        yTotalOut = []
        
        uAgnOut = []
        gAgnOut = []
        rAgnOut = []
        iAgnOut = []
        zAgnOut = []
        yAgnOut = []
        
        i=0
        for vv in varParams:
            if vv != numpy.unicode_("None"):           
                deltaMag=self.applyVariability(vv)
                uAgnOut.append(uAgn[i]+deltaMag['u'])
                gAgnOut.append(gAgn[i]+deltaMag['g'])
                rAgnOut.append(rAgn[i]+deltaMag['r'])
                iAgnOut.append(iAgn[i]+deltaMag['i'])
                zAgnOut.append(zAgn[i]+deltaMag['z'])
                yAgnOut.append(yAgn[i]+deltaMag['y'])
                
                uTotalOut.append(self.sum_magnitudes(disk = uDisk[i], bulge = uBulge[i],
                        agn = uAgnOut[i]))
                
                gTotalOut.append(self.sum_magnitudes(disk = gDisk[i], bulge = gBulge[i],
                        agn = gAgnOut[i]))
                        
                rTotalOut.append(self.sum_magnitudes(disk = rDisk[i], bulge = rBulge[i],
                        agn = rAgnOut[i]))
                
                iTotalOut.append(self.sum_magnitudes(disk = iDisk[i], bulge = iBulge[i],
                        agn = iAgnOut[i]))
                        
                zTotalOut.append(self.sum_magnitudes(disk = zDisk[i], bulge = zBulge[i],
                        agn = zAgnOut[i]))
                
                yTotalOut.append(self.sum_magnitudes(disk = yDisk[i], bulge = yBulge[i],
                        agn = yAgnOut[i]))
            
            else:
                uTotalOut.append(uTotal[i])
                gTotalOut.append(gTotal[i])
                rTotalOut.append(rTotal[i])
                iTotalOut.append(iTotal[i])
                zTotalOut.append(zTotal[i])
                yTotalOut.append(yTotal[i])
                
                uAgnOut.append(uAgn[i])
                gAgnOut.append(gAgn[i])
                rAgnOut.append(rAgn[i])
                iAgnOut.append(iAgn[i])
                zAgnOut.append(zAgn[i])
                yAgnOut.append(yAgn[i])
            
            i+=1
        
        return numpy.array([uTotalOut,gTotalOut,rTotalOut,iTotalOut,zTotalOut,yTotalOut,\
                           uAgnOut,gAgnOut,rAgnOut,iAgnOut,zAgnOut,yAgnOut])
        
    
    def applyVariability(self, varParams):
        """
        varParams will be the varParamStr column from the data base
        
        This method uses json to convert that into a machine-readable object
        
        it uses the varMethodName to select the correct variability method from the
        dict self.variabilityMethods
        
        it uses then feeds the pars array to that method, under the assumption
        that the parameters needed by the method can be found therein
        """
        if self.variabilityInitialized == False:
            self.initializeVariability()
            
        varCmd = json.loads(varParams)
        method = varCmd['varMethodName']
        params = varCmd['pars']
        #expmjd=numpy.asarray(self.obs_metadata.metadata['Opsim_expmjd'][0],dtype=float)
        expmjd=self.obs_metadata.mjd
        output = self.variabilityMethods[method](params,expmjd)
        return output
    
    def applyStdPeriodic(self, params, keymap, expmjd, inPeriod=None,
            inDays=True, interpFactory=None):
        expmjd = numpy.asarray(expmjd)
        #expmjd = numpy.asarray(self.obs_metadata.metadata['Opsim_expmjd'][0],dtype=float)
        
        filename = params[keymap['filename']]
        toff = float(params[keymap['t0']])
        epoch = expmjd - toff
        if self.variabilityLcCache.has_key(filename):
            splines = self.variabilityLcCache[filename]['splines']
            period = self.variabilityLcCache[filename]['period']
        else:
            lc = numpy.loadtxt(self.variabilityDataDir+"/"+filename, unpack=True, comments='#')
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
                if self.variabilityCache:
                    self.variabilityLcCache[filename] = {'splines':splines, 'period':period}
            else:
                splines['u'] = interp1d(lc[0], lc[1])
                splines['g'] = interp1d(lc[0], lc[2])
                splines['r'] = interp1d(lc[0], lc[3])
                splines['i'] = interp1d(lc[0], lc[4])
                splines['z'] = interp1d(lc[0], lc[5])
                splines['y'] = interp1d(lc[0], lc[6])
                if self.variabilityCache:
                    self.variabilityLcCache[filename] = {'splines':splines, 'period':period}

        phase = epoch/period - epoch//period
        magoff = {}
        for k in splines:
            magoff[k] = splines[k](phase)
        return magoff

    @variabilityRegistration
    def applyMflare(self, params, expmjd):
        
        params['lcfilename'] = "mflare/"+params['lcfilename'][:-5]+"1.dat"
        keymap = {'filename':'lcfilename', 't0':'t0'}
        magoff = self.applyStdPeriodic(params, keymap, expmjd, inPeriod=params['length'])
        for k in magoff.keys():
            magoff[k] = -magoff[k]
        return magoff

    @variabilityRegistration
    def applyRRly(self, params, expmjd):
    
        keymap = {'filename':'filename', 't0':'tStartMjd'}
        return self.applyStdPeriodic(params, keymap, expmjd,
                interpFactory=InterpolatedUnivariateSpline)
    
    @variabilityRegistration
    def applyCepheid(self, params, expmjd):
    
        keymap = {'filename':'lcfile', 't0':'t0'}
        return self.applyStdPeriodic(params, keymap, expmjd, inPeriod=params['period'], inDays=False,
                interpFactory=InterpolatedUnivariateSpline)

    @variabilityRegistration
    def applyEb(self, params, expmjd):
        keymap = {'filename':'lcfile', 't0':'t0'}
        dMags = self.applyStdPeriodic(params, keymap, expmjd,
                interpFactory=InterpolatedUnivariateSpline)
        for k in dMags.keys():
            dMags[k] = -2.5*numpy.log10(dMags[k])
        return dMags

    @variabilityRegistration
    def applyMicrolens(self, params, expmjd_in):
        #expmjd = numpy.asarray(self.obs_metadata.metadata['Opsim_expmjd'][0],dtype=float)
        expmjd = numpy.asarray(expmjd_in,dtype=float)
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


    @variabilityRegistration
    def applyAgn(self, params, expmjd_in):
        dMags = {}
        #expmjd = numpy.asarray(self.obs_metadata.metadata['Opsim_expmjd'][0],dtype=float)
        expmjd = numpy.asarray(expmjd_in,dtype=float)
        toff = numpy.float(params['t0_mjd'])
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

    @variabilityRegistration
    def applyMicrolensing(self, params, expmjd_in):
        dMags = {}
        #epochs = numpy.asarray(self.obs_metadata.metadata['Opsim_expmjd'][0],dtype=float) - params['t0']
        epochs = numpy.asarray(expmjd_in,dtype=float) - params['t0']
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

    @variabilityRegistration
    def applyAmcvn(self, params, expmjd_in):
        maxyears = 10.
        dMag = {}
        #epochs = numpy.asarray(self.obs_metadata.metadata['Opsim_expmjd'][0],dtype=float)
        epochs = numpy.asarray(expmjd_in,dtype=float)
        # get the light curve of the typical variability
        uLc   = params['amplitude']*numpy.cos((epochs - params['t0'])/params['period'])
        gLc   = uLc
        rLc   = uLc
        iLc   = uLc
        zLc   = uLc
        yLc   = uLc

        # add in the flux from any bursting
        if params['does_burst']:
            adds = numpy.zeros(epochs.size)
            for o in numpy.linspace(params['t0'] + params['burst_freq'],\
                                 params['t0'] + maxyears*365.25, \
                                 numpy.ceil(maxyears*365.25/params['burst_freq'])):
                tmp = numpy.exp( -1*(epochs - o)/params['burst_scale'])/numpy.exp(-1.)
                adds -= params['amp_burst']*tmp*(tmp < 1.0)  ## kill the contribution 
            ## add some blue excess during the outburst
            uLc += adds +  2.0*params['color_excess_during_burst']*adds/min(adds)
            gLc += adds + params['color_excess_during_burst']*adds/min(adds)
            rLc += adds + 0.5*params['color_excess_during_burst']*adds/min(adds)
            iLc += adds
            zLc += adds
            yLc += adds
                      
        dMag['u'] = uLc
        dMag['g'] = gLc
        dMag['r'] = rLc
        dMag['i'] = iLc
        dMag['z'] = zLc
        dMag['y'] = yLc
        return dMag

    @variabilityRegistration
    def applyBHMicrolens(self, params, expmjd_in):
        #expmjd = numpy.asarray(self.obs_metadata.metadata['Opsim_expmjd'][0],dtype=float)
        expmjd = numpy.asarray(expmjd_in,dtype=float)
        filename = params['filename']
        toff = float(params['t0'])
        epoch = expmjd - toff
        lc = numpy.loadtxt(self.variabilityDataDir+"/"+filename, unpack=True, comments='#')
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
        if expmjd.size == 1:
            epoch = [epoch]
        for ep in epoch:
            if ep < minage or ep > maxage:
                moff.append(1.)
            else:
                moff.append(magnification(ep))
        moff = numpy.asarray(moff)
        for k in ['u','g','r','i','z','y']:
            magoff[k] = -2.5*numpy.log(moff)
        return magoff
