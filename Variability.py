import numpy
import math
import os
from scipy.interpolate import InterpolatedUnivariateSpline

class Variability(object):
    """Variability class for adding temporal variation to the magnitudes of
    objects in the base catalog.  All methods should return a dictionary of
    magnitude offsets.
    """

    def __init__(self, baseDir=None, cache=False):
        self.rrlyLc = {}
        self.cache = cache
        if baseDir is None:
            try:
                self.datadir = os.path.join(os.environ.get("CAT_SHARE_DATA"),"data","LightCurves")
            except:
                raise("No directory specified and $CAT_SHARE_DATA is undefined")
            

    def applyRRly(self, params, expmjd):
        filename = params['filename']
        toff = params['tStartMjd']
        epoch = expmjd - toff
        if self.rrlyLc.has_key(filename):
            splines = self.rrlyLc[filename]['splines']
            period = self.rrlyLc[filename]['period']
        else:
            lc = numpy.loadtxt(self.datadir+"/"+filename, unpack=True, comments='#')
            dt = lc[0][1] - lc[0][0]
            period = lc[0][-1] + dt
            lc[0] /= period
            splines  = {}
            splines['u'] = InterpolatedUnivariateSpline(lc[0], lc[1])
            splines['g'] = InterpolatedUnivariateSpline(lc[0], lc[2])
            splines['r'] = InterpolatedUnivariateSpline(lc[0], lc[3])
            splines['i'] = InterpolatedUnivariateSpline(lc[0], lc[4])
            splines['z'] = InterpolatedUnivariateSpline(lc[0], lc[5])
            splines['y'] = InterpolatedUnivariateSpline(lc[0], lc[6])
            if self.cache:
                self.rrlyLc[filename] = {'splines':splines, 'period':period}
        phase = epoch/period - epoch//period
        magoff = {}
        for k in splines:
            magoff[k] = splines[k](phase)
        return magoff

    def applyAgn(self, params, expmjd):
        #This assumes monochromatic variability.  Should we go the extra step
        #and at least use different sf's?
        toff = params['t0_mjd']
        seed = params['seed']
        sfint = params['agn_sfr']
        tau = params['agn_tau']
        epoch = expmjd - toff

        dt = tau/100.
        nbins = math.ceil(epoch/dt)
        dt = (epoch/nbins)/tau
        sdt = math.sqrt(dt)
        s2 = math.sqrt(2.)

        es = numpy.random.normal(0., 1., nbins)
        dx = 0.
        for e in es:
            dx = -dx*dt + s2*sfint*e*sdt
        return {'u':dx, 'g':dx, 'r':dx, 'i':dx, 'z':dx, 'y':dx}



