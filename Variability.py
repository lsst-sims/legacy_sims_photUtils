import numpy
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
        self.rrlyLc = {}
        self.cache = cache
        if baseDir is None:
            try:
                self.datadir = os.path.join(os.environ.get("CAT_SHARE_DATA"),"data","LightCurves")
            except:
                raise("No directory specified and $CAT_SHARE_DATA is undefined")
            

    def applyRRly(self, params, expmjd):
        expmjd = numpy.asarray(expmjd)
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
        dMags = {}
        expmjd = numpy.asarray(expmjd)
        toff = params['t0_mjd']
        seed = params['seed']
        sfint = {}
        sfint['u'] = params['agn_sfu']
        sfint['g'] = params['agn_sfg']
        sfint['r'] = params['agn_sfr']
        sfint['i'] = params['agn_sfi']
        sfint['z'] = params['agn_sfz']
        sfint['y'] = params['agn_sfy']
        tau = params['agn_tau']
        epochs = expmjd - toff
        endepoch = epochs.max()

        dt = tau/100.
        nbins = math.ceil(endepoch/dt)
        dt = (endepoch/nbins)/tau
        sdt = math.sqrt(dt)
        s2 = math.sqrt(2.)

        es = numpy.random.normal(0., 1., nbins)
        for k in sfint.keys():
            dx = numpy.zeros(nbins+1)
            dx[0] = 0.
            for i in range(nbins):
                dx[i+1] = -dx[i]*dt + s2*sfint[k]*es[i]*sdt
            x = numpy.linspace(0, endepoch, nbins+1)
            intdx = interp1d(x, dx)
            magoff = intdx(epochs)
            dMags[k] = magoff
        return dMags



