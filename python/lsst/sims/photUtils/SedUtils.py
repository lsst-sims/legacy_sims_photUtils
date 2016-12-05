import numpy as np
from lsst.sims.photUtils import Bandpass


__all__ = ["getImsimFluxNorm"]


def getImsimFluxNorm(sed, magmatch):
    """
    Calculate the flux normalization of an SED in the imsim bandpass.

    Parameters:
    -----------
    sed is the SED to be normalized

    magmatch is the desired magnitude in the imsim bandpass

    Output
    ------
    The factor by which the flux of sed needs to be multiplied to achieve
    the desired magnitude.
    """

    # This method works based on the assumption that the imsim bandpass
    # is a delta function.  If that ever ceases to be true, the unit test
    # testSedUtils.py, which checks that the results of this method are
    # identical to calling Sed.calcFluxNorm and passing in the imsim bandpass,
    # will fail and we will know to modify this method.

    if not hasattr(getImsimFluxNorm, 'imsim_wavelen'):
        bp = Bandpass()
        bp.imsimBandpass()
        non_zero_dex = np.where(bp.sb>0.0)[0][0]
        getImsimFluxNorm.imsim_wavelen = bp.wavelen[non_zero_dex]

    if sed.fnu is None:
        sed.flambdaTofnu()

    mag = -2.5*np.log10(np.interp(getImsimFluxNorm.imsim_wavelen, sed.wavelen, sed.fnu)) - sed.zp
    dmag = magmatch - mag
    return np.power(10, (-0.4*dmag))