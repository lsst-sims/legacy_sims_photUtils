# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 14:07:03 2015

@author: Bryce Kalmbach
"""
import numpy as np

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Bandpass import Bandpass

class rgUtils():

    """
    This class is designed to provide methods that will be useful to both selectStarSED and selectGalaxySED.
    """

    def calcMagNorm(self, objectMags, sedObj, photObj, redshift = None, stepSize = 0.01, initBand = 0):

        """
        This will find the magNorm value that gives the closest match to the magnitudes of the object
        using the matched SED.

        @param [in] objectMags are the magnitude values for the object with extinction matching that of
        the SED object. In the normal case using the selectSED routines above it will be dereddened mags.

        @param [in] sedObj is an Sed class instance that is set with the wavelength and flux of the
        matched SED

        @param [in] photObj is a PhotometryBase class instance with the Bandpasses set to those
        for the magnitudes given for the catalog object

        @param [in] redshift is the redshift of the object if the magnitude is observed

        @param [in] stepSize is the accuracy you want to match your magNorm within

        @param [in] initBand is the number of the bandpass in the magnitude array that you will use for
        the first naive match guess. Since imsimbandpass uses 500nm the best option is to use that closest
        or encompassing 500 nm. Be aware, this starts at 0, but is initialized to 1 meaning second in array.

        @param [out] testMagNorm is the magnitude normalization for the given magnitudes and SED
        """

        sedTest = Sed()
        sedTest.setSED(sedObj.wavelen, flambda = sedObj.flambda)
        if redshift is not None:
            sedTest.redshiftSED(redshift)
        imSimBand = Bandpass()
        imSimBand.imsimBandpass()
        #Use the object's magnitude in the first band as a naive estimate
        testMagNorm = objectMags[initBand]
        testFluxNorm = sedTest.calcFluxNorm(testMagNorm, imSimBand)
        normedSED = Sed()
        norm_wavelen, norm_fnu = sedTest.multiplyFluxNorm(testFluxNorm, wavelen = sedTest.wavelen,
                                                          fnu = sedTest.fnu)
        normedSED.setSED(norm_wavelen, fnu = norm_fnu)
        sedMags = np.array(photObj.manyMagCalc_list(normedSED))
        diff = np.sort(objectMags - sedMags)
        diffSq = np.sum(diff**2, dtype=np.float64)
        diffSqPrev = np.sum(diff**2, dtype=np.float64)
        #Search either downward or upward along magNorm axis based upon greatest difference
        if diff[np.argmax(np.abs(diff))] < 0:
            alphaAdd = -stepSize
        else:
            alphaAdd = stepSize
        #Recursively adjust the magNorm until you reach a minimum in the sum squared error of the mags

        bestMagNorm = testMagNorm
        bestDiffSq = diffSq
        while diffSq - diffSqPrev < 1.0e-10:
            diffSqPrev = np.sum(diff**2, dtype=np.float64)
            testMagNorm += alphaAdd
            testFluxNorm = sedTest.calcFluxNorm(testMagNorm, imSimBand)
            norm_wavelen, norm_fnu = sedTest.multiplyFluxNorm(testFluxNorm, wavelen = sedTest.wavelen,
                                                              fnu = sedTest.fnu)
            normedSED.setSED(norm_wavelen, fnu = norm_fnu)
            sedMags = np.array(photObj.manyMagCalc_list(normedSED))
            diff = np.sort(objectMags - sedMags)
            diffSq = np.sum(diff**2, dtype=np.float64)
            if diffSq < bestDiffSq:
                bestMagNorm = testMagNorm
                bestDiffSq = diffSq

        return bestMagNorm
