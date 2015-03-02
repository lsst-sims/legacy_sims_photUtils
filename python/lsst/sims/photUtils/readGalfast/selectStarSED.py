import os
import numpy as np
import eups

from lsst.sims.photUtils.readGalfast.rgUtils import rgStar
from lsst.sims.photUtils.Photometry import PhotometryBase as phot
from lsst.sims.photUtils.EBV import EBVbase as ebv

__all__ = ["selectStarSED"]

class selectStarSED(rgStar):

    """
    This class provides a way to match incoming star colors to those of the approriate SED.
    """

    def findSED(self, sedList, catMags, catRA, catDec, reddening = True, magNormAcc = 2, bandpassList = None, 
                colors = None, extCoeffs = (4.239, 3.303, 2.285, 1.698, 1.263)):

        """
        This will find the closest match to an input SED. sEDDict must have 'kurucz', 'mlt', 'wdH', or 'wdHE'
        depending on which types you are inputting and will determine which it is based upon its 
        galfast component value

        @param [in] sEDDict is a dictionary from the load(type) routines that are a part of this class

        @param [in] magU is the U-band magnitude from galfast output

        @param [in] magG is the G-band magnitude from galfast output

        @param [in] magR is the R-band magnitude from galfast output

        @param [in] magI is the I-band magnitude from galfast output

        @param [in] magZ is the Z-band magnitude form galfast output

        @param [in] am is the extinction parameter for the object from galfast output

        @param [in] comp is the component number from the galfast output. Used to determine the type
        of SED to use. See galfast documentation for more info.

        @param [in] reddening indicates whether there is extinction and reddening included
        in the galfast output
        
        @param [in] magNormAcc is the number of decimal places within the magNorm result will be accurate.

        @param [in] coeffs is the set of coefficients from galfast's photometry.conf file that scale
        the extinction in each band, usually calibrated to 1.0 in R-band        

        @param [out] sEDName is the name of the SED file that most closely matches the input mags
        accounting for type of star
        
        @param [out] magNorm is the magnitude normalization for the given magnitudes and SED
        
        """
        
        starPhot = phot()
        if bandpassList is None:
            starPhot.loadBandPassesFromFiles(['u','g','r','i','z'], 
                                            bandPassDir = os.path.join(eups.productDir('throughputs'),
                                                                       'sdss'),
                                            bandPassRoot = 'sdss_')
        else:
            starPhot.bandPassList = bandpassList
        starPhot.setupPhiArray_dict()
        
        if colors is None:
            modelColors = self.calcBasicColors(sedList, starPhot)
        else:
            modelColors = colors
        #Transpose so that all values for one color are in one row as needed for the matching loop below
        modelColors = np.transpose(modelColors)

        if reddening == True:
            calcEBV = ebv()
            raDec = np.array((catRA,catDec))
            #If only matching one object need to reshape for calculateEbv
            if len(raDec.shape) == 1:
                raDec = raDec.reshape((2,1))
            ebvVals = calcEBV.calculateEbv(equatorialCoordinates = raDec)
            objMags = self.deReddenMags(ebvVals, catMags, extCoeffs)
        else:
            objMags = catMags
            
        matchColors = []

        for filtNum in range(0, len(extCoeffs)-1):
            matchColors.append(np.transpose(objMags)[filtNum] - np.transpose(objMags)[filtNum+1])

        matchColors = np.transpose(matchColors)

        numCatMags = len(catMags)
        numOn = 0
        sedMatches = []
        magNormMatches = []

        for catObject in matchColors:
            distanceArray = np.zeros(len(sedList))
            for filtNum in range(0, len(starPhot.bandPassList)-1):
                distanceArray += np.power((modelColors[filtNum] - catObject[filtNum]),2)
            matchedSEDNum = np.nanargmin(distanceArray)
            sedMatches.append(sedList[matchedSEDNum].name)
            magNorm = self.calcMagNorm(objMags[numOn], sedList[matchedSEDNum], 
                                       starPhot, stepSize = np.power(10, -float(magNormAcc)))
            magNormMatches.append(magNorm)
            numOn += 1
            if numOn % 10000 == 0:
                print 'Matched %i of %i catalog objects to SEDs' % (numOn, numCatMags)
        
        if numCatMags > 1:        
            print 'Done Matching. Matched %i catalog objects to SEDs' % (numCatMags)

        return sedMatches, magNormMatches
