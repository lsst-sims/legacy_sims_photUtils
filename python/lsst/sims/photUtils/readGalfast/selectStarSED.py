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
        This will find the closest match to the magnitudes of a galaxy catalog if those magnitudes are in
        the observed frame and can correct for extinction from within the milky way as well if needed.
        In order to make things faster it first calculates colors for all model SEDs at redshifts between
        the minimum and maximum redshifts of the catalog objects provided with a grid spacing in redshift
        defined by the parameter dzAcc.

        @param [in] sedList is the set of spectral objects from the models SEDs provided by loadBC03
        or other custom loader routine.
        
        @param [in] catMags is an array of the magnitudes of catalog objects to be matched with a model SED.
        It should be organized so that there is one object's magnitudes along each row.

        @param [in] catRA is an array of the RA positions for each catalog object.

        @param [in] catDec is an array of the Dec position for each catalog object.

        @param [in] reddening is a boolean that determines whether to correct catalog magnitudes for 
        dust in the milky way. This uses calculateEBV from EBV.py to find an EBV value for the object's
        ra and dec coordinates and then uses the coefficients provided by extCoeffs which should come
        from Schlafly and Finkbeiner (2011) for the correct filters and in the same order as provided
        in bandpassList.
        
        @param [in] magNormAcc is the number of decimal places within the magNorm result will be accurate.

        @param [in] bandpassList is a list of bandpass objects with which to calculate magnitudes. If left
        equal to None it will by default load the SDSS [u,g,r,i,z] bandpasses and therefore agree with 
        default extCoeffs.
        
        @param [in] colors is None if you are just providing a list of SED objects to match, but is the 
        array holding the colors of those SED models (each row should be the colors for one model in the 
        same order as sedList) if you have already calculated the colors.

        @param [in] extCoeffs are the Schlafly and Finkbeiner (2011) coefficients for the given filters
        from bandpassList and need to be in the same order as bandpassList. The default given are the SDSS
        [u,g,r,i,z] values.

        @param [out] sedMatches is a list with the name of a model SED that matches most closely to each
        object in the catalog.

        @param [out] magNormMatches are the magnitude normalizations for the given magnitudes and 
        matched SED.
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
