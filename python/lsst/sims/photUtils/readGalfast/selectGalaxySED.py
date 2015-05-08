import os
import numpy as np
import eups

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.readGalfast.rgUtils import rgGalaxy
from lsst.sims.photUtils.Photometry import PhotometryBase as phot
from lsst.sims.photUtils.EBV import EBVbase as ebv

__all__ = ["selectGalaxySED"]

class selectGalaxySED(rgGalaxy):

    """
    This class provides methods to match galaxy catalog magnitudes to an SED.
    """

    def matchToRestFrame(self, sedList, catMags, mag_error = None, bandpassDict = None, makeCopy = False):

        """
        This will find the closest match to the magnitudes of a galaxy catalog if those magnitudes are in
        the rest frame. Objects without magnitudes in at least two adjacent bandpasses will return as none
        and print out a message.

        @param [in] sedList is the set of spectral objects from the models SEDs provided by loadBC03
        or other custom loader routine.

        @param [in] catMags is an array of the magnitudes of catalog objects to be matched with a model SED.
        It should be organized so that there is one object's magnitudes along each row.

        @param [in] mag_error are provided error values for magnitudes in objectMags. If none provided
        then this defaults to 1.0. This should be an array of the same size as catMags.

        @param [in] bandpassDict is an OrderedDict of bandpass objects with which to calculate magnitudes. 
        If left equal to None it will by default load the SDSS [u,g,r,i,z] bandpasses.

        @param [in] makeCopy indicates whether or not to operate on copies of the SED objects in sedList
        since this method will change the wavelength grid.

        @param [out] sedMatches is a list with the name of a model SED that matches most closely to each
        object in the catalog.

        @param [out] magNormMatches are the magnitude normalizations for the given magnitudes and
        matched SED.

        @param [out] matchErrors contains the Mean Squared Error between the colors of each object and 
        the colors of the matched SED.
        """

        #Set up photometry to calculate model Mags
        galPhot = phot()
        if bandpassDict is None:
            galPhot.loadTotalBandpassesFromFiles(['u','g','r','i','z'],
                                            bandpassDir = os.path.join(eups.productDir('throughputs'),'sdss'),
                                            bandpassRoot = 'sdss_')
        else:
            galPhot.bandpassDict = bandpassDict
        galPhot.setupPhiArray_dict()

        modelColors = []
        sedMatches = []
        magNormMatches = []

        #Find the colors for all model SEDs
        modelColors = self.calcBasicColors(sedList, galPhot, makeCopy = makeCopy)
        modelColors = np.transpose(modelColors)

        #Match the catalog colors to models
        numCatMags = len(catMags)
        numOn = 0
        notMatched = 0
        matchColors = []
        matchErrors = []

        for filtNum in range(0, len(galPhot.bandpassDict)-1):
            matchColors.append(np.transpose(catMags)[filtNum] - np.transpose(catMags)[filtNum+1])

        matchColors = np.transpose(matchColors)

        for catObject in matchColors:
            #This is done to handle objects with incomplete magnitude data
            colorRange = np.arange(0, len(galPhot.bandpassDict)-1)
            filtNums = np.arange(0, len(galPhot.bandpassDict))
            if np.isnan(np.amin(catObject))==True:
                colorRange = np.where(np.isnan(catObject)==False)[0]
                filtNums = np.unique([colorRange, colorRange+1]) #To pick out right filters in calcMagNorm
            if len(colorRange) == 0:
                print 'Could not match object #%i. No magnitudes for two adjacent bandpasses.' % (numOn)
                notMatched += 1
                sedMatches.append(None)
                magNormMatches.append(None)
                matchErrors.append(None)
            else:
                distanceArray = np.zeros(len(sedList))
                for colorNum in colorRange:
                    distanceArray += np.power((modelColors[colorNum] - catObject[colorNum]),2)
                matchedSEDNum = np.nanargmin(distanceArray)
                sedMatches.append(sedList[matchedSEDNum].name)
                magNorm = self.calcMagNorm(np.array(catMags[numOn]), sedList[matchedSEDNum],
                                           galPhot, mag_error = mag_error, filtRange = filtNums)
                magNormMatches.append(magNorm)
                matchErrors.append(distanceArray[matchedSEDNum]/len(colorRange))
            numOn += 1
            if numOn % 10000 == 0:
                print 'Matched %i of %i catalog objects to SEDs' % (numOn-notMatched, numCatMags)

        print 'Done Matching. Matched %i of %i catalog objects to SEDs' % (numCatMags-notMatched, numCatMags)
        if notMatched > 0:
            print '%i objects did not get matched' % (notMatched)

        return sedMatches, magNormMatches, matchErrors

    def matchToObserved(self, sedList, catMags, catRedshifts, catRA = None, catDec = None,
                        mag_error = None, bandpassDict = None, dzAcc = 2, reddening = True,
                        extCoeffs = (4.239, 3.303, 2.285, 1.698, 1.263)):

        """
        This will find the closest match to the magnitudes of a galaxy catalog if those magnitudes are in
        the observed frame and can correct for reddening from within the milky way as well if needed.
        In order to make things faster it first calculates colors for all model SEDs at redshifts between
        the minimum and maximum redshifts of the catalog objects provided with a grid spacing in redshift
        defined by the parameter dzAcc. Objects without magnitudes in at least two adjacent bandpasses will
        return as none and print out a message.

        @param [in] sedList is the set of spectral objects from the models SEDs provided by loadBC03
        or other custom loader routine.

        @param [in] catMags is an array of the magnitudes of catalog objects to be matched with a model SED.
        It should be organized so that there is one object's magnitudes along each row.

        @param [in] catRedshifts is an array of the redshifts of each catalog object.

        @param [in] catRA is an array of the RA positions for each catalog object.

        @param [in] catDec is an array of the Dec position for each catalog object.

        @param [in] mag_error are provided error values for magnitudes in objectMags. If none provided
        then this defaults to 1.0. This should be an array of the same size as catMags.

        @param [in] bandpassDict is an OrderedDict of bandpass objects with which to calculate magnitudes.
        If left equal to None it will by default load the SDSS [u,g,r,i,z] bandpasses and therefore agree with
        default extCoeffs.

        @param [in] dzAcc is the number of decimal places you want to use when building the redshift grid.
        For example, dzAcc = 2 will create a grid between the minimum and maximum redshifts with colors
        calculated at every 0.01 change in redshift.

        @param [in] reddening is a boolean that determines whether to correct catalog magnitudes for
        dust in the milky way. By default, it is True.
        If true, this uses calculateEBV from EBV.py to find an EBV value for the object's
        ra and dec coordinates and then uses the coefficients provided by extCoeffs which should come
        from Schlafly and Finkbeiner (2011) for the correct filters and in the same order as provided
        in bandpassDict.
        If false, this means it will not run the dereddening procedure.

        @param [in] extCoeffs are the Schlafly and Finkbeiner (2011) (ApJ, 737, 103) coefficients for the
        given filters from bandpassDict and need to be in the same order as bandpassDict. The default given
        are the SDSS [u,g,r,i,z] values.

        @param [out] sedMatches is a list with the name of a model SED that matches most closely to each
        object in the catalog.

        @param [out] magNormMatches are the magnitude normalizations for the given magnitudes and
        matched SED.

        @param [out] matchErrors contains the Mean Squared Error between the colors of each object and 
        the colors of the matched SED.
        """

        #Set up photometry to calculate model Mags
        galPhot = phot()
        if bandpassDict is None:
            galPhot.loadTotalBandpassesFromFiles(['u','g','r','i','z'],
                                            bandpassDir = os.path.join(eups.productDir('throughputs'),'sdss'),
                                            bandpassRoot = 'sdss_')
        else:
            galPhot.bandpassDict = bandpassDict
        galPhot.setupPhiArray_dict()

        #Calculate ebv from ra, dec coordinates if needed
        if reddening == True:
            #Check that catRA and catDec are included
            if catRA is None or catDec is None:
                raise RuntimeError("Reddening is True, but catRA and catDec are not included.")
            calcEBV = ebv()
            raDec = np.array((catRA,catDec))
            #If only matching one object need to reshape for calculateEbv
            if len(raDec.shape) == 1:
                raDec = raDec.reshape((2,1))
            ebvVals = calcEBV.calculateEbv(equatorialCoordinates = raDec)
            objMags = self.deReddenMags(ebvVals, catMags, extCoeffs)
        else:
            objMags = catMags

        minRedshift = np.round(np.min(catRedshifts), dzAcc)
        maxRedshift = np.round(np.max(catRedshifts), dzAcc)
        dz = np.power(10., (-1*dzAcc))

        redshiftRange = np.round(np.arange(minRedshift - dz, maxRedshift + (2*dz), dz), dzAcc)
        numRedshifted = 0
        sedMatches = [None] * len(catRedshifts)
        magNormMatches = [None] * len(catRedshifts)
        matchErrors = [None] * len(catRedshifts)
        redshiftIndex = np.argsort(catRedshifts)

        numOn = 0
        notMatched = 0
        lastRedshift = -100
        print 'Starting Matching. Arranged by redshift value.'
        for redshift in redshiftRange:

            if numRedshifted % 10 == 0:
                print '%i out of %i redshifts gone through' % (numRedshifted, len(redshiftRange))
            numRedshifted += 1

            colorSet = []
            for galSpec in sedList:
                sedColors = []
                fileSED = Sed()
                fileSED.setSED(wavelen = galSpec.wavelen, flambda = galSpec.flambda)
                fileSED.redshiftSED(redshift)
                sedColors = self.calcBasicColors([fileSED], galPhot, makeCopy = True)
                colorSet.append(sedColors)
            colorSet = np.transpose(colorSet)
            for currentIndex in redshiftIndex[numOn:]:
                matchMags = objMags[currentIndex]
                if lastRedshift < np.round(catRedshifts[currentIndex],dzAcc) <= redshift:
                    colorRange = np.arange(0, len(galPhot.bandpassDict)-1)
                    matchColors = []
                    for colorNum in colorRange:
                        matchColors.append(matchMags[colorNum] - matchMags[colorNum+1])
                    #This is done to handle objects with incomplete magnitude data
                    filtNums = np.arange(0, len(galPhot.bandpassDict))
                    if np.isnan(np.amin(matchColors))==True:
                        colorRange = np.where(np.isnan(matchColors)==False)[0]
                        filtNums = np.unique([colorRange, colorRange+1]) #Pick right filters in calcMagNorm
                    if len(colorRange) == 0:
                        print 'Could not match object #%i. No magnitudes for two adjacent bandpasses.' \
                              % (currentIndex)
                        notMatched += 1
                        #Don't need to assign 'None' here in result array, b/c 'None' is default value
                    else:
                        distanceArray = [np.zeros(len(sedList))]
                        for colorNum in colorRange:
                            distanceArray += np.power((colorSet[colorNum] - matchColors[colorNum]),2)
                        matchedSEDNum = np.nanargmin(distanceArray)
                        sedMatches[currentIndex] = sedList[matchedSEDNum].name
                        magNormVal = self.calcMagNorm(np.array(matchMags), sedList[matchedSEDNum], 
                                                      galPhot, mag_error = mag_error,
                                                      redshift = catRedshifts[currentIndex],
                                                      filtRange = filtNums)
                        magNormMatches[currentIndex] = magNormVal
                        matchErrors[currentIndex] = (distanceArray[0,matchedSEDNum]/len(colorRange))
                    numOn += 1
                else:
                    break
            lastRedshift = redshift

        print 'Done Matching. Matched %i of %i catalog objects to SEDs' % (len(catMags)-notMatched, 
                                                                           len(catMags))
        if notMatched > 0:
            print '%i objects did not get matched.' % (notMatched)

        return sedMatches, magNormMatches, matchErrors
