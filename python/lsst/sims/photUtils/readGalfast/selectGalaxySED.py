import os
import gzip
import pyfits
import numpy as np
import re

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.photUtils import Photometry as phot
from lsst.sims.photUtils.EBV import EBVbase as ebv
from lsst.sims.photUtils.readGalfast.rgUtils import Spectrum
from lsst.sims.catalogs.measures.instance.fileMaps import SpecMap

__all__ = ["selectGalaxySED"]

class selectGalaxySED():

    def __init__(self, galDir = None):
        
        """
        @param [in] galDir is the directory where the galaxy SEDs are stored
        """

        if galDir == None:
            #Use SpecMap to pull in directory's location in LSST Stack
            specMap = SpecMap()
            specFileStart = 'Exp' #Start of sample BC03 name in sims_sed_library
            for key, val in sorted(specMap.subdir_map.iteritems()):
                if re.match(key, specFileStart):
                    galSpecDir = str(val)
            self.galDir = str(os.environ['SIMS_SED_LIBRARY_DIR'] + '/' + galSpecDir)            
        else:
            self.galDir = galDir

    def loadBC03(self, subset = None):

        """
        This loads the Bruzual and Charlot SEDs that are currently in the SIMS_SED_LIBRARY.
        If the user wants to use different SEDs another loading method can be created and used in place
        of this.
        
        @param [in] subset is the list of the subset of files in the galDir that the user
        can specify if using all the SEDs in the directory is not desired.

        @param [out] sedList is the set of model SED spectra objects to be passed onto the matching routines.
        """

        files = []

        if subset is None:
            for fileName in os.listdir(self.galDir):
                files.append(fileName)
        else:
            for fileName in subset:
                files.append(fileName)

        numFiles = len(files)
        numOn = 0

        sedList = []

        for fileName in files:
            if numOn % 100 == 0:
                print 'Loading %i of %i: BC Galaxy SEDs' % (numOn, numFiles)
 
            try:
                spec = Spectrum(str(self.galDir + '/' + fileName))
                spec.name = fileName
                spec.type = fileName.split('.')[0]
                spec.age = float(fileName.split('.')[1])
                metallicity = fileName.split('.')[2].split('Z')[0]
                #Final form is z/zSun
                spec.metallicity = float(metallicity) * (10 ** ((len(metallicity)-1)*-1))

            except:
                continue

            sedList.append(spec)

            numOn += 1

        return sedList

    def matchToRestFrame(self, sedList, catMags, filterList = ('u', 'g', 'r', 'i', 'z'),
                         throughputDir = os.getenv('SDSS_THROUGHPUTS'), filterRoot = 'sdss_'):

        """
        This will find the closest match to the magnitudes of a galaxy catalog if those magnitudes are in
        the rest frame.

        @param [in] sedList is the set of spectral objects from the models SEDs provided by loadBC03
        or other custom loader routine.

        @param [in] catMags is an array of the magnitudes of catalog objects to be matched with a model SED.
        It should be organized so that there is one object's magnitudes along each row.

        @param [in] filterList is the set of filters corresponding to the object magnitudes in the order
        they are organized in each row of the array.

        @param [in] throughputDir is the directory where the filter throughputs are stored.

        @param [in] filterroot is the root name of the throughputs in the directory that you want to use.

        @param [out] sedMatches is a list with the name of a model SED that matches most closely to each
        object in the catalog.
        """

        #Set up photometry to calculate model Mags
        galPhot = phot()
        bandpassDict = galPhot.loadBandpasses(filterlist = filterList,
                                                   dataDir = throughputDir,
                                                   filterroot = filterRoot)
        phiArray, wavelenstep = galPhot.setupPhiArray_dict(bandpassDict, filterList)

        modelColors = []
        sedMatches = []

        #Find the colors for all model SEDs
        for galSpec in sedList:
            fileSED = Sed()
            fileSED.setSED(galSpec.wave, flambda = galSpec.flux)
            sEDMagDict = galPhot.manyMagCalc_dict(fileSED, phiArray, wavelenstep, bandpassDict, filterList)
            colorInfo = []
            for filtNum in range(0, len(filterList)-1):
                colorInfo.append(sEDMagDict[filterList[filtNum]] - sEDMagDict[filterList[filtNum+1]])
            modelColors.append(colorInfo)
        modelColors = np.transpose(modelColors)

        #Match the catalog colors to models
        numCatMags = len(catMags)
        numOn = 0
        matchColors = []

        for filtNum in range(0, len(filterList)-1):
            matchColors.append(np.transpose(catMags)[filtNum] - np.transpose(catMags)[filtNum+1])

        matchColors = np.transpose(matchColors)

        for catObject in matchColors:
            if numOn % 10000 == 0:
                print 'Matched %i of %i catalog objects to SEDs' % (numOn, numCatMags)
            distanceArray = np.zeros(len(sedList))
            for filtNum in range(0, len(filterList)-1):
                distanceArray += np.power((modelColors[filtNum] - catObject[filtNum]),2)
            sedMatches.append(sedList[np.nanargmin(distanceArray)].name)
            numOn += 1
            
        return sedMatches

    def matchToObserved(self, sedList, catRA, catDec, catRedshifts, catMags, 
                        filterList = ('u','g','r','i','z'), throughputDir = os.getenv('SDSS_THROUGHPUTS'), 
                        filterRoot = 'sdss_', dzAcc = 2, extinction = True, 
                        extCoeffs = (4.239, 3.303, 2.285, 1.698, 1.263)):

        """
        This will find the closest match to the magnitudes of a galaxy catalog if those magnitudes are in
        the observed frame and can correct for extinction from within the milky way as well if needed.
        In order to make things faster it first calculates colors for all model SEDs at redshifts between
        the minimum and maximum redshifts of the catalog objects provided with a grid spacing in redshift
        defined by the parameter dzAcc.

        @param [in] sedList is the set of spectral objects from the models SEDs provided by loadBC03
        or other custom loader routine

        @param [in] catRA is an array of the RA positions for each catalog object.

        @param [in] catDec is an array of the Dec position for each catalog object.

        @param [in] catRedshifts is an array of the redshifts of each catalog object.

        @param [in] catMags is an array of the magnitudes of catalog objects to be matched with a model SED.
        It should be organized so that there is one object's magnitudes along each row.

        @param [in] filterList is the set of filters corresponding to the object magnitudes in the order
        they are organized in each row of the array.

        @param [in] throughputDir is the directory where the filter throughputs are stored.

        @param [in] filterroot is the root name of the throughputs in the directory that you want to use.

        @param [in] dzAcc is the number of decimal places you want to use when building the redshift grid.
        For example, dzAcc = 2 will create a grid between the minimum and maximum redshifts with colors
        calculated at every 0.01 change in redshift.

        @param [in] extinction is a boolean that determines whether to correct catalog magnitudes for 
        dust in the milky way. This uses calculateEBV from EBV.py to find an EBV value for the object's
        ra and dec coordinates and then uses the coefficients provided by extCoeffs which should come
        from Schlafly and Finkbeiner (2011) for the correct filters and in the same order as provided
        in filterList.

        @param [in] extCoords are the Schlafly and Finkbeiner (2011) coefficients for the given filters
        from filterList and need to be in the same order as filterList.

        @param [out] sedMatches is a list with the name of a model SED that matches most closely to each
        object in the catalog.
        """

        #Set up photometry to calculate model Mags
        galPhot = phot()
        bandpassDict = galPhot.loadBandpasses(filterlist = filterList,
                                              dataDir = throughputDir,
                                              filterroot = filterRoot)
        phiArray, wavelenstep = galPhot.setupPhiArray_dict(bandpassDict, filterList)
        
        minRedshift = np.round(np.min(catRedshifts), dzAcc)
        maxRedshift = np.round(np.max(catRedshifts), dzAcc)
        dz = np.power(10., (-1*dzAcc))

        redshiftColors = {}
        redshiftRange = np.arange(minRedshift - dz, maxRedshift + dz, dz)
        numRedshifted = 0
        print 'Building Redshifted Color Set.'
        for redshift in redshiftRange:
            colorSet = []
            for galSpec in sedList:
                sedColors = []
                fileSED = Sed()
                fileSED.setSED(galSpec.wave, flambda = galSpec.flux)
                fileSED.redshiftSED(redshift)
                sEDMagDict = galPhot.manyMagCalc_dict(fileSED, phiArray, wavelenstep,
                                                      bandpassDict, filterList)
                for filtNum in range(0, len(filterList)-1):
                    sedColors.append(sEDMagDict[filterList[filtNum]] - sEDMagDict[filterList[filtNum+1]])
                colorSet.append(sedColors)
            colorSet = np.transpose(colorSet)
            redshiftColors[str(np.round(redshift,dzAcc))] = colorSet
            if numRedshifted % 10 == 0:
                print '%i out of %i redshifts gone through' % (numRedshifted, len(redshiftRange))
            numRedshifted += 1

        print 'Done Building Set. Starting Matching.'

        sedMatches = []
        numCatMags = len(catMags)
        numOn = 0

        #Calculate ebv from ra, dec coordinates if needed
        if extinction == True:
            calcEBV = ebv()
            raDec = np.array((catRA,catDec))
            #If only matching one object need to reshape for calculateEbv
            if len(raDec.shape) == 1:
                raDec = raDec.reshape((2,1))
            ebvVals = calcEBV.calculateEbv(equatorialCoordinates = raDec)
        #If extinction is false is won't be used anyway so just set it to ones for the loop
        else:
            ebvVals = np.ones(len(catMags))

        for ebvValue, matchMags, matchRedshift in zip(ebvVals, catMags, catRedshifts):
            if numOn % 10000 == 0:
                print 'Matched %i of %i catalog objects to SEDs' % (numOn, numCatMags)
            
            distanceArray = np.zeros(len(sedList))
            modelColors = redshiftColors[str(np.round(matchRedshift, dzAcc))]
            if extinction == True:
                for filtNum in range(0, len(filterList)):
                    matchMags[filtNum] = matchMags[filtNum] - (extCoeffs[filtNum]*ebvValue)
            for filtNum in range(0, len(filterList)-1):
                matchColor = matchMags[filtNum] - matchMags[filtNum+1]
                distanceArray += np.power((modelColors[filtNum] - matchColor),2)
            sedMatches.append(sedList[np.nanargmin(distanceArray)].name)
            numOn += 1

        return sedMatches
