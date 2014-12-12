import os
import gzip
import pyfits
import numpy as np
import re
from collections import defaultdict

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.photUtils import Photometry as phot
from lsst.sims.catalogs.measures.instance.fileMaps import SpecMap

__all__ = ["selectGalaxySED"]

class Spectrum(object):
    def __init__(self, fileName):
        self.wave, self.flux = self.readSpectrum(fileName)

    def readSpectrum(self, fileName):
        wave, flux = np.genfromtxt(fileName, unpack=True)
        return wave, flux

class selectGalaxySED():

    def __init__(self, bcDir = None):
        
        if bcDir == None:
            #Use SpecMap to pull in directory's location in LSST Stack
            specMap = SpecMap()
            specFileStart = 'Exp' #Start of sample BC03 name in sims_sed_library
            for key, val in sorted(specMap.subdir_map.iteritems()):
                if re.match(key, specFileStart):
                    galSpecDir = str(val)
            self.bcDir = str(os.environ['SIMS_SED_LIBRARY_DIR'] + '/' + galSpecDir)            
        else:
            self.bcDir = bcDir

    def loadBC03(self, subset = None):

        files = []

        if subset is None:
            for fileName in os.listdir(self.bcDir):
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
                spec = Spectrum(str(self.bcDir + '/' + fileName))
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

    def matchToObserved(self, sedList, catRedshifts, catMags, filterList = ('u','g','r','i','z'),
                        throughputDir = os.getenv('SDSS_THROUGHPUTS'), filterRoot = 'sdss_'):

        #Set up photometry to calculate model Mags
        galPhot = phot()
        bandpassDict = galPhot.loadBandpasses(filterlist = filterList,
                                              dataDir = throughputDir,
                                              filterroot = filterRoot)
        phiArray, wavelenstep = galPhot.setupPhiArray_dict(bandpassDict, filterList)
        
        sedMatches = []
        numCatMags = len(catMags)
        numOn = 0
        for matchMags, matchRedshift in zip(catMags, catRedshifts):
            if numOn % 100 == 0:
                print 'Matched %i of %i catalog objects to SEDs' % (numOn, numCatMags)
            matchColors = []
            for filtNum in range(0, len(matchMags)-1):
                matchColors.append(matchMags[filtNum]-matchMags[filtNum+1])
            matchColors = np.transpose(matchColors)
                
            modelColors = []

            for galSpec in sedList:
                fileSED = Sed()
                fileSED.setSED(galSpec.wave, flambda = galSpec.flux)
                fileSED.redshiftSED(matchRedshift)
                sEDMagDict = galPhot.manyMagCalc_dict(fileSED, phiArray, wavelenstep, 
                                                      bandpassDict, filterList)
                colorInfo = []
                for filtNum in range(0, len(filterList)-1):
                    colorInfo.append(sEDMagDict[filterList[filtNum]] - sEDMagDict[filterList[filtNum+1]])
                modelColors.append(colorInfo)
            modelColors = np.transpose(modelColors)
            
            distanceArray = np.zeros(len(sedList))
            for filtNum in range(0, len(filterList)-1):
                distanceArray += np.power((modelColors[filtNum] - matchColors[filtNum]),2)
            sedMatches.append(sedList[np.argmin(distanceArray)].name)
            print sedList[np.nanargmin(distanceArray)].name
            numOn += 1

        return sedMatches
