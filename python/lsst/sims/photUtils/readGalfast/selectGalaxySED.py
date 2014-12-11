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

        for fileName in files[0:10]:
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

#    def loadCatalog(self, catFile, redshiftCol = None, filters = ('u', 'g', 'r', 'i', 'z'), filterCols = (1,2,3,4,5)):
        
#        if os.path.isfile(catFile) == False:
#            raise RuntimeError, '***File does not exist'

#        if filename.endswith(('.txt', '.gz', '.fits')):

    def matchRestFrameSEDs(self, sedList, catMags, filterList = ('u', 'g', 'r', 'i', 'z'),
                           throughputDir = os.getenv('SDSS_THROUGHPUTS'), filterRoot = 'sdss_'):

        #Set up photometry pieces for SDSS Mags
        galPhot = phot()
        bandpassDict = galPhot.loadBandpasses(filterlist = filterList,
                                                   dataDir = throughputDir,
                                                   filterroot = filterRoot)
        phiArray, wavelenstep = galPhot.setupPhiArray_dict(bandpassDict, filterList)

        colorName = []
        specNum = 0
        sedMatches = []

        for galSpec in sedList:
            fileSED = Sed()
            fileSED.setSED(galSpec.wave, flambda = galSpec.flux)
            sEDMagDict = galPhot.manyMagCalc_dict(fileSED, phiArray, wavelenstep, bandpassDict, filterList)
            colorInfo = []
            for filtNum in range(0, len(filterList)-1):
                colorInfo.append(sEDMagDict[filterList[filtNum]] - sEDMagDict[filterList[filtNum+1]])
            colorInfo.append(galSpec.name)
            colorName.append(colorInfo)
            specNum += 1

        for matchMags in catMags:
            matchColors = []
            for filtNum in range(0, len(matchMags)-1):
                matchColors.append(matchMags[filtNum]-matchMags[filtNum+1])
            distanceArray = []
            for modelColor in colorName:
                distance = 0.
                for filtNum in range(0, len(filterList)-1):
                    distance += np.power((modelColor[filtNum] - matchColors[filtNum]),2)
                distanceArray.append(distance)
            sedMatch = sedList[np.argmin(distanceArray)].name
            sedMatches.append(sedMatch)
        return sedMatches
