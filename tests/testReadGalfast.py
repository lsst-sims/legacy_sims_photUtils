import unittest
import os
import shutil
import numpy as np
import gzip
import pyfits
import re
import lsst.utils.tests as utilsTests
from lsst.sims.photUtils.readGalfast.selectStarSED import selectStarSED
from lsst.sims.photUtils.readGalfast.readGalfast import readGalfast
from lsst.sims.photUtils.readGalfast.selectGalaxySED import selectGalaxySED
from lsst.sims.photUtils.EBV import EBVbase as ebv
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Photometry import PhotometryBase as phot
from lsst.sims.catalogs.measures.instance.fileMaps import SpecMap

class TestSelectGalaxySED(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        
        #because SCons only knows about a limited subset of the LSST environment variables
        os.environ['LSST_THROUGHPUTS_DEFAULT'] = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline')
        os.environ['SDSS_THROUGHPUTS'] = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'sdss')
    
        specMap = SpecMap()
        specFileStart = 'Exp'
        for key, val in sorted(specMap.subdir_map.iteritems()):
            if re.match(key, specFileStart):
                galSpecDir = str(val)
        cls.galDir = str(os.environ['SIMS_SED_LIBRARY_DIR'] + '/' + galSpecDir + '/')

        #Set up Test Spectra Directory
        cls.testSpecDir = 'testGalaxySEDSpectrum/'
        
        if os.path.exists(cls.testSpecDir):
            shutil.rmtree(cls.testSpecDir)

        os.mkdir(cls.testSpecDir)

        galList = os.listdir(cls.galDir)[0:20]

        for galFile in galList:
            shutil.copy(str(cls.galDir + galFile), str(cls.testSpecDir + galFile))

        cls.filterList = ('u', 'g', 'r', 'i', 'z')

    def testLoadBC03(self):

        loadTestBC03 = selectGalaxySED(galDir = self.testSpecDir)
        testSEDs = loadTestBC03.loadBC03()

        #Read in a list of the SEDs in the test galaxy sed directory
        testGalList = os.listdir(self.testSpecDir)

        #Make sure the names of seds in folder and set that was read in are the same
        #This also tests that the name attribute is assigned to each Spectrum object correctly
        testNames = []
        for testSED in testSEDs:
            testNames.append(testSED.name)
        self.assertItemsEqual(testGalList, testNames)

        #Test same condition if a subset is provided
        testSubsetList = testGalList[0:2]
        testSEDsubset = loadTestBC03.loadBC03(subset = testSubsetList)
        testSubsetNames = []
        for testSED in testSEDsubset:
            testSubsetNames.append(testSED.name)
        self.assertItemsEqual(testSubsetList, testSubsetNames)

        #Test that attributes have been assigned
        for testSED in testSEDsubset:
            self.assertIsNotNone(testSED.name)
            self.assertIsNotNone(testSED.type)
            self.assertIsNotNone(testSED.age)
            self.assertIsNotNone(testSED.metallicity)

    def testMatchToRestFrame(self):

        galPhot = phot()
        galPhot.loadBandPassesFromFiles(self.filterList, bandPassDir = os.getenv('SDSS_THROUGHPUTS'), 
                               bandPassRoot = 'sdss_')
        galPhot.setupPhiArray_dict()

        testMatching = selectGalaxySED(galDir = self.testSpecDir)
        testSEDList = testMatching.loadBC03()

        testSEDNames = []
        testMags = []

        for testSED in testSEDList:
            
            getSEDMags = Sed()
            testSEDNames.append(testSED.name)
            getSEDMags.setSED(wavelen = testSED.wavelen, flambda = testSED.flambda)
            testMags.append(galPhot.manyMagCalc_list(getSEDMags))

        testMatchingResults = testMatching.matchToRestFrame(testSEDList, testMags, 
                                                            bandpassDir = os.getenv('SDSS_THROUGHPUTS'))

        self.assertEqual(testSEDNames, testMatchingResults)

    def testMatchToObserved(self):
        
        galPhot = phot()
        galPhot.loadBandPassesFromFiles(self.filterList, bandPassDir = os.getenv('SDSS_THROUGHPUTS'), 
                               bandPassRoot = 'sdss_')
        galPhot.setupPhiArray_dict()

        testMatching = selectGalaxySED(galDir = self.testSpecDir)
        testSEDList = testMatching.loadBC03()

        testSEDNames = []
        testRA = []
        testDec = []
        testRedshifts = []
        extCoeffs = np.arange(1, 1.5, 0.1)
        testMags = []
        testMagsRedshift = []
        testMagsExt = []

        for testSED in testSEDList:

            #As a check make sure that it matches when no extinction and no redshift are present
            getSEDMags = Sed()
            testSEDNames.append(testSED.name)
            getSEDMags.setSED(wavelen = testSED.wavelen, flambda = testSED.flambda)
            testMags.append(galPhot.manyMagCalc_list(getSEDMags))

            #Check Extinction corrections
            sedRA = np.random.uniform(10,170)
            sedDec = np.random.uniform(10,80)
            testRA.append(sedRA)
            testDec.append(sedDec)
            raDec = np.array((sedRA, sedDec)).reshape((2,1))
            ebvVal = ebv().calculateEbv(equatorialCoordinates = raDec)
            extVal = ebvVal*extCoeffs
            testMagsExt.append(galPhot.manyMagCalc_list(getSEDMags) + extVal)

            #Setup magnitudes for testing matching to redshifted values
            getRedshiftMags = Sed()
            testZ = np.round(np.random.uniform(1.1,1.2), 3)
            testRedshifts.append(testZ)
            getRedshiftMags.setSED(wavelen = testSED.wavelen, flambda = testSED.flambda)
            getRedshiftMags.redshiftSED(testZ)
            testMagsRedshift.append(galPhot.manyMagCalc_list(getRedshiftMags))
            
        testNoExtNoRedshift = testMatching.matchToObserved(testSEDList, testRA, testDec, np.zeros(20), 
                                                           testMags, 
                                                           bandpassDir = os.getenv('SDSS_THROUGHPUTS'),
                                                           extinction = False)
        testMatchingExtVals = testMatching.matchToObserved(testSEDList, testRA, testDec, np.zeros(20), 
                                                           testMagsExt, 
                                                           bandpassDir = os.getenv('SDSS_THROUGHPUTS'),
                                                           extinction = True, extCoeffs = extCoeffs)
        testMatchingRedshift = testMatching.matchToObserved(testSEDList, testRA, testDec, testRedshifts,
                                                            testMagsRedshift, 
                                                            bandpassDir = os.getenv('SDSS_THROUGHPUTS'),
                                                            dzAcc = 3, extinction = False)

        self.assertEqual(testSEDNames, testNoExtNoRedshift)
        self.assertEqual(testSEDNames, testMatchingExtVals)
        self.assertEqual(testSEDNames, testMatchingRedshift)

    @classmethod
    def tearDownClass(cls):
        del cls.galDir
        del cls.filterList

        shutil.rmtree(cls.testSpecDir)

class TestSelectStarSED(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        #because SCons only knows about a limited subset of the LSST environment variables
        os.environ['LSST_THROUGHPUTS_DEFAULT'] = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline')
        os.environ['SDSS_THROUGHPUTS'] = os.path.join(os.getenv('THROUGHPUTS_DIR'),'sdss')

        #Left this in after removing loading SEDs so that we can make sure that if the structure of
        #sims_sed_library changes in a way that affects testReadGalfast we can detect it.
        specMap = SpecMap()
        cls._specMapDict = {}
        specFileStart = ['kp', 'burrows', 'bergeron'] #The beginning of filenames of different SED types
        specFileTypes = ['kurucz', 'mlt','wd']
        for specStart, specKey in zip(specFileStart, specFileTypes):
            for key, val in sorted(specMap.subdir_map.iteritems()):
                if re.match(key, specStart):
                    cls._specMapDict[specKey] = str(val)

        cls.kmTestName = 'km99_9999.fits_g99_9999'
        cls.mTestName = 'm99.99Full.dat'

        #Set up Test Spectra Directory
        cls.testSpecDir = 'testReadGalfastSpectra'
        cls.testKDir = str(cls.testSpecDir + '/starSED/kurucz/')
        cls.testMLTDir = str(cls.testSpecDir + '/starSED/mlt/')
        cls.testWDDir = str(cls.testSpecDir + '/starSED/wDs/')

        if os.path.exists(cls.testSpecDir):
            shutil.rmtree(cls.testSpecDir)

        os.makedirs(cls.testKDir)
        os.mkdir(cls.testMLTDir)
        os.mkdir(cls.testWDDir)
        cls.kDir = os.environ['SIMS_SED_LIBRARY_DIR'] + '/' + cls._specMapDict['kurucz'] + '/'
        cls.mltDir = os.environ['SIMS_SED_LIBRARY_DIR'] + '/' + cls._specMapDict['mlt'] + '/'
        cls.wdDir = os.environ['SIMS_SED_LIBRARY_DIR'] + '/' + cls._specMapDict['wd'] + '/'
        kList = os.listdir(cls.kDir)[0:20]
        #Use particular indices to get different types of seds within mlt and wds
        for kFile, mltFile, wdFile in zip(kList, 
                                          np.array(os.listdir(cls.mltDir))[np.arange(-10,11)], 
                                          np.array(os.listdir(cls.wdDir))[np.arange(-10,11)]):
            shutil.copyfile(str(cls.kDir + kFile), str(cls.testKDir + kFile))
            shutil.copyfile(str(cls.mltDir + mltFile), str(cls.testMLTDir + mltFile))
            shutil.copyfile(str(cls.wdDir + wdFile), str(cls.testWDDir + wdFile))
        #Load in extra kurucz to test Logz Readout
        if 'km01_7000.fits_g40_7140.gz' not in kList:
            shutil.copyfile(str(cls.kDir + 'km01_7000.fits_g40_7140.gz'), 
                            str(cls.testKDir + 'km01_7000.fits_g40_7140.gz'))
        if 'kp01_7000.fits_g40_7240.gz' not in kList:
            shutil.copyfile(str(cls.kDir + 'kp01_7000.fits_g40_7240.gz'), 
                            str(cls.testKDir + 'kp01_7000.fits_g40_7240.gz'))        

    @classmethod
    def tearDownClass(cls):
        del cls._specMapDict
        del cls.kDir
        del cls.mltDir
        del cls.wdDir

        if os.path.exists(cls.kmTestName):
            os.unlink(cls.kmTestName)

        if os.path.exists(cls.mTestName):
            os.unlink(cls.mTestName)

        del cls.kmTestName
        del cls.mTestName

        shutil.rmtree(cls.testSpecDir)

    def testDefaults(self):
        """Make sure that if there are Nones for the init that they load the correct dirs"""
        loadTest = selectStarSED()
        self.assertEqual(loadTest.kuruczDir, self.kDir)
        self.assertEqual(loadTest.mltDir, self.mltDir)
        self.assertEqual(loadTest.wdDir, self.wdDir)

    def testLoadKurucz(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for """
        #Test Matching to Kurucz SEDs
        loadTestKurucz = selectStarSED(kuruczDir = self.testKDir)
        testSEDs = loadTestKurucz.loadKuruczSEDs()

        #Read in a list of the SEDs in the kurucz test sed directory
        testKuruczList = os.listdir(self.testKDir)

        #First make sure that all SEDs are correctly accounted for if no subset provided
        self.assertEqual(testSEDs['sEDName'], testKuruczList)

        #Then test to make sure no values are null in any of the lists within the dictionary
        for key in testSEDs.keys():
            for item in testSEDs[key]:
                self.assertIsNotNone(item)

        #Test same condition if subset is provided
        testSubsetList = ['km01_7000.fits_g40_7140.gz', 'kp01_7000.fits_g40_7240.gz']
        testSEDsSubset = loadTestKurucz.loadKuruczSEDs(subset = testSubsetList)

        #Next make sure that correct subset loads if subset is provided
        self.assertEqual(testSEDsSubset['sEDName'], testSubsetList)

        #Finally make sure that values are being entered into dictionaries correctly
        self.assertEqual(testSEDsSubset['logZ'], [-0.1, 0.1]) #Test both pos. and neg. get in right
        self.assertEqual(testSEDsSubset['logg'], [4.0, 4.0])
        self.assertEqual(testSEDsSubset['temp'], [7140, 7240])

        #Now test colors are being calculated correctly by loading a flat sed into it
        testSED = Sed()
        testSEDsColors = selectStarSED()

        testSED.setFlatSED(wavelen_min = 280.0, wavelen_max = 1170.0)
        testSED.multiplyFluxNorm(testSED.calcFluxNorm(10, testSEDsColors.starPhot.bandPassList[2]))
        #Give kurucz like filename just so it works
        testSED.writeSED(self.kmTestName)

        testSEDsColors.kuruczDir = str(os.getcwd() + '/')
        testColors = testSEDsColors.loadKuruczSEDs(subset = [self.kmTestName])
        self.assertAlmostEqual(testColors['umg'][0], 0.0)
        self.assertAlmostEqual(testColors['gmr'][0], 0.0)
        self.assertAlmostEqual(testColors['rmi'][0], 0.0)
        self.assertAlmostEqual(testColors['imz'][0], 0.0)

    def testLoadMLT(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for"""
        #Test Matching to mlt SEDs
        loadTestMLT = selectStarSED(mltDir = self.testMLTDir)
        testSEDs = loadTestMLT.loadmltSEDs()

        #Read in a list of the SEDs in the mlt test sed directory
        testMLTList = os.listdir(self.testMLTDir)

        #First make sure that all SEDs are correctly accounted for if no subset provided
        self.assertItemsEqual(testSEDs['sEDName'], testMLTList)

        #Then test to make sure no values are null in any of the lists within the dictionary
        for key in testSEDs.keys():
            for item in testSEDs[key]:
                self.assertIsNotNone(item)

        #Test same condition if subset is provided
        testSubsetList = ['m0.0Full.dat.gz']
        testSEDsSubset = selectStarSED().loadmltSEDs(subset = testSubsetList)

        #Next make sure that correct subset loads if subset is provided
        self.assertItemsEqual(testSEDsSubset['sEDName'], testSubsetList)

        #Now test colors are being calculated correctly by loading a flat sed into it
        testSED = Sed()
        testSEDsColors = selectStarSED()

        testSED.setFlatSED(wavelen_min = 280.0, wavelen_max = 1170.0)
        testSED.multiplyFluxNorm(testSED.calcFluxNorm(10, testSEDsColors.starPhot.bandPassList[2]))
        #Give mlt like filename just so it works
        testSED.writeSED(self.mTestName)

        testSEDsColors.mltDir = str(os.getcwd() + '/')
        testColors = testSEDsColors.loadmltSEDs(subset = [self.mTestName])
        self.assertAlmostEqual(testColors['umg'][0], 0.0)
        self.assertAlmostEqual(testColors['gmr'][0], 0.0)
        self.assertAlmostEqual(testColors['rmi'][0], 0.0)
        self.assertAlmostEqual(testColors['imz'][0], 0.0)

    def testLoadWD(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for and
        values for each property have been calculated."""
        #Test Matching to WD SEDs
        loadTestWD = selectStarSED(wdDir = self.testWDDir)
        testSEDs = loadTestWD.loadwdSEDs()

        #Add extra step because WD SEDs are separated into helium and hydrogen
        testSEDNamesLists = []
        testSEDNamesLists.append(testSEDs['H']['sEDName'])
        testSEDNamesLists.append(testSEDs['HE']['sEDName'])
        testSEDNames = [name for nameList in testSEDNamesLists for name in nameList]

        #Read in a list of the SEDs in the wd test sed directory
        testWDList = os.listdir(self.testWDDir)

        #First make sure that all SEDs are correctly accounted for if no subset provided
        self.assertItemsEqual(testSEDNames, testWDList)

        #Then test to make sure no values are null in any of the lists within the dictionary
        for wdType in ['H', 'HE']:
            for key in testSEDs[wdType].keys():
                for item in testSEDs[wdType][key]:
                    self.assertIsNotNone(item)

        #Test same condition if subset is provided
        testSubsetList = ['bergeron_10000_75.dat_10100.gz', 'bergeron_He_9000_80.dat_9400.gz']

        testSEDsSubset = selectStarSED().loadwdSEDs(subset = testSubsetList)

        #Add extra step because WD SEDs are separated into helium and hydrogen
        testSEDNamesSubsetLists = []
        testSEDNamesSubsetLists.append(testSEDsSubset['H']['sEDName'])
        testSEDNamesSubsetLists.append(testSEDsSubset['HE']['sEDName'])
        testSEDNamesSubset = [name for nameList in testSEDNamesSubsetLists for name in nameList]

        #Next make sure that correct subset loads if subset is provided
        self.assertItemsEqual(testSEDNamesSubset, testSubsetList)

        #Make sure that the names get separated into correct wd type
        self.assertItemsEqual(testSEDsSubset['H']['sEDName'], [testSubsetList[0]])
        self.assertItemsEqual(testSEDsSubset['HE']['sEDName'], [testSubsetList[1]])

        #Now test colors are being calculated correctly by loading a flat sed into it
        testSED = Sed()
        testSEDsColors = selectStarSED()

        testSED.setFlatSED(wavelen_min = 280.0, wavelen_max = 1170.0)
        testSED.multiplyFluxNorm(testSED.calcFluxNorm(10, testSEDsColors.starPhot.bandPassList[2]))
        #Give WD like filename just so it works
        testName = 'bergeron_9999_99.dat_9999'
        testNameHe = 'bergeron_He_9999_99.dat_9999'
        testSED.writeSED(testName)
        testSED.writeSED(testNameHe)

        testSEDsColors.wdDir = str(os.getcwd() + '/')
        testColors = testSEDsColors.loadwdSEDs(subset = [testName, testNameHe])
        self.assertAlmostEqual(testColors['H']['umg'][0], 0.0)
        self.assertAlmostEqual(testColors['H']['gmr'][0], 0.0)
        self.assertAlmostEqual(testColors['H']['rmi'][0], 0.0)
        self.assertAlmostEqual(testColors['H']['imz'][0], 0.0)
        self.assertAlmostEqual(testColors['HE']['umg'][0], 0.0)
        self.assertAlmostEqual(testColors['HE']['gmr'][0], 0.0)
        self.assertAlmostEqual(testColors['HE']['rmi'][0], 0.0)
        self.assertAlmostEqual(testColors['HE']['imz'][0], 0.0)

        os.unlink(testName)
        os.unlink(testNameHe)

    def testDeReddenGalfast(self):

        """Test that consistent numbers come out of deReddening procedure"""

        am = 0.5
        coeffs = np.ones(5)
        mags = np.arange(2,-3,-1)

        testDeRedColors = selectStarSED().deReddenGalfast(am, mags[0], mags[1], mags[2], mags[3],
                                                          mags[4], coeffs)

        #Test Output
        expectedDeRed = np.ones(4)
        np.testing.assert_equal(testDeRedColors, expectedDeRed)

    def testFindSED(self):

        """Pull SEDs from each type and make sure that each SED gets matched to itself."""

        testMatching = selectStarSED(sEDDir = self.testSpecDir, kuruczDir = self.testKDir,
                                     mltDir = self.testMLTDir, wdDir = self.testWDDir)

        testSubsetList = []
        testSubsetType = []
        testSubsetComp = []
        testOutputName = []
        testOutputNameReddened = []
        #Populate Lists
        for fileName in os.listdir(testMatching.kuruczDir)[0:10]:
            testSubsetList.append(fileName)
            testSubsetType.append('kurucz')
            testSubsetComp.append(1)
        for fileName in os.listdir(testMatching.mltDir)[0:10]:
            testSubsetList.append(fileName)
            testSubsetType.append('mlt')
            testSubsetComp.append(1)
        #This is to get both H and HE WDs
        wdLists = [os.listdir(testMatching.wdDir)[-10:], os.listdir(testMatching.wdDir)[:10]]
        wdNames = [wdName for wdList in wdLists for wdName in wdList]
        for fileName in wdNames:
            testSubsetList.append(fileName)
            testSubsetType.append('wDs')
            if 'He' in fileName:
                testSubsetComp.append(15)
            else:
                testSubsetComp.append(10)

        testMatchingDict = {}
        testMatchingDict['kurucz'] = testMatching.loadKuruczSEDs()
        testMatchingDict['mlt'] = testMatching.loadmltSEDs()
        testMatchingWDs = testMatching.loadwdSEDs()
        testMatchingDict['wdH'] = testMatchingWDs['H']
        testMatchingDict['wdHE'] = testMatchingWDs['HE']

        specMap = SpecMap()

        for testSEDNum in range(0, len(testSubsetList)):

            if testSEDNum % 10 == 0:
                print 'Calculating test colors for SED %i of %i' % (testSEDNum, len(testSubsetList))

            getSEDColors = Sed()
            testSEDName = testSubsetList[testSEDNum].strip('.gz')
            getSEDColors.readSED_flambda(str(testMatching.sEDDir + '/' +
                                             str(specMap.__getitem__(testSEDName))))
            getSEDColors.multiplyFluxNorm(getSEDColors.calcFluxNorm(10, testMatching.starPhot.bandPassList[2]))
            testMagList = testMatching.starPhot.manyMagCalc_list(getSEDColors)

            #First test without reddening
            testOutputName.append(testMatching.findSED(testMatchingDict, testMagList[0], testMagList[1],
                                                       testMagList[2], testMagList[3], testMagList[4],
                                                       0, testSubsetComp[testSEDNum], reddening = False))
            #Next test with reddening and custom coeffs
            am = 0.5
            reddenCoeffs = np.array([1.2, 1.1, 1.0, 0.9, 0.8])
            testReddening = am * reddenCoeffs
            testReddenedMagList = testMagList + testReddening
            testOutputNameReddened.append(testMatching.findSED(testMatchingDict, testReddenedMagList[0],
                                                               testReddenedMagList[1],
                                                               testReddenedMagList[2],
                                                               testReddenedMagList[3],
                                                               testReddenedMagList[4], am,
                                                               testSubsetComp[testSEDNum], True,
                                                               reddenCoeffs))
        self.assertEqual(testOutputName, testSubsetList)
        self.assertEqual(testOutputNameReddened, testSubsetList)

class TestReadGalfast(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        #because SCons only knows about a limited subset of the LSST environment variables
        os.environ['LSST_THROUGHPUTS_DEFAULT'] = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline')
        os.environ['SDSS_THROUGHPUTS'] = os.path.join(os.getenv('THROUGHPUTS_DIR'),'sdss')

        #Left this in after removing loading SEDs so that we can make sure that if the structure of
        #sims_sed_library changes in a way that affects testReadGalfast we can detect it.
        specMap= SpecMap()
        cls._specMapDict = {}
        specFileStart = ['kp', 'burrows', 'bergeron'] #The beginning of filenames of different SED types
        specFileTypes = ['kurucz', 'mlt','wd']
        for specStart, specKey in zip(specFileStart, specFileTypes):
            for key, val in sorted(specMap.subdir_map.iteritems()):
                if re.match(key, specStart):
                    cls._specMapDict[specKey] = str(val)

        #Set up Test Spectra Directory
        cls.testSpecDir = 'testReadGalfastSpectra'
        cls.testKDir = str(cls.testSpecDir + '/starSED/kurucz/')
        cls.testMLTDir = str(cls.testSpecDir + '/starSED/mlt/')
        cls.testWDDir = str(cls.testSpecDir + '/starSED/wDs/')

        if os.path.exists(cls.testSpecDir):
            shutil.rmtree(cls.testSpecDir)

        os.makedirs(cls.testKDir)
        os.mkdir(cls.testMLTDir)
        os.mkdir(cls.testWDDir)
        cls.kDir = os.environ['SIMS_SED_LIBRARY_DIR'] + '/' + cls._specMapDict['kurucz'] + '/'
        cls.mltDir = os.environ['SIMS_SED_LIBRARY_DIR'] + '/' + cls._specMapDict['mlt'] + '/'
        cls.wdDir = os.environ['SIMS_SED_LIBRARY_DIR'] + '/' + cls._specMapDict['wd'] + '/'
        #Use particular indices to get different types of seds within mlt and wds
        for kFile, mltFile, wdFile in zip(os.listdir(cls.kDir)[0:20], 
                                          np.array(os.listdir(cls.mltDir))[np.arange(-10,11)], 
                                          np.array(os.listdir(cls.wdDir))[np.arange(-10,11)]):
            shutil.copyfile(str(cls.kDir + kFile), str(cls.testKDir + kFile))
            shutil.copyfile(str(cls.mltDir + mltFile), str(cls.testMLTDir + mltFile))
            shutil.copyfile(str(cls.wdDir + wdFile), str(cls.testWDDir + wdFile))

    @classmethod
    def tearDownClass(cls):
        del cls._specMapDict

        if os.path.exists('exampleOutput.txt'):
            os.unlink('exampleOutput.txt')

        if os.path.exists('exampleOutputGzip.txt'):
            os.unlink('exampleOutputGzip.txt')

        if os.path.exists('exampleOutputFits.txt'):
            os.unlink('exampleOutputFits.txt')

        if os.path.exists('example.txt'):
             os.unlink('example.txt')

        if os.path.exists('gzipExample.txt.gz'):
            os.unlink('gzipExample.txt.gz')

        if os.path.exists('exampleFits.fits'):
             os.unlink('exampleFits.fits')

        shutil.rmtree(cls.testSpecDir)

    def testParseGalfast(self):

        """Test Read-in of Galfast Header"""

        testRG = readGalfast()
        #First test that exception is raised when an invalid header label is used
        testInvalidHeader = 'lb[2] header'
        self.assertRaises(RuntimeError, testRG.parseGalfast, testInvalidHeader)

        #Next test that '#' is ignored
        testSymbolHeader = '# lb[2]'
        testSymbolDict = testRG.parseGalfast(testSymbolHeader)
        self.assertEqual(testSymbolDict, {'l':0, 'b':1})

        #Test that new line marker at end of line is ignored
        testNewLineHeader = '# lb[2] \n '
        testNewLineDict = testRG.parseGalfast(testNewLineHeader)
        self.assertEqual(testNewLineDict, {'l':0, 'b':1})

        #Next test that extra spaces are ignored
        testSpaceHeader = 'lb[2]     radec[2]'
        testSpaceDict = testRG.parseGalfast(testSpaceHeader)
        self.assertEqual(testSpaceDict, {'l':0, 'b':1, 'ra':2, 'dec':3})

        #Test that every header value gets put in the right place using different order than usually come out
        testFullHeader = '# lb[2] XYZ[3] radec[2] absSDSSr{alias=M1;alias=absmag;band=SDSSr;} ' +\
                         'DM comp FeH vcyl[3] pmlb[3] pmradec[3] Am AmInf SDSSugriz[5]{class=magnitude;' +\
                         'fieldNames=0:SDSSu,1:SDSSg,2:SDSSr,3:SDSSi,4:SDSSz;} ' +\
                         'SDSSugrizPhotoFlags{class=flags;}'
        actualFullHeaderDict = {'l':0, 'b':1, 'X':2, 'Y':3, 'Z':4, 'ra':5, 'dec':6, 'absSDSSr':7, 'DM':8,
                                'comp':9, 'FeH':10, 'Vr':11, 'Vphi':12, 'Vz':13, 'pml':14, 'pmb':15,
                                'vRadlb':16, 'pmra':17, 'pmdec':18, 'vRad':19, 'Am':20, 'AmInf':21,
                                'SDSSu':22, 'SDSSg':23, 'SDSSr':24, 'SDSSi':25, 'SDSSz':26,
                                'SDSSPhotoFlags':27}
        testFullHeaderDict = testRG.parseGalfast(testFullHeader)
        self.assertEqual(testFullHeaderDict, actualFullHeaderDict)

    def testFindLSSTMags(self):

        """Make sure that the Magnitudes you get from spectra are correct"""

        testRG = readGalfast()
        sEDParams = selectStarSED()
        absLSSTr = 10.
        DM = 10.
        am = 0.5
        lsstExtCoords = [1.25, 1.15, 1.05, 0.95, 0.85, 0.75]
        sEDDict = {}

        #Test a kurucz, an mlt, and a wd
        testSpectra = ['bergeron_10000_75.dat_10100.gz', 'm0.0Full.dat.gz', 'km01_7000.fits_g40_7140.gz']
        testTypes = ['wdH', 'mlt', 'kurucz']
        testDirs = ['wDs', 'mlt', 'kurucz']
        #Build SedDict
        for testSpectrum, testType, testDir in zip(testSpectra, testTypes, testDirs):
            sEDObj = Sed()
            sEDObj.readSED_flambda(sEDParams.sEDDir + '/starSED/' + testDir + '/' + testSpectrum)
            rMagDict = {}
            testTypeDict = {}
            lsstgmrDict = {}; lsstumgDict = {}; lsstrmiDict = {}; lsstimzDict = {}; lsstzmyDict = {};
            lsstrMagsDict = {}
            rMagDict[testSpectrum] = sEDObj.calcMag(sEDParams.lsstPhot.bandPassList[2])
            testTypeDict['rMags'] = rMagDict
            lsstTestMagList =  sEDParams.lsstPhot.manyMagCalc_list(sEDObj)
            lsstgmrDict[testSpectrum] = lsstTestMagList[1] - lsstTestMagList[2]
            lsstumgDict[testSpectrum] = lsstTestMagList[0] - lsstTestMagList[1]
            lsstrmiDict[testSpectrum] = lsstTestMagList[2] - lsstTestMagList[3]
            lsstimzDict[testSpectrum] = lsstTestMagList[3] - lsstTestMagList[4]
            lsstzmyDict[testSpectrum] = lsstTestMagList[4] - lsstTestMagList[5]
            lsstrMagsDict[testSpectrum] = lsstTestMagList[2]
            testTypeDict['lsstgmr'] = lsstgmrDict
            testTypeDict['lsstumg'] = lsstumgDict
            testTypeDict['lsstrmi'] = lsstrmiDict
            testTypeDict['lsstimz'] = lsstimzDict
            testTypeDict['lsstzmy'] = lsstzmyDict
            testTypeDict['lsstrMags'] = lsstrMagsDict
            sEDDict[testType] = testTypeDict

        for testSpectrum, testType, testDir in zip(testSpectra, testTypes, testDirs):
            sEDObj = Sed()
            sEDObj.readSED_flambda(sEDParams.sEDDir + '/starSED/' + testDir + '/' + testSpectrum)
            #Calculate the SDSSr if the SED's LSSTr is set to absLSSTr above
            #Then Make sure with this input you get the same LSSTr out of findLSSTMags
            lsstFluxNorm = sEDObj.calcFluxNorm(absLSSTr, sEDParams.lsstPhot.bandPassList[2])
            sEDObj.multiplyFluxNorm(lsstFluxNorm)
            absSDSSr = sEDObj.calcMag(sEDParams.lsstPhot.bandPassList[2])
            testFluxNorm, testMagDict = testRG.findLSSTMags(testSpectrum, sEDDict, absSDSSr,
                                                            DM, am, lsstExtCoords)
            self.assertAlmostEqual(testFluxNorm, lsstFluxNorm)
            self.assertAlmostEqual(testMagDict['r'] - DM - (am * lsstExtCoords[2]), absLSSTr)

    def testConvDMtoKpc(self):

        """Make sure Distance Modulus get correctly converted to distance in kpc"""

        testRG = readGalfast()
        dm = 10.
        actualKpcDistance = 1.
        testKpcDistance = testRG.convDMtoKpc(dm)
        self.assertEqual(testKpcDistance, actualKpcDistance)

    def testLoadGalfast(self):

        """Make sure all desired input file types are correctly processed and return correct output files"""

        testRG = readGalfast()
        #First test that it makes sure file exists
        self.assertRaises(RuntimeError, testRG.loadGalfast, ['notarealfile.txt'], ['noOutput.txt'])

        #Next test that if an unknown file format is entered it exits
        self.assertRaises(RuntimeError, testRG.loadGalfast, ['notarealfile.dat'], ['noOutput.txt'])

        #Write example files and then load in and make sure example output files are created
        #First .txt
        exampleIn = open('example.txt', 'w')
        inHeader = '# lb[2] XYZ[3] radec[2] absSDSSr{alias=M1;alias=absmag;band=SDSSr;} DM comp FeH ' +\
                   'vcyl[3] pmlb[3] pmradec[3] Am AmInf SDSSugriz[5]{class=magnitude;fieldNames=0:SDSSu,' +\
                   '1:SDSSg,2:SDSSr,3:SDSSi,4:SDSSz;} SDSSugrizPhotoFlags{class=flags;} \n'
        testComment = '# Comment\n'
        inData = '   1.79371816  -89.02816704   11.92064832  -27.62775082       7.15       0.22   ' +\
                 '-421.87   8.126   4.366   0 -0.095    13.7  -183.4    -6.2   -20.58   -12.60    ' +\
                 '13.02    21.34   -11.26    13.02  0.037  0.037  14.350  12.949  12.529  12.381  12.358 0'
        exampleIn.write(inHeader)
        exampleIn.write(testComment)
        exampleIn.write(inData)
        exampleIn.close()

        #Then gzipped
        exampleGzipIn = gzip.open('gzipExample.txt.gz', 'w')
        exampleGzipIn.write(inHeader)
        exampleGzipIn.write(testComment)
        exampleGzipIn.write(inData)
        exampleGzipIn.close()

        #Finally a fits file, but first make sure to remove pre-existing file
        if os.path.isfile('exampleFits.fits'):
            os.remove('exampleFits.fits')
        columnNames = ['lb', 'XYZ', 'radec', 'absSDSSr', 'DM', 'comp', 'FeH', 'vcyl', 'pmlb', 'pmradec',
                       'Am', 'AmInf', 'SDSSugriz', 'SDSSugrizPhotoFlags']
        columnArrays = [[[1.79371816, -89.02816704]],  [[7.15, 0.22, -421.87]],
                        [[11.92064832, -27.62775082]], [[8.126]], [[4.366]], [[0]], [[-0.095]],
                        [[13.7,  -183.4, -6.2]], [[-20.58, -12.60, 13.02]], [[21.34, -11.26, 13.02]],
                        [[0.037]], [[0.037]],  [[14.350, 12.949, 12.529, 12.381, 12.358]], [[0]]]
        columnFormats = ['2E', '3E', '2E', 'E', 'E', 'E', 'E', '3E', '3E', '3E', 'E', 'E', '5E', 'E']
        cols = pyfits.ColDefs([pyfits.Column(name = columnNames[0], format = columnFormats[0],
                                             array = columnArrays[0])])
        for colName, colArray, colFormat in zip(columnNames[1:], columnArrays[1:], columnFormats[1:]):
            cols.add_col(pyfits.Column(name = colName, format = colFormat, array = colArray))
        exampleTable = pyfits.new_table(cols)
        exampleTable.writeto('exampleFits.fits')
        testRG.loadGalfast(['example.txt', 'gzipExample.txt.gz', 'exampleFits.fits'],
                           ['exampleOutput.txt', 'exampleOutputGzip.txt', 'exampleOutputFits.txt'],
                           kuruczPath = self.testKDir,
                           mltPath = self.testMLTDir,
                           wdPath = self.testWDDir)
        self.assertTrue(os.path.isfile('exampleOutput.txt'))
        self.assertTrue(os.path.isfile('exampleOutputGzip.txt'))
        self.assertTrue(os.path.isfile('exampleOutputFits.txt'))

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(TestSelectGalaxySED)
    suites += unittest.makeSuite(TestSelectStarSED)
    suites += unittest.makeSuite(TestReadGalfast)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)
if __name__ == "__main__":
    run(True)
