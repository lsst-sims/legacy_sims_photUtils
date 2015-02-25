import unittest
import os
import shutil
import numpy as np
import gzip
import pyfits
import re
import eups
import lsst.utils.tests as utilsTests
from lsst.sims.photUtils.readGalfast.selectStarSED import selectStarSED
from lsst.sims.photUtils.readGalfast.readGalfast import readGalfast
from lsst.sims.photUtils.readGalfast.selectGalaxySED import selectGalaxySED
from lsst.sims.photUtils.readGalfast.rgUtils import rgBase
from lsst.sims.photUtils.readGalfast.rgUtils import rgStar
from lsst.sims.photUtils.readGalfast.rgUtils import rgGalaxy
from lsst.sims.photUtils.EBV import EBVbase as ebv
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Photometry import PhotometryBase as phot
from lsst.sims.catalogs.measures.instance.fileMaps import SpecMap

class TestRGBase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        specMap = SpecMap()
        specFileStart = 'Exp'
        for key, val in sorted(specMap.subdir_map.iteritems()):
            if re.match(key, specFileStart):
                galSpecDir = str(val)
        cls.galDir = str(eups.productDir('sims_sed_library') + '/' + galSpecDir + '/')
        cls.filterList = ('u', 'g', 'r', 'i', 'z')

    @classmethod
    def tearDownClass(cls):
        del cls.galDir
        del cls.filterList

    def testCalcMagNorm(self):

        """Tests the calculation of magnitude normalization for an SED with the given magnitudes
        in the given bandpasses."""

        testUtils = rgBase()        
        testPhot = phot()
        testPhot.loadBandPassesFromFiles(self.filterList, 
                                        bandPassDir = os.path.join(eups.productDir('throughputs'),'sdss'),
                                        bandPassRoot = 'sdss_')
        testPhot.setupPhiArray_dict()

        unChangedSED = Sed()
        unChangedSED.readSED_flambda(str(self.galDir + os.listdir(self.galDir)[0]))

        imSimBand = Bandpass()
        imSimBand.imsimBandpass()
        testSED = Sed()
        testSED.setSED(unChangedSED.wavelen, flambda = unChangedSED.flambda)
        magNorm = 20.0
        redVal = 0.1
        testSED.redshiftSED(redVal)
        fluxNorm = testSED.calcFluxNorm(magNorm, imSimBand)
        testSED.multiplyFluxNorm(fluxNorm)
        sedMags = testPhot.manyMagCalc_list(testSED)
        stepSize = 0.001
        testMagNorm = testUtils.calcMagNorm(sedMags, unChangedSED, testPhot, 
                                            redshift = redVal, stepSize = stepSize)
        
        self.assertAlmostEqual(magNorm, testMagNorm, delta = stepSize)

    def testCalcBasicColors(self):

        """Tests the calculation of the colors of an SED in given bandpasses."""

        testUtils = rgBase()        
        testSED = Sed()
        testPhot = phot()
        testPhot.loadBandPassesFromFiles(self.filterList, 
                                        bandPassDir = os.path.join(eups.productDir('throughputs'),'sdss'),
                                        bandPassRoot = 'sdss_')
        testPhot.setupPhiArray_dict()

        testSED.readSED_flambda(str(self.galDir + os.listdir(self.galDir)[0]))
        testMags = testPhot.manyMagCalc_list(testSED)
        testColors = []
        for filtNum in range(0, len(self.filterList)-1):
            testColors.append(testMags[filtNum] - testMags[filtNum+1])

        testOutput = testUtils.calcBasicColors([testSED], testPhot)
        np.testing.assert_equal([testColors], testOutput)

    def testDeReddenMags(self):

        """Test that consistent numbers come out of deReddening procedure"""

        am = 0.5
        coeffs = np.ones(5)
        mags = np.arange(2,-3,-1)

        testDeRed = rgBase().deReddenMags(am, mags, coeffs)

        #Test Output
        np.testing.assert_equal(testDeRed,[ mags-(am*coeffs)])

class TestRGStar(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

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
        cls.kDir = eups.productDir('sims_sed_library') + '/' + cls._specMapDict['kurucz'] + '/'
        cls.mltDir = eups.productDir('sims_sed_library') + '/' + cls._specMapDict['mlt'] + '/'
        cls.wdDir = eups.productDir('sims_sed_library') + '/' + cls._specMapDict['wd'] + '/'
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

    def testDefaults(self):
        """Make sure that if there are Nones for the init that they load the correct dirs"""
        loadTest = rgStar()
        self.assertEqual(loadTest.kuruczDir, self.kDir)
        self.assertEqual(loadTest.mltDir, self.mltDir)
        self.assertEqual(loadTest.wdDir, self.wdDir)

    def testLoadKurucz(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for """
        #Test Matching to Kurucz SEDs
        loadTestKurucz = rgStar(kuruczDir = self.testKDir)
        testSEDs = loadTestKurucz.loadKuruczSEDs()

        #Read in a list of the SEDs in the kurucz test sed directory
        testKuruczList = os.listdir(self.testKDir)

        #First make sure that all SEDs are correctly accounted for if no subset provided
        testNames = []
        for testSED in testSEDs:
            testNames.append(testSED.name)
        self.assertItemsEqual(testKuruczList, testNames)

        #Test same condition if subset is provided
        testSubsetList = ['km01_7000.fits_g40_7140.gz', 'kp01_7000.fits_g40_7240.gz']
        testSEDsSubset = loadTestKurucz.loadKuruczSEDs(subset = testSubsetList)

        #Next make sure that correct subset loads if subset is provided
        testSubsetNames = []
        testSubsetLogZ = []
        testSubsetLogG = []
        testSubsetTemp = []
        for testSED in testSEDsSubset:
            testSubsetNames.append(testSED.name)
            testSubsetLogZ.append(testSED.logZ)
            testSubsetLogG.append(testSED.logg)
            testSubsetTemp.append(testSED.temp)
        self.assertItemsEqual(testSubsetList, testSubsetNames)
        self.assertEqual(testSubsetLogZ, [-0.1, 0.1]) #Test both pos. and neg. get in right
        self.assertEqual(testSubsetLogG, [4.0, 4.0]) #Test storage of logg and temp
        self.assertEqual(testSubsetTemp, [7140, 7240])

        #Test that attributes have been assigned
        for testSED in testSEDsSubset:
            self.assertIsNotNone(testSED.name)
            self.assertIsNotNone(testSED.logZ)
            self.assertIsNotNone(testSED.logg)
            self.assertIsNotNone(testSED.temp)

    def testLoadMLT(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for"""
        #Test Matching to mlt SEDs
        loadTestMLT = rgStar(mltDir = self.testMLTDir)
        testSEDs = loadTestMLT.loadmltSEDs()

        #Read in a list of the SEDs in the mlt test sed directory
        testMLTList = os.listdir(self.testMLTDir)

        #First make sure that all SEDs are correctly accounted for if no subset provided
        testNames = []
        for testSED in testSEDs:
            testNames.append(testSED.name)
        self.assertItemsEqual(testMLTList, testNames)

        #Next make sure that correct subset loads if subset is provided
        testSubsetList = testMLTList[0:2]
        testSEDsubset = loadTestMLT.loadmltSEDs(subset = testSubsetList)
        testSubsetNames = []
        for testSED in testSEDsubset:
            testSubsetNames.append(testSED.name)
        self.assertItemsEqual(testSubsetList, testSubsetNames)

        #Test that attributes have been assigned
        for testSED in testSEDsubset:
            self.assertIsNotNone(testSED.name)

    def testLoadWD(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for and
        that there are separate lists for H and HE."""
        #Test Matching to WD SEDs
        loadTestWD = rgStar(wdDir = self.testWDDir)
        testSEDsH, testSEDsHE = loadTestWD.loadwdSEDs()

        #Add extra step because WD SEDs are separated into helium and hydrogen
        testNames = []
        for testH in testSEDsH:
            testNames.append(testH.name)
        for testHE in testSEDsHE:
            testNames.append(testHE.name)

        #Read in a list of the SEDs in the wd test sed directory
        testWDList = os.listdir(self.testWDDir)

        #First make sure that all SEDs are correctly accounted for if no subset provided
        self.assertItemsEqual(testNames, testWDList)

        #Test same condition if subset is provided
        testSubsetList = ['bergeron_10000_75.dat_10100.gz', 'bergeron_He_9000_80.dat_9400.gz']

        testSEDsSubsetH, testSEDsSubsetHE = selectStarSED().loadwdSEDs(subset = testSubsetList)

        testNamesSubset = []
        for testH in testSEDsSubsetH:
            testNamesSubset.append(testH.name)
        for testHE in testSEDsSubsetHE:
            testNamesSubset.append(testHE.name)

        #Next make sure that correct subset loads if subset is provided
        self.assertItemsEqual(testNamesSubset, testSubsetList)

        #Make sure that the names get separated into correct wd type
        self.assertEqual(testSEDsSubsetH[0].name, testSubsetList[0])
        self.assertEqual(testSEDsSubsetHE[0].name, testSubsetList[1])

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

class TestRGGalaxy(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):

        specMap = SpecMap()
        specFileStart = 'Exp'
        for key, val in sorted(specMap.subdir_map.iteritems()):
            if re.match(key, specFileStart):
                galSpecDir = str(val)
        cls.galDir = str(eups.productDir('sims_sed_library') + '/' + galSpecDir + '/')

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
        """Test Loader for Bruzual and Charlot Galaxies"""
        loadTestBC03 = rgGalaxy(galDir = self.testSpecDir)
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

    @classmethod
    def tearDownClass(cls):
        del cls.galDir
        del cls.filterList

        shutil.rmtree(cls.testSpecDir)

class TestSelectGalaxySED(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        specMap = SpecMap()
        specFileStart = 'Exp'
        for key, val in sorted(specMap.subdir_map.iteritems()):
            if re.match(key, specFileStart):
                galSpecDir = str(val)
        cls.galDir = str(eups.productDir('sims_sed_library') + '/' + galSpecDir + '/')

        #Set up Test Spectra Directory
        cls.testSpecDir = 'testGalaxySEDSpectrum/'
        
        if os.path.exists(cls.testSpecDir):
            shutil.rmtree(cls.testSpecDir)

        os.mkdir(cls.testSpecDir)

        galList = os.listdir(cls.galDir)[0:20]

        for galFile in galList:
            shutil.copy(str(cls.galDir + galFile), str(cls.testSpecDir + galFile))

    def testMatchToRestFrame(self):
        """Test that Galaxies with no effects added into catalog mags are matched correctly."""
        galPhot = phot()
        galPhot.loadTotalBandpassesFromFiles()
        galPhot.setupPhiArray_dict()
        
        imSimBand = Bandpass()
        imSimBand.imsimBandpass()

        testMatching = selectGalaxySED(galDir = self.testSpecDir)
        testSEDList = testMatching.loadBC03()

        testSEDNames = []
        testMags = []
        testMagNormList = []
        magNormStep = 1

        for testSED in testSEDList:

            getSEDMags = Sed()
            testSEDNames.append(testSED.name)
            getSEDMags.setSED(wavelen = testSED.wavelen, flambda = testSED.flambda)
            testMagNorm = np.round(np.random.uniform(20.0,22.0),magNormStep)
            testMagNormList.append(testMagNorm)
            fluxNorm = getSEDMags.calcFluxNorm(testMagNorm, imSimBand)
            getSEDMags.multiplyFluxNorm(fluxNorm)
            testMags.append(galPhot.manyMagCalc_list(getSEDMags))

        #Also testing to make sure passing in non-default bandpasses works
        testMatchingResults = testMatching.matchToRestFrame(testSEDList, testMags, magNormAcc = magNormStep,
                                                            bandpassList = galPhot.bandPassList)

        self.assertEqual(testSEDNames, testMatchingResults[0])
        np.testing.assert_almost_equal(testMagNormList, testMatchingResults[1], decimal = magNormStep)

    def testMatchToObserved(self):
        """Test that Galaxy SEDs with extinction or redshift are matched correctly"""
        galPhot = phot()
        galPhot.loadTotalBandpassesFromFiles()
        galPhot.setupPhiArray_dict()
        
        imSimBand = Bandpass()
        imSimBand.imsimBandpass()

        testMatching = selectGalaxySED(galDir = self.testSpecDir)
        testSEDList = testMatching.loadBC03()

        testSEDNames = []
        testRA = []
        testDec = []
        testRedshifts = []
        testMagNormList = []
        magNormStep = 1
        extCoeffs = [1.8140, 1.4166, 0.9947, 0.7370, 0.5790, 0.4761]
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
            testZ = np.round(np.random.uniform(1.1,1.3),3)
            testRedshifts.append(testZ)
            testMagNorm = np.round(np.random.uniform(20.0,22.0),magNormStep)
            testMagNormList.append(testMagNorm)
            getRedshiftMags.setSED(wavelen = testSED.wavelen, flambda = testSED.flambda)
            getRedshiftMags.redshiftSED(testZ)
            fluxNorm = getRedshiftMags.calcFluxNorm(testMagNorm, imSimBand)
            getRedshiftMags.multiplyFluxNorm(fluxNorm)
            testMagsRedshift.append(galPhot.manyMagCalc_list(getRedshiftMags))
            
        #Will also test in passing of non-default bandpass
        testNoExtNoRedshift = testMatching.matchToObserved(testSEDList, testRA, testDec, np.zeros(20), 
                                                           testMags, reddening = False,
                                                           bandpassList = galPhot.bandpassDict.values())
        testMatchingEbvVals = testMatching.matchToObserved(testSEDList, testRA, testDec, np.zeros(20), 
                                                           testMagsExt,
                                                           reddening = True, extCoeffs = extCoeffs,
                                                           bandpassList = galPhot.bandPassList)
        testMatchingRedshift = testMatching.matchToObserved(testSEDList, testRA, testDec, testRedshifts,
                                                            testMagsRedshift, dzAcc = 3,
                                                            magNormAcc = magNormStep, reddening = False,
                                                            bandpassList = galPhot.bandPassList)

        self.assertEqual(testSEDNames, testNoExtNoRedshift[0])
        self.assertEqual(testSEDNames, testMatchingEbvVals[0])
        self.assertEqual(testSEDNames, testMatchingRedshift[0])
        np.testing.assert_almost_equal(testMagNormList, testMatchingRedshift[1], 
                                       decimal = magNormStep)

    @classmethod
    def tearDownClass(cls):
        del cls.galDir
        shutil.rmtree(cls.testSpecDir)

class TestSelectStarSED(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

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
        cls.kDir = eups.productDir('sims_sed_library') + '/' + cls._specMapDict['kurucz'] + '/'
        cls.mltDir = eups.productDir('sims_sed_library') + '/' + cls._specMapDict['mlt'] + '/'
        cls.wdDir = eups.productDir('sims_sed_library') + '/' + cls._specMapDict['wd'] + '/'
        kList = os.listdir(cls.kDir)[0:3]
        #Use particular indices to get different types of seds within mlt and wds
        for kFile, mltFile, wdFile in zip(kList, 
                                          np.array(os.listdir(cls.mltDir))[np.arange(-3,0)], 
                                          np.array(os.listdir(cls.wdDir))[np.arange(-1,2)]):
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

    def testFindSED(self):
        """Pull SEDs from each type and make sure that each SED gets matched to itself.
        Includes testing with extinction and passing in only colors."""
        starPhot = phot()
        starPhot.loadBandPassesFromFiles(('u','g','r','i','z'), 
                                        bandPassDir = os.path.join(eups.productDir('throughputs'),'sdss'),
                                        bandPassRoot = 'sdss_')
        starPhot.setupPhiArray_dict()
        
        imSimBand = Bandpass()
        imSimBand.imsimBandpass()

        testMatching = selectStarSED(sEDDir = self.testSpecDir, kuruczDir = self.testKDir,
                                     mltDir = self.testMLTDir, wdDir = self.testWDDir)
        testSEDList = []
        testSEDList.append(testMatching.loadKuruczSEDs())
        testSEDList.append(testMatching.loadmltSEDs())
        testSEDListH, testSEDListHE = testMatching.loadwdSEDs()
        testSEDList.append(testSEDListH)
        testSEDList.append(testSEDListHE)

        testSEDNames = []
        testMags = []
        testMagNormList = []
        magNormStep = 1

        for typeList in testSEDList:
            if len(typeList) != 0:
                typeSEDNames = []
                typeMags = []
                typeMagNorms = []
                for testSED in typeList:
                    getSEDMags = Sed()
                    typeSEDNames.append(testSED.name)
                    getSEDMags.setSED(wavelen = testSED.wavelen, flambda = testSED.flambda)
                    testMagNorm = np.round(np.random.uniform(20.0,22.0),magNormStep)
                    typeMagNorms.append(testMagNorm)
                    fluxNorm = getSEDMags.calcFluxNorm(testMagNorm, imSimBand)
                    getSEDMags.multiplyFluxNorm(fluxNorm)
                    typeMags.append(starPhot.manyMagCalc_list(getSEDMags))
                testSEDNames.append(typeSEDNames)
                testMags.append(typeMags)
                testMagNormList.append(typeMagNorms)
            
        fakeRA = np.ones(len(testSEDList[0]))
        fakeDec = np.ones(len(testSEDList[0]))

        #Since default bandPassList should be SDSS ugrizy shouldn't need to specify it
        for typeList, names, mags, magNorms in zip(testSEDList, testSEDNames, testMags, testMagNormList):
            testMatchingResults = testMatching.findSED(typeList, mags, fakeRA, fakeDec,
                                                       magNormAcc = magNormStep, reddening = False)
            self.assertEqual(names, testMatchingResults[0])
            np.testing.assert_almost_equal(magNorms, testMatchingResults[1], decimal = magNormStep)

        #Now test what happens if we pass in a bandPassList
        testMatchingResultsNoDefault = testMatching.findSED(testSEDList[0], testMags[0], fakeRA, fakeDec,
                                                            bandpassList = starPhot.bandPassList,
                                                            reddening = False)
        self.assertEqual(testSEDNames[0], testMatchingResultsNoDefault[0])
        np.testing.assert_almost_equal(testMagNormList[0], testMatchingResultsNoDefault[1], 
                                       decimal = magNormStep)
        
        #Test Reddening
        testRA = np.random.uniform(10,170,len(testSEDList[0]))
        testDec = np.random.uniform(10,80,len(testSEDList[0]))
        extFactor = .5
        raDec = np.array((testRA, testDec))
        ebvVals = ebv().calculateEbv(equatorialCoordinates = raDec)
        extVals = ebvVals*extFactor
        testRedMags = []
        for extVal, testMagSet in zip(extVals, testMags[0]):
            testRedMags.append(testMagSet + extVal)
        testMatchingResultsRed = testMatching.findSED(testSEDList[0], testRedMags, testRA, testDec,
                                                      extCoeffs = np.ones(5)*extFactor)
        self.assertEqual(testSEDNames[0], testMatchingResultsRed[0])
        np.testing.assert_almost_equal(testMagNormList[0], testMatchingResultsRed[1], 
                                       decimal = magNormStep)
        
        #Finally, test color input
        testColors = []
        for testMagSet in testMags[0]:
            testColorSet = []
            for filtNum in range(0, len(starPhot.bandPassList)-1):
                testColorSet.append(testMagSet[filtNum] - testMagSet[filtNum+1])
            testColors.append(testColorSet)
        testMatchingColorsInput = testMatching.findSED(testSEDList[0], testMags[0], testRA, testDec,
                                                       reddening = False, colors = testColors)
        self.assertEqual(testSEDNames[0], testMatchingColorsInput[0])
        np.testing.assert_almost_equal(testMagNormList[0], testMatchingColorsInput[1], 
                                       decimal = magNormStep)

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

class TestReadGalfast(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

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
        cls.kDir = eups.productDir('sims_sed_library') + '/' + cls._specMapDict['kurucz'] + '/'
        cls.mltDir = eups.productDir('sims_sed_library') + '/' + cls._specMapDict['mlt'] + '/'
        cls.wdDir = eups.productDir('sims_sed_library') + '/' + cls._specMapDict['wd'] + '/'
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

        if os.path.exists('exampleOutputGzip.txt.gz'):
            os.unlink('exampleOutputGzip.txt.gz')

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
        inHeader = '# lb[2] radec[2] XYZ[3] absSDSSr{alias=M1;alias=absmag;band=SDSSr;} DM comp FeH ' +\
                   'vcyl[3] pmlb[3] pmradec[3] Am AmInf SDSSugriz[5]{class=magnitude;fieldNames=0:SDSSu,' +\
                   '1:SDSSg,2:SDSSr,3:SDSSi,4:SDSSz;} SDSSugrizPhotoFlags{class=flags;} \n'
        testComment = '# Comment\n'
        inData = '   1.79371816  -89.02816704   11.92064832  -27.62775082       7.15       0.22   ' +\
                 '-421.87   8.126   4.366   0 -0.095    13.7  -183.4    -6.2   -20.58   -12.60    ' +\
                 '13.02    21.34   -11.26    13.02  0.037  0.037  14.350  12.949  12.529  12.381  12.358 0\n'
        exampleIn.write(inHeader)
        exampleIn.write(testComment)
        exampleIn.write(inData)
        exampleIn.close()

        #Then gzipped. Also testing multiple lines in catalog.
        exampleGzipIn = gzip.open('gzipExample.txt.gz', 'w')
        exampleGzipIn.write(inHeader)
        exampleGzipIn.write(testComment)
        exampleGzipIn.write(inData)
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
                           ['exampleOutput.txt', 'exampleOutputGzip.txt.gz', 'exampleOutputFits.txt'],
                           kuruczPath = self.testKDir,
                           mltPath = self.testMLTDir,
                           wdPath = self.testWDDir)
        self.assertTrue(os.path.isfile('exampleOutput.txt'))
        self.assertTrue(os.path.isfile('exampleOutputGzip.txt.gz'))
        self.assertTrue(os.path.isfile('exampleOutputFits.txt'))

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(TestRGBase)
    suites += unittest.makeSuite(TestRGStar)
    suites += unittest.makeSuite(TestRGGalaxy)
    suites += unittest.makeSuite(TestSelectGalaxySED)
    suites += unittest.makeSuite(TestSelectStarSED)
    suites += unittest.makeSuite(TestReadGalfast)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)
if __name__ == "__main__":
    run(True)
