import unittest
import os
import shutil
import numpy as np
import re
import lsst.utils
import lsst.utils.tests as utilsTests
from lsst.sims.photUtils.selectStarSED import selectStarSED
from lsst.sims.photUtils.selectGalaxySED import selectGalaxySED
from lsst.sims.photUtils.matchUtils import matchBase
from lsst.sims.photUtils.matchUtils import matchStar
from lsst.sims.photUtils.matchUtils import matchGalaxy
from lsst.sims.photUtils.EBV import EBVbase as ebv
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils import BandpassDict
from lsst.sims.utils import SpecMap

class TestMatchBase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        specMap = SpecMap()
        specFileStart = 'Exp'
        for key, val in sorted(specMap.subdir_map.iteritems()):
            if re.match(key, specFileStart):
                galSpecDir = str(val)
        cls.galDir = str(lsst.utils.getPackageDir('sims_sed_library') + '/' + galSpecDir + '/')
        cls.filterList = ('u', 'g', 'r', 'i', 'z')

    @classmethod
    def tearDownClass(cls):
        del cls.galDir
        del cls.filterList

    def testCalcMagNorm(self):

        """Tests the calculation of magnitude normalization for an SED with the given magnitudes
        in the given bandpasses."""

        testUtils = matchBase()
        testPhot = BandpassDict.loadTotalBandpassesFromFiles(self.filterList,
                                        bandpassDir = os.path.join(lsst.utils.getPackageDir('throughputs'),'sdss'),
                                        bandpassRoot = 'sdss_')

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
        sedMags = testPhot.magListForSed(testSED)
        stepSize = 0.001
        testMagNorm = testUtils.calcMagNorm(sedMags, unChangedSED, testPhot, redshift = redVal)
        #Test adding in mag_errors. If an array of np.ones is passed in we should get same result
        testMagNormWithErr = testUtils.calcMagNorm(sedMags, unChangedSED, testPhot,
                                                   mag_error = np.ones(len(sedMags)), redshift = redVal)
        #Also need to add in test for filtRange
        sedMagsIncomp = sedMags
        sedMagsIncomp[1] = None
        filtRangeTest = [0, 2, 3, 4]
        testMagNormFiltRange = testUtils.calcMagNorm(sedMagsIncomp, unChangedSED, testPhot,
                                                     redshift = redVal, filtRange = filtRangeTest)
        self.assertAlmostEqual(magNorm, testMagNorm, delta = stepSize)
        self.assertAlmostEqual(magNorm, testMagNormWithErr, delta = stepSize)
        self.assertAlmostEqual(magNorm, testMagNormFiltRange, delta = stepSize)

    def testCalcBasicColors(self):

        """Tests the calculation of the colors of an SED in given bandpasses."""

        testUtils = matchBase()
        testSED = Sed()
        testPhot = BandpassDict.loadTotalBandpassesFromFiles(self.filterList,
                                        bandpassDir = os.path.join(lsst.utils.getPackageDir('throughputs'),'sdss'),
                                        bandpassRoot = 'sdss_')

        testSED.readSED_flambda(str(self.galDir + os.listdir(self.galDir)[0]))
        testMags = testPhot.magListForSed(testSED)
        testColors = []
        for filtNum in range(0, len(self.filterList)-1):
            testColors.append(testMags[filtNum] - testMags[filtNum+1])

        testOutput = testUtils.calcBasicColors([testSED], testPhot)
        np.testing.assert_equal([testColors], testOutput)

    def testSEDCopyBasicColors(self):

        """Tests that when makeCopy=True in calcBasicColors the SED object is unchanged after calling
        and that colors are still accurately calculated"""

        testUtils = matchBase()
        testSED = Sed()
        copyTest = Sed()
        testPhot = BandpassDict.loadTotalBandpassesFromFiles(self.filterList,
                                        bandpassDir = os.path.join(lsst.utils.getPackageDir('throughputs'),'sdss'),
                                        bandpassRoot = 'sdss_')
        testSED.readSED_flambda(str(self.galDir + os.listdir(self.galDir)[0]))
        copyTest.setSED(wavelen = testSED.wavelen, flambda = testSED.flambda)
        testLambda = copyTest.wavelen[0]
        testMags = testPhot.magListForSed(testSED)
        testColors = []
        for filtNum in range(0, len(self.filterList)-1):
            testColors.append(testMags[filtNum] - testMags[filtNum+1])
        testOutput = testUtils.calcBasicColors([copyTest], testPhot, makeCopy=True)

        self.assertEqual(testLambda, copyTest.wavelen[0])
        np.testing.assert_equal([testColors], testOutput)

    def testDeReddenMags(self):

        """Test that consistent numbers come out of deReddening procedure"""

        am = 0.5
        coeffs = np.ones(5)
        mags = np.arange(2,-3,-1)

        testDeRed = matchBase().deReddenMags(am, mags, coeffs)

        #Test Output
        np.testing.assert_equal(testDeRed,[ mags-(am*coeffs)])

class TestMatchStar(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        #Left this in after removing loading SEDs so that we can make sure that if the structure of
        #sims_sed_library changes in a way that affects testMatchSEDs we can detect it.
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
        cls.testSpecDir = 'testMatchingSpectra'
        cls.testKDir = str(cls.testSpecDir + '/starSED/kurucz/')
        cls.testMLTDir = str(cls.testSpecDir + '/starSED/mlt/')
        cls.testWDDir = str(cls.testSpecDir + '/starSED/wDs/')

        if os.path.exists(cls.testSpecDir):
            shutil.rmtree(cls.testSpecDir)

        os.makedirs(cls.testKDir)
        os.mkdir(cls.testMLTDir)
        os.mkdir(cls.testWDDir)
        cls.kDir = lsst.utils.getPackageDir('sims_sed_library') + '/' + cls._specMapDict['kurucz'] + '/'
        cls.mltDir = lsst.utils.getPackageDir('sims_sed_library') + '/' + cls._specMapDict['mlt'] + '/'
        cls.wdDir = lsst.utils.getPackageDir('sims_sed_library') + '/' + cls._specMapDict['wd'] + '/'
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
        loadTest = matchStar()
        self.assertEqual(loadTest.kuruczDir, self.kDir)
        self.assertEqual(loadTest.mltDir, self.mltDir)
        self.assertEqual(loadTest.wdDir, self.wdDir)

    def testLoadKurucz(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for """
        #Test Matching to Kurucz SEDs
        loadTestKurucz = matchStar(kuruczDir = self.testKDir)
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
        loadTestMLT = matchStar(mltDir = self.testMLTDir)
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
        loadTestWD = matchStar(wdDir = self.testWDDir)
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

class TestMatchGalaxy(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        specMap = SpecMap()
        specFileStart = 'Exp'
        for key, val in sorted(specMap.subdir_map.iteritems()):
            if re.match(key, specFileStart):
                galSpecDir = str(val)
        cls.galDir = str(lsst.utils.getPackageDir('sims_sed_library') + '/' + galSpecDir + '/')

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
        loadTestBC03 = matchGalaxy(galDir = self.testSpecDir)
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
        cls.galDir = str(lsst.utils.getPackageDir('sims_sed_library') + '/' + galSpecDir + '/')

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
        np.random.seed(42)
        galPhot = BandpassDict.loadTotalBandpassesFromFiles()

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
            testMags.append(galPhot.magListForSed(getSEDMags))

        #Also testing to make sure passing in non-default bandpasses works
        #Substitute in nan values to simulate incomplete data.
        testMags[0][1] = np.nan
        testMags[0][2] = np.nan
        testMags[0][4] = np.nan
        testMags[1][1] = np.nan
        testMatchingResults = testMatching.matchToRestFrame(testSEDList, testMags,
                                                            bandpassDict = galPhot)
        self.assertEqual(None, testMatchingResults[0][0])
        self.assertEqual(testSEDNames[1:], testMatchingResults[0][1:])
        self.assertEqual(None, testMatchingResults[1][0])
        np.testing.assert_almost_equal(testMagNormList[1:], testMatchingResults[1][1:], decimal = magNormStep)

        #Test Match Errors
        errMags = np.array((testMags[2], testMags[2], testMags[2], testMags[2]))
        errMags[1,1] += 1. #Total MSE will be 2/(5 colors) = 0.4
        errMags[2, 0:2] = np.nan
        errMags[2, 3] += 1. #Total MSE will be 2/(3 colors) = 0.667
        errMags[3, :] = None
        errSED = testSEDList[2]
        testMatchingResultsErrors = testMatching.matchToRestFrame([errSED], errMags,
                                                                  bandpassDict = galPhot)
        np.testing.assert_almost_equal(np.array((0.0, 0.4, 2./3.)), testMatchingResultsErrors[2][0:3],
                                       decimal = 3)
        self.assertEqual(None, testMatchingResultsErrors[2][3])

    def testReddeningException(self):
        """Test that if reddening=True in matchToObserved CatRA & CatDec are defined or exception is raised"""
        testException = selectGalaxySED(galDir = self.testSpecDir)
        testSEDList = testException.loadBC03()
        magnitudes = [[1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0]]
        redshifts = [1.0, 1.0]
        self.assertRaises(RuntimeError, testException.matchToObserved, testSEDList, magnitudes, redshifts,
                          reddening = True)

    def testMatchToObserved(self):
        """Test that Galaxy SEDs with extinction or redshift are matched correctly"""
        np.random.seed(42)
        galPhot = BandpassDict.loadTotalBandpassesFromFiles()

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
            testMags.append(galPhot.magListForSed(getSEDMags))

            #Check Extinction corrections
            sedRA = np.random.uniform(10,170)
            sedDec = np.random.uniform(10,80)
            testRA.append(sedRA)
            testDec.append(sedDec)
            raDec = np.array((sedRA, sedDec)).reshape((2,1))
            ebvVal = ebv().calculateEbv(equatorialCoordinates = raDec)
            extVal = ebvVal*extCoeffs
            testMagsExt.append(galPhot.magListForSed(getSEDMags) + extVal)

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
            testMagsRedshift.append(galPhot.magListForSed(getRedshiftMags))

        #Will also test in passing of non-default bandpass
        testNoExtNoRedshift = testMatching.matchToObserved(testSEDList, testMags, np.zeros(20),
                                                           reddening = False,
                                                           bandpassDict = galPhot)
        testMatchingEbvVals = testMatching.matchToObserved(testSEDList, testMagsExt, np.zeros(20),
                                                           catRA = testRA, catDec = testDec,
                                                           reddening = True, extCoeffs = extCoeffs,
                                                           bandpassDict = galPhot)
        #Substitute in nan values to simulate incomplete data and make sure magnorm works too.
        testMagsRedshift[0][1] = np.nan
        testMagsRedshift[0][3] = np.nan
        testMagsRedshift[0][4] = np.nan
        testMagsRedshift[1][1] = np.nan
        testMatchingRedshift = testMatching.matchToObserved(testSEDList, testMagsRedshift, testRedshifts,
                                                            dzAcc = 3, reddening = False,
                                                            bandpassDict = galPhot)

        self.assertEqual(testSEDNames, testNoExtNoRedshift[0])
        self.assertEqual(testSEDNames, testMatchingEbvVals[0])
        self.assertEqual(None, testMatchingRedshift[0][0])
        self.assertEqual(testSEDNames[1:], testMatchingRedshift[0][1:])
        self.assertEqual(None, testMatchingRedshift[1][0])
        np.testing.assert_almost_equal(testMagNormList[1:], testMatchingRedshift[1][1:],
                                       decimal = magNormStep)

        #Test Match Errors
        errMag = testMagsRedshift[2]
        errRedshift = testRedshifts[2]
        errMags = np.array((errMag, errMag, errMag, errMag))
        errRedshifts = np.array((errRedshift, errRedshift, errRedshift, errRedshift))
        errMags[1,1] += 1. #Total MSE will be 2/(5 colors) = 0.4
        errMags[2, 0:2] = np.nan
        errMags[2, 3] += 1. #Total MSE will be 2/(3 colors) = 0.667
        errMags[3, :] = None
        errSED = testSEDList[2]
        testMatchingResultsErrors = testMatching.matchToObserved([errSED], errMags, errRedshifts,
                                                                 reddening = False,
                                                                 bandpassDict = galPhot)
        print testMatchingResultsErrors
        np.testing.assert_almost_equal(np.array((0.0, 0.4, 2./3.)), testMatchingResultsErrors[2][0:3],
                                       decimal = 1) #Give a little more leeway due to redshifting effects
        self.assertEqual(None, testMatchingResultsErrors[2][3])

    @classmethod
    def tearDownClass(cls):
        del cls.galDir
        shutil.rmtree(cls.testSpecDir)

class TestSelectStarSED(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        #Left this in after removing loading SEDs so that we can make sure that if the structure of
        #sims_sed_library changes in a way that affects testMatchSEDs we can detect it.
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
        cls.testSpecDir = 'testMatchingSpectra'
        cls.testKDir = str(cls.testSpecDir + '/starSED/kurucz/')
        cls.testMLTDir = str(cls.testSpecDir + '/starSED/mlt/')
        cls.testWDDir = str(cls.testSpecDir + '/starSED/wDs/')

        if os.path.exists(cls.testSpecDir):
            shutil.rmtree(cls.testSpecDir)

        os.makedirs(cls.testKDir)
        os.mkdir(cls.testMLTDir)
        os.mkdir(cls.testWDDir)
        cls.kDir = lsst.utils.getPackageDir('sims_sed_library') + '/' + cls._specMapDict['kurucz'] + '/'
        cls.mltDir = lsst.utils.getPackageDir('sims_sed_library') + '/' + cls._specMapDict['mlt'] + '/'
        cls.wdDir = lsst.utils.getPackageDir('sims_sed_library') + '/' + cls._specMapDict['wd'] + '/'
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

    def testReddeningException(self):
        """Test that if reddening=True in matchToObserved CatRA & CatDec are defined or exception is raised"""
        testException = selectStarSED(sEDDir = self.testSpecDir, kuruczDir = self.testKDir,
                                      mltDir = self.testMLTDir, wdDir = self.testWDDir)
        testSEDList = testException.loadKuruczSEDs()
        magnitudes = [[1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0]]
        self.assertRaises(RuntimeError, testException.findSED, testSEDList, magnitudes,
                          reddening = True)

    def testFindSED(self):
        """Pull SEDs from each type and make sure that each SED gets matched to itself.
        Includes testing with extinction and passing in only colors."""
        np.random.seed(42)
        starPhot = BandpassDict.loadTotalBandpassesFromFiles(('u','g','r','i','z'),
                                        bandpassDir = os.path.join(lsst.utils.getPackageDir('throughputs'),'sdss'),
                                        bandpassRoot = 'sdss_')

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
                    typeMags.append(starPhot.magListForSed(getSEDMags))
                testSEDNames.append(typeSEDNames)
                testMags.append(typeMags)
                testMagNormList.append(typeMagNorms)

        fakeRA = np.ones(len(testSEDList[0]))
        fakeDec = np.ones(len(testSEDList[0]))

        #Since default bandpassDict should be SDSS ugrizy shouldn't need to specify it
        #Substitute in nan values to simulate incomplete data.
        for typeList, names, mags, magNorms in zip(testSEDList, testSEDNames, testMags, testMagNormList):
            if len(typeList) > 2:
                nanMags = np.array(mags)
                nanMags[0][0] = np.nan
                nanMags[0][2] = np.nan
                nanMags[0][3] = np.nan
                nanMags[1][1] = np.nan
                testMatchingResults = testMatching.findSED(typeList, nanMags, reddening = False)
                self.assertEqual(None, testMatchingResults[0][0])
                self.assertEqual(names[1:], testMatchingResults[0][1:])
                self.assertEqual(None, testMatchingResults[1][0])
                np.testing.assert_almost_equal(magNorms[1:], testMatchingResults[1][1:],
                                               decimal = magNormStep)
            else:
                testMatchingResults = testMatching.findSED(typeList, mags, reddening = False)
                self.assertEqual(names, testMatchingResults[0])
                np.testing.assert_almost_equal(magNorms, testMatchingResults[1], decimal = magNormStep)

        #Test Null Values option
        nullMags = np.array(testMags[0])
        nullMags[0][0] = -99.
        nullMags[0][4] = -99.
        nullMags[1][0] = -99.
        nullMags[1][1] = -99.
        testMatchingResultsNull = testMatching.findSED(testSEDList[0], nullMags,
                                                       nullValues = -99., reddening = False)
        self.assertEqual(testSEDNames[0], testMatchingResultsNull[0])
        np.testing.assert_almost_equal(testMagNormList[0], testMatchingResultsNull[1],
                                       decimal = magNormStep)

        #Test Error Output
        errMags = np.array((testMags[0][0], testMags[0][0], testMags[0][0], testMags[0][0]))
        errMags[1,1] += 1. #Total MSE will be 2/(4 colors) = 0.5
        errMags[2, 0:2] = np.nan
        errMags[2, 3] += 1. #Total MSE will be 2/(2 colors) = 1.0
        errMags[3, :] = None
        errSED = testSEDList[0][0]
        testMatchingResultsErrors = testMatching.findSED([errSED], errMags, reddening = False)
        np.testing.assert_almost_equal(np.array((0.0, 0.5, 1.0)), testMatchingResultsErrors[2][0:3],
                                       decimal = 3)
        self.assertEqual(None, testMatchingResultsErrors[2][3])

        #Now test what happens if we pass in a bandpassDict
        testMatchingResultsNoDefault = testMatching.findSED(testSEDList[0], testMags[0],
                                                            bandpassDict = starPhot,
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
        testMatchingResultsRed = testMatching.findSED(testSEDList[0], testRedMags, catRA = testRA,
                                                      catDec = testDec, reddening = True,
                                                      extCoeffs = np.ones(5)*extFactor)
        self.assertEqual(testSEDNames[0], testMatchingResultsRed[0])
        np.testing.assert_almost_equal(testMagNormList[0], testMatchingResultsRed[1],
                                       decimal = magNormStep)

        #Finally, test color input
        testColors = []
        for testMagSet in testMags[0]:
            testColorSet = []
            for filtNum in range(0, len(starPhot)-1):
                testColorSet.append(testMagSet[filtNum] - testMagSet[filtNum+1])
            testColors.append(testColorSet)
        testMatchingColorsInput = testMatching.findSED(testSEDList[0], testMags[0],
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

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(TestMatchBase)
    suites += unittest.makeSuite(TestMatchStar)
    suites += unittest.makeSuite(TestMatchGalaxy)
    suites += unittest.makeSuite(TestSelectGalaxySED)
    suites += unittest.makeSuite(TestSelectStarSED)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)
if __name__ == "__main__":
    run(True)
