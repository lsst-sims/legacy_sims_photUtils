import numpy

import os
import unittest
import lsst.utils.tests as utilsTests

from lsst.sims.catalogs.measures.instance import InstanceCatalog, register_method, register_class, compound
from lsst.sims.catalogs.generation.db import ObservationMetaData
from lsst.sims.coordUtils import AstrometryStars, AstrometryGalaxies
from lsst.sims.catalogs.generation.utils import myTestGals, myTestStars, \
                                                makeStarTestDB, makeGalTestDB, getOneChunk

from lsst.sims.photUtils.Photometry import PhotometryGalaxies, PhotometryStars

from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.EBV import EBVbase, EBVmixin

from lsst.sims.photUtils.Variability import Variability

from lsst.sims.photUtils.utils import MyVariability, testDefaults, cartoonPhotometryStars, \
                                      cartoonPhotometryGalaxies, testCatalog, cartoonStars, \
                                      cartoonGalaxies, testStars, testGalaxies

class variabilityUnitTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create test databases
        if os.path.exists('PhotometryTestDatabase.db'):
            print "deleting database"
            os.unlink('PhotometryTestDatabase.db')
        makeStarTestDB(size=100000, seedVal=1, filename='PhotometryTestDatabase.db')
        makeGalTestDB(size=100000, seedVal=1, filename='PhotometryTestDatabase.db')

    @classmethod
    def tearDownClass(cls):
        if os.path.exists('PhotometryTestDatabase.db'):
            os.unlink('PhotometryTestDatabase.db')

    def setUp(self):
        self.obs_metadata = ObservationMetaData(mjd=52000.7, bandpassName='i',
                            boundType = 'circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                            boundLength=1.0,m5=dict(u=23.9, g=25.0, r=24.7, i=24.0, z=23.3, y=22.1))

        self.galaxy = myTestGals(address='sqlite:///PhotometryTestDatabase.db')
        self.star = myTestStars(address='sqlite:///PhotometryTestDatabase.db')

    def tearDown(self):
        del self.galaxy
        del self.star
        del self.obs_metadata

    def testGalaxyVariability(self):

        galcat = testGalaxies(self.galaxy, obs_metadata=self.obs_metadata)
        results = self.galaxy.query_columns(['varParamStr'], obs_metadata=self.obs_metadata,
                                             constraint='VarParamStr is not NULL')
        result = getOneChunk(results)
        for row in result:
            mags=galcat.applyVariability(row['varParamStr'])

    def testStarVariability(self):
        starcat = testStars(self.star, obs_metadata=self.obs_metadata)
        results = self.star.query_columns(['varParamStr'], obs_metadata=self.obs_metadata,
                                         constraint='VarParamStr is not NULL')
        result = getOneChunk(results)
        for row in result:
            mags=starcat.applyVariability(row['varParamStr'])

class photometryUnitTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create test databases
        if os.path.exists('PhotometryTestDatabase.db'):
            print "deleting database"
            os.unlink('PhotometryTestDatabase.db')
        makeStarTestDB(size=100000, seedVal=1, filename='PhotometryTestDatabase.db')
        makeGalTestDB(size=100000, seedVal=1, filename='PhotometryTestDatabase.db')

    @classmethod
    def tearDownClass(cls):
        if os.path.exists('PhotometryTestDatabase.db'):
            os.unlink('PhotometryTestDatabase.db')

    def setUp(self):
        self.obs_metadata = ObservationMetaData(mjd=52000.7, bandpassName='i',
                            boundType='circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                            boundLength=1.0, m5 = 25.0)

        self.galaxy = myTestGals(address='sqlite:///PhotometryTestDatabase.db')
        self.star = myTestStars(address='sqlite:///PhotometryTestDatabase.db')

    def tearDown(self):
        del self.galaxy
        del self.star
        del self.obs_metadata

    def testStars(self):
        test_cat=testStars(self.star, obs_metadata=self.obs_metadata)
        test_cat.write_catalog("testStarsOutput.txt")
        results = self.star.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)
        os.unlink("testStarsOutput.txt")


    def testGalaxies(self):
        test_cat=testGalaxies(self.galaxy, obs_metadata=self.obs_metadata)
        test_cat.write_catalog("testGalaxiesOutput.txt")
        results = self.galaxy.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)
        os.unlink("testGalaxiesOutput.txt")

    def testAlternateBandpassesStars(self):
        """
        This will test our ability to do photometry using non-LSST bandpasses.

        It will first calculate the magnitudes using the getters in cartoonPhotometryStars.

        It will then load the alternate bandpass files 'by hand' and re-calculate the magnitudes
        and make sure that the magnitude values agree.  This is guarding against the possibility
        that some default value did not change and the code actually ended up loading the
        LSST bandpasses.
        """

        obs_metadata_pointed=ObservationMetaData(mjd=2013.23,
                                                 boundType='circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                                                 boundLength=1.0)
        obs_metadata_pointed.metadata = {}
        obs_metadata_pointed.metadata['Opsim_filter'] = 'i'
        test_cat=cartoonStars(self.star,obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog("testStarsCartoon.txt")

        cartoonDir = os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'
        testBandPasses = {}
        keys = ['u','g','r','i','z']

        bplist = []

        for kk in keys:
            testBandPasses[kk] = Bandpass()
            testBandPasses[kk].readThroughput(os.path.join(cartoonDir,"test_bandpass_%s.dat" % kk))
            bplist.append(testBandPasses[kk])

        sedObj = Sed()
        phiArray, waveLenStep = sedObj.setupPhiArray(bplist)

        i = 0

        #since all of the SEDs in the cartoon database are the same, just test on the first
        #if we ever include more SEDs, this can be somethin like
        #for ss in test_cata.sedMasterList:
        #
        ss=test_cat.sedMasterList[0]
        ss.resampleSED(wavelen_match = bplist[0].wavelen)
        ss.flambdaTofnu()
        mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
        for j in range(len(mags)):
            self.assertAlmostEqual(mags[j],test_cat.magnitudeMasterList[i][j],10)
        i += 1

        os.unlink("testStarsCartoon.txt")

    def testAlternateBandpassesGalaxies(self):
        """
        the same as testAlternateBandpassesStars, but for galaxies
        """

        obs_metadata_pointed=ObservationMetaData(mjd=50000.0,
                               boundType='circle',unrefractedRA=0.0,unrefractedDec=0.0,
                               boundLength=10.0)
        obs_metadata_pointed.metadata = {}
        obs_metadata_pointed.metadata['Opsim_filter'] = 'i'
        test_cat=cartoonGalaxies(self.galaxy,obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog("testGalaxiesCartoon.txt")

        cartoonDir = os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'
        testBandPasses = {}
        keys = ['u','g','r','i','z']

        bplist = []

        for kk in keys:
            testBandPasses[kk] = Bandpass()
            testBandPasses[kk].readThroughput(os.path.join(cartoonDir,"test_bandpass_%s.dat" % kk))
            bplist.append(testBandPasses[kk])

        sedObj = Sed()
        phiArray, waveLenStep = sedObj.setupPhiArray(bplist)

        components = ['Bulge', 'Disk', 'Agn']

        for cc in components:
            i = 0

            for ss in test_cat.sedMasterDict[cc]:
                if ss.wavelen != None:
                    ss.resampleSED(wavelen_match = bplist[0].wavelen)
                    ss.flambdaTofnu()
                    mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
                    for j in range(len(mags)):
                        self.assertAlmostEqual(mags[j],test_cat.magnitudeMasterDict[cc][i][j],10)
                i += 1

        os.unlink("testGalaxiesCartoon.txt")

    def testEBV(self):

        ebvObject = EBVbase()
        ra = []
        dec = []
        gLat = []
        gLon = []
        for i in range(10):
            ra.append(i*2.0*numpy.pi/10.0)
            dec.append(i*numpy.pi/10.0)

            gLat.append(-0.5*numpy.pi+i*numpy.pi/10.0)
            gLon.append(i*2.0*numpy.pi/10.0)

            equatorialCoordinates=numpy.array([ra,dec])
            galacticCoordinates=numpy.array([gLon,gLat])

        ebvOutput = ebvObject.calculateEbv(equatorialCoordinates=equatorialCoordinates)
        self.assertEqual(len(ebvOutput),len(ra))

        ebvOutput = ebvObject.calculateEbv(galacticCoordinates=galacticCoordinates)
        self.assertEqual(len(ebvOutput),len(gLon))

        self.assertRaises(RuntimeError, ebvObject.calculateEbv, equatorialCoordinates=equatorialCoordinates,
        galacticCoordinates=galacticCoordinates)
        self.assertRaises(RuntimeError, ebvObject.calculateEbv, equatorialCoordinates=None, galacticCoordinates=None)
        self.assertRaises(RuntimeError, ebvObject.calculateEbv)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(variabilityUnitTest)
    suites += unittest.makeSuite(photometryUnitTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
