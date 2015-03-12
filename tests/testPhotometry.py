import numpy

import os
import unittest
import eups
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.db import ObservationMetaData
from lsst.sims.catalogs.generation.utils import myTestGals, myTestStars, \
                                                makeStarTestDB, makeGalTestDB, getOneChunk

from lsst.sims.catalogs.measures.instance import defaultSpecMap
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils import PhotometryStars, PhotometryGalaxies, PhotometryBase
from lsst.sims.photUtils import PhotometricDefaults, setM5
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

        makeStarTestDB(filename='PhotometryTestDatabase.db', size=100000, seedVal=1)
        makeGalTestDB(filename='PhotometryTestDatabase.db', size=100000, seedVal=1)

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
        ct = 0
        for row in result:
            mags=galcat.applyVariability(row['varParamStr'])
            ct += 1
        self.assertTrue(ct>0) #to make sure that the test was actually performed

    def testStarVariability(self):
        starcat = testStars(self.star, obs_metadata=self.obs_metadata)
        results = self.star.query_columns(['varParamStr'], obs_metadata=self.obs_metadata,
                                         constraint='VarParamStr is not NULL')
        result = getOneChunk(results)
        ct = 0
        for row in result:
            ct += 1
            mags=starcat.applyVariability(row['varParamStr'])
        self.assertTrue(ct>0) #to make sure that the test was actually performed

class photometryUnitTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create test databases
        if os.path.exists('PhotometryTestDatabase.db'):
            print "deleting database"
            os.unlink('PhotometryTestDatabase.db')

        makeStarTestDB(filename='PhotometryTestDatabase.db', size=100000, seedVal=1)
        makeGalTestDB(filename='PhotometryTestDatabase.db', size=100000, seedVal=1)

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
        cat = open("testStarsOutput.txt")
        lines = cat.readlines()
        self.assertTrue(len(lines)>1) #to make sure we did not write an empty catalog
        cat.close()
        results = self.star.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)
        self.assertTrue(len(result)>0) #to make sure some results are returned
        os.unlink("testStarsOutput.txt")


    def testGalaxies(self):
        test_cat=testGalaxies(self.galaxy, obs_metadata=self.obs_metadata)
        test_cat.write_catalog("testGalaxiesOutput.txt")
        cat = open("testGalaxiesOutput.txt")
        lines = cat.readlines()
        self.assertTrue(len(lines)>1) #to make sure we did not write an empty catalog
        cat.close()
        results = self.galaxy.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)
        self.assertTrue(len(result)>0) #to make sure some results are returned
        os.unlink("testGalaxiesOutput.txt")

    def testStandAloneStellarPhotometry(self):
        """
        Test that it is possible to run PhotometryStars.calculate_magnitudes
        outside of the context of an InstanceCatalog
        """
        idNames = ['1','2','3']
        sedNames = ['km20_5750.fits_g40_5790','m2.0Full.dat',
                     'bergeron_6500_85.dat_6700']
        magNorm = [28.5, 23.0, 21.0]

        dummyId = ['1', '2']
        dummySed = ['km20_5750.fits_g40_5790','m2.0Full.dat']
        dummyMagNorm = [28.5, 23.0]
        bandpassNames = ['u','g','r','i','z','y']

        phot = PhotometryStars()
        phot.loadTotalBandpassesFromFiles(bandpassNames)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes,
                          idNames=dummyId, sedNames=sedNames, magNorm=magNorm)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes,
                          idNames=idNames, sedNames=dummySed, magNorm=magNorm)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes,
                          idNames=idNames, sedNames=sedNames, magNorm=dummyMagNorm)

        magnitudes = phot.calculate_magnitudes(idNames=idNames, sedNames=sedNames,
                                               magNorm=magNorm)
        for n in idNames:
            self.assertTrue(len(magnitudes[n])==len(bandpassNames)) #to make sure we calculated all the magnitudes

    def testGalaxyPhotometryStandAlone(self):
        idNames = ['Alice', 'Bob', 'Charlie']

        diskSeds = ['Const.80E07.02Z.spec','Inst.80E07.002Z.spec','Burst.19E07.0005Z.spec']
        diskMagNorm = [24.2, 28.1, 29.0]
        diskAv = [3.1, 3.2, 2.9]

        bulgeSeds = ['Inst.80E07.002Z.spec','Const.80E07.02Z.spec','Burst.19E07.0005Z.spec']
        bulgeMagNorm = [25.0, 28.0, 27.1]
        bulgeAv = [2.8, 3.2, 3.3]

        agnSeds = ['agn.spec', 'agn.spec', 'agn.spec']
        agnMagNorm = [22.0, 23.0, 26.0]

        redshift = [0.2, 0.3, 1.1]
        cosmologicalDistanceModulus = [5.0, 3.0, 4.5]

        diskSedsDummy = ['Inst.80E07.002Z.spec','Burst.19E07.0005Z.spec']
        diskMagNormDummy = [28.1, 29.0]
        diskAvDummy = [3.2, 2.9]

        bulgeSedsDummy = ['Const.80E07.02Z.spec','Burst.19E07.0005Z.spec']
        bulgeMagNormDummy = [28.0, 27.1]
        bulgeAvDummy = [3.2, 3.3]

        agnSeds = ['agn.spec', 'agn.spec', 'agn.spec']
        agnMagNorm = [22.0, 23.0, 26.0]

        agnSedsDummy = ['agn.spec', 'agn.spec']
        agnMagNormDummy = [22.0, 26.0]

        redshiftDummy = [0.3, 1.1]
        cosmologicalDistanceModulusDummy = [3.0, 4.5]

        phot = PhotometryGalaxies()
        phot.loadTotalBandpassesFromFiles(bandpassNames=['u','g','r','i','z','y'])

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          diskNames=diskSedsDummy, diskMagNorm=diskMagNorm, diskAv=diskAv,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          diskNames=diskSeds, diskMagNorm=diskMagNormDummy, diskAv=diskAvDummy,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          diskNames=diskSeds, diskMagNorm=diskMagNorm, diskAv=diskAvDummy,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          bulgeNames=bulgeSedsDummy, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNormDummy, bulgeAv=bulgeAv,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAvDummy,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          agnNames=agnSedsDummy, agnMagNorm=agnMagNorm)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          agnNames=agnSeds, agnMagNorm=agnMagNormDummy)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv,
                          redshift=redshiftDummy)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv,
                          redshift=redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulusDummy)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, redshift=redshift)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          bulgeNames=bulgeSeds, bulgeAv=bulgeAv, redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          diskNames=diskSeds, diskMagNorm=diskMagNorm, redshift=redshift)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          diskNames=diskSeds, diskAv=diskAv, redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          agnNames=agnSeds, redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          diskNames=diskSeds, diskMagNorm=diskMagNorm, diskAv=diskAv)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes, idNames,
                          agnNames=agnSeds, agnMagNorm=agnMagNorm)

        bulgeNamePossibilities = [bulgeSeds, None]
        diskNamePossibilities = [diskSeds, None]
        agnNamePossibilities = [agnSeds, None]

        for bulgeNames in bulgeNamePossibilities:
            for diskNames in diskNamePossibilities:
                for agnNames in agnNamePossibilities:

                    magnitudes = phot.calculate_magnitudes(idNames, redshift=redshift,
                                                           bulgeNames=bulgeNames, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv,
                                                           diskNames=diskNames, diskMagNorm=diskMagNorm, diskAv=diskAv,
                                                           agnNames=agnNames, agnMagNorm=agnMagNorm)

                    for name in idNames:
                        for i in range(len(phot.bandpassDict)):
                            flux=0.0
                            if bulgeNames is None:
                                self.assertTrue(magnitudes[name]['bulge'][i] is None)
                            else:
                                self.assertTrue(magnitudes[name]['bulge'][i] is not None)
                                self.assertFalse(numpy.isnan(magnitudes[name]['bulge'][i]))
                                flux += numpy.power(10.0, -0.4*(magnitudes[name]['bulge'][i]-22.0))

                            if diskNames is None:
                                self.assertTrue(magnitudes[name]['disk'][i] is None)
                            else:
                                self.assertTrue(magnitudes[name]['disk'][i] is not None)
                                self.assertFalse(numpy.isnan(magnitudes[name]['disk'][i]))
                                flux += numpy.power(10.0, -0.4*(magnitudes[name]['disk'][i]-22.0))

                            if agnNames is None:
                                self.assertTrue(magnitudes[name]['agn'][i] is None)
                            else:
                                self.assertTrue(magnitudes[name]['agn'][i] is not None)
                                self.assertFalse(numpy.isnan(magnitudes[name]['agn'][i]))
                                flux += numpy.power(10.0, -0.4*(magnitudes[name]['agn'][i]-22.0))

                            if agnNames is None and diskNames is None and bulgeNames is None:
                                self.assertTrue(magnitudes[name]['total'][i] is None)
                            else:
                                self.assertTrue(magnitudes[name]['total'][i] is not None)
                                self.assertFalse(numpy.isnan(magnitudes[name]['total'][i]))
                                self.assertAlmostEqual(magnitudes[name]['total'][i], -2.5*numpy.log10(flux)+22.0, 10)

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
        #if we ever include more SEDs, this can be something like
        #for ss in test_cata.sedMasterList:
        #
        ss=test_cat.sedMasterList[0]
        ss.resampleSED(wavelen_match = bplist[0].wavelen)
        ss.flambdaTofnu()
        mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
        self.assertTrue(len(mags)==len(test_cat.bandpassDict))
        self.assertTrue(len(mags)>0)
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

        ct = 0
        for cc in components:
            i = 0

            for ss in test_cat.sedMasterDict[cc]:
                if ss.wavelen != None:
                    ss.resampleSED(wavelen_match = bplist[0].wavelen)
                    ss.flambdaTofnu()
                    mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
                    for j in range(len(mags)):
                        ct += 1
                        self.assertAlmostEqual(mags[j],test_cat.magnitudeMasterDict[cc][i][j],10)
                i += 1

        self.assertTrue(ct>0)
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

class uncertaintyUnitTest(unittest.TestCase):
    """
    Test the calculation of photometric uncertainties
    """

    def setUp(self):
        starName = os.path.join(eups.productDir('sims_sed_library'),defaultSpecMap['km20_5750.fits_g40_5790'])
        self.starSED = Sed()
        self.starSED.readSED_flambda(starName)
        imsimband = Bandpass()
        imsimband.imsimBandpass()
        fNorm = self.starSED.calcFluxNorm(22.0, imsimband)
        self.starSED.multiplyFluxNorm(fNorm)

        self.totalBandpasses = []
        self.hardwareBandpasses = []

        componentList = ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                         'lens1.dat', 'lens2.dat', 'lens3.dat']
        hardwareComponents = []
        for c in componentList:
            hardwareComponents.append(os.path.join(eups.productDir('throughputs'),'baseline',c))

        self.bandpasses = ['u', 'g', 'r', 'i', 'z', 'y']
        for b in self.bandpasses:
            filterName = os.path.join(eups.productDir('throughputs'),'baseline','filter_%s.dat' % b)
            components = hardwareComponents + [filterName]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.hardwareBandpasses.append(bandpassDummy)
            components = components + [os.path.join(eups.productDir('throughputs'),'baseline','atmos.dat')]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.totalBandpasses.append(bandpassDummy)



    def tearDown(self):
        del self.starSED
        del self.bandpasses
        del self.hardwareBandpasses
        del self.totalBandpasses

    def testUncertaintyExceptions(self):
        """
        Test the calculatePhotometricUncertainty raises exceptions when it needs to
        """
        phot = PhotometryBase()
        phot.loadBandpassesFromFiles()
        magnitudes = [22.0, 23.0, 24.0, 25.0, 26.0, 27.0]
        shortMagnitudes = [22.0]
        self.assertRaises(RuntimeError, phot.calculatePhotometricUncertainty, magnitudes)
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, bandpassName='g', m5=23.0)
        self.assertRaises(RuntimeError, phot.calculatePhotometricUncertainty, shortMagnitudes, obs_metadata=obs_metadata)

        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, bandpassName='g')
        self.assertRaises(ValueError, phot.calculatePhotometricUncertainty, magnitudes, obs_metadata=obs_metadata)

        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, bandpassName='g', m5={'u':22.0, 'g':24.0})
        self.assertRaises(ValueError, phot.calculatePhotometricUncertainty, magnitudes, obs_metadata=obs_metadata)

    def testRawUncertainty(self):
        """
        Test that values calculated by calculatePhotometricUncertainty agree
        with values calculated by Sed.calcSNR_psf
        """
        for ii in range(2):
            if ii == 0:
                msgroot = "m5 is a float"
                m5 = 25.0
            else:
                msgroot = "m5 is a dict"
                m5 = {'u':23.0, 'g':21.0, 'r':24.6, 'i':23.6, 'z':22.5, 'y':20.0}

            phot = PhotometryBase()
            phot.loadTotalBandpassesFromFiles()
            obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5)
            magnitudes = phot.manyMagCalc_list(self.starSED)

            skySeds = []

            for i in range(len(self.bandpasses)):
                skyDummy = Sed()
                skyDummy.readSED_flambda(os.path.join(eups.productDir('throughputs'), 'baseline', 'darksky.dat'))
                normalizedSkyDummy = setM5(obs_metadata.m5(self.bandpasses[i]), skyDummy,
                                                           self.totalBandpasses[i], self.hardwareBandpasses[i],
                                                           seeing=PhotometricDefaults.seeing[self.bandpasses[i]])
                skySeds.append(normalizedSkyDummy)

            sigma = phot.calculatePhotometricUncertainty(magnitudes, obs_metadata=obs_metadata)
            for i in range(len(self.bandpasses)):
                snr = self.starSED.calcSNR_psf(self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                                               seeing=PhotometricDefaults.seeing[self.bandpasses[i]])
                ss = 2.5*numpy.log10(1.0+1.0/snr)
                msg = '%e is not %e; failed when ' % (ss, sigma[i]) + msgroot
                self.assertAlmostEqual(ss, sigma[i], 10, msg=msg)

    def testSystematicUncertainty(self):
        """
        Test that systematic uncertainty is added correctly.
        """
        sig2sys = 0.002
        for ii in range(2):
            if ii == 0:
                msgroot = "m5 is a float"
                m5 = 25.0
            else:
                msgroot = "m5 is a dict"
                m5 = {'u':23.0, 'g':21.0, 'r':24.6, 'i':23.6, 'z':22.5, 'y':20.0}

            phot = PhotometryBase()
            phot.loadTotalBandpassesFromFiles()
            obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5)
            magnitudes = phot.manyMagCalc_list(self.starSED)

            skySeds = []

            for i in range(len(self.bandpasses)):
                skyDummy = Sed()
                skyDummy.readSED_flambda(os.path.join(eups.productDir('throughputs'), 'baseline', 'darksky.dat'))
                normalizedSkyDummy = setM5(obs_metadata.m5(self.bandpasses[i]), skyDummy,
                                                           self.totalBandpasses[i], self.hardwareBandpasses[i],
                                                           seeing=PhotometricDefaults.seeing[self.bandpasses[i]])
                skySeds.append(normalizedSkyDummy)

            sigma = phot.calculatePhotometricUncertainty(magnitudes, obs_metadata=obs_metadata, sig2sys=sig2sys)
            for i in range(len(self.bandpasses)):
                snr = self.starSED.calcSNR_psf(self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                                               seeing=PhotometricDefaults.seeing[self.bandpasses[i]])

                control = 1.0/(snr*snr) + sig2sys
                test = numpy.power(numpy.power(10.0, sigma[i]/2.5) -1.0, 2)

                msg = '%e is not %e; failed when ' % (test, control) + msgroot

                self.assertAlmostEqual(test, control, 10, msg=msg)



def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(variabilityUnitTest)
    suites += unittest.makeSuite(photometryUnitTest)
    suites += unittest.makeSuite(uncertaintyUnitTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
