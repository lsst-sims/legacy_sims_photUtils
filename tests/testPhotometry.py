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
from lsst.sims.photUtils import PhotometryStars, PhotometryGalaxies, PhotometryBase, PhotometryHardware
from lsst.sims.photUtils import LSSTdefaults, PhotometricParameters, calcSNR_gamma, calcGamma, \
                                calcM5, calcSNR_psf, expectedSkyCountsForM5
from lsst.sims.photUtils.utils import testDefaults, cartoonPhotometryStars, \
                                      cartoonPhotometryGalaxies, testCatalog, cartoonStars, \
                                      cartoonGalaxies, testStars, testGalaxies, \
                                      cartoonStarsOnlyI, cartoonStarsIZ, \
                                      cartoonGalaxiesIG, galaxiesWithHoles


def setM5(m5target, skysed, totalBandpass, hardware,
          photParams,
          seeing = None):
    """
    Take an SED representing the sky and normalize it so that
    m5 (the magnitude at which an object is detected in this
    bandpass at 5-sigma) is set to some specified value.

    The 5-sigma limiting magnitude (m5) for an observation is
    determined by a combination of the telescope and camera parameters
    (such as diameter of the mirrors and the readnoise) together with the
    sky background. This method (setM5) scales a provided sky background
    Sed so that an observation would have a target m5 value, for the
    provided hardware parameters. Using the resulting Sed in the
    'calcM5' method will return this target value for m5.

    @param [in] the desired value of m5

    @param [in] skysed is an instantiation of the Sed class representing
    sky emission

    @param [in] totalBandpass is an instantiation of the Bandpass class
    representing the total throughput of the telescope (instrumentation
    plus atmosphere)

    @param [in] hardware is an instantiation of the Bandpass class representing
    the throughput due solely to instrumentation.

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] seeing in arcseconds

    @param [out] returns an instantiation of the Sed class that is the skysed renormalized
    so that m5 has the desired value.

    Note that the returned SED will be renormalized such that calling the method
    self.calcADU(hardwareBandpass) on it will yield the number of counts per square
    arcsecond in a given bandpass.
    """

    #This is based on the LSST SNR document (v1.2, May 2010)
    #www.astro.washington.edu/users/ivezic/Astr511/LSST_SNRdoc.pdf

    if seeing is None:
        seeing = LSSTdefaults().seeing('r')

    skyCountsTarget = expectedSkyCountsForM5(m5target, totalBandpass, seeing=seeing,
                                             photParams=photParams)

    skySedOut = Sed(wavelen=numpy.copy(skysed.wavelen),
                    flambda=numpy.copy(skysed.flambda))

    skyCounts = skySedOut.calcADU(hardware, photParams=photParams) \
                    * photParams.platescale * photParams.platescale
    skySedOut.multiplyFluxNorm(skyCountsTarget/skyCounts)

    return skySedOut


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
        self.obs_metadata = ObservationMetaData(mjd=52000.7,
                            boundType = 'circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                            boundLength=1.0,
                            m5=[23.9, 25.0, 24.7, 24.0, 23.3, 22.1],
                            bandpassName=['u', 'g', 'r', 'i', 'z', 'y'])

        self.galaxy = myTestGals(database='PhotometryTestDatabase.db')
        self.star = myTestStars(database='PhotometryTestDatabase.db')

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

        self.galaxy = myTestGals(database='PhotometryTestDatabase.db')
        self.star = myTestStars(database='PhotometryTestDatabase.db')

    def tearDown(self):
        del self.galaxy
        del self.star
        del self.obs_metadata

    def testStarCatalog(self):
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


    def testGalaxyCatalog(self):
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

    def testSumMagnitudes(self):
        """
        Test that the method sum_magnitudes in PhotometryGalaxies handles
        NaNs correctly.  Test it both in the vectorized and non-vectorized form.
        """
        mm_0 = 22.0

        bulge = 15.0*numpy.ones(8)

        disk = 15.2*numpy.ones(8)

        agn = 15.4*numpy.ones(8)

        bulge[0] = numpy.NaN
        disk[1] = numpy.NaN
        agn[2] = numpy.NaN

        bulge[3] = numpy.NaN
        disk[3] = numpy.NaN

        bulge[4] = numpy.NaN
        agn[4] = numpy.NaN

        disk[5] = numpy.NaN
        agn[5] = numpy.NaN

        bulge[7] = numpy.NaN
        disk[7] = numpy.NaN
        agn[7] = numpy.NaN

        bulge_flux = numpy.power(10.0, -0.4*(bulge-mm_0))
        disk_flux = numpy.power(10.0, -0.4*(disk-mm_0))
        agn_flux = numpy.power(10.0, -0.4*(agn-mm_0))

        answer = numpy.zeros(8)
        answer[0] = -2.5*numpy.log10(disk_flux[0]+agn_flux[0]) + mm_0
        answer[1] = -2.5*numpy.log10(bulge_flux[1]+agn_flux[1]) + mm_0
        answer[2] = -2.5*numpy.log10(bulge_flux[2]+disk_flux[2]) + mm_0
        answer[3] = -2.5*numpy.log10(agn_flux[3]) + mm_0
        answer[4] = -2.5*numpy.log10(disk_flux[4]) + mm_0
        answer[5] = -2.5*numpy.log10(bulge_flux[5]) + mm_0
        answer[6] = -2.5*numpy.log10(bulge_flux[6]+disk_flux[6]+agn_flux[6]) + mm_0
        answer[7] = numpy.NaN

        phot = PhotometryGalaxies()
        test = phot.sum_magnitudes(bulge=bulge, disk=disk, agn=agn)

        numpy.testing.assert_array_almost_equal(test, answer, decimal=10)

        for ix, (bb, dd, aa, truth) in enumerate(zip(bulge, disk, agn, answer)):
            test = phot.sum_magnitudes(bulge=bb, disk=dd, agn=aa)
            if ix<7:
                self.assertAlmostEqual(test, truth, 10)
                self.assertTrue(not numpy.isnan(test))
            else:
                self.assertTrue(numpy.isnan(test))
                self.assertTrue(numpy.isnan(truth))

    def testSumMagnitudesCatalog(self):
        """
        test that sum_magnitudes handles NaNs correctly in the context
        of a catalog by outputting a catalog of galaxies with NaNs in
        different component magnitudes, reading that catalog back in,
        and then calculating the summed magnitude by hand and comparing
        """

        catName = 'galaxiesWithHoles.txt'
        obs_metadata=ObservationMetaData(mjd=50000.0,
                               boundType='circle',unrefractedRA=0.0,unrefractedDec=0.0,
                               boundLength=10.0)
        test_cat=galaxiesWithHoles(self.galaxy,obs_metadata=obs_metadata)
        test_cat.write_catalog(catName)

        dtype = numpy.dtype([
                            ('raJ2000', numpy.float),
                            ('decJ2000', numpy.float),
                            ('u', numpy.float), ('g', numpy.float), ('r', numpy.float),
                            ('i', numpy.float), ('z', numpy.float), ('y', numpy.float),
                            ('ub', numpy.float), ('gb', numpy.float), ('rb', numpy.float),
                            ('ib', numpy.float), ('zb', numpy.float), ('yb', numpy.float),
                            ('ud', numpy.float), ('gd', numpy.float), ('rd', numpy.float),
                            ('id', numpy.float), ('zd', numpy.float), ('yd', numpy.float),
                            ('ua', numpy.float), ('ga', numpy.float), ('ra', numpy.float),
                            ('ia', numpy.float), ('za', numpy.float), ('ya', numpy.float)
                            ])



        data = numpy.genfromtxt(catName, dtype=dtype, delimiter=', ')
        self.assertTrue(len(data)>16)
        phot = PhotometryGalaxies()

        test = phot.sum_magnitudes(bulge=data['ub'], disk=data['ud'], agn=data['ua'])
        numpy.testing.assert_array_almost_equal(test, data['u'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['gb'], disk=data['gd'], agn=data['ga'])
        numpy.testing.assert_array_almost_equal(test, data['g'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['rb'], disk=data['rd'], agn=data['ra'])
        numpy.testing.assert_array_almost_equal(test, data['r'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['ib'], disk=data['id'], agn=data['ia'])
        numpy.testing.assert_array_almost_equal(test, data['i'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['zb'], disk=data['zd'], agn=data['za'])
        numpy.testing.assert_array_almost_equal(test, data['z'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['yb'], disk=data['yd'], agn=data['ya'])
        numpy.testing.assert_array_almost_equal(test, data['y'], decimal=10)

        # make sure that there were some NaNs for our catalog to deal with (but that they were not
        # all NaNs
        for line in [data['u'], data['g'], data['r'], data['i'], data['z'], data['y'],
                     data['ub'], data['gb'], data['rb'], data['ib'], data['zb'], data['yb'],
                     data['ud'], data['gd'], data['rd'], data['id'], data['zd'], data['yd'],
                     data['ua'], data['ga'], data['ra'], data['ia'], data['za'], data['ya']]:

            ctNans = len(numpy.where(numpy.isnan(line))[0])
            self.assertTrue(ctNans>0)
            self.assertTrue(ctNans<len(line))

        if os.path.exists(catName):
            os.unlink(catName)

    def testStandAloneStellarPhotometry(self):
        """
        Test that it is possible to run PhotometryStars.calculate_magnitudes
        outside of the context of an InstanceCatalog
        """
        objectID = ['1','2','3']
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
                          objectID=dummyId, sedNames=sedNames, magNorm=magNorm)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes,
                          objectID=objectID, sedNames=dummySed, magNorm=magNorm)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes,
                          objectID=objectID, sedNames=sedNames, magNorm=dummyMagNorm)

        magnitudes = phot.calculate_magnitudes(objectID=objectID, sedNames=sedNames,
                                               magNorm=magNorm)
        for n in objectID:
            self.assertTrue(len(magnitudes[n])==len(bandpassNames)) #to make sure we calculated all the magnitudes

    def testGalaxyPhotometryStandAlone(self):
        objectID = ['Alice', 'Bob', 'Charlie']

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

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          diskNames=diskSedsDummy, diskMagNorm=diskMagNorm, diskAv=diskAv,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          diskNames=diskSeds, diskMagNorm=diskMagNormDummy, diskAv=diskAvDummy,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          diskNames=diskSeds, diskMagNorm=diskMagNorm, diskAv=diskAvDummy,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          bulgeNames=bulgeSedsDummy, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNormDummy, bulgeAv=bulgeAv,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAvDummy,
                          redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          agnNames=agnSedsDummy, agnMagNorm=agnMagNorm)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          agnNames=agnSeds, agnMagNorm=agnMagNormDummy)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv,
                          redshift=redshiftDummy)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv,
                          redshift=redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulusDummy)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, redshift=redshift)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          bulgeNames=bulgeSeds, bulgeAv=bulgeAv, redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          diskNames=diskSeds, diskMagNorm=diskMagNorm, redshift=redshift)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          diskNames=diskSeds, diskAv=diskAv, redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          agnNames=agnSeds, redshift=redshift)

        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          bulgeNames=bulgeSeds, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          diskNames=diskSeds, diskMagNorm=diskMagNorm, diskAv=diskAv)
        self.assertRaises(RuntimeError, phot.calculate_magnitudes, objectID,
                          agnNames=agnSeds, agnMagNorm=agnMagNorm)

        bulgeNamePossibilities = [bulgeSeds, None]
        diskNamePossibilities = [diskSeds, None]
        agnNamePossibilities = [agnSeds, None]

        for bulgeNames in bulgeNamePossibilities:
            for diskNames in diskNamePossibilities:
                for agnNames in agnNamePossibilities:

                    magnitudes = phot.calculate_magnitudes(objectID, redshift=redshift,
                                                           bulgeNames=bulgeNames, bulgeMagNorm=bulgeMagNorm, bulgeAv=bulgeAv,
                                                           diskNames=diskNames, diskMagNorm=diskMagNorm, diskAv=diskAv,
                                                           agnNames=agnNames, agnMagNorm=agnMagNorm)

                    for name in objectID:
                        for i in range(len(phot.bandpassDict)):
                            flux=0.0
                            if bulgeNames is None:
                                self.assertTrue(numpy.isnan(magnitudes[name]['bulge'][i]))
                            else:
                                self.assertTrue(magnitudes[name]['bulge'][i] is not None)
                                self.assertFalse(numpy.isnan(magnitudes[name]['bulge'][i]))
                                flux += numpy.power(10.0, -0.4*(magnitudes[name]['bulge'][i]-22.0))

                            if diskNames is None:
                                self.assertTrue(numpy.isnan(magnitudes[name]['disk'][i]))
                            else:
                                self.assertTrue(magnitudes[name]['disk'][i] is not None)
                                self.assertFalse(numpy.isnan(magnitudes[name]['disk'][i]))
                                flux += numpy.power(10.0, -0.4*(magnitudes[name]['disk'][i]-22.0))

                            if agnNames is None:
                                self.assertTrue(numpy.isnan(magnitudes[name]['agn'][i]))
                            else:
                                self.assertTrue(magnitudes[name]['agn'][i] is not None)
                                self.assertFalse(numpy.isnan(magnitudes[name]['agn'][i]))
                                flux += numpy.power(10.0, -0.4*(magnitudes[name]['agn'][i]-22.0))


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

    def testStellarPhotometryIndices(self):
        """
        A test to make sure that stellar photometry still calculates the right values
        even when it is not calculating all of the magnitudes in the getter
        """

        baselineDtype = numpy.dtype([('id',int),
                                     ('raObserved', float), ('decObserved', float),
                                     ('magNorm', float),
                                     ('cartoon_u', float), ('cartoon_g',float),
                                     ('cartoon_r', float), ('cartoon_i', float),
                                     ('cartoon_z', float)])

        baselineCatName = 'stellarBaselineCatalog.txt'

        testDtype = numpy.dtype([('id',int),
                                 ('raObserved',float), ('decObserved',float),
                                 ('cartoon_i',float)])

        testCatName = 'stellarTestCatalog.txt'


        obs_metadata_pointed=ObservationMetaData(mjd=2013.23,
                                                 boundType='circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                                                 boundLength=1.0)

        baseline_cat=cartoonStars(self.star,obs_metadata=obs_metadata_pointed)
        baseline_cat.write_catalog(baselineCatName)
        baselineData = numpy.genfromtxt(baselineCatName, dtype=baselineDtype, delimiter=',')

        test_cat=cartoonStarsOnlyI(self.star, obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(testCatName)
        testData = numpy.genfromtxt(testCatName, dtype=testDtype, delimiter=',')
        ct = 0
        for b, t in zip(baselineData, testData):
            self.assertAlmostEqual(b['cartoon_i'], t['cartoon_i'], 10)
            ct+=1
        self.assertTrue(ct>0)

        testDtype = numpy.dtype([('id',int),
                                 ('raObserved',float), ('decObserved',float),
                                 ('cartoon_i',float), ('cartoon_z',float)])


        test_cat=cartoonStarsIZ(self.star, obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(testCatName)
        testData = numpy.genfromtxt(testCatName, dtype=testDtype, delimiter=',')
        ct = 0
        for b, t in zip(baselineData, testData):
            self.assertAlmostEqual(b['cartoon_i'], t['cartoon_i'], 10)
            self.assertAlmostEqual(b['cartoon_z'], t['cartoon_z'], 10)
            ct+=1
        self.assertTrue(ct>0)

        if os.path.exists(testCatName):
            os.unlink(testCatName)
        if os.path.exists(baselineCatName):
            os.unlink(baselineCatName)

    def testGalaxyPhotometricIndices(self):
        baselineCatName = 'galaxyBaselineCatalog.txt'
        baselineDtype = numpy.dtype([('galid', int),
                                     ('raObserved', float),
                                     ('decObserved', float),
                                     ('ctotal_u', float),
                                     ('ctotal_g', float),
                                     ('ctotal_r', float),
                                     ('ctotal_i', float),
                                     ('ctotal_z', float)])

        obs_metadata_pointed=ObservationMetaData(mjd=50000.0,
                               boundType='circle',unrefractedRA=0.0,unrefractedDec=0.0,
                               boundLength=10.0)

        baseline_cat=cartoonGalaxies(self.galaxy,obs_metadata=obs_metadata_pointed)
        baseline_cat.write_catalog(baselineCatName)
        baselineData = numpy.genfromtxt(baselineCatName, dtype=baselineDtype, delimiter=',')

        testCatName = 'galaxyTestCatalog.txt'
        testDtype = numpy.dtype([('galid', int),
                                 ('raObserved', float),
                                 ('decObserved', float),
                                 ('ctotal_i', float),
                                 ('ctotal_g', float)])
        test_cat = cartoonGalaxiesIG(self.galaxy, obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(testCatName)
        testData = numpy.genfromtxt(testCatName, dtype=testDtype, delimiter=',')
        ct = 0
        for b,t in zip(baselineData, testData):
            self.assertAlmostEqual(b['ctotal_i'], t['ctotal_i'], 10)
            self.assertAlmostEqual(b['ctotal_g'], t['ctotal_g'], 10)
            ct += 1
        self.assertTrue(ct>0)

        if os.path.exists(baselineCatName):
            os.unlink(baselineCatName)

        if os.path.exists(testCatName):
            os.unlink(testCatName)

    def testPhotometricIndicesRaw(self):
        """
        Use manMagCalc_list with specified indices on an Sed.  Make sure
        that the appropriate magnitudes are or are not Nan
        """
        starName = os.path.join(eups.productDir('sims_sed_library'),defaultSpecMap['km20_5750.fits_g40_5790'])
        starPhot = PhotometryStars()
        starPhot.loadTotalBandpassesFromFiles()
        testSed = Sed()
        testSed.readSED_flambda(starName)
        indices = [1,3]
        mags = starPhot.manyMagCalc_list(testSed, indices=indices)
        self.assertTrue(numpy.isnan(mags[0]))
        self.assertFalse(numpy.isnan(mags[1]))
        self.assertTrue(numpy.isnan(mags[2]))
        self.assertFalse(numpy.isnan(mags[3]))
        self.assertTrue(numpy.isnan(mags[4]))
        self.assertTrue(numpy.isnan(mags[5]))
        self.assertTrue(len(mags)==6)


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
        Test the calculateMagnitudeUncertainty raises exceptions when it needs to
        """
        phot = PhotometryBase()
        phot.loadBandpassesFromFiles()
        magnitudes = numpy.array([22.0, 23.0, 24.0, 25.0, 26.0, 27.0])
        shortMagnitudes = numpy.array([22.0])
        self.assertRaises(RuntimeError, phot.calculateMagnitudeUncertainty, magnitudes)
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, bandpassName='g', m5=23.0)
        self.assertRaises(RuntimeError, phot.calculateMagnitudeUncertainty, shortMagnitudes, obs_metadata=obs_metadata)

        photParams = PhotometricParameters()
        shortGamma = numpy.array([1.0, 1.0])
        fluxes = numpy.power(10.0, -0.4*magnitudes)
        shortFluxes = numpy.power(10.0, -0.4*shortMagnitudes)
        self.assertRaises(RuntimeError, calcSNR_gamma, fluxes, phot.bandpassDict.values(), shortMagnitudes, photParams)
        self.assertRaises(RuntimeError, calcSNR_gamma, shortFluxes, phot.bandpassDict.values(), magnitudes, photParams)
        self.assertRaises(RuntimeError, calcSNR_gamma, fluxes, phot.bandpassDict.values(), magnitudes, photParams, gamma=shortGamma)
        snr, gg = calcSNR_gamma(fluxes, phot.bandpassDict.values(), magnitudes, photParams)


    def testSignalToNoise(self):
        """
        Test that calcSNR_gamma and calcSNR_psf give similar results
        """
        defaults = LSSTdefaults()
        photParams = PhotometricParameters()
        hardware = PhotometryHardware()
        hardware.loadBandpassesFromFiles()

        m5 = []
        for filt in hardware.bandpassDict:
            m5.append(calcM5(hardware.skySED, hardware.bandpassDict[filt],
                      hardware.hardwareBandpassDict[filt],
                      photParams, seeing=defaults.seeing(filt)))


        sedDir = eups.productDir('sims_sed_library')
        sedDir = os.path.join(sedDir, 'starSED', 'kurucz')
        fileNameList = os.listdir(sedDir)

        numpy.random.seed(42)
        offset = numpy.random.random_sample(len(fileNameList))*2.0

        for ix, name in enumerate(fileNameList):
            if ix>100:
                break
            spectrum = Sed()
            spectrum.readSED_flambda(os.path.join(sedDir, name))
            ff = spectrum.calcFluxNorm(m5[2]-offset[ix], hardware.bandpassDict.values()[2])
            spectrum.multiplyFluxNorm(ff)
            fluxList = []
            controlList = []
            magList = []
            for filt in hardware.bandpassDict:
                controlList.append(calcSNR_psf(spectrum, hardware.bandpassDict[filt],
                                               hardware.skySED,
                                               hardware.hardwareBandpassDict[filt],
                                               photParams, defaults.seeing(filt)))

                fluxList.append(spectrum.calcFlux(hardware.bandpassDict[filt]))

            testList, gammaList = calcSNR_gamma(numpy.array(fluxList),
                                        numpy.array(hardware.bandpassDict.values()),
                                        numpy.array(m5),
                                        photParams)

            for tt, cc in zip(controlList, testList):
                msg = '%e != %e ' % (tt, cc)
                self.assertTrue(numpy.abs(tt/cc - 1.0) < 0.001, msg=msg)



    def testRawUncertainty(self):
        """
        Test that values calculated by calculatePhotometricUncertainty agree
        with values calculated by calcSNR_psf
        """

        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]
        phot = PhotometryBase()
        phot.loadTotalBandpassesFromFiles()
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5, bandpassName=self.bandpasses)
        magnitudes = phot.manyMagCalc_list(self.starSED)

        skySeds = []

        for i in range(len(self.bandpasses)):
            skyDummy = Sed()
            skyDummy.readSED_flambda(os.path.join(eups.productDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSkyDummy = setM5(obs_metadata.m5[self.bandpasses[i]], skyDummy,
                                       self.totalBandpasses[i], self.hardwareBandpasses[i],
                                       seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                                       photParams=PhotometricParameters())

            skySeds.append(normalizedSkyDummy)

        sigma = phot.calculateMagnitudeUncertainty(magnitudes, obs_metadata=obs_metadata)
        for i in range(len(self.bandpasses)):
            snr = calcSNR_psf(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            ss = 2.5*numpy.log10(1.0+1.0/snr)
            msg = '%e is not %e; failed' % (ss, sigma[i])
            self.assertAlmostEqual(ss, sigma[i], 10, msg=msg)

    def testSystematicUncertainty(self):
        """
        Test that systematic uncertainty is added correctly.
        """
        sigmaSysSq = 0.002
        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]

        phot = PhotometryBase()
        phot.loadTotalBandpassesFromFiles()
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5, bandpassName=self.bandpasses)
        magnitudes = phot.manyMagCalc_list(self.starSED)

        skySeds = []

        for i in range(len(self.bandpasses)):
            skyDummy = Sed()
            skyDummy.readSED_flambda(os.path.join(eups.productDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSkyDummy = setM5(obs_metadata.m5[self.bandpasses[i]], skyDummy,
                                                       self.totalBandpasses[i], self.hardwareBandpasses[i],
                                                       seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                                                       photParams=PhotometricParameters())

            skySeds.append(normalizedSkyDummy)

        sigma = phot.calculateMagnitudeUncertainty(magnitudes, obs_metadata=obs_metadata, sigmaSysSq=sigmaSysSq)
        for i in range(len(self.bandpasses)):
            snr = calcSNR_psf(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            testSNR, gamma = calcSNR_gamma(numpy.array([Sed().fluxFromMag(magnitudes[i])]), [self.totalBandpasses[i]],
                                           numpy.array([m5[i]]), photParams=PhotometricParameters())

            self.assertAlmostEqual(snr, testSNR[0], 10, msg = 'failed on calcSNR_gamma test %e != %e ' \
                                                               % (snr, testSNR[0]))

            control = 1.0/(snr*snr) + sigmaSysSq
            test = numpy.power(numpy.power(10.0, sigma[i]/2.5) -1.0, 2)

            msg = '%e is not %e; failed' % (test, control)

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
