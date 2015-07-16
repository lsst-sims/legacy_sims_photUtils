import numpy

import os
import unittest
import eups
import lsst.utils.tests as utilsTests
from lsst.sims.utils import ObservationMetaData
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils import PhotometryStars, PhotometryGalaxies, PhotometryBase, PhotometryHardware
from lsst.sims.photUtils import LSSTdefaults, PhotometricParameters, calcSNR_m5, calcGamma, \
                                calcM5, calcSNR_sed, calcSkyCountsPerPixelForM5, magErrorFromSNR

from lsst.sims.photUtils.utils import setM5


class photometryUnitTest(unittest.TestCase):

    def setUp(self):
        self.obs_metadata = ObservationMetaData(mjd=52000.7, bandpassName='i',
                            boundType='circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                            boundLength=1.0, m5 = 25.0)

    def tearDown(self):
        del self.obs_metadata


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

    def testStandAlongGalaxyPhotometry(self):
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
        starFileName = os.path.join(eups.productDir('sims_sed_library'),'starSED')
        starFileName = os.path.join(starFileName, 'kurucz','km20_5750.fits_g40_5790.gz')
        starName = os.path.join(eups.productDir('sims_sed_library'),starFileName)
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






def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(photometryUnitTest)
    suites += unittest.makeSuite(uncertaintyUnitTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
