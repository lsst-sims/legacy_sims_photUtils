import numpy as np

import os
import unittest
import lsst.utils
import lsst.utils.tests
from lsst.sims.utils import ObservationMetaData
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils import LSSTdefaults, PhotometricParameters, calcSNR_m5, \
                                calcM5, calcSNR_sed, magErrorFromSNR, BandpassDict
from lsst.sims.photUtils.utils import setM5


def setup_module(module):
    lsst.utils.tests.init()


class photometryUnitTest(unittest.TestCase):

    def setUp(self):
        self.obs_metadata = ObservationMetaData(mjd=52000.7, bandpassName='i',
                            boundType='circle',pointingRA=200.0,pointingDec=-30.0,
                            boundLength=1.0, m5 = 25.0)

    def tearDown(self):
        del self.obs_metadata


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
                                                 boundType='circle',pointingRA=200.0,pointingDec=-30.0,
                                                 boundLength=1.0)

        bandpassDir=os.path.join(lsst.utils.getPackageDir('sims_photUtils'),'tests','cartoonSedTestData')

        cartoon_dict = BandpassDict.loadTotalBandpassesFromFiles(['u','g','r','i','z'],bandpassDir = bandpassDir,
                                                                 bandpassRoot = 'test_bandpass_')

        testBandPasses = {}
        keys = ['u','g','r','i','z']

        bplist = []

        for kk in keys:
            testBandPasses[kk] = Bandpass()
            testBandPasses[kk].readThroughput(os.path.join(bandpassDir,"test_bandpass_%s.dat" % kk))
            bplist.append(testBandPasses[kk])

        sedObj = Sed()
        phiArray, waveLenStep = sedObj.setupPhiArray(bplist)

        sedFileName = os.path.join(lsst.utils.getPackageDir('sims_photUtils'),
                                   'tests/cartoonSedTestData/starSed/')
        sedFileName = os.path.join(sedFileName,'kurucz','km20_5750.fits_g40_5790.gz')
        ss = Sed()
        ss.readSED_flambda(sedFileName)

        controlBandpass = Bandpass()
        controlBandpass.imsimBandpass()
        ff = ss.calcFluxNorm(22.0, controlBandpass)
        ss.multiplyFluxNorm(ff)

        testMags = cartoon_dict.magListForSed(ss)

        ss.resampleSED(wavelen_match = bplist[0].wavelen)
        ss.flambdaTofnu()
        mags = -2.5*np.log10(np.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
        self.assertEqual(len(mags), len(testMags))
        self.assertGreater(len(mags), 0)
        for j in range(len(mags)):
            self.assertAlmostEqual(mags[j],testMags[j],10)



    def testEBV(self):

        ebvObject = EBVbase()
        ra = []
        dec = []
        gLat = []
        gLon = []
        for i in range(10):
            ra.append(i*2.0*np.pi/10.0)
            dec.append(i*np.pi/10.0)

            gLat.append(-0.5*np.pi+i*np.pi/10.0)
            gLon.append(i*2.0*np.pi/10.0)

            equatorialCoordinates=np.array([ra,dec])
            galacticCoordinates=np.array([gLon,gLat])

        ebvOutput = ebvObject.calculateEbv(equatorialCoordinates=equatorialCoordinates)
        self.assertEqual(len(ebvOutput),len(ra))

        ebvOutput = ebvObject.calculateEbv(galacticCoordinates=galacticCoordinates)
        self.assertEqual(len(ebvOutput),len(gLon))
        self.assertGreater(len(ebvOutput), 0)

        self.assertRaises(RuntimeError, ebvObject.calculateEbv, equatorialCoordinates=equatorialCoordinates,
        galacticCoordinates=galacticCoordinates)
        self.assertRaises(RuntimeError, ebvObject.calculateEbv, equatorialCoordinates=None, galacticCoordinates=None)
        self.assertRaises(RuntimeError, ebvObject.calculateEbv)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
