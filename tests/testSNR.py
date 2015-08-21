import os
import numpy
import unittest
import lsst.utils
import lsst.utils.tests as utilsTests
from lsst.sims.utils import ObservationMetaData
import lsst.sims.photUtils.SignalToNoise as snr
from lsst.sims.photUtils import Sed, Bandpass, PhotometricParameters, \
                                LSSTdefaults
from lsst.sims.photUtils.utils import setM5

class TestSNRmethods(unittest.TestCase):


    def setUp(self):

        starFileName = os.path.join(lsst.utils.getPackageDir('sims_sed_library'),'starSED')
        starFileName = os.path.join(starFileName, 'kurucz','km20_5750.fits_g40_5790.gz')
        starName = os.path.join(lsst.utils.getPackageDir('sims_sed_library'),starFileName)
        self.starSED = Sed()
        self.starSED.readSED_flambda(starName)
        imsimband = Bandpass()
        imsimband.imsimBandpass()
        fNorm = self.starSED.calcFluxNorm(22.0, imsimband)
        self.starSED.multiplyFluxNorm(fNorm)

        hardwareDir = os.path.join(lsst.utils.getPackageDir('throughputs'),'baseline')
        componentList = ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                         'lens1.dat', 'lens2.dat', 'lens3.dat']
        self.skySed = Sed()
        self.skySed.readSED_flambda(os.path.join(hardwareDir,'darksky.dat'))

        totalNameList = ['total_u.dat', 'total_g.dat', 'total_r.dat', 'total_i.dat',
                         'total_z.dat', 'total_y.dat']

        self.bpList = []
        self.hardwareList = []
        for name in totalNameList:
            dummy = Bandpass()
            dummy.readThroughput(os.path.join(hardwareDir, name))
            self.bpList.append(dummy)

            dummy = Bandpass()
            hardwareNameList = [os.path.join(hardwareDir, name)]
            for component in componentList:
                hardwareNameList.append(os.path.join(hardwareDir, component))
            dummy.readThroughputList(hardwareNameList)
            self.hardwareList.append(dummy)

        self.filterNameList = ['u', 'g', 'r', 'i', 'z', 'y']


    def testMagError(self):
        """
        Make sure that calcMagError_sed and calcMagError_m5
        agree to within 0.001
        """
        defaults = LSSTdefaults()
        photParams = PhotometricParameters()

        #create a cartoon spectrum to test on
        spectrum = Sed()
        spectrum.setFlatSED()
        spectrum.multiplyFluxNorm(1.0e-9)

        #find the magnitudes of that spectrum in our bandpasses
        magList = []
        for total in self.bpList:
            magList.append(spectrum.calcMag(total))
        magList = numpy.array(magList)

        #try for different normalizations of the skySED
        for fNorm in numpy.arange(1.0, 5.0, 1.0):
            self.skySed.multiplyFluxNorm(fNorm)
            m5List = []
            magSed = []
            for total, hardware, filterName in \
            zip(self.bpList, self.hardwareList, self.filterNameList):

                seeing = defaults.seeing(filterName)

                m5List.append(snr.calcM5(self.skySed, total, hardware, photParams,seeing=seeing))

                magSed.append(snr.calcMagError_sed(spectrum, total, self.skySed,
                                                   hardware, photParams, seeing=seeing))

            magSed = numpy.array(magSed)

            magM5 = snr.calcMagError_m5(magList, self.bpList,
                                        numpy.array(m5List), photParams)


            numpy.testing.assert_array_almost_equal(magM5, magSed, decimal=3)


    def testVerboseSNR(self):
        """
        Make sure that calcSNR_sed has everything it needs to run in verbose mode
        """
        defaults = LSSTdefaults()
        photParams = PhotometricParameters()

        #create a cartoon spectrum to test on
        spectrum = Sed()
        spectrum.setFlatSED()
        spectrum.multiplyFluxNorm(1.0e-9)

        snr.calcSNR_sed(spectrum, self.bpList[0], self.skySed,
                        self.hardwareList[0], photParams, seeing=0.7, verbose=True)


    def testSNRexceptions(self):
        """
        test that calcSNR_m5 raises an exception when arguments are not of the right shape.
        """

        photParams = PhotometricParameters()
        shortGamma = numpy.array([1.0, 1.0])
        shortMagnitudes = numpy.array([22.0, 23.0])
        magnitudes = 22.0*numpy.ones(6)
        self.assertRaises(RuntimeError, snr.calcSNR_m5, magnitudes, self.bpList, shortMagnitudes, photParams)
        self.assertRaises(RuntimeError, snr.calcSNR_m5, shortMagnitudes, self.bpList, magnitudes, photParams)
        self.assertRaises(RuntimeError, snr.calcSNR_m5, magnitudes, self.bpList, magnitudes, photParams, gamma=shortGamma)
        signalToNoise, gg = snr.calcSNR_m5(magnitudes, self.bpList, magnitudes, photParams)


    def testSignalToNoise(self):
        """
        Test that calcSNR_m5 and calcSNR_sed give similar results
        """
        defaults = LSSTdefaults()
        photParams = PhotometricParameters()

        m5 = []
        for i in range(len(self.hardwareList)):
            m5.append(snr.calcM5(self.skySed, self.bpList[i],
                      self.hardwareList[i],
                      photParams, seeing=defaults.seeing(self.filterNameList[i])))


        sedDir = lsst.utils.getPackageDir('sims_sed_library')
        sedDir = os.path.join(sedDir, 'starSED', 'kurucz')
        fileNameList = os.listdir(sedDir)

        numpy.random.seed(42)
        offset = numpy.random.random_sample(len(fileNameList))*2.0

        for ix, name in enumerate(fileNameList):
            if ix>100:
                break
            spectrum = Sed()
            spectrum.readSED_flambda(os.path.join(sedDir, name))
            ff = spectrum.calcFluxNorm(m5[2]-offset[ix], self.bpList[2])
            spectrum.multiplyFluxNorm(ff)
            magList = []
            controlList = []
            magList = []
            for i in range(len(self.bpList)):
                controlList.append(snr.calcSNR_sed(spectrum, self.bpList[i],
                                               self.skySed,
                                               self.hardwareList[i],
                                               photParams, defaults.seeing(self.filterNameList[i])))

                magList.append(spectrum.calcMag(self.bpList[i]))

            testList, gammaList = snr.calcSNR_m5(numpy.array(magList),
                                        numpy.array(self.bpList),
                                        numpy.array(m5),
                                        photParams)

            for tt, cc in zip(controlList, testList):
                msg = '%e != %e ' % (tt, cc)
                self.assertTrue(numpy.abs(tt/cc - 1.0) < 0.001, msg=msg)



    def testSystematicUncertainty(self):
        """
        Test that systematic uncertainty is added correctly.
        """
        sigmaSys = 0.002
        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]
        photParams= PhotometricParameters(sigmaSys=sigmaSys)

        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5, bandpassName=self.filterNameList)
        magnitudes = []
        for bp in self.bpList:
            mag = self.starSED.calcMag(bp)
            magnitudes.append(mag)

        skySedList = []

        for bp, hardware, filterName in zip(self.bpList, self.hardwareList, self.filterNameList):
            skyDummy = Sed()
            skyDummy.readSED_flambda(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSkyDummy = setM5(obs_metadata.m5[filterName], skyDummy,
                                       bp, hardware,
                                       seeing=LSSTdefaults().seeing(filterName),
                                       photParams=photParams)

            skySedList.append(normalizedSkyDummy)

        sigmaList = snr.calcMagError_m5(numpy.array(magnitudes), numpy.array(self.bpList), \
                                        numpy.array(m5), photParams)

        for i in range(len(self.bpList)):
            snrat = snr.calcSNR_sed(self.starSED, self.bpList[i], skySedList[i], self.hardwareList[i],
                                  seeing=LSSTdefaults().seeing(self.filterNameList[i]),
                                  photParams=PhotometricParameters())

            testSNR, gamma = snr.calcSNR_m5(numpy.array([magnitudes[i]]), [self.bpList[i]],
                                           numpy.array([m5[i]]), photParams=PhotometricParameters(sigmaSys=0.0))

            self.assertAlmostEqual(snrat, testSNR[0], 10, msg = 'failed on calcSNR_m5 test %e != %e ' \
                                                               % (snrat, testSNR[0]))

            control = numpy.sqrt(numpy.power(snr.magErrorFromSNR(testSNR),2) + numpy.power(sigmaSys,2))

            msg = '%e is not %e; failed' % (sigmaList[i], control)

            self.assertAlmostEqual(sigmaList[i], control, 10, msg=msg)




    def testNoSystematicUncertainty(self):
        """
        Test that systematic uncertainty is handled correctly when set to None.
        """
        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]
        photParams= PhotometricParameters(sigmaSys=0.0)

        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5, bandpassName=self.filterNameList)

        magnitudes = []
        for bp in self.bpList:
            mag = self.starSED.calcMag(bp)
            magnitudes.append(mag)

        skySedList = []

        for bp, hardware, filterName in zip(self.bpList, self.hardwareList, self.filterNameList):
            skyDummy = Sed()
            skyDummy.readSED_flambda(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSkyDummy = setM5(obs_metadata.m5[filterName], skyDummy,
                                       bp, hardware,
                                       seeing=LSSTdefaults().seeing(filterName),
                                       photParams=photParams)

            skySedList.append(normalizedSkyDummy)

        sigmaList = snr.calcMagError_m5(numpy.array(magnitudes), numpy.array(self.bpList), \
                                        numpy.array(m5), photParams)


        for i in range(len(self.bpList)):
            snrat = snr.calcSNR_sed(self.starSED, self.bpList[i], skySedList[i], self.hardwareList[i],
                              seeing=LSSTdefaults().seeing(self.filterNameList[i]),
                              photParams=PhotometricParameters())

            testSNR, gamma = snr.calcSNR_m5(numpy.array([magnitudes[i]]), [self.bpList[i]],
                                           numpy.array([m5[i]]), photParams=PhotometricParameters(sigmaSys=0.0))

            self.assertAlmostEqual(snrat, testSNR[0], 10, msg = 'failed on calcSNR_m5 test %e != %e ' \
                                                               % (snrat, testSNR[0]))

            control = snr.magErrorFromSNR(testSNR)

            msg = '%e is not %e; failed' % (sigmaList[i], control)

            self.assertAlmostEqual(sigmaList[i], control, 10, msg=msg)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(TestSNRmethods)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
