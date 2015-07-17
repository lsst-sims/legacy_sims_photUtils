import numpy

import os
import unittest
import eups
import lsst.utils.tests as utilsTests
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.generation.utils import myTestGals, myTestStars, \
                                                makeStarTestDB, makeGalTestDB, getOneChunk

from lsst.sims.catalogs.measures.instance import defaultSpecMap
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils import PhotometryBase, PhotometryHardware
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

        bandpassDir=os.path.join(eups.productDir('sims_photUtils'),'tests','cartoonSedTestData')

        cartoon_phot = PhotometryBase()
        cartoon_phot.loadTotalBandpassesFromFiles(['u','g','r','i','z'],bandpassDir = bandpassDir,
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

        sedFileName = os.path.join(eups.productDir('sims_sed_library'),'starSED','kurucz')
        sedFileName = os.path.join(sedFileName,'km20_5750.fits_g40_5790.gz')
        ss = Sed()
        ss.readSED_flambda(sedFileName)
        
        controlBandpass = Bandpass()
        controlBandpass.imsimBandpass()
        ff = ss.calcFluxNorm(22.0, controlBandpass)
        ss.multiplyFluxNorm(ff)

        testMags = cartoon_phot.manyMagCalc_list(ss)

        ss.resampleSED(wavelen_match = bplist[0].wavelen)
        ss.flambdaTofnu()
        mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
        self.assertTrue(len(mags)==len(testMags))
        self.assertTrue(len(mags)>0)
        for j in range(len(mags)):
            self.assertAlmostEqual(mags[j],testMags[j],10)



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
        self.assertRaises(RuntimeError, calcSNR_m5, magnitudes, phot.bandpassDict.values(), shortMagnitudes, photParams)
        self.assertRaises(RuntimeError, calcSNR_m5, shortMagnitudes, phot.bandpassDict.values(), magnitudes, photParams)
        self.assertRaises(RuntimeError, calcSNR_m5, magnitudes, phot.bandpassDict.values(), magnitudes, photParams, gamma=shortGamma)
        snr, gg = calcSNR_m5(magnitudes, phot.bandpassDict.values(), magnitudes, photParams)


    def testSignalToNoise(self):
        """
        Test that calcSNR_m5 and calcSNR_sed give similar results
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
            magList = []
            controlList = []
            magList = []
            for filt in hardware.bandpassDict:
                controlList.append(calcSNR_sed(spectrum, hardware.bandpassDict[filt],
                                               hardware.skySED,
                                               hardware.hardwareBandpassDict[filt],
                                               photParams, defaults.seeing(filt)))

                magList.append(spectrum.calcMag(hardware.bandpassDict[filt]))

            testList, gammaList = calcSNR_m5(numpy.array(magList),
                                        numpy.array(hardware.bandpassDict.values()),
                                        numpy.array(m5),
                                        photParams)

            for tt, cc in zip(controlList, testList):
                msg = '%e != %e ' % (tt, cc)
                self.assertTrue(numpy.abs(tt/cc - 1.0) < 0.001, msg=msg)



    def testRawUncertainty(self):
        """
        Test that values calculated by calculatePhotometricUncertainty agree
        with values calculated by calcSNR_sed
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
            snr = calcSNR_sed(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            ss = 2.5*numpy.log10(1.0+1.0/snr)
            ss = numpy.sqrt(ss*ss + numpy.power(phot.photParams.sigmaSys,2))
            msg = '%e is not %e; failed' % (ss, sigma[i])
            self.assertAlmostEqual(ss, sigma[i], 10, msg=msg)

    def testSystematicUncertainty(self):
        """
        Test that systematic uncertainty is added correctly.
        """
        sigmaSys = 0.002
        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]
        photParams= PhotometricParameters(sigmaSys=sigmaSys)

        phot = PhotometryBase()
        phot.photParams = photParams

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
            snr = calcSNR_sed(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            testSNR, gamma = calcSNR_m5(numpy.array([magnitudes[i]]), [self.totalBandpasses[i]],
                                           numpy.array([m5[i]]), photParams=PhotometricParameters(sigmaSys=0.0))

            self.assertAlmostEqual(snr, testSNR[0], 10, msg = 'failed on calcSNR_m5 test %e != %e ' \
                                                               % (snr, testSNR[0]))

            control = numpy.sqrt(numpy.power(magErrorFromSNR(testSNR),2) + numpy.power(sigmaSys,2))

            msg = '%e is not %e; failed' % (sigma[i], control)

            self.assertAlmostEqual(sigma[i], control, 10, msg=msg)


    def testNoSystematicUncertainty(self):
        """
        Test that systematic uncertainty is handled correctly when set to None.
        """
        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]
        photParams= PhotometricParameters(sigmaSys=0.0)

        phot = PhotometryBase()
        phot.photParams = photParams

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
            snr = calcSNR_sed(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            testSNR, gamma = calcSNR_m5(numpy.array([magnitudes[i]]), [self.totalBandpasses[i]],
                                           numpy.array([m5[i]]), photParams=PhotometricParameters(sigmaSys=0.0))

            self.assertAlmostEqual(snr, testSNR[0], 10, msg = 'failed on calcSNR_m5 test %e != %e ' \
                                                               % (snr, testSNR[0]))

            control = magErrorFromSNR(testSNR)

            msg = '%e is not %e; failed' % (sigma[i], control)

            self.assertAlmostEqual(sigma[i], control, 10, msg=msg)


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
