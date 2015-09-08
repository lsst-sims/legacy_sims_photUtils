import numpy

import os
import unittest
import lsst.utils
import lsst.utils.tests as utilsTests
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import defaultSpecMap
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils import loadBandpassesFromFiles, \
                                loadTotalBandpassesFromFiles
from lsst.sims.photUtils import LSSTdefaults, PhotometricParameters, calcSNR_m5, \
                                calcM5, calcSNR_sed, magErrorFromSNR
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

        bandpassDir=os.path.join(lsst.utils.getPackageDir('sims_photUtils'),'tests','cartoonSedTestData')

        cartoon_dict = loadTotalBandpassesFromFiles(['u','g','r','i','z'],bandpassDir = bandpassDir,
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

        sedFileName = os.path.join(lsst.utils.getPackageDir('sims_sed_library'),'starSED','kurucz')
        sedFileName = os.path.join(sedFileName,'km20_5750.fits_g40_5790.gz')
        ss = Sed()
        ss.readSED_flambda(sedFileName)

        controlBandpass = Bandpass()
        controlBandpass.imsimBandpass()
        ff = ss.calcFluxNorm(22.0, controlBandpass)
        ss.multiplyFluxNorm(ff)

        testMags = cartoon_dict.calcMagListFromSed(ss)

        ss.resampleSED(wavelen_match = bplist[0].wavelen)
        ss.flambdaTofnu()
        mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
        self.assertTrue(len(mags)==len(testMags))
        self.assertTrue(len(mags)>0)
        for j in range(len(mags)):
            self.assertAlmostEqual(mags[j],testMags[j],10)



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
        starName = os.path.join(lsst.utils.getPackageDir('sims_sed_library'),defaultSpecMap['km20_5750.fits_g40_5790'])
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
            hardwareComponents.append(os.path.join(lsst.utils.getPackageDir('throughputs'),'baseline',c))

        self.bandpasses = ['u', 'g', 'r', 'i', 'z', 'y']
        for b in self.bandpasses:
            filterName = os.path.join(lsst.utils.getPackageDir('throughputs'),'baseline','filter_%s.dat' % b)
            components = hardwareComponents + [filterName]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.hardwareBandpasses.append(bandpassDummy)
            components = components + [os.path.join(lsst.utils.getPackageDir('throughputs'),'baseline','atmos.dat')]
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
        Test that calcSNR_m5 raises exceptions when it needs to
        """
        totalDict, hardwareDict = loadBandpassesFromFiles()
        magnitudes = numpy.array([22.0, 23.0, 24.0, 25.0, 26.0, 27.0])
        shortMagnitudes = numpy.array([22.0])
        photParams = PhotometricParameters()
        shortGamma = numpy.array([1.0, 1.0])
        self.assertRaises(RuntimeError, calcSNR_m5, magnitudes, totalDict.values(), shortMagnitudes, photParams)
        self.assertRaises(RuntimeError, calcSNR_m5, shortMagnitudes, totalDict.values(), magnitudes, photParams)
        self.assertRaises(RuntimeError, calcSNR_m5, magnitudes, totalDict.values(), magnitudes, photParams, gamma=shortGamma)
        snr, gg = calcSNR_m5(magnitudes, totalDict.values(), magnitudes, photParams)


    def testSignalToNoise(self):
        """
        Test that calcSNR_m5 and calcSNR_sed give similar results
        """
        defaults = LSSTdefaults()
        photParams = PhotometricParameters()
        totalDict, hardwareDict = loadBandpassesFromFiles()

        skySED = Sed()
        skySED.readSED_flambda(os.path.join(lsst.utils.getPackageDir('throughputs'),
                                            'baseline', 'darksky.dat'))

        m5 = []
        for filt in totalDict:
            m5.append(calcM5(skySED, totalDict[filt],
                      hardwareDict[filt],
                      photParams, seeing=defaults.seeing(filt)))


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
            ff = spectrum.calcFluxNorm(m5[2]-offset[ix], totalDict.values()[2])
            spectrum.multiplyFluxNorm(ff)
            magList = []
            controlList = []
            magList = []
            for filt in totalDict:
                controlList.append(calcSNR_sed(spectrum, totalDict[filt],
                                               skySED,
                                               hardwareDict[filt],
                                               photParams, defaults.seeing(filt)))

                magList.append(spectrum.calcMag(totalDict[filt]))

            testList, gammaList = calcSNR_m5(numpy.array(magList),
                                        numpy.array(totalDict.values()),
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

        bandpassDict = loadTotalBandpassesFromFiles()
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5, bandpassName=self.bandpasses)
        magnitudes = bandpassDict.calcMagListFromSed(self.starSED)

        skySeds = []

        for i in range(len(self.bandpasses)):
            skyDummy = Sed()
            skyDummy.readSED_flambda(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSkyDummy = setM5(obs_metadata.m5[self.bandpasses[i]], skyDummy,
                                                       self.totalBandpasses[i], self.hardwareBandpasses[i],
                                                       seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                                                       photParams=PhotometricParameters())

            skySeds.append(normalizedSkyDummy)


        for i in range(len(self.bandpasses)):
            snr = calcSNR_sed(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            testSNR, gamma = calcSNR_m5(numpy.array([magnitudes[i]]), [self.totalBandpasses[i]],
                                           numpy.array([m5[i]]), photParams=PhotometricParameters(sigmaSys=0.0))

            self.assertAlmostEqual(snr, testSNR[0], 10, msg = 'failed on calcSNR_m5 test %e != %e ' \
                                                               % (snr, testSNR[0]))

            control = numpy.sqrt(numpy.power(magErrorFromSNR(testSNR),2) + numpy.power(sigmaSys,2))



    def testNoSystematicUncertainty(self):
        """
        Test that systematic uncertainty is handled correctly when set to None.
        """
        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]
        photParams= PhotometricParameters(sigmaSys=0.0)

        bandpassDict = loadTotalBandpassesFromFiles()
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5, bandpassName=self.bandpasses)
        magnitudes = bandpassDict.calcMagListFromSed(self.starSED)

        skySeds = []

        for i in range(len(self.bandpasses)):
            skyDummy = Sed()
            skyDummy.readSED_flambda(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSkyDummy = setM5(obs_metadata.m5[self.bandpasses[i]], skyDummy,
                                                       self.totalBandpasses[i], self.hardwareBandpasses[i],
                                                       seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                                                       photParams=PhotometricParameters())

            skySeds.append(normalizedSkyDummy)

        for i in range(len(self.bandpasses)):
            snr = calcSNR_sed(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            testSNR, gamma = calcSNR_m5(numpy.array([magnitudes[i]]), [self.totalBandpasses[i]],
                                           numpy.array([m5[i]]), photParams=PhotometricParameters(sigmaSys=0.0))

            self.assertAlmostEqual(snr, testSNR[0], 10, msg = 'failed on calcSNR_m5 test %e != %e ' \
                                                               % (snr, testSNR[0]))



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
