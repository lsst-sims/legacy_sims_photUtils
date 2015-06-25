import os
import numpy
import eups
import unittest
import lsst.utils.tests as utilsTests

import lsst.sims.photUtils.SignalToNoise as snr
from lsst.sims.photUtils import Sed, PhotometryHardware, PhotometricParameters, \
                                LSSTdefaults

class TestSNRmethods(unittest.TestCase):

    def testMagError(self):
        """
        Make sure that calcMagError_sed and calcMagError_m5
        agree to within 0.001
        """
        defaults = LSSTdefaults()
        photParams = PhotometricParameters()
        phot = PhotometryHardware()
        phot.loadBandpassesFromFiles()

        #create a cartoon spectrum to test on
        spectrum = Sed()
        spectrum.setFlatSED()
        spectrum.multiplyFluxNorm(1.0e-9)

        #find the magnitudes of that spectrum in our bandpasses
        magList = []
        for total in phot.bandpassDict.values():
            magList.append(spectrum.calcMag(total))
        magList = numpy.array(magList)

        #try for different normalizations of the skySED
        for fNorm in numpy.arange(1.0, 5.0, 1.0):
            phot.skySED.multiplyFluxNorm(fNorm)
            m5List = []
            magSed = []
            for bp in phot.bandpassDict:
                total = phot.bandpassDict[bp]
                hardware = phot.hardwareBandpassDict[bp]
                seeing = defaults.seeing(bp)

                m5List.append(snr.calcM5(phot.skySED, total, hardware, photParams,seeing=seeing))

                magSed.append(snr.calcMagError_sed(spectrum, total, phot.skySED,
                                                   hardware, photParams, seeing=seeing))

            magSed = numpy.array(magSed)

            magM5 = snr.calcMagError_m5(magList, phot.bandpassDict.values(),
                                        numpy.array(m5List), photParams)


            numpy.testing.assert_array_almost_equal(magM5, magSed, decimal=3)


    def testSkyCounts(self):
        """
        Make sure that sky counts per pixel and total sky counts are properly
        related.
        """
        phot = PhotometryHardware()
        phot.loadBandpassesFromFiles()
        seeing = 0.6
        platescale = 0.1
        photParams = PhotometricParameters(platescale=platescale)
        neff = snr.calcNeff(seeing=seeing, platescale=platescale)
        m5 = 22.0
        for bp in phot.bandpassDict:
            ctsPerPixel = snr.calcSkyCountsPerPixelForM5(m5, phot.bandpassDict[bp],
                                                         photParams, seeing=seeing)

            ctsTotal = snr.calcSkyCountsForM5(m5, phot.bandpassDict[bp], photParams, seeing=seeing)

            self.assertAlmostEqual(ctsTotal/ctsPerPixel, neff, 7)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(TestSNRmethods)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
