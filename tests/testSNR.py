import os
import numpy
import eups
import unittest
import lsst.utils.tests as utilsTests

import lsst.sims.photUtils.SignalToNoise as snr
from lsst.sims.photUtils import Sed, Bandpass, PhotometricParameters, \
                                LSSTdefaults

class TestSNRmethods(unittest.TestCase):


    def setUp(self):
        hardwareDir = os.path.join(eups.productDir('throughputs'),'baseline')
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


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(TestSNRmethods)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
