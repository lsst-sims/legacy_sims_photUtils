import os
import unittest
import eups
import lsst.utils.tests as utilsTests

from lsst.sims.photUtils import Bandpass, Sed, PhotometricParameters

class ParametersUnitTest(unittest.TestCase):

    def testAssignment(self):
        """
        Test that the getters and setters of PhotometricParameters work
        properly
        """
        defaults = PhotometricParameters()
        params = ['exptime', 'nexp', 'effarea',
                  'gain', 'readnoise', 'darkcurrent',
                  'othernoise', 'platescale']

        for attribute in params:
            kwargs = {}
            kwargs[attribute] = -100.0
            testCase = PhotometricParameters(**kwargs)

            for pp in params:
                if pp != attribute:
                    self.assertEqual(defaults.__getattribute__(pp),
                                     testCase.__getattribute__(pp))
                else:
                    self.assertNotEqual(defaults.__getattribute__(pp),
                                        testCase.__getattribute__(pp))

                    self.assertEqual(testCase.__getattribute__(pp), -100.0)


    def testApplication(self):
        """
        Test that PhotometricParameters get properly propagated into
        Sed methods.  We will test this using Sed.calcADU, since the ADU
        scale linearly with the appropriate parameter.
        """

        testSed = Sed()
        testSed.setFlatSED()

        testBandpass = Bandpass()
        testBandpass.readThroughput(os.path.join(eups.productDir('throughputs'),
                                                 'baseline','total_g.dat'))

        control = testSed.calcADU(testBandpass)

        testCase = PhotometricParameters(exptime=30.0)

        test = testSed.calcADU(testBandpass, photParams=testCase)

        self.assertTrue(control>0.0)
        self.assertEqual(control, 0.5*test)



def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(ParametersUnitTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
