import os
import unittest
import eups
import lsst.utils.tests as utilsTests

from lsst.sims.photUtils import Bandpass, Sed, PhotometricParameters, \
                                PhysicalParameters

class PhotometricParametersUnitTest(unittest.TestCase):

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


class PhysicalParametersUnitTest(unittest.TestCase):

    def testAssignment(self):
        """
        Make sure it is impossible to change the values stored in
        PhysicalParameters
        """

        pp = PhysicalParameters()
        success = 0
        msg = ''

        try:
            pp.minwavelen = 2.0
            success += 1
            msg += 'was able to assign minwavelen; '
        except:
            pass

        try:
            pp.maxwavelen = 2.0
            success += 1
            msg += 'was able to assign maxwavelen; '
        except:
            pass

        try:
            pp.wavelenstep = 2.0
            success += 1
            msg += 'was able to assign wavelenstep; '
        except:
            pass

        try:
            pp.lightspeed = 2.0
            success += 1
            msg += 'was able to assign lightspeed; '
        except:
            pass

        try:
            pp.planck = 2.0
            success += 1
            msg += 'was able to assign planck; '
        except:
            pass

        try:
            pp.nm2m = 2.0
            success += 1
            msg += 'was able to assign nm2m; '
        except:
            pass

        try:
            pp.ergsetc2jansky = 2.0
            msg += 'was able to assign ergsetc2jansky; '
            success += 1
        except:
            pass

        self.assertEqual(success, 0, msg=msg)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(PhotometricParametersUnitTest)
    suites += unittest.makeSuite(PhysicalParametersUnitTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
