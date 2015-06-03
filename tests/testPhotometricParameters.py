import os
import unittest
import eups
import lsst.utils.tests as utilsTests

from lsst.sims.photUtils import Bandpass, Sed, PhotometricParameters, \
                                PhysicalParameters

class PhotometricParametersUnitTest(unittest.TestCase):

    def testInit(self):
        """
        Test that the init and getters of PhotometricParameters work
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


    def testAssignment(self):
        """
        Test that it is impossible to set PhotometricParameters on the fly
        """
        testCase = PhotometricParameters()
        controlCase = PhotometricParameters()
        success = 0

        msg = ''
        try:
            testCase.exptime = -1.0
            success += 1
            msg += 'was able to assign exptime; '
        except:
            self.assertEqual(testCase.exptime, controlCase.exptime)

        try:
            testCase.nexp = -1.0
            success += 1
            msg += 'was able to assign nexp; '
        except:
            self.assertEqual(testCase.nexp, controlCase.nexp)

        try:
            testCase.effarea = -1.0
            success += 1
            msg += 'was able to assign effarea; '
        except:
            self.assertEqual(testCase.effarea, controlCase.effarea)

        try:
            testCase.gain = -1.0
            success += 1
            msg += 'was able to assign gain; '
        except:
            self.assertEqual(testCase.gain, controlCase.gain)

        try:
            testCase.readnoise = -1.0
            success += 1
            msg += 'was able to assign readnoise; '
        except:
            self.assertEqual(testCase.readnoise, controlCase.readnoise)

        try:
            testCase.darkcurrent = -1.0
            success += 1
            msg += 'was able to assign darkcurrent; '
        except:
            self.assertEqual(testCase.darkcurrent, controlCase.darkcurrent)

        try:
            testCase.othernoise = -1.0
            success += 1
            msg += 'was able to assign othernoise; '
        except:
            self.assertEqual(testCase.othernoise, controlCase.othernoise)

        try:
            testCase.platescale = -1.0
            success += 1
            msg += 'was able to assign platescale; '
        except:
            self.assertEqual(testCase.platescale, controlCase.platescale)

        self.assertEqual(success,0)


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

        control = testSed.calcADU(testBandpass,
                                  photParams=PhotometricParameters())

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
        control = PhysicalParameters()
        success = 0
        msg = ''

        try:
            pp.minwavelen = 2.0
            success += 1
            msg += 'was able to assign minwavelen; '
        except:
            self.assertEqual(pp.minwavelen, control.minwavelen)

        try:
            pp.maxwavelen = 2.0
            success += 1
            msg += 'was able to assign maxwavelen; '
        except:
            self.assertEqual(pp.maxwavelen, control.maxwavelen)

        try:
            pp.wavelenstep = 2.0
            success += 1
            msg += 'was able to assign wavelenstep; '
        except:
            self.assertEqual(pp.wavelenstep, control.wavelenstep)

        try:
            pp.lightspeed = 2.0
            success += 1
            msg += 'was able to assign lightspeed; '
        except:
            self.assertEqual(pp.lightspeed, control.lightspeed)

        try:
            pp.planck = 2.0
            success += 1
            msg += 'was able to assign planck; '
        except:
            self.assertEqual(pp.planck, control.planck)

        try:
            pp.nm2m = 2.0
            success += 1
            msg += 'was able to assign nm2m; '
        except:
            self.assertEqual(pp.nm2m, control.nm2m)

        try:
            pp.ergsetc2jansky = 2.0
            msg += 'was able to assign ergsetc2jansky; '
            success += 1
        except:
            self.assertEqual(pp.ergsetc2jansky, control.ergsetc2jansky)

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
