import numpy as np
import warnings
import unittest
import os
import eups
import lsst.utils.tests as utilsTests
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass

class TestSedWavelenLimits(unittest.TestCase):
    def setUp(self):
        self.wmin = 500
        self.wmax = 1500
        self.bandpasswavelen = np.arange(self.wmin, self.wmax+.5, 1)
        self.bandpasssb = np.ones(len(self.bandpasswavelen))
        self.testbandpass = Bandpass(wavelen=self.bandpasswavelen, sb=self.bandpasssb)

    def tearDown(self):
        del self.bandpasswavelen
        del self.bandpasssb
        del self.testbandpass
        del self.wmin
        del self.wmax

    def testSedWavelenRange(self):
        """Test setting sed with wavelength range different from standard values works properly."""
        sedwavelen = self.bandpasswavelen * 1.0
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda, name='TestSed')
        np.testing.assert_equal(testsed.wavelen, sedwavelen)
        np.testing.assert_equal(testsed.flambda, sedflambda)
        self.assertEqual(testsed.name, 'TestSed')

    def testSedBandpassMatch(self):
        """Test errors when bandpass and sed do not completely overlap in wavelength range."""
        # Test case where they do match (no error message)
        sedwavelen = np.arange(self.wmin, self.wmax+.5, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        print ''
        # Test that no warning is made.
        with warnings.catch_warnings(record=True) as wa:
            w,f = testsed.resampleSED(wavelen_match=self.testbandpass.wavelen, wavelen=testsed.wavelen, flux=testsed.flambda)
            self.assertEqual(len(wa), 0)
        np.testing.assert_equal(w, testsed.wavelen)
        np.testing.assert_equal(f, testsed.flambda)
        # Test that warning is given for non-overlap at either top or bottom end of wavelength range.
        sedwavelen = np.arange(self.wmin, self.wmax - 50, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        with warnings.catch_warnings(record=True) as wa:
            testsed.resampleSED(wavelen_match=self.testbandpass.wavelen)
            self.assertEqual(len(wa), 1)
            self.assertTrue('non-overlap' in str(wa[-1].message))
        self.assertTrue(np.isnan(testsed.flambda[-1:]))
        sedwavelen = np.arange(self.wmin+50, self.wmax, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        with warnings.catch_warnings(record=True) as wa:
            testsed.resampleSED(wavelen_match=self.testbandpass.wavelen)
            self.assertEqual(len(wa), 1)
            self.assertTrue('non-overlap' in str(wa[-1].message))
        self.assertTrue(np.isnan(testsed.flambda[0]))
        self.assertTrue(np.isnan(testsed.flambda[49]))

    def testSedMagErrors(self):
        """Test error handling at mag and adu calculation levels of sed."""
        sedwavelen = np.arange(self.wmin+50, self.wmax, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        # Test handling in calcMag
        with warnings.catch_warnings(record=True) as w:
            mag = testsed.calcMag(self.testbandpass)
            self.assertEqual(len(w), 1)
            self.assertTrue("non-overlap" in str(w[-1].message))
        self.assertTrue(np.isnan(mag))
        # Test handling in calcADU
        with warnings.catch_warnings(record=True) as w:
            adu = testsed.calcADU(self.testbandpass)
            self.assertEqual(len(w), 1)
            self.assertTrue("non-overlap" in str(w[-1].message))
        self.assertTrue(np.isnan(adu))
        # Test handling in calcFlux
        with warnings.catch_warnings(record=True) as w:
            flux = testsed.calcFlux(self.testbandpass)
            self.assertEqual(len(w), 1)
            self.assertTrue("non-overlap" in str(w[-1].message))
        self.assertTrue(np.isnan(flux))


class TestSedName(unittest.TestCase):
    def setUp(self):
        self.wmin = 500
        self.wmax = 1500
        self.wavelen = np.arange(self.wmin, self.wmax+.5, 1)
        self.flambda = np.ones(len(self.wavelen))
        self.name = 'TestSed'
        self.testsed = Sed(self.wavelen, self.flambda, name=self.name)

    def tearDown(self):
        del self.wmin, self.wmax, self.wavelen, self.flambda
        del self.name
        del self.testsed


    def testSetName(self):
        self.assertEqual(self.testsed.name, self.name)

    def testRedshiftName(self):
        testsed = Sed(self.testsed.wavelen, self.testsed.flambda, name=self.testsed.name)
        redshift = .2
        testsed.redshiftSED(redshift=redshift)
        newname = testsed.name + '_Z' + '%.2f' %(redshift)
        testsed.name = newname
        self.assertEqual(testsed.name, newname)

class TestM5(unittest.TestCase):

    def testSettingM5(self):
        """
        Test that the setM5 method in Bandpass actually does set m5 to the
        desired value
        """
        np.random.seed(32)
        m5Targets = 15.0 + np.random.sample(10)*7.0
        thrDir = eups.productDir('throughputs')
        skysedName = os.path.join(thrDir,'baseline','darksky.dat')
        components = ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                      'lens1.dat', 'lens2.dat', 'lens3.dat']
        bandpasses = ['u','g','r','i','z','y']
        for b in bandpasses:
            totalBandpass = Bandpass()
            totalBandpass.readThroughput(os.path.join(thrDir,'baseline','total_%s.dat'%b))

            bplist = []
            for c in components:
                bplist.append(os.path.join(thrDir,'baseline',c))
            bplist += [os.path.join(thrDir,'baseline','filter_%s.dat'%b)]

            hardwareBandpass = Bandpass()
            hardwareBandpass.readThroughputList(bplist)

            for m in m5Targets:
                skysed = Sed()
                skysed.readSED_flambda(skysedName)

                totalBandpass.setM5(m, skysed, hardwareBandpass)
                m5Result = totalBandpass.calcM5(skysed, hardwareBandpass)
                self.assertAlmostEqual(m/m5Result, 1.0, 6)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(TestSedWavelenLimits)
    suites += unittest.makeSuite(TestSedName)
    suites += unittest.makeSuite(TestM5)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)

