import numpy as np
import warnings
import unittest
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

    def testSedWavelenRange(self):
        """Test setting sed with wavelength range different from standard values works properly."""
        sedwavelen = self.bandpasswavelen * 1.0
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        np.testing.assert_equal(testsed.wavelen, sedwavelen)
        np.testing.assert_equal(testsed.flambda, sedflambda)

    def testSedBandpassMatch(self):
        """Test errors when bandpass and sed do not completely overlap in wavelength range."""
        # Test case where they do match (no error message)
        sedwavelen = np.arange(self.wmin, self.wmax+.5, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        print ''
        # Test that no warning or exception is made.
        w,f = testsed.resampleSED(wavelen_match=self.testbandpass.wavelen, wavelen=testsed.wavelen, flux=testsed.flambda)
        np.testing.assert_equal(w, testsed.wavelen)
        np.testing.assert_equal(f, testsed.flambda)
        # Test that exception is made for non-overlap at either top or bottom end of wavelength range.
        sedwavelen = np.arange(self.wmin, self.wmax - 50, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        self.assertRaises(ValueError, testsed.resampleSED, wavelen_match=self.testbandpass.wavelen)
        sedwavelen = np.arange(self.wmin+50, self.wmax, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        self.assertRaises(ValueError, testsed.resampleSED, wavelen_match=self.testbandpass.wavelen)

    def testSedMagErrors(self):
        """Test error handling at mag and adu calculation levels of sed."""
        sedwavelen = np.arange(self.wmin+50, self.wmax, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        # Test handling in calcMag
        with warnings.catch_warnings(record=True) as w:
            mag = testsed.calcMag(self.testbandpass)
            self.assertEqual(len(w), 1)
            self.assertTrue("wavelength" in str(w[-1].message))
            self.assertTrue("non-overlap" in str(w[-1].message))
        self.assertTrue(np.isnan(mag))
        # Test handling in calcADU
        with warnings.catch_warnings(record=True) as w:
            adu = testsed.calcADU(self.testbandpass)
            self.assertEqual(len(w), 1)
            self.assertTrue("wavelength" in str(w[-1].message))
            self.assertTrue("non-overlap" in str(w[-1].message))
        self.assertTrue(np.isnan(adu))
        # Test handling in calcFlux
        with warnings.catch_warnings(record=True) as w:
            flux = testsed.calcFlux(self.testbandpass)
            self.assertEqual(len(w), 1)
            self.assertTrue("wavelength" in str(w[-1].message))
            self.assertTrue("non-overlap" in str(w[-1].message))
        self.assertTrue(np.isnan(flux))
        
if __name__ == "__main__":
    suitelist = []
    suitelist.append(unittest.TestLoader().loadTestsFromTestCase(TestSedWavelenLimits))
    suite = unittest.TestSuite(suitelist)
    unittest.TextTestRunner(verbosity=2).run(suite)                            
