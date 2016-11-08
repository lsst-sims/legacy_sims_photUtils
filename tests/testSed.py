from __future__ import with_statement
import numpy as np
import warnings
import unittest
import gzip
import os

import lsst.utils.tests
from  lsst.utils import getPackageDir
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.sims.photUtils import PhotometricParameters


def setup_module(module):
    lsst.utils.tests.init()


class TestSedWavelenLimits(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter('always')
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
            w, f = testsed.resampleSED(wavelen_match=self.testbandpass.wavelen,
                                       wavelen=testsed.wavelen, flux=testsed.flambda)
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
            self.assertIn('non-overlap', str(wa[-1].message))
        np.testing.assert_equal(testsed.flambda[-1:], np.NaN)
        sedwavelen = np.arange(self.wmin+50, self.wmax, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        with warnings.catch_warnings(record=True) as wa:
            testsed.resampleSED(wavelen_match=self.testbandpass.wavelen)
            self.assertEqual(len(wa), 1)
            self.assertIn('non-overlap', str(wa[-1].message))
        np.testing.assert_equal(testsed.flambda[0], np.NaN)
        np.testing.assert_equal(testsed.flambda[49], np.NaN)

    def testSedMagErrors(self):
        """Test error handling at mag and adu calculation levels of sed."""
        sedwavelen = np.arange(self.wmin+50, self.wmax, 1)
        sedflambda = np.ones(len(sedwavelen))
        testsed = Sed(wavelen=sedwavelen, flambda=sedflambda)
        # Test handling in calcMag
        with warnings.catch_warnings(record=True) as w:
            mag = testsed.calcMag(self.testbandpass)
            self.assertEqual(len(w), 1)
            self.assertIn("non-overlap", str(w[-1].message))
        np.testing.assert_equal(mag, np.NaN)
        # Test handling in calcADU
        with warnings.catch_warnings(record=True) as w:
            adu = testsed.calcADU(self.testbandpass,
                                  photParams=PhotometricParameters())
            self.assertEqual(len(w), 1)
            self.assertIn("non-overlap", str(w[-1].message))
        np.testing.assert_equal(adu, np.NaN)
        # Test handling in calcFlux
        with warnings.catch_warnings(record=True) as w:
            flux = testsed.calcFlux(self.testbandpass)
            self.assertEqual(len(w), 1)
            self.assertIn("non-overlap", str(w[-1].message))
        np.testing.assert_equal(flux, np.NaN)


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
        newname = testsed.name + '_Z' + '%.2f' % (redshift)
        testsed.name = newname
        self.assertEqual(testsed.name, newname)


class SedBasicFunctionsTestCase(unittest.TestCase):

    longMessage = True

    def test_read_sed_flambda(self):
        """
        Test how readSED_flambda handles the reading of SED filenames
        when we fail to correctly specify their gzipped state.
        """
        scratch_dir = os.path.join(getPackageDir("sims_photUtils"),
                                   "tests", "scratchSpace")

        rng = np.random.RandomState(88)
        zipped_name = os.path.join(scratch_dir, "zipped_sed.txt.gz")
        unzipped_name = os.path.join(scratch_dir, "unzipped_sed.txt")
        if os.path.exists(zipped_name):
            os.unlink(zipped_name)
        if os.path.exists(unzipped_name):
            os.unlink(unzipped_name)
        wv = np.arange(100.0, 1000.0, 10.0)
        flux = rng.random_sample(len(wv))
        with gzip.open(zipped_name, "w") as output_file:
            for ww, ff in zip(wv, flux):
                output_file.write("%e %e\n" % (ww, ff))
        with open(unzipped_name, "w") as output_file:
            for ww, ff in zip(wv, flux):
                output_file.write("%e %e\n" % (ww, ff))

        ss = Sed()
        ss.readSED_flambda(zipped_name)
        ss.readSED_flambda(zipped_name[:-3])
        ss.readSED_flambda(unzipped_name)
        ss.readSED_flambda(unzipped_name+'.gz')

        # make sure an error is raised when you try to read
        # a file that does not exist
        with self.assertRaises(IOError) as context:
            ss.readSED_flambda(os.path.join(scratch_dir, "nonsense.txt"))
        self.assertIn("sed file", context.exception.message)

        if os.path.exists(zipped_name):
            os.unlink(zipped_name)
        if os.path.exists(unzipped_name):
            os.unlink(unzipped_name)

    def test_eq(self):
        """
        Test that __eq__ in Sed works correctly
        """
        sed_dir = os.path.join(getPackageDir('sims_sed_library'), 'starSED', 'kurucz')
        list_of_seds = os.listdir(sed_dir)
        sedname1 = os.path.join(sed_dir, list_of_seds[0])
        sedname2 = os.path.join(sed_dir, list_of_seds[1])
        ss1 = Sed()
        ss1.readSED_flambda(sedname1)
        ss2 = Sed()
        ss2.readSED_flambda(sedname2)
        ss3 = Sed()
        ss3.readSED_flambda(sedname1)

        self.assertFalse(ss1 == ss2)
        self.assertTrue(ss1 != ss2)
        self.assertTrue(ss1 == ss3)
        self.assertFalse(ss1 != ss3)

        ss3.flambdaTofnu()

        self.assertFalse(ss1 == ss3)
        self.assertTrue(ss1 != ss3)

    def test_cache(self):
        """
        Verify that loading an SED from the cache gives identical
        results to loading the same SED from ASCII (since we are
        not calling cache_LSST_seds(), as soon as we load an SED
        with readSED_flambda, it should get stored in the
        _global_misc_sed_cache)
        """
        sed_dir = os.path.join(getPackageDir('sims_sed_library'),
                               'starSED', 'kurucz')

        dtype = np.dtype([('wavelen', float), ('flambda', float)])

        sed_name_list = os.listdir(sed_dir)
        msg = ('An SED loaded from the cache is not '
               'identical to the same SED loaded from disk')
        for ix in range(5):
            full_name = os.path.join(sed_dir, sed_name_list[ix])
            from_np = np.genfromtxt(full_name, dtype=dtype)
            ss_uncache = Sed()
            ss_uncache.readSED_flambda(full_name)
            ss_cache  = Sed()
            ss_cache.readSED_flambda(full_name)

            self.assertEqual(ss_cache, ss_uncache, msg=msg)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()

