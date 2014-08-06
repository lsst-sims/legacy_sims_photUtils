import numpy as np
import unittest
import warnings
import os
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.applyIGM import applyIGM

class TestApplyIGM(unittest.TestCase):

    def testApplyIGM(self):

        """Test application of IGM from Lookup Tables to SED objects"""

        #Test than a warning comes up if input redshift is out of range and that no changes occurs to SED
        testSed = Sed()
        testSed.readSED_flambda(os.environ['SIMS_SED_LIBRARY_DIR'] + '/galaxySED/Inst.80E09.25Z.spec.gz')
        testFlambda = []
        for fVal in testSed.flambda:
            testFlambda.append(fVal)
        testIGM = applyIGM()
        with warnings.catch_warnings(record=True) as wa:
            testIGM.applyIGM(1.1, testSed)
            self.assertEqual(len(wa), 1)
            self.assertTrue('IGM Lookup tables' in str(wa[-1].message))
        np.testing.assert_equal(testFlambda, testSed.flambda)

        #Test that lookup table is read in correctly
        testTable15 = np.genfromtxt(str(os.environ['SIMS_PHOTUTILS_DIR'] + '/python/lsst/sims/photUtils/' +
                                        'IGMLookupTables/MeanLookupTable_zSource1.5.tbl'))
        np.testing.assert_equal(testTable15, testIGM.meanLookups['1.5'])

        #Test output by making sure that an incoming sed with flambda = 1.0 everywhere will return the
        #transmission values of the lookup table as its flambda output
        testSed.setSED(testSed.wavelen, flambda = np.ones(len(testSed.wavelen)))
        testIGM.applyIGM(1.5, testSed)
        testTable15Above300 = testTable15[np.where(testTable15[:,0] >= 300.0)]
        testSed.resampleSED(wavelen_match = testTable15Above300[:,0])
        np.testing.assert_allclose(testTable15Above300[:,1], testSed.flambda, 1e-4)

if __name__ == "__main__":

    suite = unittest.TestLoader().loadTestsFromTestCase(TestApplyIGM)
    unittest.TextTestRunner(verbosity=2).run(suite)
        
