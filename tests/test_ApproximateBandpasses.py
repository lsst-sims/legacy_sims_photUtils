"""
Module to test that Approximate Bandpasses are in sync with the official LSST
Bandpasses from SYSENG

Note: While the LSST bandpasses list throughput values corresponding to
wavelengths in the range of 300.0-1150.0 nm, the `approximate_baseline`
directory of throughputs is created by a script manually. It is thus possible
for this directory to fall out of sync with the SYSENG values in `baseline`.
This module is intended to test whether this is happening.
"""
from __future__ import with_statement
import os
import unittest
import numpy as np
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.photUtils import BandpassDict

def setup_module(module):
    lsst.utils.tests.init()

class ApproximateBandPassTest(unittest.TestCase):
    """
    Tests for the approximate Bandpasses in the throughputs directory
    """
    
    def setUp(self):
        """
        setup before tests
        """
        throughputsDir = getPackageDir('throughputs')
        self.approxBandPassDir = os.path.join(throughputsDir, 'approximate_baseline')
        self.refBandPassDir = os.path.join(throughputsDir, 'baseline')
        self.refBandPassDict = BandpassDict.loadTotalBandpassesFromFiles()
        self.approxBandPassDict = \
            BandpassDict.loadTotalBandpassesFromFiles(bandpassDir=self.approxBandPassDir)
        self.errorMsg = "The failure of this test indicates that the"
        " approximate bandpasses in the lsst throughputs directory do not"
        "sync up with the baseline bandpasses is throughputs. This may require running"
        " the script : throughputs.approximate_baseline/approximateBandpasses.py"

    def test_BandPassIntegrals(self):
        """
        Test that the ratio of the quantity
        \int d\lambda T(\lambda) = band flux for a SED proportional to $\lambda$
        for the approximate bandpasses to the SYSENG band passes is 1.0 to an
        absolute tolerance hard coded to be 1.0e-14

        """
        for bn in 'ugrizy':
            refBandPass = self.refBandPassDict[bn]
            approxBandPass = self.approxBandPassDict[bn]
            refStep = np.diff(refBandPass.wavelen)
            approxStep = np.diff(approxBandPass.wavelen)

            # Currently we have uniform sampling, but the end points have
            # very slightly different steps. This accounts for 3 possible values
            # If there are more, then the steps are non-uniform
            if len(np.unique(approxStep)) > 3 :
                raise ValueError('The step sizes in {} seem to be unequal', 'Approximate Baseline')
            if len(np.unique(refStep)) > 3 :
                raise ValueError('The step sizes in {} seem to be unequal', 'Baseline')

            ratio = approxStep[1] * approxBandPass.wavelen.sum() / refStep[1]/ refBandPass.wavelen.sum()
            self.assertAlmostEqual(ratio, 1.0, delta=1.0e-14, msg=self.errorMsg)

        def tearDown(self):
            pass

class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
