import numpy as np
import unittest
import warnings
import os
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.applyIGM import ApplyIGM

class TestApplyIGM(unittest.TestCase):

    def testInitializeIGM(self):

        "Test Initialization Method"

        #Make sure that if we initialize IGM with new inputs that it is initializing with them
        testIGM = ApplyIGM()
        testSed = Sed()
        testSed.readSED_flambda(os.environ['SIMS_SED_LIBRARY_DIR'] + '/galaxySED/Inst.80E09.25Z.spec.gz')
        testIGM.applyIGM(1.8, testSed)
        testZmin = 1.8
        testZmax = 2.2
        #Want new values for testing,
        #so make sure we are not just putting in the same values as are already there
        self.assertNotEqual(testZmin, testIGM.zMin)
        self.assertNotEqual(testZmax, testIGM.zMax)
        testIGM.initializeIGM(zMin = testZmin, zMax = testZmax)
        self.assertEqual(testZmin, testIGM.zMin)
        self.assertEqual(testZmax, testIGM.zMax)

    def testLoadTables(self):

        "Test Readin of IGM Lookup Tables"

        tableDirectory = str(os.environ['SIMS_SED_LIBRARY_DIR'] + '/igm')
        #First make sure that if variance option is turned on but there are no variance files that
        #the correct error is raised
        testIGM = ApplyIGM()
        testIGM.initializeIGM(zMax = 1.5)
        testMeanLookupTable = open('MeanLookupTable_zSource1.5.tbl', 'w')
        testMeanLookupTable.write('300.0        0.9999')
        testMeanLookupTable.close()
        self.assertRaisesRegexp(IOError, "Cannot find variance tables.", testIGM.loadTables(os.getcwd()))
        os.remove('MeanLookupTable_zSource1.5.tbl')

        #Then make sure that the mean lookup tables and var lookup tables all get loaded into proper dicts
        testIGMDicts = ApplyIGM()
        testIGMDicts.initializeIGM()
        testIGMDicts.loadTables(tableDirectory)
        redshiftValues = ['1.5', '1.6', '1.7', '1.8', '1.9', '2.0', '2.1', '2.2', '2.3', '2.4', '2.5',
                          '2.6', '2.7', '2.8', '2.9']
        self.assertItemsEqual(testIGMDicts.meanLookups.keys(), redshiftValues)
        self.assertItemsEqual(testIGMDicts.varLookups.keys(), redshiftValues)

        #Finally make sure that if Variance Boolean is false that nothing is passed in to varLookups
        testIGMVar = ApplyIGM()
        testIGMVar.initializeIGM()
        testIGMVar.loadTables(tableDirectory, varianceTbl = False)
        self.assertEqual(testIGMVar.varLookups, {})


    def testApplyIGM(self):

        """Test application of IGM from Lookup Tables to SED objects"""

        #Test that a warning comes up if input redshift is out of range and that no changes occurs to SED
        testSed = Sed()
        testSed.readSED_flambda(os.environ['SIMS_SED_LIBRARY_DIR'] + '/galaxySED/Inst.80E09.25Z.spec.gz')
        testFlambda = []
        for fVal in testSed.flambda:
            testFlambda.append(fVal)
        testIGM = ApplyIGM()
        testIGM.initializeIGM()
        with warnings.catch_warnings(record=True) as wa:
            testIGM.applyIGM(1.1, testSed)
            self.assertEqual(len(wa), 1)
            self.assertTrue('IGM Lookup tables' in str(wa[-1].message))
        np.testing.assert_equal(testFlambda, testSed.flambda)

        #Test that lookup table is read in correctly
        testTable15 = np.genfromtxt(str(os.environ['SIMS_SED_LIBRARY_DIR'] + '/igm/' +
                                        'MeanLookupTable_zSource1.5.tbl'))
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

    unittest.main()
