from __future__ import with_statement
import os
import numpy
import unittest
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.utils import makeStarTestDB, myTestStars
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.photUtils import PhotometryStars, PhotometryGalaxies

class FakeStellarVariabilityMixin(object):

    @compound('delta_test_u', 'delta_test_g', 'delta_test_r')
    def get_variability(self):
        ra = self.column_by_name('raJ2000')
        return numpy.array([4.0*numpy.ones(len(ra)),
                            10.0*numpy.ones(len(ra)),
                            20.0*numpy.ones(len(ra))])


class StellarBaselineCatalogClass(InstanceCatalog, PhotometryStars):

    star_seds = ['km20_5750.fits_g40_5790',
                 'kp10_9250.fits_g40_9250',
                 'bergeron_6500_85.dat_6700']

    def get_sedFilename(self):
        ra = self.column_by_name('raJ2000')
        return numpy.array([self.star_seds[i%3] for i in range(len(ra))])

    @compound('test_u', 'test_g', 'test_r', 'test_i', 'test_z', 'test_y')
    def get_test_mags(self):
        objectID = self.column_by_name('id')

        columnNames = [name for name in self.get_test_mags._colnames]

        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandpassesFromFiles()

        indices = [ii for ii, name in enumerate(self.get_test_mags._colnames) \
                   if name in self.all_calculated_columns]

        if len(indices) == 6:
            indices = None

        return self.meta_magnitudes_getter(objectID, columnNames, indices=indices)

class StellarVariabilityCatalogClass(StellarBaselineCatalogClass, FakeStellarVariabilityMixin):
    pass

class VariabilitySetupTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.starDbName = 'VariabilityInfrastructureTestDB.db'
        if os.path.exists(cls.starDbName):
            os.unlink(cls.starDbName)

        makeStarTestDB(filename=cls.starDbName)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.starDbName):
            os.unlink(cls.starDbName)

    def setUp(self):
        self.starDB = myTestStars(address='sqlite:///'+self.starDbName)

    def tearDown(self):
        del self.starDB


    def testStellarVariabilityInfrastructure(self):

        outputs = ['id', 'test_u', 'test_g', 'test_r', 'test_i', 'test_z', 'test_y']
        baseline = StellarBaselineCatalogClass(self.starDB, column_outputs=outputs)
        variable = StellarVariabilityCatalogClass(self.starDB, column_outputs=outputs)

        for bb, vv in zip(baseline.iter_catalog(), variable.iter_catalog()):
            self.assertEqual(bb[0], vv[0])
            self.assertAlmostEqual(vv[1]-bb[1], 4.0, 10)
            self.assertAlmostEqual(vv[2]-bb[2], 10.0, 10)
            self.assertAlmostEqual(vv[3]-bb[3], 20.0, 10)
            self.assertAlmostEqual(vv[4], bb[4], 10)
            self.assertAlmostEqual(vv[5], bb[5], 10)
            self.assertAlmostEqual(vv[6], bb[6], 10)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(VariabilitySetupTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
