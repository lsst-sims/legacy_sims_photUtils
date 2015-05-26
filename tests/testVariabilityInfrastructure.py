from __future__ import with_statement
import os
import numpy
import unittest
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.utils import makeStarTestDB, myTestStars
from lsst.sims.catalogs.generation.utils import makeGalTestDB, myTestGals
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

    def get_sedFilename(self):

        star_seds = ['km20_5750.fits_g40_5790',
                     'kp10_9250.fits_g40_9250',
                     'bergeron_6500_85.dat_6700']

        ra = self.column_by_name('raJ2000')
        return numpy.array([star_seds[i%3] for i in range(len(ra))])

    @compound('test_u', 'test_g', 'test_r', 'test_i', 'test_z', 'test_y')
    def get_test_mags(self):
        objectID = self.column_by_name('id')

        columnNames = [name for name in self.get_test_mags._colnames]

        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandpassesFromFiles()

        indices = [ii for ii, name in enumerate(self.get_test_mags._colnames) \
                   if name in self._actually_calculated_columns]

        if len(indices) == 6:
            indices = None

        return self.meta_magnitudes_getter(objectID, columnNames, indices=indices)

class StellarVariabilityCatalogClass(StellarBaselineCatalogClass, FakeStellarVariabilityMixin):
    pass


class FakeGalaxyVariabilityMixin(object):

    @compound('delta_test_agn_u', 'delta_test_agn_g', 'delta_test_agn_r',
              'delta_test_bulge_u', 'delta_test_disk_i')
    def get_variability(self):
        ra = self.column_by_name('raJ2000')
        return numpy.array([1.0*numpy.ones(len(ra)),
                            2.0*numpy.ones(len(ra)),
                            3.0*numpy.ones(len(ra)),
                            4.0*numpy.ones(len(ra)),
                            5.0*numpy.ones(len(ra))])


class GalaxyBaselineCatalogClass(InstanceCatalog, PhotometryGalaxies):

    @compound('internalAvBulge', 'internalAvDisk')
    def get_internalAv(self):
        ra = self.column_by_name('raJ2000')
        return numpy.array([2.5*numpy.ones(len(ra)), 2.5*numpy.ones(len(ra))])

    @compound('sedFilenameBulge', 'sedFilenameDisk', 'sedFilenameAgn')
    def get_filenames(self):
        ra = self.column_by_name('raJ2000')

        galaxy_seds = ['Const.80E07.02Z.spec','Inst.80E07.002Z.spec','Burst.19E07.0005Z.spec']
        agn_sed = 'agn.spec'

        agnSeds = []
        for ii in range(len(ra)):
            agnSeds.append(agn_sed)

        bulgeSeds = [galaxy_seds[(ii+1)%3] for ii in range(len(ra))]
        diskSeds = [galaxy_seds[ii%3] for ii in range(len(ra))]

        return numpy.array([bulgeSeds, diskSeds, agnSeds])


    @compound('test_u', 'test_g', 'test_r', 'test_i', 'test_z', 'test_y',
              'test_bulge_u', 'test_bulge_g', 'test_bulge_r',
              'test_bulge_i', 'test_bulge_z', 'test_bulge_y',
              'test_disk_u', 'test_disk_g', 'test_disk_r',
              'test_disk_i', 'test_disk_z', 'test_disk_y',
              'test_agn_u', 'test_agn_g', 'test_agn_r',
              'test_agn_i', 'test_agn_z', 'test_agn_y')
    def get_test_mags(self):
        """
        Getter for test galaxy magnitudes

        """
        objectID = self.column_by_name('id')

        columnNames = [name for name in self.get_test_mags._colnames]

        """
        Here is where we need some code to load a list of bandpass objects
        into self.bandpassDict so that the bandpasses are available to the
        mixin.  Ideally, we would only do this once for the whole catalog
        """
        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandpassesFromFiles()

        indices = numpy.unique([ii % 6 for ii, name in enumerate(self.get_test_mags._colnames) \
                               if name in self._actually_calculated_columns])

        if len(indices)==6:
            indices=None

        return self.meta_magnitudes_getter(objectID, columnNames, indices=indices)


class GalaxyVariabilityCatalogClass(GalaxyBaselineCatalogClass, FakeGalaxyVariabilityMixin):
    pass

class VariabilityDesignTest(unittest.TestCase):
    """
    This unit test case will test that the general
    variability design was correclty implemented.
    It will not test an particular model of variability.
    It will merely test that, given a mean and delta magnitude,
    the InstanceCatalog class can correctly calculate the varying
    magnitude of an object.
    """

    @classmethod
    def setUpClass(cls):
        cls.starDbName = 'VariabilityInfrastructureTestStarDB.db'
        if os.path.exists(cls.starDbName):
            os.unlink(cls.starDbName)

        makeStarTestDB(filename=cls.starDbName)

        cls.galaxyDbName = 'VariabilityInfrastructureTestGalaxyDB.db'
        if os.path.exists(cls.galaxyDbName):
            os.unlink(cls.galaxyDbname)

        makeGalTestDB(filename=cls.galaxyDbName)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.starDbName):
            os.unlink(cls.starDbName)

        if os.path.exists(cls.galaxyDbName):
            os.unlink(cls.galaxyDbName)

    def setUp(self):
        self.starDB = myTestStars(address='sqlite:///'+self.starDbName)
        self.galaxyDB = myTestGals(address='sqlite:///'+self.galaxyDbName)

    def tearDown(self):
        del self.starDB
        del self.galaxyDB


    def testStellarVariabilityInfrastructure(self):
        """
        Test that the variability design was correctly implemented
        in the case of stars
        """

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


    def testGalaxyVariabilityInfrastructure(self):
        """
        Test that the variability design was correctly implemented in
        the case of galaxies
        """

        outputs = ['id',
                   'test_u', 'test_g', 'test_r', 'test_i', 'test_z', 'test_y',
                   'test_bulge_u', 'test_bulge_g', 'test_bulge_r', 'test_bulge_i',
                   'test_bulge_z', 'test_bulge_y',
                   'test_disk_u', 'test_disk_g', 'test_disk_r', 'test_disk_i',
                   'test_disk_z', 'test_disk_y',
                   'test_agn_u', 'test_agn_g', 'test_agn_r',
                   'test_agn_i', 'test_agn_z','test_agn_y']

        baseline = GalaxyBaselineCatalogClass(self.galaxyDB, column_outputs=outputs)
        variable = GalaxyVariabilityCatalogClass(self.galaxyDB, column_outputs=outputs)

        phot = PhotometryGalaxies()

        variable_indices = [19, 20, 21, 7, 16]

        for bb, vv in zip(baseline.iter_catalog(), variable.iter_catalog()):
            self.assertEqual(bb[0], vv[0])

            #test that the variable components are altered
            #the way they ought to be
            self.assertAlmostEqual(bb[19]+1.0, vv[19], 10)
            self.assertAlmostEqual(bb[20]+2.0, vv[20], 10)
            self.assertAlmostEqual(bb[21]+3.0, vv[21], 10)
            self.assertAlmostEqual(bb[7]+4.0, vv[7], 10)
            self.assertAlmostEqual(bb[16]+5.0, vv[16], 10)

            #test that the components which do not vary are equal
            for ix in range(7,25):
                if ix not in variable_indices:
                    self.assertAlmostEqual(bb[ix], vv[ix], 10)

            #test that the total magnitudes are correctly calculated
            for ix in range(6):

                self.assertAlmostEqual(bb[ix+1],
                                       phot.sum_magnitudes(bulge=bb[7+ix],
                                                           disk=bb[13+ix],
                                                           agn=bb[19+ix]),
                                       10)

                self.assertAlmostEqual(vv[ix+1],
                                       phot.sum_magnitudes(bulge=vv[7+ix],
                                                           disk=vv[13+ix],
                                                           agn=vv[19+ix]),
                                       10)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(VariabilityDesignTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
