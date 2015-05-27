import numpy

import os
import unittest
import eups
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.db import ObservationMetaData, CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.coordUtils import AstrometryStars, AstrometryGalaxies
from lsst.sims.photUtils import PhotometryStars, PhotometryGalaxies, setupPhotometryCatalog
from lsst.sims.photUtils.utils import makeStarDatabase, makeGalaxyDatabase

class testStarCatalog(InstanceCatalog, AstrometryStars, PhotometryStars):
    """
    A class with no photometry columns.  Meant to be passed to setupPhotometryCatalog
    where it will be given photometry columns
    """
    column_outputs = ['raObserved', 'decObserved']
    default_formats = {'f':'%.12e'}

class baselineStarCatalog(InstanceCatalog, AstrometryStars, PhotometryStars):
    """
    Baseline photometry catalog against which to compare testStarCatalog
    """
    column_outputs = ['raObserved', 'decObserved',
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r', 'sigma_lsst_i',
                      'sigma_lsst_z', 'sigma_lsst_y']
    default_formats = {'f':'%.12e'}

class testGalaxyCatalog(InstanceCatalog, AstrometryGalaxies, PhotometryGalaxies):
    """
    A class with no photometry columns.  Meant to be passed to setupPhotometryCatalog
    where it will be given photometry columns
    """
    column_outputs = ['raObserved', 'decObserved']
    default_formats = {'f':'%.12e'}

class baselineGalaxyCatalog(InstanceCatalog, AstrometryGalaxies, PhotometryGalaxies):
    """
    Baseline photometry catalog against which to compare testGalaxyCatalog
    """
    column_outputs = ['raObserved', 'decObserved',
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r', 'sigma_lsst_i',
                      'sigma_lsst_z', 'sigma_lsst_y']
    default_formats = {'f':'%.12e'}

class testStarDBObject(CatalogDBObject):
    """
    CatalogDBObject to map our test database of stars
    """
    tableid = 'starsALL_forceseek'
    idColKey = 'id'
    raColName = 'ra'
    decColName = 'decl'
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('magNorm', None),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 40)]

class testGalaxyDBObject(CatalogDBObject):

    #: This is the base table for the galaxies
    #tableid = 'final_clone_db'
    tableid = 'galaxy'
    idColKey = 'galtileid'
    raColName = '((CAST(ra AS NUMERIC(9,6))%360.)+360.)%360.'
    decColName = 'dec'

    columns = [('galtileid', None, numpy.int64),
            ('galid', None, str, 30),
            ('raJ2000', 'ra*PI()/180.'),
            ('decJ2000', 'dec*PI()/180.'),
            ('raJ2000Bulge', 'bra*PI()/180.'),
            ('decJ2000Bulge', 'bdec*PI()/180.'),
            ('raJ2000Disk', 'dra*PI()/180.'),
            ('decJ2000Disk', 'ddec*PI()/180.'),
            ('raJ2000Agn', 'agnra*PI()/180.'),
            ('decJ2000Agn', 'agndec*PI()/180.'),
            ('magNormBulge', 'magnorm_bulge'),
            ('magNormDisk', 'magnorm_disk'),
            ('magNormAgn', 'magnorm_agn'),
            ('sedFilenameBulge', 'sedname_bulge', unicode, 40),
            ('sedFilenameDisk', 'sedname_disk', unicode, 40),
            ('sedFilenameAgn', 'sedname_agn', unicode, 40),
            ('majorAxisBulge', 'a_b*PI()/648000.'),
            ('minorAxisBulge', 'b_b*PI()/648000.'),
            ('positionAngleBulge', 'pa_bulge*PI()/180.'),
            ('sindexBulge', 'bulge_n', int),
            ('majorAxisDisk', 'a_d*PI()/648000.'),
            ('minorAxisDisk', 'b_d*PI()/648000.'),
            ('positionAngleDisk', 'pa_disk*PI()/180.'),
            ('sindexDisk', 'disk_n', int),
            ('internalExtinctionModelBulge', 'ext_model_b', str, 3),
            ('internalAvBulge', 'av_b'),
            ('internalRvBulge', 'rv_b'),
            ('internalExtinctionModelDisk', 'ext_model_d', str, 3),
            ('internalAvDisk', 'av_d'),
            ('internalRvDisk', 'rv_d'),
            ('lsst_u', 'u_ab'),
            ('lsst_g', 'g_ab'),
            ('lsst_r', 'r_ab'),
            ('lsst_i', 'i_ab'),
            ('lsst_z', 'z_ab'),
            ('lsst_y', 'y_ab')]

class InstanceCatalogSetupUnittest(unittest.TestCase):

    def setUp(self):
        self.driver = 'sqlite'
        self.StarDBName = 'setupTestStars.db'
        self.GalaxyDBName = 'setupTestGalaxies.db'
        if os.path.exists(self.StarDBName):
            os.unlink(self.StarDBName)

        if os.path.exists(self.GalaxyDBName):
            os.unlink(self.GalaxyDBName)

        self.unrefractedRA = 50.0
        self.unrefractedDec = -5.0
        self.radius = 1.0
        makeStarDatabase(filename=self.StarDBName, size=100,
                         unrefractedRA=self.unrefractedRA,
                         unrefractedDec=self.unrefractedDec,
                         radius=self.radius)


        makeGalaxyDatabase(filename=self.GalaxyDBName, size=100,
                           unrefractedRA=self.unrefractedRA,
                           unrefractedDec=self.unrefractedDec,
                           radius=self.radius)

        self.starDBObj = testStarDBObject(driver=self.driver, database= self.StarDBName)
        self.galaxyDBObj = testGalaxyDBObject(driver=self.driver, database=self.GalaxyDBName)

        self.obs_metadata = ObservationMetaData(unrefractedRA=self.unrefractedRA,
                                                unrefractedDec=self.unrefractedDec,
                                                boundType='circle', boundLength=self.radius,
                                                bandpassName='g', mjd=57000.0,
                                                m5=24.5)

        self.obs_metadata_compound = ObservationMetaData(unrefractedRA=self.unrefractedRA,
                                                         unrefractedDec=self.unrefractedDec,
                                                         boundType='circle', boundLength=self.radius,
                                                         bandpassName=['g','i'], mjd=57000.0,
                                                         m5=[24.5, 17.5])

    def tearDown(self):
        if os.path.exists(self.StarDBName):
            os.unlink(self.StarDBName)

        if os.path.exists(self.GalaxyDBName):
            os.unlink(self.GalaxyDBName)

        del self.starDBObj
        del self.galaxyDBObj
        del self.StarDBName
        del self.GalaxyDBName
        del self.unrefractedRA
        del self.unrefractedDec
        del self.radius
        del self.obs_metadata

    def testExceptions(self):
        """
        Make sure that setupPhotometryCatalog throws errors when it is supposed to
        """

        class dummyClass(object):
            def __init__(self):
                pass

        xx = dummyClass()
        self.assertRaises(RuntimeError, setupPhotometryCatalog, obs_metadata=xx,
                          dbConnection=self.starDBObj, catalogClass=testStarCatalog)

        self.assertRaises(RuntimeError, setupPhotometryCatalog, obs_metadata=self.obs_metadata,
                          dbConnection=xx, catalogClass=testStarCatalog)

        self.assertRaises(RuntimeError, setupPhotometryCatalog, obs_metadata=self.obs_metadata,
                          dbConnection=self.starDBObj, catalogClass=dummyClass)


    def testSetupPhotometry(self):
        """
        Make sure that catalogs instantiated by setupPhotometryCatalog contain the
        correct columns.
        """

        #test case with a single bandpass
        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata, dbConnection=self.starDBObj,
                                     catalogClass=testStarCatalog)

        self.assertTrue('lsst_g' in cat.iter_column_names())
        self.assertFalse('lsst_u' in cat.iter_column_names())
        self.assertFalse('lsst_r' in cat.iter_column_names())
        self.assertFalse('lsst_i' in cat.iter_column_names())
        self.assertFalse('lsst_z' in cat.iter_column_names())
        self.assertFalse('lsst_y' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_g' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_u' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_r' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_i' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_z' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_y' in cat.iter_column_names())

        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata, dbConnection=self.starDBObj,
                                     catalogClass=testStarCatalog, uncertainty=True)

        self.assertTrue('lsst_g' in cat.iter_column_names())
        self.assertFalse('lsst_u' in cat.iter_column_names())
        self.assertFalse('lsst_r' in cat.iter_column_names())
        self.assertFalse('lsst_i' in cat.iter_column_names())
        self.assertFalse('lsst_z' in cat.iter_column_names())
        self.assertFalse('lsst_y' in cat.iter_column_names())
        self.assertTrue('sigma_lsst_g' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_u' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_r' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_i' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_z' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_y' in cat.iter_column_names())

        #test case with two bandpasses
        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                     dbConnection=self.starDBObj, catalogClass=testStarCatalog)

        self.assertTrue('lsst_g' in cat.iter_column_names())
        self.assertTrue('lsst_i' in cat.iter_column_names())
        self.assertFalse('lsst_u' in cat.iter_column_names())
        self.assertFalse('lsst_r' in cat.iter_column_names())
        self.assertFalse('lsst_z' in cat.iter_column_names())
        self.assertFalse('lsst_y' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_g' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_u' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_r' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_i' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_z' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_y' in cat.iter_column_names())

        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                     dbConnection=self.starDBObj, catalogClass=testStarCatalog,
                                     uncertainty=True)

        self.assertTrue('lsst_g' in cat.iter_column_names())
        self.assertTrue('lsst_i' in cat.iter_column_names())
        self.assertFalse('lsst_u' in cat.iter_column_names())
        self.assertFalse('lsst_r' in cat.iter_column_names())
        self.assertFalse('lsst_z' in cat.iter_column_names())
        self.assertFalse('lsst_y' in cat.iter_column_names())
        self.assertTrue('sigma_lsst_g' in cat.iter_column_names())
        self.assertTrue('sigma_lsst_i' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_u' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_r' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_z' in cat.iter_column_names())
        self.assertFalse('sigma_lsst_y' in cat.iter_column_names())

        #make sure that class default columns did not get overwritten
        cat = testStarCatalog(self.starDBObj, obs_metadata=self.obs_metadata)

        self.assertFalse('lsst_u' in cat.iter_column_names())
        self.assertFalse('lsst_g' in cat.iter_column_names())
        self.assertFalse('lsst_r' in cat.iter_column_names())
        self.assertFalse('lsst_i' in cat.iter_column_names())
        self.assertFalse('lsst_z' in cat.iter_column_names())
        self.assertFalse('lsst_y' in cat.iter_column_names())



    def testActualCatalog(self):
        """
        Make sure that the values written to catalogs that are instantiated using
        setupPhotometryCatalog are correct
        """
        msgroot = ['failed on stars; ', 'failed on galaxies; ']

        testCatClasses = [testStarCatalog, testGalaxyCatalog]
        testCatDBs = [self.starDBObj, self.galaxyDBObj]
        baselineCats = []
        baselineCats.append(baselineStarCatalog(self.starDBObj, obs_metadata=self.obs_metadata))
        baselineCats.append(baselineGalaxyCatalog(self.galaxyDBObj, obs_metadata=self.obs_metadata))

        testName = 'testSetupCat.txt'
        baseName = 'baseSetupCat.txt'


        basedtype = numpy.dtype([('raObserved', numpy.float), ('decObserved', numpy.float),
                                 ('lsst_u', numpy.float), ('lsst_g', numpy.float),
                                 ('lsst_r', numpy.float), ('lsst_i', numpy.float),
                                 ('lsst_z', numpy.float), ('lsst_y', numpy.float),
                                 ('sigma_lsst_u', numpy.float), ('sigma_lsst_g',numpy.float),
                                 ('sigma_lsst_r', numpy.float), ('sigma_lsst_i', numpy.float),
                                 ('sigma_lsst_z', numpy.float), ('sigma_lsst_y', numpy.float)])

        for (testCatClass, dbo, baselineCat, msgr) in zip(testCatClasses, testCatDBs, baselineCats, msgroot):

            testdtype = numpy.dtype([('raObserved', numpy.float), ('decObserved', numpy.float),
                                 ('lsst_g', numpy.float)])


            testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata,
                                              dbConnection=dbo,
                                              catalogClass=testCatClass)

            testCat.write_catalog(testName)
            baselineCat.write_catalog(baseName)

            testData = numpy.genfromtxt(testName, dtype=testdtype, delimiter=',')
            baseData = numpy.genfromtxt(baseName, dtype=basedtype, delimiter=',')

            ct = 0
            for b, t in zip(baseData, testData):
                self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12,
                                       msg = '%s single column; %.12e != %.12e' % (msgr, b['lsst_g'], t['lsst_g']))
                ct +=1

            self.assertTrue(ct>0)

            testdtype = numpy.dtype([('raObserved', numpy.float), ('decObserved', numpy.float),
                                     ('lsst_g', numpy.float), ('lsst_i', numpy.float)])

            testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                             dbConnection=dbo,
                                             catalogClass=testCatClass)
            testCat.write_catalog(testName)
            testData = numpy.genfromtxt(testName, dtype=testdtype, delimiter=',')
            ct = 0
            for b, t in zip(baseData, testData):
                self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['lsst_g'], t['lsst_g']))
                self.assertAlmostEqual(b['lsst_i'], t['lsst_i'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['lsst_i'], t['lsst_i']))
                ct += 1

            self.assertTrue(ct>0)

            if os.path.exists(testName):
                os.unlink(testName)
            if os.path.exists(baseName):
                os.unlink(baseName)

    def testActualCatalogWithUncertainty(self):
        """
        Make sure that the values written to catalogs that are instantiated using
        setupPhotometryCatalog are correct (include photometric uncertainty)
        """

        msgroot = ['failed on stars; ', 'failed on galaxies; ']

        testCatClasses = [testStarCatalog, testGalaxyCatalog]
        testCatDBs = [self.starDBObj, self.galaxyDBObj]
        baselineCats = []

        #need to set up the baseline catalogs with the compound obs_metadata so that they get the
        #correct m5 values for both magnitudes (otherwise, they will use LSST defaults, which
        #disagree with our cartoon test case)
        baselineCats.append(baselineStarCatalog(self.starDBObj, obs_metadata=self.obs_metadata_compound))
        baselineCats.append(baselineGalaxyCatalog(self.galaxyDBObj, obs_metadata=self.obs_metadata_compound))

        testName = 'testSetupCatUncertainty.txt'
        baseName = 'baseSetupCatUncertainty.txt'

        basedtype = numpy.dtype([('raObserved', numpy.float), ('decObserved', numpy.float),
                                 ('lsst_u', numpy.float), ('lsst_g', numpy.float),
                                 ('lsst_r', numpy.float), ('lsst_i', numpy.float),
                                 ('lsst_z', numpy.float), ('lsst_y', numpy.float),
                                 ('sigma_lsst_u', numpy.float), ('sigma_lsst_g',numpy.float),
                                 ('sigma_lsst_r', numpy.float), ('sigma_lsst_i', numpy.float),
                                 ('sigma_lsst_z', numpy.float), ('sigma_lsst_y', numpy.float)])

        for (testCatClass, dbo, baselineCat, msgr) in zip(testCatClasses, testCatDBs, baselineCats, msgroot):

            testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata,
                                             dbConnection=dbo,
                                             catalogClass=testCatClass,
                                             uncertainty=True)

            testdtype = numpy.dtype([('raObserved', numpy.float), ('decObserved', numpy.float),
                                     ('lsst_g', numpy.float), ('sigma_lsst_g', numpy.float)])


            testCat.write_catalog(testName)
            baselineCat.write_catalog(baseName)

            testData = numpy.genfromtxt(testName, dtype=testdtype, delimiter=',')
            baseData = numpy.genfromtxt(baseName, dtype=basedtype, delimiter=',')

            ct = 0
            for b, t in zip(baseData, testData):
                self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12,
                                       msg = '%s single column; %.12e != %.12e ' % (msgr, b['lsst_g'], t['lsst_g']))
                self.assertAlmostEqual(b['sigma_lsst_g'], t['sigma_lsst_g'], 12,
                                       msg = '%s sigle column; %.12e != %.12e ' % (msgr, b['sigma_lsst_i'], t['sigma_lsst_g']))
                ct +=1

            self.assertTrue(ct>0)

            testdtype = numpy.dtype([('raObserved', numpy.float), ('decObserved', numpy.float),
                                     ('lsst_g', numpy.float), ('sigma_lsst_g', numpy.float),
                                     ('lsst_i', numpy.float), ('sigma_lsst_i', numpy.float)])

            testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                             dbConnection=dbo,
                                             catalogClass=testCatClass,
                                             uncertainty=True)
            testCat.write_catalog(testName)
            testData = numpy.genfromtxt(testName, dtype=testdtype, delimiter=',')
            ct = 0
            for b, t in zip(baseData, testData):
                self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['lsst_g'], t['lsst_g']))
                self.assertAlmostEqual(b['lsst_i'], t['lsst_i'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['lsst_i'], t['lsst_i']))
                self.assertAlmostEqual(b['sigma_lsst_g'], t['sigma_lsst_g'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['sigma_lsst_g'], t['lsst_g']))
                self.assertAlmostEqual(b['sigma_lsst_i'], t['sigma_lsst_i'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['sigma_lsst_i'], t['sigma_lsst_i']))
                ct +=1

            self.assertTrue(ct>0)

            if os.path.exists(testName):
                os.unlink(testName)
            if os.path.exists(baseName):
                os.unlink(baseName)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(InstanceCatalogSetupUnittest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
