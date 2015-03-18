import numpy

import os
import unittest
import eups
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.db import ObservationMetaData, CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.coordUtils import AstrometryStars
from lsst.sims.photUtils import PhotometryStars, setupPhotometryCatalog
from lsst.sims.photUtils.utils import makeStarDatabase

class testCatalog(InstanceCatalog, AstrometryStars, PhotometryStars):
    column_outputs = ['raObserved', 'decObserved']
    default_formats = {'f':'%.12e'}

class baselineCatalog(InstanceCatalog, AstrometryStars, PhotometryStars):
    column_outputs = ['raObserved', 'decObserved',
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y']
    default_formats = {'f':'%.12e'}

class testDBObject(CatalogDBObject):
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

class InstanceCatalogSetupUnittest(unittest.TestCase):

    def setUp(self):
        self.dbName = 'setupTestStars.db'
        if os.path.exists(self.dbName):
            os.unlink(self.dbName)

        self.unrefractedRA = 50.0
        self.unrefractedDec = -5.0
        self.radius = 1.0
        obs_metadata = makeStarDatabase(filename=self.dbName, size=100,
                                        unrefractedRA=self.unrefractedRA,
                                        unrefractedDec=self.unrefractedDec,
                                        radius=self.radius)
        self.dbObj = testDBObject(address='sqlite:///' + self.dbName)

        self.obs_metadata = ObservationMetaData(unrefractedRA=self.unrefractedRA,
                                                unrefractedDec=self.unrefractedDec,
                                                boundType='circle', boundLength=self.radius,
                                                bandpassName='g', mjd=57000.0)

        self.obs_metadata_compound = ObservationMetaData(unrefractedRA=self.unrefractedRA,
                                                         unrefractedDec=self.unrefractedDec,
                                                         boundType='circle', boundLength=self.radius,
                                                         bandpassName=['g','i'], mjd=57000.0)

    def tearDown(self):
        if os.path.exists(self.dbName):
            os.unlink(self.dbName)

        del self.dbObj
        del self.dbName
        del self.unrefractedRA
        del self.unrefractedDec
        del self.radius
        del self.obs_metadata

    def testExceptions(self):

        class dummyClass(object):
            def __init__(self):
                pass

        xx = dummyClass()
        self.assertRaises(RuntimeError, setupPhotometryCatalog, obs_metadata=xx,
                          dbConnection=self.dbObj, catalogClass=testCatalog)

        self.assertRaises(RuntimeError, setupPhotometryCatalog, obs_metadata=self.obs_metadata,
                          dbConnection=xx, catalogClass=testCatalog)

        self.assertRaises(RuntimeError, setupPhotometryCatalog, obs_metadata=self.obs_metadata,
                          dbConnection=self.dbObj, catalogClass=dummyClass)


    def testSetupPhotometry(self):

        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata, dbConnection=self.dbObj,
                                     catalogClass=testCatalog)

        self.assertTrue('lsst_g' in cat.iter_column_names())
        self.assertFalse('lsst_u' in cat.iter_column_names())
        self.assertFalse('lsst_r' in cat.iter_column_names())
        self.assertFalse('lsst_i' in cat.iter_column_names())
        self.assertFalse('lsst_z' in cat.iter_column_names())
        self.assertFalse('lsst_y' in cat.iter_column_names())

        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                     dbConnection=self.dbObj, catalogClass=testCatalog)

        self.assertTrue('lsst_g' in cat.iter_column_names())
        self.assertTrue('lsst_i' in cat.iter_column_names())
        self.assertFalse('lsst_u' in cat.iter_column_names())
        self.assertFalse('lsst_r' in cat.iter_column_names())
        self.assertFalse('lsst_z' in cat.iter_column_names())
        self.assertFalse('lsst_y' in cat.iter_column_names())

        cat = testCatalog(self.dbObj, obs_metadata=self.obs_metadata)

        self.assertFalse('lsst_u' in cat.iter_column_names())
        self.assertFalse('lsst_g' in cat.iter_column_names())
        self.assertFalse('lsst_r' in cat.iter_column_names())
        self.assertFalse('lsst_i' in cat.iter_column_names())
        self.assertFalse('lsst_z' in cat.iter_column_names())
        self.assertFalse('lsst_y' in cat.iter_column_names())



    def testActualCatalog(self):


        testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata,
                                         dbConnection=self.dbObj,
                                         catalogClass=testCatalog)

        baselineCat = baselineCatalog(self.dbObj, obs_metadata=self.obs_metadata)

        testdtype = numpy.dtype([('raObserved', numpy.float), ('decObserved', numpy.float),
                                 ('lsst_g', numpy.float)])

        basedtype = numpy.dtype([('raObserved', numpy.float), ('decObserved', numpy.float),
                                 ('lsst_u', numpy.float), ('lsst_g', numpy.float),
                                 ('lsst_r', numpy.float), ('lsst_i', numpy.float),
                                 ('lsst_z', numpy.float), ('lsst_y', numpy.float)])

        testName = 'testCat.txt'
        baseName = 'baseCat.txt'
        testCat.write_catalog(testName)
        baselineCat.write_catalog(baseName)

        testData = numpy.genfromtxt(testName, dtype=testdtype, delimiter=',')
        baseData = numpy.genfromtxt(baseName, dtype=basedtype, delimiter=',')

        ct = 0
        for b, t in zip(baseData, testData):
            self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12)
            ct +=1

        self.assertTrue(ct>0)

        testdtype = numpy.dtype([('raObserved', numpy.float), ('decObserved', numpy.float),
                                 ('lsst_g', numpy.float), ('lsst_i', numpy.float)])

        testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                         dbConnection=self.dbObj,
                                         catalogClass=testCatalog)
        testCat.write_catalog(testName)
        testData = numpy.genfromtxt(testName, dtype=testdtype, delimiter=',')
        ct = 0
        for b, t in zip(baseData, testData):
            self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12)
            self.assertAlmostEqual(b['lsst_i'], t['lsst_i'], 12)
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
