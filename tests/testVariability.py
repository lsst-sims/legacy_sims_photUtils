from __future__ import with_statement
import os
import unittest
import lsst.utils.tests as utilsTests
import numpy
import sqlalchemy
import sqlite3
import json
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.photUtils.Photometry import PhotometryStars
from lsst.sims.photUtils.Variability import Variability


def makeMflareTable(size=10, **kwargs):
    """
    Make a test database to serve information to the mflareTest object
    """
    sedFiles = ['m2.0Full.dat','m5.1Full.dat', 'm4.9Full.dat']
    lcFiles = ['flare_lc_bin3_4.dat', 'flare_lc_bin1_4.dat', 'flare_lc_bin3_3.dat']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE mFlare
                     (varsimobjid int, varParamStr text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")
    
    for i in xrange(size):
        sedFile = sedFiles[numpy.random.randint(0,3)]
        varParam = {'varMethodName':'applyMflare', 
           'pars':{'t0':48000.0, 'lcfilename':lcFiles[numpy.random.randint(0,3)], 'dt':0.00069444418, 'length': 1825}}
        paramStr = json.dumps(varParam)
        
        qstr = '''INSERT INTO mFlare VALUES (%i, '%s', '%s')''' % (i, paramStr,sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

class mflareDB(CatalogDBObject):
    objid = 'mflareTest'
    tableid = 'mFlare'
    dbAddress = 'sqlite:///VariabilityTestDatabase.db'
    idColKey = 'varsimobjid'
    columns = [('id','varsimobjid',int),
               ('sedFilename','sedfilename',str,40)]

class mflareCatalog(InstanceCatalog,PhotometryStars,Variability):
    catalog_type = 'mflareCatalog'
    column_outputs = ['varsimobjid','sedFilename','lsst_u_var']
    
    default_columns=[('magNorm',14.0,float)]
    

class VariabilityTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if os.path.exists('VariabilityTestDatabase.db'):
            os.unlink('VariabilityTestDatabase.db')

    @classmethod
    def tearDownClass(cls):
        if os.path.exists('VariabilityTestDatabase.db'):
            os.unlink('VariabilityTestDatabase.db')

    def setUp(self):
        self.obs_metadata = ObservationMetaData(mjd=52000.0)

    def tearDown(self):
        del self.obs_metadata

    def testMflares(self):
        makeMflareTable()
        myDB = CatalogDBObject.from_objid('mflareTest')
        myCatalog = myDB.getCatalog('mflareCatalog',obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('mFlareTestCatalog.dat',chunk_size=1)

        if os.path.exists('mFlareTestCatalog.dat'):
            os.unlink('mFlareTestCatalog.dat')

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(VariabilityTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)

