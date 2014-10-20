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
    Make a test database to serve information to the flare test
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
           'pars':{'t0':48000.0, 'lcfilename':lcFiles[numpy.random.randint(0,len(lcFiles))], 'dt':0.00069444418, 'length': 1825}}
        paramStr = json.dumps(varParam)
        
        qstr = '''INSERT INTO mFlare VALUES (%i, '%s', '%s')''' % (i, paramStr,sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeRRlyTable(size=100, **kwargs):
    """
    Make a test database to serve information to the rrlyrae test
    """
    sedFiles = ['kp10_8750.fits_g35_8950','kp03_10500.fits_g45_10600','km50_6750.fits_g20_6750']
    lcFiles = ['rrly_lc/RRc/959802_per.txt','rrly_lc/RRc/1078860_per.txt','rrly_lc/RRab/98874_per.txt',
               'rrly_lc/RRab/3879827_per.txt']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE RRly
                     (varsimobjid int, varParamStr text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    numpy.random.seed(32)
    mjDisplacement = (numpy.random.sample(size)-50.0)*50.0
    for i in xrange(size):
        sedFile = sedFiles[numpy.random.randint(0,3)]
        varParam = {'varMethodName':'applyRRly', 
           'pars':{'tStartMjd':48000.0+mjDisplacement[i], 'filename':lcFiles[numpy.random.randint(0,len(lcFiles))]}}
        paramStr = json.dumps(varParam)
        
        qstr = '''INSERT INTO RRly VALUES (%i, '%s', '%s')''' % (i, paramStr,sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeCepheidTable(size=100, **kwargs):
    """
    Make a test database to serve information to the cepheid test
    """
    sedFiles = ['kp10_8750.fits_g35_8950','kp03_10500.fits_g45_10600','km50_6750.fits_g20_6750']
    lcFiles = ['cepheid_lc/classical_longPer_specfile', 'cepheid_lc/classical_medPer_specfile',
               'cepheid_lc/classical_shortPer_specfile', 'cepheid_lc/classical_shortPer_specfile',
               'cepheid_lc/popII_longPer_specfile', 'cepheid_lc/popII_shortPer_specfile']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE cepheid
                     (varsimobjid int, varParamStr text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    numpy.random.seed(32)
    periods = numpy.random.sample(size)*50.0
    mjDisplacement = (numpy.random.sample(size)-0.5)*50.0
    for i in xrange(size):
        sedFile = sedFiles[numpy.random.randint(0,3)]
        varParam = {'varMethodName':'applyCepheid', 
           'pars':{'period':periods[i], 'lcfile':lcFiles[i%len(lcFiles)], 't0':48000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)
        
        qstr = '''INSERT INTO cepheid VALUES (%i, '%s', '%s')''' % (i, paramStr,sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()


class variabilityDB(CatalogDBObject):
    dbAddress = 'sqlite:///VariabilityTestDatabase.db'
    idColKey = 'varsimobjid'
    columns = [('id','varsimobjid',int),
               ('sedFilename','sedfilename',str,40)]

class mflareDB(variabilityDB):
    objid = 'mflareTest'
    tableid = 'mFlare'

class rrlyDB(variabilityDB):
    objid = 'rrlyTest'
    tableid = 'RRly'

class cepheidDB(variabilityDB):
    objid = 'cepheidTest'
    tableid = 'cepheid'

class variabilityCatalog(InstanceCatalog,PhotometryStars,Variability):
    catalog_type = 'variabilityCatalog'
    column_outputs = ['varsimobjid','sedFilename','lsst_u_var','lsst_u']
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
    
    @unittest.skip("until final test; this is still very slow")
    def testMflares(self):
        makeMflareTable()
        myDB = CatalogDBObject.from_objid('mflareTest')
        myCatalog = myDB.getCatalog('variabilityCatalog',obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('mFlareTestCatalog.dat',chunk_size=1000)

        if os.path.exists('mFlareTestCatalog.dat'):
            os.unlink('mFlareTestCatalog.dat')

    def testRRlyrae(self):
        makeRRlyTable()
        myDB = CatalogDBObject.from_objid('rrlyTest')
        myCatalog = myDB.getCatalog('variabilityCatalog',obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('rrlyTestCatalog.dat',chunk_size=1000)

        if os.path.exists('rrlyTestCatalog.dat'):
            os.unlink('rrlyTestCatalog.dat')

    def testCepheids(self):
        makeCepheidTable()
        myDB = CatalogDBObject.from_objid('cepheidTest')
        myCatalog = myDB.getCatalog('variabilityCatalog',obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('cepheidTestCatalog.dat',chunk_size=1000)

        #if os.path.exists('cepheidTestCatalog.dat'):
        #    os.unlink('cepheidTestCatalog.dat')

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

