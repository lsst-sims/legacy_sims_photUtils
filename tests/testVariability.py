from __future__ import with_statement
import os
import unittest
import lsst.utils.tests as utilsTests
import numpy
import sqlalchemy
import sqlite3
import json
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.photUtils.Photometry import PhotometryStars, PhotometryGalaxies
from lsst.sims.photUtils.Variability import VariabilityStars, VariabilityGalaxies

def makeMflareTable(size=10, **kwargs):
    """
    Make a test database to serve information to the flare test
    """

    #a haphazard sample of mdwarf SEDs
    sedFiles = ['m2.0Full.dat', 'm5.1Full.dat', 'm4.9Full.dat']

    #a haphazard sample of mflare light curves
    lcFiles = ['flare_lc_bin3_4.dat', 'flare_lc_bin1_4.dat', 'flare_lc_bin3_3.dat']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE mFlare
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    for i in xrange(size):
        sedFile = sedFiles[numpy.random.randint(0,len(sedFiles))]
        varParam = {'varMethodName':'applyMflare',
           'pars':{'t0':48000.0, 'lcfilename':lcFiles[numpy.random.randint(0,len(lcFiles))], 'dt':0.00069444418, 'length': 1825}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO mFlare VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeRRlyTable(size=100, **kwargs):
    """
    Make a test database to serve information to the rrlyrae test
    """

    #a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    #a haphazard sample of RRLyrae light curves
    lcFiles = ['rrly_lc/RRc/959802_per.txt', 'rrly_lc/RRc/1078860_per.txt', 'rrly_lc/RRab/98874_per.txt',
               'rrly_lc/RRab/3879827_per.txt']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE RRly
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    numpy.random.seed(32)
    mjDisplacement = (numpy.random.sample(size)-50.0)*50.0
    for i in xrange(size):
        sedFile = sedFiles[numpy.random.randint(0,len(sedFiles))]
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

    #a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    #a haphazard sample of cepheid light curves
    lcFiles = ['cepheid_lc/classical_longPer_specfile', 'cepheid_lc/classical_medPer_specfile',
               'cepheid_lc/classical_shortPer_specfile', 'cepheid_lc/classical_shortPer_specfile',
               'cepheid_lc/popII_longPer_specfile', 'cepheid_lc/popII_shortPer_specfile']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE cepheid
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    numpy.random.seed(32)
    periods = numpy.random.sample(size)*50.0
    mjDisplacement = (numpy.random.sample(size)-0.5)*50.0
    for i in xrange(size):
        sedFile = sedFiles[numpy.random.randint(0,len(sedFiles))]
        varParam = {'varMethodName':'applyCepheid',
           'pars':{'period':periods[i], 'lcfile':lcFiles[numpy.random.randint(0,len(lcFiles))], 't0':48000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO cepheid VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeEbTable(size=100, **kwargs):
    """
    Make a test database to serve information to the Eb test
    """

    #a haphazard sample of eclipsing binary light curves
    lcFiles = ['eb_lc/EB.2294.inp', 'eb_lc/EB.1540.inp', 'eb_lc/EB.2801.inp']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE eb
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    numpy.random.seed(32)
    periods = numpy.random.sample(size)*50.0
    mjDisplacement = (numpy.random.sample(size)-0.5)*50.0
    for i in xrange(size):
        sedFile = 'sed_flat_norm.txt'
        varParam = {'varMethodName':'applyEb',
           'pars':{'period':periods[i], 'lcfile':lcFiles[numpy.random.randint(0,len(lcFiles))], 't0':48000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO eb VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeMicrolensingTable(size=100, **kwargs):
    """
    Make a test database to serve information to the microlensing test
    """

    #a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    #there are two microlensing methods; they should be equivalent
    method = ['applyMicrolensing', 'applyMicrolens']
    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE microlensing
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    numpy.random.seed(32)
    that = numpy.random.sample(size)*40.0+40.0
    umin = numpy.random.sample(size)
    mjDisplacement = numpy.random.sample(size)*50.0
    for i in xrange(size):
        sedFile = sedFiles[0]
        varParam = {'varMethodName':method[i%len(method)],
           'pars':{'that':that[i], 'umin':umin[i], 't0':52000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO microlensing VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeBHMicrolensingTable(size=100, **kwargs):
    """
    Make a test database to serve information to the BHmicrolensing test
    """

    #a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    #a sample of black hole microlensing light curves that do not repeat time steps
    #(repeating time steps causes the scipy spline interpolation routine to return Nan)
    lcFiles = ['microlens/bh_binary_source/lc_14_25_75_8000_0_0.05_316',
               'microlens/bh_binary_source/lc_14_25_4000_8000_0_phi1.09_0.005_100',
               'microlens/bh_binary_source/lc_14_25_75_8000_0_tets2.09_0.005_316']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE bhmicrolensing
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    numpy.random.seed(32)
    mjDisplacement = numpy.random.sample(size)*5.0*365.25
    for i in xrange(size):
        sedFile = sedFiles[numpy.random.randint(0,len(sedFiles))]
        varParam = {'varMethodName':'applyBHMicrolens',
           'pars':{'filename':lcFiles[numpy.random.randint(0,len(lcFiles))], 't0':52000.0-mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO bhmicrolensing VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeAmcvnTable(size=100, **kwargs):
    """
    Make a test database to serve information to the AMCVN test
    """

    #a haphazard sample of white dwarf SEDs
    sedFiles = ['bergeron_He_4750_70.dat_4950', 'bergeron_50000_85.dat_54000']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE amcvn
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    numpy.random.seed(32)
    doesBurst = numpy.random.randint(0,1,size=size)
    burst_freq = numpy.random.randint(10,150,size=size)
    burst_scale = 115
    amp_burst = numpy.random.sample(size)*8.0
    color_excess_during_burst = numpy.random.sample(size)*0.2-0.4
    amplitude = numpy.random.sample(size)*0.2
    period = numpy.random.sample(size)*200.0
    mjDisplacement = numpy.random.sample(size)*50.0
    for i in xrange(size):
        sedFile = sedFiles[numpy.random.randint(0,len(sedFiles))]
        varParam = {'varMethodName':'applyAmcvn',
           'pars':{'does_burst':int(doesBurst[i]), #have to cast to int from numpy.int for json
                   'burst_freq':int(burst_freq[i]),
                   'burst_scale':burst_scale,
                   'amp_burst':amp_burst[i],
                   'color_excess_during_burst':color_excess_during_burst[i],
                   'amplitude':amplitude[i],
                   'period':period[i],
                   't0':52000.0+mjDisplacement[i]}}

        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO amcvn VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeAgnTable(size=100, **kwargs):
    """
    Make a test database to serve information to the microlensing test
    """

    #a haphazard sample of galaxy SEDs
    sedFiles = ['Exp.31E06.0005Z.spec', 'Inst.79E06.1Z.spec', 'Const.50E07.0005Z.spec']
    method = ['applyAgn']
    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE agn
                     (galid int, varsimobjid int,
                      internalAvBulge real, internalAvDisk real, redshift real,
                      variability text,
                      sedFilenameBulge text, sedFilenameDisk text, sedFilenameAgn text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    numpy.random.seed(32)
    agn_tau = numpy.random.sample(size)*100.0+100.0
    agn_sfu = numpy.random.sample(size)*2.0
    agn_sfg = numpy.random.sample(size)*2.0
    agn_sfr = numpy.random.sample(size)*2.0
    agn_sfi = numpy.random.sample(size)*2.0
    agn_sfz = numpy.random.sample(size)*2.0
    agn_sfy = numpy.random.sample(size)*2.0
    mjDisplacement = numpy.random.sample(size)*5.0
    avBulge = numpy.random.sample(size)*0.5+2.6
    avDisk = numpy.random.sample(size)*0.5+2.6
    redshift = numpy.random.sample(size)*0.5
    for i in xrange(size):
        varParam = {'varMethodName':'applyAgn',
           'pars':{'agn_tau':agn_tau[i], 'agn_sfu':agn_sfu[i], 'agn_sfg':agn_sfg[i],
                    'agn_sfr':agn_sfr[i], 'agn_sfi':agn_sfi[i], 'agn_sfz':agn_sfz[i],
                    'agn_sfy':agn_sfy[i], 't0_mjd':48000.0+mjDisplacement[i],
                    'seed':numpy.random.randint(0,200000)}}

        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO agn VALUES (%i, %i, %f, %f, %f, '%s', '%s', '%s', '%s')''' % \
               (i, i, avBulge[i], avDisk[i], redshift[i],
               paramStr,
               sedFiles[numpy.random.randint(0,len(sedFiles))],
               sedFiles[numpy.random.randint(0,len(sedFiles))],
               'agn.spec')

        c.execute(qstr)
    conn.commit()
    conn.close()

class variabilityDB(CatalogDBObject):
    driver = 'sqlite'
    database = 'VariabilityTestDatabase.db'
    idColKey = 'varsimobjid'
    columns = [('id', 'varsimobjid', int),
               ('sedFilename', 'sedfilename', str, 40),
               ('varParamStr', 'variability', str, 600)]

class mflareDB(variabilityDB):
    objid = 'mflareTest'
    tableid = 'mFlare'

class rrlyDB(variabilityDB):
    objid = 'rrlyTest'
    tableid = 'RRly'

class cepheidDB(variabilityDB):
    objid = 'cepheidTest'
    tableid = 'cepheid'

class ebDB(variabilityDB):
    objid = 'ebTest'
    tableid = 'eb'

class microlensDB(variabilityDB):
    objid = 'microlensTest'
    tableid = 'microlensing'

class BHmicrolensDB(variabilityDB):
    objid = 'bhmicrolensTest'
    tableid = 'bhmicrolensing'

class amcvnDB(variabilityDB):
    objid = 'amcvnTest'
    tableid = 'amcvn'

class agnDB(variabilityDB):
    objid = 'agnTest'
    tableid = 'agn'

class StellarVariabilityCatalog(InstanceCatalog, PhotometryStars, VariabilityStars):
    catalog_type = 'stellarVariabilityCatalog'
    column_outputs = ['varsimobjid', 'sedFilename', 'delta_lsst_u']
    default_columns=[('magNorm', 14.0, float)]


class GalaxyVariabilityCatalog(InstanceCatalog, PhotometryGalaxies, VariabilityGalaxies):
    catalog_type = 'galaxyVariabilityCatalog'
    column_outputs = ['varsimobjid', 'sedFilenameAgn', 'lsstUdiff', 'delta_uAgn']
    default_columns=[('magNormAgn', 14.0, float),
                     ('magNormDisk', 14.0, float),
                     ('magNormBulge', 14.0, float)]

    def get_lsstUdiff(self):
        lsstUvar = self.column_by_name('lsst_u')

        bulge = self.column_by_name('uBulge')
        disk = self.column_by_name('uDisk')
        agn = self.column_by_name('uAgn') - self.column_by_name('delta_uAgn')
        lsstU = self.sum_magnitudes(bulge=bulge, disk=disk, agn=agn)

        return lsstUvar - lsstU

    def get_agnUdiff(self):
        lsstU = self.column_by_name('uAgn')
        lsstUvar = self.column_by_name('uAgn_var')
        return lsstUvar - lsstU

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
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('mFlareTestCatalog.dat', chunk_size=1000)

        if os.path.exists('mFlareTestCatalog.dat'):
            os.unlink('mFlareTestCatalog.dat')

    def testRRlyrae(self):
        makeRRlyTable()
        myDB = CatalogDBObject.from_objid('rrlyTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('rrlyTestCatalog.dat', chunk_size=1000)

        if os.path.exists('rrlyTestCatalog.dat'):
            os.unlink('rrlyTestCatalog.dat')

    def testCepheids(self):
        makeCepheidTable()
        myDB = CatalogDBObject.from_objid('cepheidTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('cepheidTestCatalog.dat', chunk_size=1000)

        if os.path.exists('cepheidTestCatalog.dat'):
            os.unlink('cepheidTestCatalog.dat')

    def testEb(self):
        makeEbTable()
        myDB = CatalogDBObject.from_objid('ebTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('ebTestCatalog.dat', chunk_size=1000)

        if os.path.exists('ebTestCatalog.dat'):
            os.unlink('ebTestCatalog.dat')

    def testMicrolensing(self):
        #Note: this test assumes that the parameters for the microlensing variability
        #model occur in a standard varParamStr column in the database.
        #Actually, the current database of microlensing events simply store the variability
        #parameters as independent columns in the database.
        #The varParamStr formalism is how the applyMicrolensing methods are written, however,
        #so that is what we will test.

        makeMicrolensingTable()
        myDB = CatalogDBObject.from_objid('microlensTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('microlensTestCatalog.dat', chunk_size=1000)

        if os.path.exists('microlensTestCatalog.dat'):
            os.unlink('microlensTestCatalog.dat')

    def testBHMicrolensing(self):
        #Note: this test assumes that the parameters for the BHmicrolensing variability
        #model occur in a standard varParamStr column in the database.
        #Actually, the current database of BHmicrolensing events simply store the variability
        #parameters as independent columns in the database.
        #The varParamStr formalism is how the applyBHMicrolens method is written, however,
        #so that is what we will test.

        makeBHMicrolensingTable()
        myDB = CatalogDBObject.from_objid('bhmicrolensTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('bhmicrolensTestCatalog.dat', chunk_size=1000)

        if os.path.exists('bhmicrolensTestCatalog.dat'):
            os.unlink('bhmicrolensTestCatalog.dat')

    def testAmcvn(self):
        #Note: this test assumes that the parameters for the Amcvn variability
        #model occur in a standard varParamStr column in the database.
        #Actually, the current database of Amcvn events simply store the variability
        #parameters as independent columns in the database.
        #The varParamStr formalism is how the applyAmcvn method is written, however,
        #so that is what we will test.

        makeAmcvnTable()
        myDB = CatalogDBObject.from_objid('amcvnTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('amcvnTestCatalog.dat', chunk_size=1000)

        if os.path.exists('amcvnTestCatalog.dat'):
            os.unlink('amcvnTestCatalog.dat')

    def testAgn(self):

        makeAgnTable()
        myDB = CatalogDBObject.from_objid('agnTest')
        myCatalog = myDB.getCatalog('galaxyVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog('agnTestCatalog.dat', chunk_size=1000)

        if os.path.exists('agnTestCatalog.dat'):
            os.unlink('agnTestCatalog.dat')

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

