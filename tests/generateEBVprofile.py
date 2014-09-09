#this will run a bunch of fake galaxies through get_EBV so that
#we can determine how long it is taking

import numpy, os
from lsst.sims.catalogs.generation.utils import makeStarTestDB
from lsst.sims.catalogs.generation.db import DBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.photUtils import EBVmixin
from lsst.sims.coordUtils import AstrometryStars

class myCatalogClass(InstanceCatalog,EBVmixin,AstrometryStars):
    column_outputs = ['id','raJ2000','decJ2000','EBV']
    catalog_type = 'myCatalog'

class myTestStars(DBObject):
    objid = 'teststars'
    tableid = 'stars'
    idColKey = 'id'
    #Make this implausibly large?  
    appendint = 1023
    dbAddress = 'sqlite:///testEBVdatabase.db'
    raColName = 'ra'
    decColName = 'decl'
    columns = [('id', None, int),
               ('raJ2000', 'ra*%f'%(numpy.pi/180.)),
               ('decJ2000', 'decl*%f'%(numpy.pi/180.)),
               ('umag', None),
               ('gmag', None),
               ('rmag', None),
               ('imag', None),
               ('zmag', None),
               ('ymag', None),
               ('magNorm', 'mag_norm', float)]


#if os.path.exists('testEBVdatabase.db'):
#    os.unlink('testEBVdatabase.db')
#makeStarTestDB(filename='testEBVdatabase.db',size=10000)

myDB=DBObject.from_objid('teststars')
myCatalog=myDB.getCatalog('myCatalog')
myCatalog.write_catalog('test_catalog.sav')


