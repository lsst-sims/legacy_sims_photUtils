#this will run a bunch of fake galaxies through get_EBV so that
#we can determine how long it is taking

import numpy, os
from lsst.sims.catalogs.generation.utils import makeStarTestDB, myTestStars
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.photUtils import EBVmixin

class myCatalogClass(InstanceCatalog,EBVmixin):
    column_outputs = ['id','raJ2000','decJ2000','EBV']

if os.path.exists('testEBVdatabase.db'):
    os.unlink('testEBVdatabase.db')
makeStarTestDB(filename='testEBVdatabase.db',size=10000)



