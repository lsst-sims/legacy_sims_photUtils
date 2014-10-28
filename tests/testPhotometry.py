import numpy

import os
import unittest
import lsst.utils.tests as utilsTests

from lsst.sims.catalogs.measures.instance import InstanceCatalog, register_method, register_class, compound
from lsst.sims.catalogs.generation.db import ObservationMetaData
from lsst.sims.coordUtils import AstrometryStars, AstrometryGalaxies
from lsst.sims.catalogs.generation.utils import myTestGals, myTestStars, \
                                                makeStarTestDB, makeGalTestDB, getOneChunk
                                                
from lsst.sims.photUtils.Photometry import PhotometryGalaxies, PhotometryStars

from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.EBV import EBVbase, EBVmixin

from lsst.sims.photUtils.Variability import Variability

@register_class
class MyVariability(Variability):
    @register_method('testVar')
    def applySineVar(self, varParams, expmjd):
        period = varParams['period']
        amplitude = varParams['amplitude']
        phase = expmjd%period
        magoff = amplitude*numpy.sin(2*numpy.pi*phase)
        return {'u':magoff, 'g':magoff, 'r':magoff, 'i':magoff, 'z':magoff, 'y':magoff}

class testDefaults(object):
    """
    This class just provides default values for quantities that
    the astrometry mixins require in order to run
    """
    
    def get_proper_motion_ra(self):
        ra=self.column_by_name('raJ2000')
        out=numpy.zeros(len(ra))
        for i in range(len(ra)):
            out[i]=0.0
        
        return out
  
    
    def get_proper_motion_dec(self):
        ra=self.column_by_name('raJ2000')
        out=numpy.zeros(len(ra))
        for i in range(len(ra)):
            out[i]=0.0
        
        return out
    
    def get_parallax(self):
        ra=self.column_by_name('raJ2000')
        out=numpy.zeros(len(ra))
        for i in range(len(ra)):
            out[i]=1.2
        
        return out
    
    def get_radial_velocity(self):
        ra=self.column_by_name('raJ2000')
        out=numpy.zeros(len(ra))
        for i in range(len(ra)):
            out[i]=0.0
        
        return out

class cartoonPhotometryStars(PhotometryStars):
    """
    This is a class to support loading cartoon bandpasses into photometry so that we can be sure
    that the photometry mixin is loading the right files and calculating the right magnitudes.
    """

    @compound('cartoon_u','cartoon_g','cartoon_r','cartoon_i','cartoon_z')
    def get_magnitudes(self):
        """
        Example photometry getter for alternative (i.e. non-LSST) bandpasses
        """
        
        idNames = self.column_by_name('id')
        bandPassNames=['u','g','r','i','z']
        bandPassDir=os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'
        
        if self.bandPassList is None or self.phiArray is None:
            self.loadBandPassesFromFiles(bandPassNames,bandPassDir = bandPassDir, 
                    bandPassRoot = 'test_bandpass_')
                    
            self.setupPhiArray_dict()
        
        output = self.meta_magnitudes_getter(idNames)
        
        #############################################################################
        #Everything below this comment exists solely for the purposes of the unit test
        #if you need to write a customized getter for photometry that uses non-LSST
        #bandpasses, you only need to emulate the code above this comment.
       
       
        magNormList = self.column_by_name('magNorm')
        sedNames = self.column_by_name('sedFilename')
        
        #the two variables below will allow us to get at the SED and magnitude
        #data from within the unit test class, so that we can be sure
        #that the mixin loaded the correct bandPasses
        sublist = self.loadSeds(sedNames,magNorm = magNormList)
        for ss in sublist:
            self.sedMasterList.append(ss)
        
        if len(output) > 0:
            for i in range(len(output[0])):
                subList = []
                for j in range(len(output)):
                    subList.append(output[j][i])
            
                self.magnitudeMasterList.append(subList)
       
        return output


class cartoonPhotometryGalaxies(PhotometryGalaxies):
    """
    This is a class to support loading cartoon bandpasses into photometry so that we can be sure
    that the photometry mixin is loading the right files and calculating the right magnitudes.
    """

    @compound('ctotal_u','ctotal_g','ctotal_r','ctotal_i','ctotal_z',
              'cbulge_u','cbulge_g','cbulge_r','cbulge_i','cbulge_z',
              'cdisk_u','cdisk_g','cdisk_r','cdisk_i','cdisk_z',
              'cagn_u','cagn_g','cagn_r','cagn_i','cagn_z')
    def get_magnitudes(self):
        """
        getter for photometry of galaxies using non-LSST bandpasses
        """
        
        idNames = self.column_by_name('galid')
        bandPassNames=['u','g','r','i','z']
        bandPassDir=os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'
        
        if self.bandPassList is None or self.phiArray is None:
            self.loadBandPassesFromFiles(bandPassNames,bandPassDir = bandPassDir, 
                      bandPassRoot = 'test_bandpass_')
                      
            self.setupPhiArray_dict()
        
        output = self.meta_magnitudes_getter(idNames)
        
        ##########################################################################
        #Everything below this comment exists only for the purposes of the unittest.
        #If you need to write your own customized getter for photometry using 
        #non-LSST bandpasses, you only need to emulate the code above this comment
        
        if len(output) > 0:
            for i in range(len(output[0])):
                j = 5
                subList = []
                while j < 10:
                    subList.append(output[j][i])
                    j += 1
                self.magnitudeMasterDict['Bulge'].append(subList)
                
                subList = []
                while j < 15:
                    subList.append(output[j][i])
                    j += 1
                self.magnitudeMasterDict['Disk'].append(subList)
                
                subList = []
                while j < 20:
                    subList.append(output[j][i])
                    j += 1
                self.magnitudeMasterDict['Agn'].append(subList)
                
                
        
        
        componentNames = ['Bulge','Disk','Agn']
        
        for cc in componentNames:
            magName = "magNorm" + cc
            magNormList = self.column_by_name(magName)
            sName = "sedFilename" + cc
            sedNames = self.column_by_name(sName)
            
            if cc == 'Bulge' or cc == 'Disk':
                AvName = "internalAv"+cc
                Av = self.column_by_name(AvName)
            else:
                Av = None
        
            
            redshift = self.column_by_name("redshift")
            
            sublist = self.loadSeds(sedNames, magNorm = magNormList)
            self.applyAvAndRedshift(sublist, internalAv = Av, redshift = redshift)
            
            for ss in sublist:
                self.sedMasterDict[cc].append(ss)
       
        return output




class testCatalog(InstanceCatalog,AstrometryStars,Variability,testDefaults):
    catalog_type = 'MISC'
    default_columns=[('expmjd',5000.0,float)]
    
    def db_required_columns(self):
        return ['raJ2000'],['varParamStr']


class cartoonStars(InstanceCatalog,AstrometryStars,EBVmixin,Variability,cartoonPhotometryStars,testDefaults):
    catalog_type = 'cartoonStars'
    column_outputs=['id','raObserved','decObserved','magNorm',\
    'cartoon_u','cartoon_g','cartoon_r','cartoon_i','cartoon_z']
    
    #the lists below will contain the SED objects and the magnitudes
    #in a form that unittest can access and validate
    
    sedMasterList = []
    magnitudeMasterList = []
    
    #I need to give it the name of an actual SED file that spans the expected wavelength range
    defSedName = 'km30_5250.fits_g00_5370'
    default_columns = [('sedFilename', defSedName, (str,len(defSedName))), ('glon', 180., float), 
                       ('glat', 30., float)]

   

class cartoonGalaxies(InstanceCatalog,AstrometryGalaxies,EBVmixin,Variability,cartoonPhotometryGalaxies,testDefaults):
    catalog_type = 'cartoonGalaxies'
    column_outputs=['galid','raObserved','decObserved',\
    'ctotal_u','ctotal_g','ctotal_r','ctotal_i','ctotal_z']
    
    #I need to give it the name of an actual SED file that spans the expected wavelength range
    defSedName = "Inst.80E09.25Z.spec"
    default_columns = [('sedFilename', defSedName, (str, len(defSedName))) ,
                       ('sedFilenameAgn', defSedName, (str, len(defSedName))),
                       ('sedFilenameBulge', defSedName, (str, len(defSedName))),
                       ('sedFilenameDisk', defSedName, (str, len(defSedName))),
                       ('glon', 210., float),
                       ('glat', 70., float),
                       ('internalAvBulge',3.1,float),
                       ('internalAvDisk',3.1,float)]

    def get_galid(self):
        return self.column_by_name('id')
    
    #the dicts below will contain the SED objects and the magnitudes
    #in a form that unittest can access and validate
    
    sedMasterDict = {}
    sedMasterDict["Bulge"] = []
    sedMasterDict["Disk"] = []
    sedMasterDict["Agn"] = []
    
    magnitudeMasterDict = {}
    magnitudeMasterDict["Bulge"] = []
    magnitudeMasterDict["Disk"] = []
    magnitudeMasterDict["Agn"] = []
 
        
class testStars(InstanceCatalog, EBVmixin,MyVariability,PhotometryStars,testDefaults):
    catalog_type = 'test_stars'
    column_outputs=['id','raJ2000','decJ2000','magNorm',\
    'stellar_magNorm_var', \
    'lsst_u','sigma_lsst_u','lsst_u_var','sigma_lsst_u_var',
    'lsst_g','sigma_lsst_g','lsst_g_var','sigma_lsst_g_var',\
    'lsst_r','sigma_lsst_r','lsst_r_var','sigma_lsst_r_var',\
    'lsst_i','sigma_lsst_i','lsst_i_var','sigma_lsst_i_var',\
    'lsst_z','sigma_lsst_z','lsst_z_var','sigma_lsst_z_var',\
    'lsst_y','sigma_lsst_y','lsst_y_var','sigma_lsst_y_var',\
    'EBV','varParamStr']
    defSedName = 'sed_flat.txt'
    default_columns = [('sedFilename', defSedName, (str,len(defSedName))), ('glon', 180., float), 
                       ('glat', 30., float)]

class testGalaxies(InstanceCatalog,EBVmixin,MyVariability,PhotometryGalaxies,testDefaults):
    catalog_type = 'test_galaxies'
    column_outputs=['galid','raJ2000','decJ2000',\
        'redshift',
        'magNorm_Recalc_var', 'magNormAgn', 'magNormBulge', 'magNormDisk', \
        'uRecalc', 'sigma_uRecalc', 'uRecalc_var','sigma_uRecalc_var',\
        'gRecalc', 'sigma_gRecalc', 'gRecalc_var','sigma_gRecalc_var',\
        'rRecalc', 'sigma_rRecalc', 'rRecalc_var', 'sigma_rRecalc_var',\
         'iRecalc', 'sigma_iRecalc', 'iRecalc_var','sigma_iRecalc_var',\
         'zRecalc', 'sigma_zRecalc', 'zRecalc_var', 'sigma_zRecalc_var',\
         'yRecalc', 'sigma_yRecalc', 'yRecalc_var', 'sigma_yRecalc_var',\
        'sedFilenameBulge','uBulge', 'sigma_uBulge', 'gBulge', 'sigma_gBulge', \
        'rBulge', 'sigma_rBulge', 'iBulge', 'sigma_iBulge', 'zBulge', 'sigma_zBulge',\
         'yBulge', 'sigma_yBulge', \
        'sedFilenameDisk','uDisk', 'sigma_uDisk', 'gDisk', 'sigma_gDisk', 'rDisk', 'sigma_rDisk', \
        'iDisk', 'sigma_iDisk', 'zDisk', 'sigma_zDisk', 'yDisk', 'sigma_yDisk', \
        'sedFilenameAgn',\
        'uAgn', 'sigma_uAgn', 'uAgn_var', 'sigma_uAgn_var',\
        'gAgn', 'sigma_gAgn', 'gAgn_var', 'sigma_gAgn_var',\
        'rAgn', 'sigma_rAgn', 'rAgn_var', 'sigma_rAgn_var',\
        'iAgn', 'sigma_iAgn', 'iAgn_var', 'sigma_iAgn_var',\
        'zAgn', 'sigma_zAgn', 'zAgn_var', 'sigma_zAgn_var',\
        'yAgn', 'sigma_yAgn', 'yAgn_var', 'sigma_yAgn_var', 'varParamStr']
    defSedName = "sed_flat.txt"
    default_columns = [('sedFilename', defSedName, (str, len(defSedName))) ,
                       ('sedFilenameAgn', defSedName, (str, len(defSedName))),
                       ('sedFilenameBulge', defSedName, (str, len(defSedName))),
                       ('sedFilenameDisk', defSedName, (str, len(defSedName))),
                       ('glon', 210., float),
                       ('glat', 70., float),
                      ]

    def get_internalAvDisk(self):
        return numpy.ones(len(self._current_chunk))*0.1

    def get_internalAvBulge(self):
        return numpy.ones(len(self._current_chunk))*0.1

    def get_galid(self):
        return self.column_by_name('id')

class variabilityUnitTest(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # Create test databases
        if os.path.exists('PhotometryTestDatabase.db'):
            print "deleting database"
            os.unlink('PhotometryTestDatabase.db')
        makeStarTestDB(size=100000, seedVal=1, filename='PhotometryTestDatabase.db')
        makeGalTestDB(size=100000, seedVal=1, filename='PhotometryTestDatabase.db')

    @classmethod
    def tearDownClass(cls):
        if os.path.exists('PhotometryTestDatabase.db'):
            os.unlink('PhotometryTestDatabase.db')
    
    def setUp(self):
        self.obs_metadata = ObservationMetaData(mjd=52000.7, bandpassName='i', 
                            boundType = 'circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                            boundLength=1.0,m5=dict(u=23.9, g=25.0, r=24.7, i=24.0, z=23.3, y=22.1))

        self.galaxy = myTestGals(address='sqlite:///PhotometryTestDatabase.db')
        self.star = myTestStars(address='sqlite:///PhotometryTestDatabase.db')

    def tearDown(self):
        del self.galaxy
        del self.star
        del self.obs_metadata

    def testGalaxyVariability(self):

        galcat = testGalaxies(self.galaxy, obs_metadata=self.obs_metadata)
        results = self.galaxy.query_columns(['varParamStr'], obs_metadata=self.obs_metadata, 
                                             constraint='VarParamStr is not NULL')
        result = getOneChunk(results)
        for row in result:
            mags=galcat.applyVariability(row['varParamStr'])

    def testStarVariability(self):
        starcat = testStars(self.star, obs_metadata=self.obs_metadata)
        results = self.star.query_columns(['varParamStr'], obs_metadata=self.obs_metadata, 
                                         constraint='VarParamStr is not NULL')
        result = getOneChunk(results)
        for row in result:
            mags=starcat.applyVariability(row['varParamStr'])

class photometryUnitTest(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # Create test databases
        if os.path.exists('PhotometryTestDatabase.db'):
            print "deleting database"
            os.unlink('PhotometryTestDatabase.db')
        makeStarTestDB(size=100000, seedVal=1, filename='PhotometryTestDatabase.db')
        makeGalTestDB(size=100000, seedVal=1, filename='PhotometryTestDatabase.db')

    @classmethod
    def tearDownClass(cls):
        if os.path.exists('PhotometryTestDatabase.db'):
            os.unlink('PhotometryTestDatabase.db')

    def setUp(self):
        self.obs_metadata = ObservationMetaData(mjd=52000.7, bandpassName='i', 
                            boundType='circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                            boundLength=1.0, m5 = 25.0)

        self.galaxy = myTestGals(address='sqlite:///PhotometryTestDatabase.db')
        self.star = myTestStars(address='sqlite:///PhotometryTestDatabase.db')

    def tearDown(self):
        del self.galaxy
        del self.star
        del self.obs_metadata

    def testStars(self):
        test_cat=testStars(self.star, obs_metadata=self.obs_metadata)
        test_cat.write_catalog("testStarsOutput.txt")
        results = self.star.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)


    def testGalaxies(self):
        test_cat=testGalaxies(self.galaxy, obs_metadata=self.obs_metadata)
        test_cat.write_catalog("testGalaxiesOutput.txt")
        results = self.galaxy.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)

    def testAlternateBandpassesStars(self):
        """
        This will test our ability to do photometry using non-LSST bandpasses.
        
        It will first calculate the magnitudes using the getters in cartoonPhotometryStars.
        
        It will then load the alternate bandpass files 'by hand' and re-calculate the magnitudes
        and make sure that the magnitude values agree.  This is guarding against the possibility
        that some default value did not change and the code actually ended up loading the
        LSST bandpasses.
        """
        
        obs_metadata_pointed=ObservationMetaData(mjd=2013.23,
                                                 boundType='circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                                                 boundLength=1.0)
        obs_metadata_pointed.metadata = {}
        obs_metadata_pointed.metadata['Opsim_filter'] = 'i'
        test_cat=cartoonStars(self.star,obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog("testStarsCartoon.txt")
        
        cartoonDir = os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'
        testBandPasses = {}
        keys = ['u','g','r','i','z']
        
        bplist = []

        for kk in keys:
            testBandPasses[kk] = Bandpass()
            testBandPasses[kk].readThroughput(os.path.join(cartoonDir,"test_bandpass_%s.dat" % kk))
            bplist.append(testBandPasses[kk])
        
        sedObj = Sed()
        phiArray, waveLenStep = sedObj.setupPhiArray(bplist)

        i = 0
        
        #since all of the SEDs in the cartoon database are the same, just test on the first
        #if we ever include more SEDs, this can be somethin like
        #for ss in test_cata.sedMasterList:
        #
        ss=test_cat.sedMasterList[0]
        ss.resampleSED(wavelen_match = bplist[0].wavelen)
        ss.flambdaTofnu()
        mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
        for j in range(len(mags)):
            self.assertAlmostEqual(mags[j],test_cat.magnitudeMasterList[i][j],10)
        i += 1
    
    def testAlternateBandpassesGalaxies(self):
        """
        the same as testAlternateBandpassesStars, but for galaxies
        """
        
        obs_metadata_pointed=ObservationMetaData(mjd=50000.0,
                               boundType='circle',unrefractedRA=0.0,unrefractedDec=0.0,
                               boundLength=10.0)
        obs_metadata_pointed.metadata = {}
        obs_metadata_pointed.metadata['Opsim_filter'] = 'i'
        test_cat=cartoonGalaxies(self.galaxy,obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog("testGalaxiesCartoon.txt")
        
        cartoonDir = os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'
        testBandPasses = {}
        keys = ['u','g','r','i','z']
        
        bplist = []

        for kk in keys:
            testBandPasses[kk] = Bandpass()
            testBandPasses[kk].readThroughput(os.path.join(cartoonDir,"test_bandpass_%s.dat" % kk))
            bplist.append(testBandPasses[kk])
        
        sedObj = Sed()
        phiArray, waveLenStep = sedObj.setupPhiArray(bplist)
        
        components = ['Bulge', 'Disk', 'Agn']
        
        for cc in components:
            i = 0
            
            for ss in test_cat.sedMasterDict[cc]:
                if ss.wavelen != None:
                    ss.resampleSED(wavelen_match = bplist[0].wavelen)
                    ss.flambdaTofnu()
                    mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
                    for j in range(len(mags)):
                        self.assertAlmostEqual(mags[j],test_cat.magnitudeMasterDict[cc][i][j],10)
                i += 1
 
    
    def testEBV(self):
        
        ebvObject = EBVbase()
        ra = []
        dec = []
        gLat = []
        gLon = []
        for i in range(10):
            ra.append(i*2.0*numpy.pi/10.0)
            dec.append(i*numpy.pi/10.0)
            
            gLat.append(-0.5*numpy.pi+i*numpy.pi/10.0)
            gLon.append(i*2.0*numpy.pi/10.0)
            
            equatorialCoordinates=numpy.array([ra,dec])
            galacticCoordinates=numpy.array([gLon,gLat])
            
        ebvOutput = ebvObject.calculateEbv(equatorialCoordinates=equatorialCoordinates)
        self.assertEqual(len(ebvOutput),len(ra))
        
        ebvOutput = ebvObject.calculateEbv(galacticCoordinates=galacticCoordinates)
        self.assertEqual(len(ebvOutput),len(gLon))
        
        self.assertRaises(RuntimeError, ebvObject.calculateEbv, equatorialCoordinates=equatorialCoordinates, 
        galacticCoordinates=galacticCoordinates)
        self.assertRaises(RuntimeError, ebvObject.calculateEbv, equatorialCoordinates=None, galacticCoordinates=None)
        self.assertRaises(RuntimeError, ebvObject.calculateEbv)
     
def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(variabilityUnitTest)
    suites += unittest.makeSuite(photometryUnitTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
