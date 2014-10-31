"""
This file defines some test catalog and DBObject classes for use with unit tests.

To date (30 October 2014) testPhotometry.py and testCosmology.py import from this module
"""

import numpy
import os
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

__all__ = ["MyVariability", "testDefaults", "cartoonPhotometryStars",
           "cartoonPhotometryGalaxies", "testCatalog", "cartoonStars",
           "cartoonGalaxies", "testStars", "testGalaxies",
           "comovingDistanceIntegrand", "cosmologicalOmega"]

@register_class
class MyVariability(Variability):
    """
    This is a mixin which provides a dummy variability method for use in unit tests
    """
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

    In addition to creating a catalog, when the get_magnitude method below is called, it will
    add sedMasterList and magnitudeMasterList to the catalog instantiation.  These are lists
    containing the magnitudes output to the catalog and the SEDs used to calculate them.

    Having these variables allows the unittest to verify the output of the catalog
    (see testAlternateBandpassesStars in testPhotometry.py to see how this works)
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

    In addition to writing the catalog, when the get_magnitudes method below is called, the
    variables sedMasterDict and mangitudeMasterDict are added to the catalog instantiation.
    These store the magnitudes calculated for the catalog and the SEDs used to find them.

    This allows the unittest to verify the contents of the catalog
    (see testAlternateBandpassesGalaxies in testPhotometry.py to see how this works)
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
    """
    A catalog of stars relying on the cartoon photometry methods (which use non-LSST bandpasses
    and output extra data for use by unit tests)
    """
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
    """
    A catalog of galaxies relying on the cartoon photometry methods (which use non-LSST bandpasses
    and output extra data for use by unit tests)
    """
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
    """
    A generic catalog of stars
    """
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
    """
    A generic catalog of galaxies
    """
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

def cosmologicalOmega(redshift, H0, Om0, Ode0 = None, Og0=0.0, Onu0=0.0, w0=-1.0, wa=0.0):
    """
    A method to compute the evolution of the Hubble and density parameters
    with redshift (as a baseline against which to test the cosmology unittest)
    """

    if Ode0 is None:
        Ode0 = 1.0 - Om0 - Og0 - Onu0

    Ok0 = 1.0 - Om0 - Ode0 - Og0 - Onu0

    aa = 1.0/(1.0+redshift)
    Omz = Om0 * numpy.power(1.0+redshift, 3)
    Ogz = Og0 * numpy.power(1.0+redshift, 4)
    Onuz = Onu0 * numpy.power(1.0+redshift, 4)
    Okz = Ok0 * numpy.power(1.0+redshift, 2)
    Odez = Ode0 * numpy.exp(-3.0*(numpy.log(aa)*(w0 + wa +1.0) - wa*(aa - 1.0)))

    Ototal = Omz + Ogz + Onuz + Odez + Okz

    return H0*numpy.sqrt(Ototal), Omz/Ototal, Odez/Ototal, Ogz/Ototal, Onuz/Ototal, Okz/Ototal

def comovingDistanceIntegrand(redshift, H0, Om0, Ode0, Og0, Onu0, w0, wa):
    """
    The integrand of comoving distance (as a baseline for cosmology unittest)
    """
    hh, mm, de, gg, nn, kk = cosmologicalOmega(redshift, H0, Om0, Ode0=Ode0,
                                          Og0=Og0, Onu0=Onu0, w0=w0, wa=wa)
    return 1.0/hh
