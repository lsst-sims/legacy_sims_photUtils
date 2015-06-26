"""
This file defines some test catalog and DBObject classes for use with unit tests.

To date (30 October 2014) testPhotometry.py and testCosmology.py import from this module
"""

import numpy
import os
import sqlite3
import json
from lsst.sims.catalogs.measures.instance import InstanceCatalog, register_method, register_class, compound
from lsst.sims.catalogs.generation.db import ObservationMetaData
from lsst.sims.coordUtils import AstrometryStars, AstrometryGalaxies
from lsst.sims.catalogs.generation.utils import myTestGals, myTestStars, \
                                                makeStarTestDB, makeGalTestDB, getOneChunk

from lsst.sims.photUtils.Photometry import PhotometryGalaxies, PhotometryStars
from lsst.sims.photUtils.SignalToNoise import calcSkyCountsPerPixelForM5

from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.EBV import EBVbase, EBVmixin

from lsst.sims.photUtils.Variability import Variability, VariabilityStars, VariabilityGalaxies

__all__ = ["setM5",
           "makeStarDatabase", "makeGalaxyDatabase",
           "TestVariabilityMixin", "testDefaults", "cartoonPhotometryStars",
           "cartoonPhotometryGalaxies", "testCatalog", "cartoonStars",
           "cartoonStarsOnlyI", "cartoonStarsIZ",
           "cartoonGalaxies", "cartoonGalaxiesIG", "testStars", "testGalaxies",
           "comovingDistanceIntegrand", "cosmologicalOmega",
           "galaxiesWithHoles"]

def setM5(m5target, skysed, totalBandpass, hardware,
          photParams,
          seeing = None):
    """
    Take an SED representing the sky and normalize it so that
    m5 (the magnitude at which an object is detected in this
    bandpass at 5-sigma) is set to some specified value.

    The 5-sigma limiting magnitude (m5) for an observation is
    determined by a combination of the telescope and camera parameters
    (such as diameter of the mirrors and the readnoise) together with the
    sky background. This method (setM5) scales a provided sky background
    Sed so that an observation would have a target m5 value, for the
    provided hardware parameters. Using the resulting Sed in the
    'calcM5' method will return this target value for m5.

    @param [in] the desired value of m5

    @param [in] skysed is an instantiation of the Sed class representing
    sky emission

    @param [in] totalBandpass is an instantiation of the Bandpass class
    representing the total throughput of the telescope (instrumentation
    plus atmosphere)

    @param [in] hardware is an instantiation of the Bandpass class representing
    the throughput due solely to instrumentation.

    @param [in] photParams is an instantiation of the
    PhotometricParameters class that carries details about the
    photometric response of the telescope.

    @param [in] seeing in arcseconds

    @param [out] returns an instantiation of the Sed class that is the skysed renormalized
    so that m5 has the desired value.

    Note that the returned SED will be renormalized such that calling the method
    self.calcADU(hardwareBandpass) on it will yield the number of counts per square
    arcsecond in a given bandpass.
    """

    #This is based on the LSST SNR document (v1.2, May 2010)
    #www.astro.washington.edu/users/ivezic/Astr511/LSST_SNRdoc.pdf

    if seeing is None:
        seeing = LSSTdefaults().seeing('r')

    skyCountsTarget = calcSkyCountsPerPixelForM5(m5target, totalBandpass, seeing=seeing,
                                             photParams=photParams)

    skySedOut = Sed(wavelen=numpy.copy(skysed.wavelen),
                    flambda=numpy.copy(skysed.flambda))

    skyCounts = skySedOut.calcADU(hardware, photParams=photParams) \
                    * photParams.platescale * photParams.platescale
    skySedOut.multiplyFluxNorm(skyCountsTarget/skyCounts)

    return skySedOut


def makeStarDatabase(filename='StellarPhotometryDB.db', size=1000, seedVal=32,
                     radius=1.0, unrefractedRA=50.0, unrefractedDec=-10.0):

    star_seds = ['km20_5750.fits_g40_5790','m2.0Full.dat','bergeron_6500_85.dat_6700']

    #Now begin building the database.
    #First create the tables.
    conn = sqlite3.connect(filename)
    c = conn.cursor()

    numpy.random.seed(seedVal)

    rr = numpy.random.sample(size)*radius
    theta = numpy.random.sample(size)*2.0*numpy.pi

    try:
        c.execute('''CREATE TABLE starsALL_forceseek
                  (simobjid int, ra real, decl real, magNorm real,
                  mudecl real, mura real, galacticAv real, vrad real, varParamStar text, sedFilename text, parallax real)''')
    except:
        raise RuntimeError("Error creating starsALL_forceseek table.")

    magnormStar = numpy.random.sample(size)*5.0+17.0
    magnormStar = numpy.random.sample(size)*4.0 + 17.0
    mudecl = numpy.random.sample(size)*0.0001
    mura = numpy.random.sample(size)*0.0001
    galacticAv = numpy.random.sample(size)*0.05*3.1
    vrad = numpy.random.sample(size)*1.0
    parallax = 0.00045+numpy.random.sample(size)*0.00001

    for i in range(size):
        raStar = unrefractedRA + rr[i]*numpy.cos(theta[i])
        decStar = unrefractedDec + rr[i]*numpy.sin(theta[i])

        cmd = '''INSERT INTO starsALL_forceseek VALUES (%i, %f, %f, %f, %f, %f, %f, %f, %s, '%s', %f)''' %\
                  (i, raStar, decStar, magnormStar[i], mudecl[i], mura[i],
                  galacticAv[i], vrad[i], 'NULL', star_seds[i%len(star_seds)], parallax[i])

        c.execute(cmd)

    conn.commit()
    conn.close()

def makeGalaxyDatabase(filename='GalaxyPhotometryDB.db', size=1000, seedVal=32,
                       radius=1.0, unrefractedRA=50.0, unrefractedDec=-10.0):

    galaxy_seds = ['Const.80E07.02Z.spec','Inst.80E07.002Z.spec','Burst.19E07.0005Z.spec']
    agn_sed = 'agn.spec'

    #Now begin building the database.
    #First create the tables.
    conn = sqlite3.connect(filename)
    c = conn.cursor()

    try:
        c.execute('''CREATE TABLE galaxy
                     (galtileid int, galid int, ra real, dec real,
                      bra real, bdec real, dra real, ddec real,
                      agnra real, agndec real,
                      magnorm_bulge, magnorm_disk, magnorm_agn,
                      sedname_bulge text, sedname_disk text, sedname_agn text,
                      varParamStr text,
                      a_b real, b_b real, pa_bulge real, bulge_n int,
                      a_d real, b_d real, pa_disk real, disk_n int,
                      ext_model_b text, av_b real, rv_b real,
                      ext_model_d text, av_d real, rv_d real,
                      u_ab real, g_ab real, r_ab real, i_ab real,
                      z_ab real, y_ab real,
                      redshift real, BulgeHalfLightRadius real, DiskHalfLightRadius real)''')

        conn.commit()
    except:
        raise RuntimeError("Error creating galaxy table.")

    mjd = 52000.0

    numpy.random.seed(seedVal)

    rr = numpy.random.sample(size)*radius
    theta = numpy.random.sample(size)*2.0*numpy.pi

    ra = unrefractedRA + rr*numpy.cos(theta)
    dec = unrefractedDec + rr*numpy.sin(theta)

    bra = numpy.radians(ra+numpy.random.sample(size)*0.01*radius)
    bdec = numpy.radians(dec+numpy.random.sample(size)*0.01*radius)
    dra = numpy.radians(ra + numpy.random.sample(size)*0.01*radius)
    ddec = numpy.radians(dec + numpy.random.sample(size)*0.01*radius)
    agnra = numpy.radians(ra + numpy.random.sample(size)*0.01*radius)
    agndec = numpy.radians(dec + numpy.random.sample(size)*0.01*radius)

    magnorm_bulge = numpy.random.sample(size)*4.0 + 17.0
    magnorm_disk = numpy.random.sample(size)*5.0 + 17.0
    magnorm_agn = numpy.random.sample(size)*5.0 + 17.0
    b_b = numpy.random.sample(size)*0.2
    a_b = b_b+numpy.random.sample(size)*0.05
    b_d = numpy.random.sample(size)*0.5
    a_d = b_d+numpy.random.sample(size)*0.1

    BulgeHalfLightRadius = numpy.random.sample(size)*0.2
    DiskHalfLightRadius = numpy.random.sample(size)*0.5

    pa_bulge = numpy.random.sample(size)*360.0
    pa_disk = numpy.random.sample(size)*360.0

    av_b = numpy.random.sample(size)*0.4
    av_d = numpy.random.sample(size)*0.4
    rv_b = numpy.random.sample(size)*0.1 + 3.0
    rv_d = numpy.random.sample(size)*0.1 + 3.0

    u_ab = numpy.random.sample(size)*4.0 + 17.0
    g_ab = numpy.random.sample(size)*4.0 + 17.0
    r_ab = numpy.random.sample(size)*4.0 + 17.0
    i_ab = numpy.random.sample(size)*4.0 + 17.0
    z_ab = numpy.random.sample(size)*4.0 + 17.0
    y_ab = numpy.random.sample(size)*4.0 +17.0
    redshift = numpy.random.sample(size)*2.0

    t0_mjd = numpy.random.sample(size)*10.0+mjd
    agn_tau = numpy.random.sample(size)*1000.0 + 1000.0
    agnSeed = numpy.random.random_integers(low=2, high=4000, size=size)
    agn_sfu = numpy.random.sample(size)
    agn_sfg = numpy.random.sample(size)
    agn_sfr = numpy.random.sample(size)
    agn_sfi = numpy.random.sample(size)
    agn_sfz = numpy.random.sample(size)
    agn_sfy = numpy.random.sample(size)

    for i in range(size):
        varParam = {'varMethodName':'applyAgn',
                    'pars':{'agn_tau':agn_tau[i], 't0_mjd':t0_mjd[i],
                    'agn_sfu':agn_sfu[i], 'agn_sfg':agn_sfg[i], 'agn_sfr':agn_sfr[i],
                    'agn_sfi':agn_sfi[i], 'agn_sfz':agn_sfz[i], 'agn_sfy':agn_sfy[i],
                    'seed':int(agnSeed[i])}}

        paramStr = json.dumps(varParam)

        cmd = '''INSERT INTO galaxy VALUES (%i, %i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,
                                            '%s', '%s', '%s', '%s',
                                            %f, %f, %f, %i,
                                            %f, %f, %f, %i,
                                            '%s', %f, %f,
                                            '%s', %f, %f,
                                            %f, %f, %f, %f, %f, %f,
                                            %f, %f, %f)''' %\
                     (i, i, ra[i], dec[i], bra[i], bdec[i], dra[i], ddec[i], agnra[i], agndec[i],
                     magnorm_bulge[i], magnorm_disk[i], magnorm_agn[i],
                     galaxy_seds[(i+1)%len(galaxy_seds)], galaxy_seds[i%len(galaxy_seds)], agn_sed,
                     paramStr,
                     a_b[i], b_b[i], pa_bulge[i], 4,
                     a_d[i], b_d[i], pa_disk[i], 1,
                     'CCM', av_b[i], rv_b[i],
                     'CCM', av_d[i], rv_d[i],
                     u_ab[i], g_ab[i], r_ab[i], i_ab[i], z_ab[i], y_ab[i], redshift[i],
                     BulgeHalfLightRadius[i], DiskHalfLightRadius[i])
        c.execute(cmd)

    conn.commit()
    conn.close()

@register_class
class TestVariabilityMixin(Variability):
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
        columnNames = [name for name in self.get_magnitudes._colnames]
        bandpassNames=['u','g','r','i','z']
        bandpassDir=os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'

        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandpassesFromFiles(bandpassNames,bandpassDir = bandpassDir,
                    bandpassRoot = 'test_bandpass_')

        output = self.meta_magnitudes_getter(idNames, columnNames)

        #############################################################################
        #Everything below this comment exists solely for the purposes of the unit test
        #if you need to write a customized getter for photometry that uses non-LSST
        #bandpasses, you only need to emulate the code above this comment.


        magNormList = self.column_by_name('magNorm')
        sedNames = self.column_by_name('sedFilename')

        #the two variables below will allow us to get at the SED and magnitude
        #data from within the unit test class, so that we can be sure
        #that the mixin loaded the correct bandpasses
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

        columnNames = [name for name in self.get_magnitudes._colnames]

        bandpassNames=['u','g','r','i','z']
        bandpassDir=os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'

        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandpassesFromFiles(bandpassNames,bandpassDir = bandpassDir,
                      bandpassRoot = 'test_bandpass_')

        output = self.meta_magnitudes_getter(idNames, columnNames)

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
            if Av is not None:
                self.applyAv(sublist, Av)

            self.applyRedshift(sublist, redshift)

            for ss in sublist:
                self.sedMasterDict[cc].append(ss)

        return output

class testCatalog(InstanceCatalog,AstrometryStars,VariabilityStars,testDefaults):
    catalog_type = 'MISC'
    default_columns=[('expmjd',5000.0,float)]

    def db_required_columns(self):
        return ['raJ2000'],['varParamStr']


class cartoonStars(InstanceCatalog,AstrometryStars,EBVmixin,VariabilityStars,cartoonPhotometryStars,testDefaults):
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

class cartoonStarsOnlyI(InstanceCatalog, AstrometryStars ,EBVmixin, VariabilityStars, PhotometryStars):
    catalog_type = 'cartoonStarsOnlyI'
    column_outputs = ['id','raObserved','decObserved','cartoon_i']

    #I need to give it the name of an actual SED file that spans the expected wavelength range
    defSedName = 'km30_5250.fits_g00_5370'
    default_columns = [('sedFilename', defSedName, (str,len(defSedName))), ('glon', 180., float),
                       ('glat', 30., float)]

    @compound('cartoon_u','cartoon_g','cartoon_r','cartoon_i','cartoon_z')
    def get_magnitudes(self):
        """
        Example photometry getter for alternative (i.e. non-LSST) bandpasses
        """

        idNames = self.column_by_name('id')
        columnNames = [name for name in self.get_magnitudes._colnames]
        bandpassNames=['u','g','r','i','z']
        bandpassDir=os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'

        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandpassesFromFiles(bandpassNames,bandpassDir = bandpassDir,
                    bandpassRoot = 'test_bandpass_')

        output = self.meta_magnitudes_getter(idNames, columnNames)
        return output

class cartoonStarsIZ(cartoonStarsOnlyI):
    catalog_type = 'cartoonStarsIR'
    column_outputs = ['id', 'raObserved', 'decObserved', 'cartoon_i', 'cartoon_z']

class cartoonGalaxies(InstanceCatalog, AstrometryGalaxies, EBVmixin, VariabilityGalaxies, cartoonPhotometryGalaxies, testDefaults):
    """
    A catalog of galaxies relying on the cartoon photometry methods (which use non-LSST bandpasses
    and output extra data for use by unit tests)
    """
    catalog_type = 'cartoonGalaxies'
    column_outputs=['galid','raObserved','decObserved',\
    'ctotal_u','ctotal_g','ctotal_r','ctotal_i','ctotal_z']

    #I need to give it the name of an actual SED file that spans the expected wavelength range
    defSedName = "Inst.80E09.25Z.spec"
    default_columns = [('sedFilenameBulge', defSedName, (str,len(defSedName))),
                       ('sedFilenameDisk', defSedName, (str,len(defSedName))),
                       ('sedFilenameAgn', defSedName, (str,len(defSedName))),
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


class cartoonGalaxiesIG(InstanceCatalog,AstrometryGalaxies,EBVmixin,VariabilityGalaxies,PhotometryGalaxies):

    catalog_type = 'cartoonGalaxiesIG'
    column_outputs=['galid','raObserved','decObserved','ctotal_i','ctotal_g']

    #I need to give it the name of an actual SED file that spans the expected wavelength range
    defSedName = "Inst.80E09.25Z.spec"
    default_columns = [('sedFilenameBulge', defSedName, (str,len(defSedName))),
                       ('sedFilenameDisk', defSedName, (str,len(defSedName))),
                       ('sedFilenameAgn', defSedName, (str,len(defSedName))),
                       ('glon', 210., float),
                       ('glat', 70., float),
                       ('internalAvBulge',3.1,float),
                       ('internalAvDisk',3.1,float)]

    def get_galid(self):
        return self.column_by_name('id')

    @compound('ctotal_u','ctotal_g','ctotal_r','ctotal_i','ctotal_z',
              'cbulge_u','cbulge_g','cbulge_r','cbulge_i','cbulge_z',
              'cdisk_u','cdisk_g','cdisk_r','cdisk_i','cdisk_z',
              'cagn_u','cagn_g','cagn_r','cagn_i','cagn_z')
    def get_magnitudes(self):
        """
        getter for photometry of galaxies using non-LSST bandpasses
        """

        idNames = self.column_by_name('galid')
        columnNames = [name for name in self.get_magnitudes._colnames]

        bandpassNames=['u','g','r','i','z']
        bandpassDir=os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'

        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandpassesFromFiles(bandpassNames,bandpassDir = bandpassDir,
                      bandpassRoot = 'test_bandpass_')

        output = self.meta_magnitudes_getter(idNames, columnNames)
        return output


class galaxiesWithHoles(InstanceCatalog, PhotometryGalaxies):
    """
    This is an InstanceCatalog of galaxies that sets some of the
    component Seds to 'None' so that we can test how sum_magnitudes
    handles NaN's in the context of a catalog.
    """
    column_outputs = ['raJ2000','decJ2000',
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
                      'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
                      'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn']

    default_formats = {'f':'%.12f'}

    defSedName = "Inst.80E09.25Z.spec"
    default_columns = [('glon', 210., float),
                       ('glat', 70., float),
                       ('internalAvBulge',3.1,float),
                       ('internalAvDisk',3.1,float)]


    def get_galid(self):
        return self.column_by_name('id')

    @compound('sedFilenameBulge', 'sedFilenameDisk', 'sedFilenameAgn')
    def get_sedNames(self):
        ra = self.column_by_name('raJ2000')
        elements = len(ra)
        bulge = []
        disk = []
        agn = []
        for ix in range(elements):
            bulge.append(self.defSedName)
            disk.append(self.defSedName)
            agn.append(self.defSedName)


        for ix in range(elements/8):
            ibase = ix*8
            if ibase+1<elements:
                bulge[ibase+1] = 'None'
            if ibase+2<elements:
                disk[ibase+2] = 'None'
            if ibase+3<elements:
                agn[ibase+3] = 'None'
            if ibase+4<elements:
                bulge[ibase+4] = 'None'
                disk[ibase+4] = 'None'
            if ibase+5<elements:
                bulge[ibase+5] = 'None'
                agn[ibase+5] = 'None'
            if ibase+6<elements:
                disk[ibase+6] = 'None'
                agn[ibase+6] = 'None'
            if ibase+7<elements:
                bulge[ibase+7] = 'None'
                disk[ibase+7] = 'None'
                agn[ibase+7] = 'None'


        return numpy.array([bulge, disk, agn])

class testStars(InstanceCatalog, EBVmixin, VariabilityStars, TestVariabilityMixin, PhotometryStars,testDefaults):
    """
    A generic catalog of stars
    """
    catalog_type = 'test_stars'
    column_outputs=['id','raJ2000','decJ2000','magNorm',\
    'lsst_u','sigma_lsst_u',
    'lsst_g','sigma_lsst_g',\
    'lsst_r','sigma_lsst_r',\
    'lsst_i','sigma_lsst_i',\
    'lsst_z','sigma_lsst_z',\
    'lsst_y','sigma_lsst_y',\
    'EBV','varParamStr']
    defSedName = 'sed_flat.txt'
    default_columns = [('sedFilename', defSedName, (str,len(defSedName))), ('glon', 180., float),
                       ('glat', 30., float)]

class testGalaxies(InstanceCatalog,EBVmixin,VariabilityGalaxies,TestVariabilityMixin,PhotometryGalaxies,testDefaults):
    """
    A generic catalog of galaxies
    """
    catalog_type = 'test_galaxies'
    column_outputs=['galid','raJ2000','decJ2000',\
        'redshift',
        'magNormAgn', 'magNormBulge', 'magNormDisk', \
        'lsst_u', 'sigma_lsst_u',\
        'lsst_g', 'sigma_lsst_g',\
        'lsst_r', 'sigma_lsst_r',\
         'lsst_i', 'sigma_lsst_i',\
         'lsst_z', 'sigma_lsst_z',\
         'lsst_y', 'sigma_lsst_y',\
        'sedFilenameBulge','uBulge', 'sigma_uBulge', 'gBulge', 'sigma_gBulge', \
        'rBulge', 'sigma_rBulge', 'iBulge', 'sigma_iBulge', 'zBulge', 'sigma_zBulge',\
         'yBulge', 'sigma_yBulge', \
        'sedFilenameDisk','uDisk', 'sigma_uDisk', 'gDisk', 'sigma_gDisk', 'rDisk', 'sigma_rDisk', \
        'iDisk', 'sigma_iDisk', 'zDisk', 'sigma_zDisk', 'yDisk', 'sigma_yDisk', \
        'sedFilenameAgn',\
        'uAgn', 'sigma_uAgn',\
        'gAgn', 'sigma_gAgn',\
        'rAgn', 'sigma_rAgn',\
        'iAgn', 'sigma_iAgn',\
        'zAgn', 'sigma_zAgn',\
        'yAgn', 'sigma_yAgn', 'varParamStr']
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

    @param [in] redshift is the redshift at which the output is desired

    @param [in] H0 is the Hubble parameter at the present epoch in km/s/Mpc

    @param [in] Om0 is the density parameter (fraction of critical) for matter at the
    present epoch

    @param [in] Ode0 is the density parameter for Dark Energy at the present epoch.
    If left as None, will be set to 1.0-Om0-Og0-Onu0 (i.e. a flat universe)

    @param [in] Og0 is the density parameter for photons at the present epoch

    @param [in] Onu0 is the density parameter for neutrinos at the present epoch
    (assume massless neutrinos)

    @param [in] w0 is a parameter for calculating the equation of state for Dark Energy
    w = w0 + wa * z/(1 + z)

    @param [in] wa is the other parameter for calculating the equation of state for Dark
    Energy

    @returns Hubble parameter at desired redshift (in km/s/Mpc)

    @returns matter density paramter at desired redshift

    @returns Dark Energy density parameter at desired redshift

    @returns photon density parameter at desired redshift

    @returns neutrino density parameter at desired redshift

    @returns curvature density parameter at desired redshift
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

    @param [in] redshift is the redshift at which to evaluate the integrand

    @param [in] H0 is the Hubble parameter at the present epoch in km/s/Mpc

    @param [in] Om0 is the density parameter (fraction of critical) for matter at the
    present epoch

    @param [in] Ode0 is the density parameter for Dark Energy at the present epoch.

    @param [in] Og0 is the density parameter for photons at the present epoch

    @param [in] Onu0 is the density parameter for neutrinos at the present epoch
    (assume massless neutrinos)

    @param [in] w0 is a parameter for calculating the equation of state for Dark Energy
    w = w0 + wa * z/(1 + z)

    @param [in] wa is the other parameter for calculating the equation of state for Dark
    Energy

    @returns 1/(Hubble parameter at desired redshift in km/s/Mpc)

    """
    hh, mm, de, gg, nn, kk = cosmologicalOmega(redshift, H0, Om0, Ode0=Ode0,
                                          Og0=Og0, Onu0=Onu0, w0=w0, wa=wa)
    return 1.0/hh
