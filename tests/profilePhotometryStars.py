import numpy

import sqlite3
from sqlite3 import dbapi2 as sqlite

import os
import warnings
import sys
import lsst.utils.tests as utilsTests

from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import DBObject, ObservationMetaData
from lsst.sims.coordUtils.Astrometry import AstrometryGalaxies, AstrometryStars
from lsst.sims.photUtils.Photometry import PhotometryGalaxies, PhotometryStars
from lsst.sims.photUtils.EBV import EBVmixin

from lsst.sims.photUtils.Variability import Variability

class testDefaults(object):

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

class testStars(InstanceCatalog,AstrometryStars,EBVmixin,Variability,PhotometryStars,testDefaults):
    catalog_type = 'test_stars'
    column_outputs=['id','ra_corr','dec_corr','magNorm',\
    'stellar_magNorm_var', \
    'lsst_u','sigma_lsst_u','lsst_u_var','sigma_lsst_u_var',
    'lsst_g','sigma_lsst_g','lsst_g_var','sigma_lsst_g_var',\
    'lsst_r','sigma_lsst_r','lsst_r_var','sigma_lsst_r_var',\
    'lsst_i','sigma_lsst_i','lsst_i_var','sigma_lsst_i_var',\
    'lsst_z','sigma_lsst_z','lsst_z_var','sigma_lsst_z_var',\
    'lsst_y','sigma_lsst_y','lsst_y_var','sigma_lsst_y_var',\
    'EBV','varParamStr']


rrlyDB=DBObject.from_objid('rrly')
obs_metadata_pointed=ObservationMetaData(mjd=2013.23, circ_bounds=dict(ra=200., dec=-30, radius=9.))
obs_metadata_pointed.metadata = {}
obs_metadata_pointed.metadata['Opsim_filter'] = 'i'
test_rrly=testStars(rrlyDB,obs_metadata=obs_metadata_pointed)
test_rrly.write_catalog("test_rrly_output.txt")


msDB=DBObject.from_objid('msstars')
obs_metadata_ms=ObservationMetaData(mjd=2013.23, circ_bounds=dict(ra=200., dec=-30, radius=0.1))
obs_metadata_ms.metadata = {}
obs_metadata_ms.metadata['Opsim_filter'] = 'i'
test_ms=testStars(msDB,obs_metadata=obs_metadata_ms)
test_ms.write_catalog("test_ms_output.txt")



wdDB=DBObject.from_objid('wdstars')

obs_metadata_pointed=ObservationMetaData(mjd=2013.23, circ_bounds=dict(ra=200., dec=-30, radius=0.5))
obs_metadata_pointed.metadata = {}
obs_metadata_pointed.metadata['Opsim_filter'] = 'i'

test_wd=testStars(wdDB,obs_metadata=obs_metadata_pointed)
test_wd.write_catalog("test_wd_output.txt")

