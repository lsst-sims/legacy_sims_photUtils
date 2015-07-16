from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog

__all__ = ["setupPhotometryCatalog"]

def setupPhotometryCatalog(obs_metadata=None, dbConnection=None, catalogClass=None,
                           photometryNameRoot='lsst', uncertainty=False):
    """
    This method will read in an InstanceCatalog class (not an instantiation of that class;
    a class itself), and instantiation of ObservationMetaData and an instantiation of
    a CatalogDBObject class, instantiate the InstanceCatalog class according to the
    ObservationMetaData and CatalogDBObject and return the instantiated instanceCatalog.

    Note that this class will tell the InstanceCatalog to only output photometric columns
    defined by the ObservationMetaData (i.e. if the ObservationMetaData.bandpass is 'u',
    then the instantiated InstanceCatalog will only write the u-based column of photometry.

    @param [in] obs_metadata is an instantiation of an ObservationMetaData

    @param [in] dbConnection is an instantiation of a CatalogDBObject daughter class

    @param [in] catalogClass is a daughter class of InstanceCatalog (not an instantiation of
    that class; just the class itself)

    @param [in] photometryNameRoot is a string indicating the naming convention of the
    photometry columns.  This method will assume that columns are named

    photometryNameRoot_bandpassName

    Note: it is your responsibility to make sure that the catalogClass is able to calculate
    the columns specified this way (i.e. if photometryNameRoot is 'sdss', you should make
    sure that your catalogClass has getters for columns sdss_[u,g,r,i,z])

    @param [out] an instantiation of the provided InstanceCatalog daughter class that is
    consistent with the other inputs
    """

    if not isinstance(obs_metadata, ObservationMetaData):
        raise RuntimeError('obs_metadata needs to be an instantiation of ObservationMetaData')

    if not isinstance(dbConnection, CatalogDBObject):
        raise RuntimeError('dbConnection needs to be an instantiation of a CatalogDBObject daughter class')

    if not issubclass(catalogClass, InstanceCatalog):
        raise RuntimeError('catalogClass needs to be a daughter class of InstanceCatalog')

    column_outputs = None #this will be how we specify the photometric columns for the
                          #catalogClass

    if obs_metadata.bandpass is not None:
        if hasattr(obs_metadata.bandpass, '__iter__'):
            for b in obs_metadata.bandpass:
                if column_outputs is None:
                    column_outputs = [photometryNameRoot+'_'+b]
                else:
                    column_outputs.append(photometryNameRoot+'_'+b)

                if uncertainty:
                    column_outputs.append('sigma_'+photometryNameRoot+'_'+b)
        else:
            column_outputs = [photometryNameRoot+'_'+obs_metadata.bandpass]
            if uncertainty:
                column_outputs.append('sigma_'+photometryNameRoot+'_'+obs_metadata.bandpass)

    instantiation = catalogClass(dbConnection, obs_metadata=obs_metadata, column_outputs=column_outputs)
    return instantiation
