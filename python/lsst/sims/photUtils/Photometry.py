"""
photUtils -


ljones@astro.washington.edu  (and ajc@astro.washington.edu)

and now (2014 March 28): scott.f.daniel@gmail.com

Collection of utilities to aid usage of Sed and Bandpass with dictionaries.

"""

import os
import numpy
from collections import OrderedDict
import lsst.utils
from lsst.sims.photUtils import Sed, Bandpass, LSSTdefaults, calcGamma, \
                                calcMagError_m5, PhotometricParameters, \
                                magErrorFromSNR, CatSimBandpassDict
from lsst.sims.utils import defaultSpecMap

__all__ = ["loadBandpassesFromFiles", "loadTotalBandpassesFromFiles", "PhotometryBase"]


def loadBandpassesFromFiles(bandpassNames=['u', 'g', 'r', 'i', 'z', 'y'],
                            filedir = os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline'),
                            bandpassRoot = 'filter_',
                            componentList = ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                                             'lens1.dat', 'lens2.dat', 'lens3.dat'],
                            atmoTransmission=os.path.join(lsst.utils.getPackageDir('throughputs'),
                                                          'baseline','atmos.dat')):
    """
    Load bandpass information from files into CatSimBandpassDicts.
    This method will separate the bandpasses into contributions due to instrumentations
    and contributions due to the atmosphere.

    @param [in] bandpassNames is a list of strings labeling the bandpasses
    (e.g. ['u', 'g', 'r', 'i', 'z', 'y'])

    @param [in] filedir is a string indicating the name of the directory containing the
    bandpass files

    @param [in] bandpassRoot is the root of the names of the files associated with the
    bandpasses.  This method assumes that bandpasses are stored in
    filedir/bandpassRoot_bandpassNames[i].dat

    @param [in] componentList lists the files associated with bandpasses representing
    hardware components shared by all filters
    (defaults to ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat', 'lens1.dat',
                  'lense2.dat', 'lenst3.dat']
    for LSST).  These files are also expected to be stored in filedir

    @param [in] atmoTransmission is the absolute path to the file representing the
    transmissivity of the atmosphere (defaults to baseline/atmos.dat in the LSST
    'throughputs' package).

    @param [out] bandpassDict is a CatSimBandpassDict containing the total
    throughput (instrumentation + atmosphere)

    @param [out] hardwareBandpassDict is a CatSimBandpassDict containing
    the throughput due to instrumentation only
    """

    commonComponents = []
    for cc in componentList:
        commonComponents.append(os.path.join(filedir,cc))

    bandpassList = []
    hardwareBandpassList = []

    for w in bandpassNames:
        components = commonComponents + [os.path.join(filedir,"%s.dat" % (bandpassRoot +w))]
        bandpassDummy = Bandpass()
        bandpassDummy.readThroughputList(components)
        hardwareBandpassList.append(bandpassDummy)

        components += [atmoTransmission]
        bandpassDummy = Bandpass()
        bandpassDummy.readThroughputList(components)
        bandpassList.append(bandpassDummy)


    bandpassDict = CatSimBandpassDict(bandpassList, bandpassNames)
    hardwareBandpassDict = CatSimBandpassDict(hardwareBandpassList, bandpassNames)

    return bandpassDict, hardwareBandpassDict


def loadTotalBandpassesFromFiles(bandpassNames=['u', 'g', 'r', 'i', 'z', 'y'],
                                bandpassDir = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline'),
                                bandpassRoot = 'total_'):
    """
    This will take the list of band passes named by bandpassNames and load them into
    a CatSimBandpassDict

    The bandpasses loaded this way are total bandpasses: they account for instrumental
    and atmospheric transmission.

    @param [in] bandpassNames is a list of names identifying each filter.
    Defaults to ['u', 'g', 'r', 'i', 'z', 'y']

    @param [in] bandpassDir is the name of the directory where the bandpass files are stored

    @param [in] bandpassRoot contains the first part of the bandpass file name, i.e., it is assumed
    that the bandpasses are stored in files of the type

    bandpassDir/bandpassRoot_bandpassNames[i].dat

    if we want to load bandpasses for a telescope other than LSST, we would do so
    by altering bandpassDir and bandpassRoot

    @param [out] bandpassDict is a CatSimBandpassDict containing the loaded throughputs
    """

    bandpassList = []

    for w in bandpassNames:
        bandpassDummy = Bandpass()
        bandpassDummy.readThroughput(os.path.join(bandpassDir,"%s.dat" % (bandpassRoot + w)))
        bandpassList.append(bandpassDummy)

    return CatSimBandpassDict(bandpassList, bandpassNames)


class PhotometryBase(object):
    """
    This class provides the basic infrastructure for photometry.
    It can read in SEDs and bandpasses, apply extinction and redshift, and, given
    an SED object it can calculate magnitudes.

    In order to avoid duplication of work, the bandpasses, wavelength array, and phi array
    are stored as instance variables once they are read in by either self.loadTotalBandpassesFromFiles()
    or self.loadBandpassesFromFiles()

    See the documentation of PhotometryHardware (from which this class inherits) for instructions how
    to read in alternative bandpasses.

    Once self.loadBandPassesFromFiles() as been called, self.loadSeds() can be used to return an array
    of SED objects.  These objects can be passed to self.manyMagCalc_list() which will calculate
    the magnitudes of the the SEDs, integrated over the loaded bandpasses, and return them as a
    dict keeyed to the array of bandpass keys stored in self.bandpassKey
    """

    #an object carrying around photometric parameters like readnoise, effective area, plate scale, etc.
    #defaults to LSST values
    photParams = PhotometricParameters()

    def calculateMagnitudeUncertainty(self, magnitudes, bandpassDict, obs_metadata=None):
        """
        Calculate the uncertainty in magnitudes using the model from equation (5) of arXiv:0805.2366

        @param [in] magnitudes is a numpy array containing the object magnitudes.  Every row corresponds to
        a bandpass, which is to say that magnitudes[i][j] is the magnitude of the jth object in the
        bandpass characterized by self.bandpassDict.values()[i]

        @param [in] bandpassDict is a CatSimBandpassDict characterizing the bandpasses being used

        @param [in] obs_metadata is the metadata of this observation (mostly desired because
        it will contain information about m5, the magnitude at which objects are detected
        at the 5-sigma limit in each bandpass)

        @param [out] a 1-dimensional numpy array containing the magntiude uncertainty
        corresponding to the magnitudes passed into this method.
        """

        if obs_metadata is None:
            raise RuntimeError("Need to pass an ObservationMetaData into calculatePhotometricUncertainty")

        if magnitudes.shape[0] != len(bandpassDict):
            raise RuntimeError("Passed %d magnitudes to " % magnitudes.shape[0] + \
                                " PhotometryBase.calculatePhotometricUncertainty; " + \
                                "needed %d " % len(bandpassDict))

        #if we have not run this method before, calculate and cache the m5 and gamma parameter
        #values (see eqn 5 of arXiv:0805.2366) for future use
        m5Defaults = None

        if not hasattr(self, '_gammaList') or \
        len(self._gammaList) != bandpassDict.nBandpasses:

            mm = []
            gg = []
            for b in bandpassDict.keys():
                if b in obs_metadata.m5:
                    mm.append(obs_metadata.m5[b])
                    gg.append(calcGamma(bandpassDict[b], obs_metadata.m5[b],
                              photParams=self.photParams))
                else:
                    if m5Defaults is None:
                        m5Defaults = LSSTdefaults()

                    try:
                        mm.append(m5Defaults.m5(b))
                        gg.append(m5Defaults.gamma(b))
                    except:
                        raise RuntimeError("No way to calculate gamma or m5 for filter %s " % b)

            self._m5List = numpy.array(mm)
            self._gammaList = numpy.array(gg)

        error = calcMagError_m5(magnitudes, bandpassDict.values(), self._m5List,
                                gamma=self._gammaList, photParams=self.photParams)

        return error
