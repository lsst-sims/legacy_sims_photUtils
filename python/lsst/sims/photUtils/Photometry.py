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

__all__ = ["loadBandpassesFromFiles", "loadTotalBandpassesFromFiles"]


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


