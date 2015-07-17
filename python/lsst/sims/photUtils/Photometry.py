"""
photUtils -


ljones@astro.washington.edu  (and ajc@astro.washington.edu)

and now (2014 March 28): scott.f.daniel@gmail.com

Collection of utilities to aid usage of Sed and Bandpass with dictionaries.

"""

import os
import numpy
import eups
from collections import OrderedDict
from lsst.sims.photUtils import Sed, Bandpass, LSSTdefaults, calcGamma, \
                                calcMagError_m5, PhotometricParameters, magErrorFromSNR
from lsst.sims.utils import defaultSpecMap

__all__ = ["PhotometryHardware", "PhotometryBase"]

class PhotometryHardware(object):
    """
    This mixin provides the basic infrastructure for reading in bandpasses.

    To initiailize a different set of bandpasses, call self.loadBandPassesFromFiles() with a different
    set of arguments.

    Once self.loadBandPassesFromFiles() as been called, self.loadSeds() can be used to return an array
    of SED objects.  These objects can be passed to self.manyMagCalc_list() which will calculate
    the magnitudes of the the SEDs, integrated over the loaded bandpasses, and return them as a
    dict keeyed to the array of bandpass keys stored in self.bandpassKey
    """

    bandpassDict = None #total (hardware plus atmosphere) bandpasses loaded in this particular catalog;
                        #will be an OrderedDict
    hardwareBandpassDict = None #bandpasses associated with just the telescope hardware
    atmoTransmission = None #atmospheric transmissivity (will be an instantiation of the Bandpass class)
    nBandpasses = 0 #the number of bandpasses loaded (will need to be zero for InstanceCatalog's dry run)

    skySED = None #the emission spectrum of the atmosphere

    phiArray = None #the response curves for the bandpasses
    waveLenStep = None

    def loadBandpassesFromFiles(self, bandpassNames=['u', 'g', 'r', 'i', 'z', 'y'],
                                filedir = os.path.join(eups.productDir('throughputs'), 'baseline'),
                                bandpassRoot = 'filter_',
                                componentList = ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                                                 'lens1.dat', 'lens2.dat', 'lens3.dat'],
                                atmoTransmission='atmos.dat', skySED='darksky.dat'):
        """
        Load bandpass information from files.  This method will separate the bandpasses
        into contributions due to instrumentations and contributions due to the atmosphere.
        It will also load an SED associated with sky emission.

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

        @param [in] atmoTransmission is the name of the file representing the
        transmissivity of the atmosphere (also assumed to be in filedir)

        @param [in] skySED is the name of the file representing the emission
        spectrum of the atmosphere (also assumed to be in filedir)

        This method loads these files and stores:
        total bandpasses in self.bandpassDict
        hardware bandpasses in self.hardwareBandpassDict
        transmissivity of the atmosphere in self.atmoTransmission
        emission spectrum of the atmosphere in self.skySED

        This method will also setup the variables
        self.phiArray
        self.waveLenStep
        for purposes of photometric calculations
        """

        #This method isn't really used now.
        #It exists in anticipation of a time when we will will load
        #a skySED and normalize it based on the sky brightness, rather
        #than m5.

        self.skySED = Sed()
        self.skySED.readSED_flambda(os.path.join(filedir, skySED))

        self.atmoTransmission = Bandpass()
        self.atmoTransmission.readThroughput(os.path.join(filedir, atmoTransmission))

        commonComponents = []
        for cc in componentList:
            commonComponents.append(os.path.join(filedir,cc))

        self.bandpassDict = OrderedDict()
        self.hardwareBandpassDict = OrderedDict()
        self.gammaDict = None
        for w in bandpassNames:
            components = commonComponents + [os.path.join(filedir,"%s.dat" % (bandpassRoot +w))]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.hardwareBandpassDict[w] = bandpassDummy
            components += [os.path.join(filedir, atmoTransmission)]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.bandpassDict[w] = bandpassDummy

        self.phiArray = None
        self.waveLenStep = None
        self.setupPhiArray_dict()


    def loadTotalBandpassesFromFiles(self,bandpassNames=['u', 'g', 'r', 'i', 'z', 'y'],
                                bandpassDir = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline'),
                                bandpassRoot = 'total_'):
        """
        This will take the list of band passes named by bandpassNames and use them to set up
        self.bandpassDict (which is being cached so that
        it does not have to be loaded again unless we change which bandpasses we want)

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

        This method sets up the variables
        self.phiArray
        self.waveLenStep
        for purposes of photometric calculations
        """

        self.bandpassDict = OrderedDict()
        self.gammaDict = None

        for w in bandpassNames:
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughput(os.path.join(bandpassDir,"%s.dat" % (bandpassRoot + w)))
            self.bandpassDict[w] = bandpassDummy

        self.phiArray = None
        self.waveLenStep = None
        self.setupPhiArray_dict()

    def setupPhiArray_dict(self):
        """
        Generate 2-dimensional numpy array for Phi values associated with the bandpasses in
        self.bandpasses

        The results from this calculation will be stored in the instance variables
        self.phiArray and self.waveLenStep for future use by self.manyMagCalc_list()
        """

        sedobj = Sed()
        if self.bandpassDict is not None:
            self.phiArray, self.waveLenStep = sedobj.setupPhiArray(self.bandpassDict.values())
            self.nBandpasses = len(self.bandpassDict)


class PhotometryBase(PhotometryHardware):
    """
    This mixin provides the basic infrastructure for photometry.
    It can read in SEDs and bandpasses, apply extinction and redshift, and, given
    an SED object it can calculate magnitudes.

    In order to avoid duplication of work, the bandpasses, wavelength array, and phi array
    are stored as instance variables once they are read in by self.loadTotalBandpassesFromFiles()

    To initiailize a different set of bandpasses, call self.loadBandPassesFromFiles() with a different
    set of arguments.

    Once self.loadBandPassesFromFiles() as been called, self.loadSeds() can be used to return an array
    of SED objects.  These objects can be passed to self.manyMagCalc_list() which will calculate
    the magnitudes of the the SEDs, integrated over the loaded bandpasses, and return them as a
    dict keeyed to the array of bandpass keys stored in self.bandpassKey
    """

    #an object carrying around photometric parameters like readnoise, effective area, plate scale, etc.
    #defaults to LSST values
    photParams = PhotometricParameters()

    # cache member variables so we are not constantly
    # re-initializing the dust model for each Sed
    _av_wavelen = None
    _a_int = None
    _b_int = None


    # Handy routines for handling Sed/Bandpass routines with sets of dictionaries.
    def loadSeds(self, sedList, magNorm=None, resample_same=False, specFileMap=None):
        """
        Takes the list of filename sedList and returns an array of SED objects.

        This code will load identical SEDs twice because it is possible for
        (astronomical) objects to have the same SEDs but different magNorms

        @param [in] sedList is a list of file names containing Seds

        @param [in] magNorm is the magnitude normalization

        @param [in] resample_same governs whether or not to resample the Seds
        so that they are all on the same wavelength grid

        @param [in] specFileMap is a mapping from the names in sedList to the absolute
        path to the SED files.  It is an instantiation of the class defined in

        sims_catalogs_measures/python/lsst/sims/catalogs_measures/instance/fileMaps.py

        If not provided, a default will be instantiated.

        @param [out] sedOut is a list of Sed objects

        """

        if specFileMap is None:
            if hasattr(self, 'specFileMap'):
                specFileMap = self.specFileMap
            else:
                specFileMap = defaultSpecMap

        dataDir=os.getenv('SIMS_SED_LIBRARY_DIR')

        #initialize a delta function bandpass for use in applying magNorm
        imsimband = Bandpass()
        imsimband.imsimBandpass()

        sedOut=[]

        #uniqueSedDict will store all of the unique SED files that have been
        #loaded.  If an object requires an SED that has already been loaded,
        #it will just copy it from the dict.
        uniqueSedDict={}

        firstsed = True
        uniqueSedDict["None"] = Sed()
        for i in range(len(sedList)):
            sedName = sedList[i]

            if sedName not in uniqueSedDict:
                sed = Sed()
                sed.readSED_flambda(os.path.join(dataDir, specFileMap[sedName]))

                if resample_same:
                    if firstsed:
                        wavelen_same = sed.wavelen
                        firstsed = False
                    else:
                        sed.resampleSED(wavelen_same)

                uniqueSedDict[sedName]=sed

        #now that we have loaded and copied all of the necessary SEDs,
        #we can apply magNorms
        for i in range(len(sedList)):

            ss = uniqueSedDict[sedList[i]]

            sed=Sed(wavelen=ss.wavelen,flambda=ss.flambda,fnu=ss.fnu, name=ss.name)

            if sedList[i] != "None":
                fNorm = sed.calcFluxNorm(magNorm[i], imsimband)
                sed.multiplyFluxNorm(fNorm)

            sedOut.append(sed)

        return sedOut


    def applyAv(self, sedList, internalAvList):
        """
        Take the array of SED objects sedList and apply the arrays of extinction and redshift
        (internalAV and redshift)

        This method does not return anything.  It makes the necessary changes
        to the Seds in SedList in situ.

        @param [in] sedList is a list of Sed objects

        @param [in] internalAvList is a list of Av extinction values internal to each object

        This method will apply the O'Donnel 1994 dust model as encoded in the Sed class
        to each Sed in sedList.
        """

        if internalAvList is None:
            return

        for sedobj, Av in zip(sedList, internalAvList):
            if sedobj.wavelen is not None:
                #setupCCMab only depends on the wavelen array
                #because this is supposed to be the same for every
                #SED object in sedList, it is only called once for
                #each invocation of applyAv
                if self._av_wavelen is None or (sedobj.wavelen!=self._av_wavelen).any():
                    self._a_int, self._b_int = sedobj.setupCCMab()
                    self._av_wavelen=sedobj.wavelen

                sedobj.addCCMDust(self._a_int, self._b_int, A_v=Av)


    def applyRedshift(self, sedList, redshiftList, dimming=True):
        """
        Take the array of SED objects sedList and apply the arrays of extinction and redshift
        (internalAV and redshift)

        This method does not return anything.  It makes the necessary changes
        to the Seds in SedList in situ.

        @param [in] sedList is a list of Sed objects

        @param [in] redshiftList is a list of redshift values

        @param [in] dimming is a boolean that turns on or off cosmological
        dimming (the extra (1+z)^-1 factor in flux).  Defaults to True.

        This method will redshift each Sed object in sedList
        """

        if redshiftList is None:
            return

        for sedobj, redshift in zip(sedList, redshiftList):
            if sedobj.wavelen is not None:
                sedobj.redshiftSED(redshift, dimming=dimming)
                sedobj.name = sedobj.name + '_Z' + '%.2f' %(redshift)
                sedobj.resampleSED(wavelen_match=self.bandpassDict.values()[0].wavelen)


    def manyMagCalc_list(self, sedobj, indices=None):
        """
        Return a list of magnitudes for a single Sed object.

        Bandpass information is taken from the instance variables self.bandpassDict,
        self.phiArray, and self.waveLenStep

        @param [in] sedobj is an Sed object

        @param [in] indices is an optional list of indices indicating which bandpasses to actually
        calculate magnitudes for.  Other magnitudes will be listed as 'None' (i.e. this method will
        return as many magnitudes as were loaded with the loadBandpassesFromFiles methods; it will
        just return nonsense for magnitudes you did not actually ask for)

        @param [out] magList is a list of magnitudes in the bandpasses stored in self.bandpassDict
        """
        # Set up the SED for using manyMagCalc - note that this CHANGES sedobj
        # Have to check that the wavelength range for sedobj matches bandpass - this is why the dictionary is passed in.

        magList = []
        if sedobj.wavelen is not None:
            sedobj.resampleSED(wavelen_match=self.bandpassDict.values()[0].wavelen)

            #for some reason, moving this call to flambdaTofnu()
            #to a point earlier in the
            #process results in some SEDs having 'None' for fnu.
            #
            #I looked more carefully at the documentation in Sed.py
            #Any time you update flambda in any way, fnu gets set to 'None'
            #This is to prevent the two arrays from getting out synch
            #(e.g. renormalizing flambda but forgettint to renormalize fnu)
            #
            sedobj.flambdaTofnu()

            if indices is not None:
                for i in range(self.nBandpasses):
                    magList.append(numpy.NaN)

                magArray = sedobj.manyMagCalc(self.phiArray, self.waveLenStep, observedBandPassInd=indices)
                for i,ix in enumerate(indices):
                    magList[ix] = magArray[i]
            else:
                magList = sedobj.manyMagCalc(self.phiArray, self.waveLenStep)

        else:
            for i in range(self.nBandpasses):
                magList.append(numpy.NaN)

        return magList

    def calculateMagnitudeUncertainty(self, magnitudes, obs_metadata=None):
        """
        Calculate the uncertainty in magnitudes using the model from equation (5) of arXiv:0805.2366

        @param [in] magnitudes is a numpy array containing the object magnitudes.  Every row corresponds to
        a bandpass, which is to say that magnitudes[i][j] is the magnitude of the jth object in the
        bandpass characterized by self.bandpassDict.values()[i]

        @param [in] obs_metadata is the metadata of this observation (mostly desired because
        it will contain information about m5, the magnitude at which objects are detected
        at the 5-sigma limit in each bandpass)

        @param [out] a 1-dimensional numpy array containing the magntiude uncertainty
        corresponding to the magnitudes passed into this method.
        """

        if obs_metadata is None:
            raise RuntimeError("Need to pass an ObservationMetaData into calculatePhotometricUncertainty")

        if magnitudes.shape[0] != self.nBandpasses:
            raise RuntimeError("Passed %d magnitudes to " % magnitudes.shape[0] + \
                                " PhotometryBase.calculatePhotometricUncertainty; " + \
                                "needed %d " % self.nBandpasses)

        #if we have not run this method before, calculate and cache the m5 and gamma parameter
        #values (see eqn 5 of arXiv:0805.2366) for future use
        m5Defaults = None

        if not hasattr(self, '_gammaList') or len(self._gammaList) != self.nBandpasses:
            mm = []
            gg = []
            for b in self.bandpassDict.keys():
                if b in obs_metadata.m5:
                    mm.append(obs_metadata.m5[b])
                    gg.append(calcGamma(self.bandpassDict[b], obs_metadata.m5[b],
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

        error = calcMagError_m5(magnitudes, self.bandpassDict.values(), self._m5List,
                                gamma=self._gammaList, photParams=self.photParams)

        return error
