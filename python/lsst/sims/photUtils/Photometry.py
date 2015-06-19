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
                                calcSNR_gamma, PhotometricParameters
from lsst.sims.catalogs.measures.instance import defaultSpecMap
from lsst.sims.catalogs.measures.instance import compound

__all__ = ["PhotometryHardware", "PhotometryBase", "PhotometryGalaxies", "PhotometryStars"]

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

    sigmaSysSq = None #square of systematic noise associated with this catalog

    #an object carrying around photometric parameters like readnoise, effective area, plate scale, etc.
    #defaults to LSST values
    photParams = PhotometricParameters()

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

    def applyAvAndRedshift(self,sedList, internalAv=None, redshift=None):
        """
        Take the array of SED objects sedList and apply the arrays of extinction and redshift
        (internalAV and redshift)

        This method does not return anything.  It makes the necessary changes
        to the Seds in SedList in situ.

        @param [in] sedList is a list of Sed objects

        @param [in] internalAv is the Av extinction internal to the object

        @param [in] redshift

        """

        wavelen_sampled=[]

        for i in range(len(sedList)):
            if sedList[i].wavelen is not None:
                if internalAv is not None:
                    #setupCCMab only depends on the wavelen array
                    #because this is supposed to be the same for every
                    #SED object in sedList, it is only called once for
                    #each invocation of applyAvAndRedshift
                    if wavelen_sampled == [] or (sedList[i].wavelen!=wavelen_sampled).any():
                        a_int, b_int = sedList[i].setupCCMab()
                        wavelen_sampled=sedList[i].wavelen

                    sedList[i].addCCMDust(a_int, b_int, A_v=internalAv[i])
                if redshift is not None:
                    #17 November 2014
                    #We do not apply cosmological dimming here because that is
                    #a wavelength-independent process and the galaxy's
                    #magNorm presumably accounts for that (i.e., we are assuming
                    #that the magNorm is calibrated to the galaxy's apparent
                    #magnitude in the imsimband, rather than some absolute
                    #magnitude).
                    sedList[i].redshiftSED(redshift[i], dimming=False)
                    sedList[i].name = sedList[i].name + '_Z' + '%.2f' %(redshift[i])
                    sedList[i].resampleSED(wavelen_match=self.bandpassDict.values()[0].wavelen)

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

    def calculateFluxSignalToNoise(self, fluxes, obs_metadata=None, sigmaSysSq=None):
        """
        Calculate the signal to noise in flux using the model from equation (5) of arXiv:0805.2366

        @param [in] fluxes is a list containing the object fluxes.  Every row corresponds to
        a bandpass, which is to say that fluxes[i][j] is the magnitude of the jth object in the
        bandpass characterized by self.bandpassDict.values()[i]

        @param [in] obs_metadata is the metadata of this observation (mostly desired because
        it will contain information about m5, the magnitude at which objects are detected
        at the 5-sigma limit in each bandpass)

        @param [in] sigmaSysSq the square of the systematic noise in flux

        @param [out] snr is a 1-dimensional numpy array containing the signal-to-noise in flux
        corresponding to the fluxes passed into this method.
        """

        if obs_metadata is None:
            raise RuntimeError("Need to pass an ObservationMetaData into calculatePhotometricUncertainty")

        if fluxes.shape[0] != self.nBandpasses:
            raise RuntimeError("Passed %d magnitudes to " % fluxes.shape[0] + \
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

        snr, gamma = calcSNR_gamma(fluxes, self.bandpassDict.values(), self._m5List, gamma=self._gammaList,
                                   sigmaSysSq=sigmaSysSq, photParams=self.photParams)

        return snr


    def calculateMagnitudeUncertainty(self, magnitudes, obs_metadata=None, sigmaSysSq=None):
        """
        Calculate the uncertainty in magnitude

        @param [in] magnitudes is a numpy array in which magnitude[i][j] is the magnitude
        of the jth astronomical object in the bandpass self.bandpassDict.values()[i]

        @param [in] obs_metadata is an instantiation of the ObservationMetaData class

        @param [in] sigmaSysSq is the square of the systematic signal to noise in flux

        @param [out] an array of magnitude uncertainties corresponding to the input magnitudes
        """

        snr = self.calculateFluxSignalToNoise(Sed().fluxFromMag(magnitudes),
                                              obs_metadata=obs_metadata, sigmaSysSq=sigmaSysSq)

        #see www.ucolick.org/~bolte/AY257/s_n.pdf section 3.1
        return 2.5*numpy.log10(1.0+1.0/snr)


class PhotometryGalaxies(PhotometryBase):
    """
    This mixin provides the code necessary for calculating the component magnitudes associated with
    galaxies.  It assumes that we want LSST filters.
    """

    def calculate_component_magnitudes(self,objectID, componentNames, \
                                       magNorm = None, internalAv = None, redshift = None,
                                       cosmologicalDistanceModulus = None, specFileMap=None,
                                       indices=None):

        """
        Calculate the magnitudes for different components (disk, bulge, agn, etc) of galaxies.
        This method is designed to be used such that you feed it all of the disk Seds from your data
        base and it returns the associated magnitudes.  Then you feed it all of the bulge Seds, etc.

        @param [in] objectID is the name of the galaxies (the whole galaxies)

        @param [in] componentNames gives the name of the SED filenames

        @param [in] magNorm is the normalizing magnitude

        @param [in] internalAv is the internal Av extinction

        @param [in] redshift is pretty self-explanatory

        @param [in] cosmologicalDistanceModulus is the effective distance modulus due to cosmological
        expansion (if that has not already been accounted for in magNorm).  This is optional.

        @param [in] specFileMap is a mapping between the filenames in diskNames, bulgeNames, and agnNames
        and the absolute locations of the corresponding files.  It is an instantiation of the class defined
        in

        sims_catalogs_measures/python/lsst/sims/catalogs/measures/instance/fileMaps.py

        If not provided, a default will be instantiated.

        @param [in] indices is an optional list of indices indicating which bandpasses to actually
        calculate magnitudes for

        @param [out] componentMags is a dict of lists such that
        magnitude["objectname"][i] will return the magnitude in the ith
        for the associated component Sed

        """


        componentMags = {}

        if componentNames != [] and componentNames is not None:
            componentSed = self.loadSeds(componentNames, magNorm = magNorm, specFileMap=specFileMap)
            self.applyAvAndRedshift(componentSed, internalAv = internalAv, redshift = redshift)

            for i in range(len(objectID)):
                subList = self.manyMagCalc_list(componentSed[i], indices=indices)

                if isinstance(cosmologicalDistanceModulus, numpy.ndarray):
                    for j in range(len(subList)):
                        subList[j] += cosmologicalDistanceModulus[i]

                componentMags[objectID[i]] = subList

        else:
            subList=[]
            for i in range(self.nBandpasses):
                subList.append(numpy.NaN)
            for i in range(len(objectID)):
                componentMags[objectID[i]]=subList

        return componentMags

    def sum_magnitudes(self, disk = None, bulge = None, agn = None):
        """
        Sum the component magnitudes of a galaxy and return the answer

        @param [in] disk is the disk magnitude must be a numpy array or a float

        @param [in] bulge is the bulge magnitude must be a numpy array or a float

        @param [in] agn is the agn magnitude must be a numpy array or a float

        @param [out] outMag is the total magnitude of the galaxy
        """

        baselineType = type(None)
        if not isinstance(disk, type(None)):
            baselineType = type(disk)
            if baselineType == numpy.ndarray:
                elements=len(disk)

        if not isinstance(bulge, type(None)):
            if baselineType == type(None):
                baselineType = type(bulge)
                if baselineType == numpy.ndarray:
                    elements = len(bulge)
            elif not isinstance(bulge, baselineType):
                raise RuntimeError("All non-None arguments of sum_magnitudes need to be " +
                                   "of the same type (float or numpy array)")

        elif not isinstance(agn, type(None)):
            if baseLineType == type(None):
                baselineType = type(agn)
                if baselineType == numpy.ndarray:
                    elements = len(agn)
            elif not isinstance(agn, baselineType):
                raise RuntimeError("All non-None arguments of sum_magnitudes need to be " +
                                   "of the same type (float or numpy array)")

        if baselineType is not float and \
           baselineType is not numpy.ndarray and \
           baselineType is not numpy.float and \
           baselineType is not numpy.float64:

            raise RuntimeError("Arguments of sum_magnitudes need to be " +
                               "either floats or numpy arrays; you appear to have passed %s " % baselineType)

        mm_0 = 22.
        tol = 1.0e-30

        if baselineType == numpy.ndarray:
            nn = numpy.zeros(elements)
        else:
            nn = 0.0

        if disk is not None:
            nn += numpy.where(numpy.isnan(disk), 0.0, numpy.power(10, -0.4*(disk - mm_0)))

        if bulge is not None:
            nn += numpy.where(numpy.isnan(bulge), 0.0, numpy.power(10, -0.4*(bulge - mm_0)))

        if agn is not None:
            nn += numpy.where(numpy.isnan(agn), 0.0, numpy.power(10, -0.4*(agn - mm_0)))

        if baselineType == numpy.ndarray:
            # according to this link
            # http://stackoverflow.com/questions/25087769/runtimewarning-divide-by-zero-error-how-to-avoid-python-numpy
            # we will still get a divide by zero error from log10, but numpy.where will be
            # circumventing the offending value, so it is probably okay
            return numpy.where(nn>tol, -2.5*numpy.log10(nn) + mm_0, numpy.NaN)
        else:
            if nn>tol:
                return -2.5*numpy.log10(nn) + mm_0
            else:
                return numpy.NaN



    def calculate_magnitudes(self, objectID, diskNames=None, diskMagNorm=None, diskAv=None,
                             bulgeNames=None, bulgeMagNorm=None, bulgeAv=None,
                             agnNames=None, agnMagNorm=None,
                             redshift=None, cosmologicalDistanceModulus=None, specFileMap=None,
                             indices=None):
        """
        Take the array of bandpasses in self.bandpassDict and the array of galaxy
        names objectID ane return a dict of dicts of lists of magnitudes

        the first level key is galid (the name of the galaxy)

        the second level key is "bulge", "disk", or "agn"

        this yields a list of magnitudes corresponding to the bandpasses in self.bandpassDict

        We need to index the galaxies by some unique identifier, such as galid
        because it is possible for galaxies to have the same sed filenames but
        different normalizations

        @param [in] objectID is a list of names uniquely identifying the objects whose magnitudes
        are being calculated

        @param [in] diskNames is a list of the names of the files containing disk SEDs

        @param [in] diskMagNorm is a list of magnitude normalizations for disk SEDs

        @param [in] diskAv is a list of extinction Av due to dust internal to the disk of the galaxy

        @param [in] bulgeNames is a list of the names of the files containing bulge SEDs

        @param [in] bulgeMagNorm is a list of the magnitude normalizations of the bulge SEDs

        @param [in] bulgeAv is a ist of extinction Av due to dust internal to the bulge of the galaxy

        @param [in] agnNames is a list of the names of the files containing AGN SEDs

        @param [in] agnMagNorm is a list of the magnitude normalizations of the AGN SEDs

        @param [in] redshift is a list of the redshifts of the galaxies

        @param [in] cosmologicalDistanceModulus is a list of the distance modulii due to
        cosmological expansion (assuming that has not been accounted for by the magNorms).
        This is optional.

        @param [in] specFileMap is a mapping between the filenames in diskNames, bulgeNames, and agnNames
        and the absolute locations of the corresponding files.  It is an instantiation of the class defined
        in

        sims_catalogs_measures/python/lsst/sims/catalogs/measures/instance/fileMaps.py

        If not provided, a default will be instantiated.

        @param [in] indices is an optional list of indices indicating which bandpasses to actually
        calculate magnitudes for

        @param [out] masterDict is a dict of magnitudes such that
        masterDict['AAA']['BBB'][i] is the magnitude in the ith bandpass of component BBB of galaxy AAA


        """

        if specFileMap is None:
            if hasattr(self, 'specFileMap'):
                specFileMap = self.specFileMap
            else:
                specFileMap = defaultSpecMap

        if diskNames is not None:
            if diskAv is None:
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need diskAv')

            if diskMagNorm is None:
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need diskMagNorm')

            if len(diskNames) != len(objectID):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d diskNames'
                                   % (len(diskNames), len(objectID)))
            if len(diskNames) != len(diskAv) or len(diskNames) != len(diskMagNorm) or len(diskMagNorm) != len(diskAv):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d diskNames, %d diskAvs, and %d diskMagNorms'
                                   % (len(diskNames), len(diskAv), len(diskMagNorm)))

        if bulgeNames is not None:
            if bulgeAv is None:
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need bulgeAv')

            if bulgeMagNorm is None:
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need bulgeMagNorm')

            if len(bulgeNames) != len(objectID):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d bulgeNames'
                                   % (len(bulgeNames), len(objectID)))
            if len(bulgeNames) != len(bulgeAv) or len(bulgeNames) != len(bulgeMagNorm) or len(bulgeMagNorm) != len(bulgeAv):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d bulgeNames, %d bulgeAvs, and %d bulgeMagNorms'
                                   % (len(bulgeNames), len(bulgeAv), len(bulgeMagNorm)))

        if agnNames is not None:
            if agnMagNorm is None:
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need agnMagNorm')

            if len(agnNames) != len(objectID):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d agnNames'
                                   % (len(agnNames), len(objectID)))
            if len(agnNames) != len(agnMagNorm):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d agnNames and %d agnMagNorms'
                                   % (len(agnNames), len(agnMagNorm)))

        if redshift is None:
            raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need redshift')

        if len(objectID) != len(redshift):
            raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d redshifts'
                               % (len(objectID), len(redshift)))


        if cosmologicalDistanceModulus is not None and len(objectID) != len(cosmologicalDistanceModulus):
            raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d cosmologicalDistanceModuli'
                               % (len(objectID), len(cosmologicalDistanceModulus)))

        diskMags = self.calculate_component_magnitudes(objectID,diskNames,magNorm = diskMagNorm, \
                        internalAv = diskAv, redshift = redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulus,
                        specFileMap=specFileMap, indices=indices)

        bulgeMags = self.calculate_component_magnitudes(objectID,bulgeNames,magNorm = bulgeMagNorm, \
                        internalAv = bulgeAv, redshift = redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulus,
                        specFileMap=specFileMap, indices=indices)

        agnMags = self.calculate_component_magnitudes(objectID,agnNames,magNorm = agnMagNorm, \
                        redshift = redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulus,
                        specFileMap=specFileMap, indices=indices)

        masterDict = {}

        for i in range(len(objectID)):
            total_mags=[]

            subDict={}
            subDict["bulge"] = bulgeMags[objectID[i]]
            subDict["disk"] = diskMags[objectID[i]]
            subDict["agn"] = agnMags[objectID[i]]

            masterDict[objectID[i]] = subDict


        return masterDict


    def meta_magnitudes_getter(self, objectID, columnNameList, indices=None):
        """
        This method will return the magnitudes for galaxies in the bandpasses stored in self.bandpassDict

        @param [in] objectID is a list of object IDs

        @param [in] columnNameList is a list of the name of the columns this method will ultimately
        be returning.  It exists so that it knows how to search for variability methods associated
        with those magnitudes.

        @param [in] indices is an optional list of indices indicating which bandpasses to actually
        calculate magnitudes for.  Note: even columns associated with bandpasses not included
        in indices should appear in columnNames above.
        """

        diskNames=self.column_by_name('sedFilenameDisk')
        bulgeNames=self.column_by_name('sedFilenameBulge')
        agnNames=self.column_by_name('sedFilenameAgn')

        diskmn = self.column_by_name('magNormDisk')
        bulgemn = self.column_by_name('magNormBulge')
        agnmn = self.column_by_name('magNormAgn')

        bulgeAv = self.column_by_name('internalAvBulge')
        diskAv = self.column_by_name('internalAvDisk')

        redshift = self.column_by_name('redshift')

        if 'cosmologicalDistanceModulus' in self.iter_column_names():
            cosmologicalDistanceModulus = self.column_by_name("cosmologicalDistanceModulus")
        else:
            cosmologicalDistanceModulus = None

        magDict=self.calculate_magnitudes(objectID,
                                          diskNames=diskNames, diskMagNorm=diskmn, diskAv=diskAv,
                                          bulgeNames=bulgeNames, bulgeMagNorm=bulgemn, bulgeAv=bulgeAv,
                                          agnNames=agnNames, agnMagNorm=agnmn,
                                          redshift=redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulus,
                                          specFileMap=self.specFileMap, indices=indices)

        failure = None

        outputBulge = None
        outputDisk = None
        outputAgn = None

        for i in range(self.nBandpasses):
            rowDisk = []
            rowBulge = []
            rowAgn = []

            for name in objectID:

                if magDict[name]["bulge"] is not None:
                    rowBulge.append(magDict[name]["bulge"][i])
                else:
                    rowBulge.append(failure)

                if magDict[name]["disk"] is not None:
                    rowDisk.append(magDict[name]["disk"][i])
                else:
                    rowDisk.append(failure)

                if magDict[name]["agn"] is not None:
                    rowAgn.append(magDict[name]["agn"][i])
                else:
                    rowAgn.append(failure)

            if outputBulge is None:
                outputBulge = numpy.array(rowBulge)
                outputDisk = numpy.array(rowDisk)
                outputAgn = numpy.array(rowAgn)
            else:
                outputBulge = numpy.vstack([outputBulge,rowBulge])
                outputDisk = numpy.vstack([outputDisk,rowDisk])
                outputAgn = numpy.vstack([outputAgn,rowAgn])


        #Add variability to the bulge components (if any)
        for ix, (columnName, columnData) in \
        enumerate(zip(columnNameList[self.nBandpasses:2*self.nBandpasses], outputBulge)):

            bandpassDex = ix % self.nBandpasses
            if indices is None or bandpassDex in indices:
                variabilityName = 'delta_' + columnName
                if variabilityName in self._all_available_columns:
                    delta = self.column_by_name(variabilityName)
                    columnData += delta

        #Add variability to the disk components (if any)
        for ix, (columnName, columnData) in \
        enumerate(zip(columnNameList[2*self.nBandpasses:3*self.nBandpasses], outputDisk)):

            bandpassDex = ix % self.nBandpasses
            if indices is None or bandpassDex in indices:
                variabilityName = 'delta_' + columnName
                if variabilityName in self._all_available_columns:
                    delta = self.column_by_name(variabilityName)
                    columnData += delta

        #Add variability to the agn components (if any)
        for ix, (columnName, columnData) in \
        enumerate(zip(columnNameList[3*self.nBandpasses:4*self.nBandpasses], outputAgn)):

            bandpassDex = ix % self.nBandpasses
            if indices is None or bandpassDex in indices:
                variabilityName = 'delta_' + columnName
                if variabilityName in self._all_available_columns:
                    delta = self.column_by_name(variabilityName)
                    columnData += delta


        #Calculate the total magnitude of the galaxy.
        #We do this here so that the variability models added above
        #have an influence on the total magnitude.
        outputTotal = None
        for ib in range(self.nBandpasses):
            if outputTotal is None:
                outputTotal = self.sum_magnitudes(bulge=outputBulge[ib],
                                                  disk=outputDisk[ib],
                                                  agn=outputAgn[ib])
            else:
                outputTotal = numpy.vstack([outputTotal,
                                            self.sum_magnitudes(bulge=outputBulge[ib],
                                                                disk=outputDisk[ib],
                                                                agn=outputAgn[ib])])


        #Add variability to the total components (if any).
        #This would be in the case that the catalog class is
        #only worried about total galaxy fluxes and thus only
        #adds variability to the whole galaxy, without worrying about
        #dividing it among the galaxy's components.
        #Adding variability to the components above and then adding variability
        #here is probably unphysical.
        for ix, (columnName, columnData) in \
        enumerate(zip(columnNameList[:self.nBandpasses], outputTotal)):

            bandpassDex = ix % self.nBandpasses
            if indices is None or bandpassDex in indices:
                variabilityName = 'delta_' + columnName
                if variabilityName in self._all_available_columns:
                    delta = self.column_by_name(variabilityName)
                    columnData += delta

        return numpy.vstack([outputTotal, outputBulge, outputDisk, outputAgn])





    @compound('sigma_lsst_u','sigma_lsst_g','sigma_lsst_r',
              'sigma_lsst_i','sigma_lsst_z','sigma_lsst_y')
    def get_photometric_uncertainties_total(self):
        """
        Getter for total photometric uncertainties associated with galaxies
        """
        magnitudes = numpy.array([self.column_by_name('lsst_u'),
                                  self.column_by_name('lsst_g'),
                                  self.column_by_name('lsst_r'),
                                  self.column_by_name('lsst_i'),
                                  self.column_by_name('lsst_z'),
                                  self.column_by_name('lsst_y')])

        return self.calculateMagnitudeUncertainty(magnitudes, obs_metadata=self.obs_metadata,
                                                    sigmaSysSq=self.sigmaSysSq)

    @compound('sigma_uBulge', 'sigma_gBulge', 'sigma_rBulge',
              'sigma_iBulge', 'sigma_zBulge', 'sigma_yBulge')
    def get_photometric_uncertainties_bulge(self):
        """
        Getter for photometric uncertainties associated with galaxy bulges
        """
        magnitudes = numpy.array([self.column_by_name('uBulge'),
                                  self.column_by_name('gBulge'),
                                  self.column_by_name('rBulge'),
                                  self.column_by_name('iBulge'),
                                  self.column_by_name('zBulge'),
                                  self.column_by_name('yBulge')])

        return self.calculateMagnitudeUncertainty(magnitudes, obs_metadata=self.obs_metadata,
                                                  sigmaSysSq=self.sigmaSysSq)

    @compound('sigma_uDisk', 'sigma_gDisk', 'sigma_rDisk',
              'sigma_iDisk', 'sigma_zDisk', 'sigma_yDisk')
    def get_photometric_uncertainties_disk(self):
        """
        Getter for photometeric uncertainties associated with galaxy disks
        """
        magnitudes = numpy.array([self.column_by_name('uDisk'),
                                  self.column_by_name('gDisk'),
                                  self.column_by_name('rDisk'),
                                  self.column_by_name('iDisk'),
                                  self.column_by_name('zDisk'),
                                  self.column_by_name('yDisk')])

        return self.calculateMagnitudeUncertainty(magnitudes, obs_metadata=self.obs_metadata,
                                                  sigmaSysSq=self.sigmaSysSq)

    @compound('sigma_uAgn', 'sigma_gAgn', 'sigma_rAgn',
              'sigma_iAgn', 'sigma_zAgn', 'sigma_yAgn')
    def get_photometric_uncertainties_agn(self):
        """
        Getter for photometric uncertainties associated with Agn
        """
        magnitudes = numpy.array([self.column_by_name('uAgn'),
                                  self.column_by_name('gAgn'),
                                  self.column_by_name('rAgn'),
                                  self.column_by_name('iAgn'),
                                  self.column_by_name('zAgn'),
                                  self.column_by_name('yAgn')])

        return self.calculateMagnitudeUncertainty(magnitudes, obs_metadata=self.obs_metadata,
                                                  sigmaSysSq=self.sigmaSysSq)

    @compound('lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
              'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
              'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
              'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn')
    def get_all_mags(self):
        """
        Getter for LSST galaxy magnitudes

        """
        objectID = self.column_by_name('galid')

        columnNames = [name for name in self.get_all_mags._colnames]

        """
        Here is where we need some code to load a list of bandpass objects
        into self.bandpassDict so that the bandpasses are available to the
        mixin.  Ideally, we would only do this once for the whole catalog
        """
        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandpassesFromFiles()

        indices = numpy.unique([ii % 6 for ii, name in enumerate(self.get_all_mags._colnames) \
                               if name in self._actually_calculated_columns])

        if len(indices)==6:
            indices=None

        return self.meta_magnitudes_getter(objectID, columnNames, indices=indices)



class PhotometryStars(PhotometryBase):
    """
    This mixin provides the infrastructure for doing photometry on stars

    It assumes that we want LSST filters.
    """

    def calculate_magnitudes(self, objectID, magNorm, sedNames, indices=None, specFileMap=None):
        """
        Take the bandpasses in bandpassDict and the array of
        star names objectID and return a dict of lists of magnitudes

        The first level key will be the name of the star (idName)

        This will give you a list of magnitudes corresponding to self.bandpassDict

        As with galaxies, it is important that we identify stars by a unique
        identifier, rather than their sedFilename, because different stars
        can have identical SEDs but different magnitudes.

        @param [in] objectID is a list of names uniquely identifying the objects being considered

        @param [in] magNorm is a list of magnitude normalizations

        @param [in] sedNames is a list of sed file names

        @param [in] indices is an optional list of indices indicating which
        bandpasses to actually calculate magnitudes for

        @param [in] specFileMap is a class which maps between sedNames and the absolute path to
        the SED files.  It is an instantiation of the class defined in

        sims_catalogs_measures/python/lsst/sims/catalogs/measures/instance/fileMaps.py

        if not provided, a default will be instantiated

        @param [out] magDict['AAA'][i] is the magnitude in the ith bandpass for object AAA

        """

        if specFileMap is None:
            if hasattr(self, 'specFileMap'):
                specFileMap=self.specFileMap
            else:
                specFileMap = defaultSpecMap

        if len(objectID) != len(magNorm) or len(objectID) != len(sedNames) or len(sedNames) != len(magNorm):
            raise RuntimeError('In PhotometryStars.calculate_magnitudes, had %d objectID, %d magNorms, and %d sedNames '
                                % (len(objectID), len(magNorm), len(sedNames)))

        sedList = self.loadSeds(sedNames, magNorm=magNorm, specFileMap=specFileMap)

        magDict = {}
        for (name,sed) in zip(objectID,sedList):
            subList = self.manyMagCalc_list(sed, indices=indices)
            magDict[name] = subList

        return magDict


    def meta_magnitudes_getter(self, objectID, columnNameList, indices=None):
        """
        This method does most of the work for stellar magnitude getters

        @param [in] objectID is a list of object names

        @param [in] columnNameList is a list of the names of the columns
        this method will ultimately be returning.

        @param [in] indices is an optional list indicating the indices of the
        bandpasses to calculate magnitudes for.  Note: even if a bandpass does
        not appear in indices, its columns should be listed in columnNames.

        @param [out] output is a 2d numpy array in which the rows are the bandpasses
        from bandpassDict and the columns are the objects from objectID

        """

        magNorm = self.column_by_name('magNorm')
        sedNames = self.column_by_name('sedFilename')
        magDict = self.calculate_magnitudes(objectID, magNorm=magNorm, sedNames=sedNames, indices=indices)
        output = None

        for i in range(self.nBandpasses):
            row = []
            for name in objectID:
                row.append(magDict[name][i])

            if output is None:
                output = numpy.array(row)
            else:
                output=numpy.vstack([output,row])

        for ix, (columnName, columnData) in enumerate(zip(columnNameList, output)):
            if indices is None or ix%self.nBandpasses in indices:
                deltaName = 'delta_' + columnName
                if deltaName in self._all_available_columns:
                    delta = self.column_by_name(deltaName)
                    columnData += delta

        return output

    @compound('sigma_lsst_u','sigma_lsst_g','sigma_lsst_r','sigma_lsst_i',
              'sigma_lsst_z','sigma_lsst_y')
    def get_photometric_uncertainties(self):
        """
        Getter for photometric uncertainties associated with stellar
        magnitudes
        """

        magnitudes = numpy.array([self.column_by_name('lsst_u'),
                                  self.column_by_name('lsst_g'),
                                  self.column_by_name('lsst_r'),
                                  self.column_by_name('lsst_i'),
                                  self.column_by_name('lsst_z'),
                                  self.column_by_name('lsst_y')])

        return self.calculateMagnitudeUncertainty(magnitudes, obs_metadata=self.obs_metadata,
                                                  sigmaSysSq=self.sigmaSysSq)


    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y')
    def get_magnitudes(self):
        """
        getter for LSST stellar magnitudes

        """
        objectID = self.column_by_name('id')

        columnNames = [name for name in self.get_magnitudes._colnames]

        """
        Here is where we need some code to load a list of bandpass objects
        into self.bandpassDict so that the bandpasses are available to the
        mixin.  Ideally, we would only do this once for the whole catalog
        """
        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandpassesFromFiles()

        indices = [ii for ii, name in enumerate(self.get_magnitudes._colnames) \
                   if name in self._actually_calculated_columns]

        if len(indices) == 6:
            indices = None

        return self.meta_magnitudes_getter(objectID, columnNames, indices=indices)

