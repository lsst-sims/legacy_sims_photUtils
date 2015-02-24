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
from lsst.sims.photUtils import Sed, Bandpass, PhotometricDefaults
from lsst.sims.catalogs.measures.instance import defaultSpecMap
from lsst.sims.catalogs.measures.instance import compound

__all__ = ["PhotometryBase", "PhotometryGalaxies", "PhotometryStars"]

class PhotometryBase(object):
    """
    This mixin provides the basic infrastructure for photometry.
    It can read in SEDs and bandpasses, apply extinction and redshift, and, given
    an SED object it can calculate magnitudes.

    In order to avoid duplication of work, the bandPasses, wavelength array, and phi array
    are stored as instance variables once they are read in by self.loadTotalBandPassesFromFiles()

    To initiailize a different set of bandPasses, call self.loadBandPassesFromFiles() with a different
    set of arguments.

    Once self.loadBandPassesFromFiles() as been called, self.loadSeds() can be used to return an array
    of SED objects.  These objects can be passed to self.manyMagCalc_list() which will calculate
    the magnitudes of the the SEDs, integrated over the loaded bandPasses, and return them as a
    dict keeyed to the array of bandpass keys stored in self.bandPassKey
    """

    bandpassDict = None #total (hardware plus atmosphere) bandpasses loaded in this particular catalog;
                        #will be an OrderedDict
    hardwareBandpassDict = None #bandpasses associated with just the telescope hardware
    atmosphereBandpass = None #bandpass of the atmosphere
    nBandpasses = 0 #the number of bandpasses loaded (will need to be zero for InstanceCatalog's dry run)
    phiArray = None #the response curves for the bandpasses
    waveLenStep = None

    gammaDict = None
    sig2sys = None
    sedList = None

    atmosphereSED = None #the emission spectrum of the atmosphere
    atmosphereSEDnormalized = None #this will be a tuple of (bandpassName, magnitude) indicating how
                                   #atmosphereSED was normalized

    def setupPhiArray_dict(self):
        """
        Generate 2-dimensional numpy array for Phi values associated with the bandpasses in
        self.bandPasses

        The results from this calculation will be stored in the instance variables
        self.phiArray and self.waveLenStep for future use by self.manyMagCalc_list()
        """

        sedobj = Sed()
        if self.bandpassDict is not None:
            self.phiArray, self.waveLenStep = sedobj.setupPhiArray(self.bandpassDict.values())

    def loadBandpassesFromFiles(self, bandpassNames=['u', 'g', 'r', 'i', 'z', 'y'],
                                filedir = os.path.join(eups.productDir('throughputs'), 'baseline'),
                                bandpassRoot = 'filter_',
                                componentList = ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                                                 'lens1.dat', 'lens2.dat', 'lens3.dat'],
                                atmosphereBandpass='atmos.dat', atmosphereSED='darksky.dat'):

        self.atmosphereSED = Sed()
        self.atmosphereSED.readSED_flambda(os.path.join(filedir, atmosphereSED))
        self.atmosphereSEDnormalized = None

        self.atmosphereBandpass = Bandpass()
        self.atmosphereBandpass.readThroughput(os.path.join(filedir, atmosphereBandpass))

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
            components += [os.path.join(filedir, atmosphereBandpass)]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.bandpassDict[w] = bandpassDummy

        self.phiArray = None
        self.waveLenStep = None
        self.nBandpasses = len(self.bandpassDict)
        self.setupPhiArray_dict()


    def loadTotalBandPassesFromFiles(self,bandPassNames=['u', 'g', 'r', 'i', 'z', 'y'],
                                bandPassDir = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline'),
                                bandPassRoot = 'total_'):
        """
        This will take the list of band passes named by bandPassNames and use them to set up
        self.bandpassDict (which is being cached so that
        it does not have to be loaded again unless we change which bandpasses we want)

        The bandpasses loaded this way are total bandpasses: they account for instrumental
        and atmospheric transmission.

        @param [in] bandPassNames is a list of names identifying each filter.
        Defaults to ['u', 'g', 'r', 'i', 'z', 'y']

        @param [in] bandPassDir is the name of the directory where the bandpass files are stored

        @param [in] bandPassRoot contains the first part of the bandpass file name, i.e., it is assumed
        that the bandPasses are stored in files of the type

        bandPassDir/bandPassRoot_bandPassNames[i].dat

        if we want to load bandpasses for a telescope other than LSST, we would do so
        by altering bandPassDir and bandPassRoot
        """

        self.bandpassDict = OrderedDict()
        self.gammaDict = None

        for w in bandPassNames:
            bandPassDummy = Bandpass()
            bandPassDummy.readThroughput(os.path.join(bandPassDir,"%s.dat" % (bandPassRoot + w)))
            self.bandpassDict[w] = bandPassDummy

        self.nBandpasses = len(self.bandpassDict)
        self.phiArray = None
        self.waveLenStep = None
        self.setupPhiArray_dict()

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

        self.sedList = sedOut
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

    def manyMagCalc_list(self, sedobj):
        """
        Return a list of magnitudes for a single Sed object.

        Bandpass information is taken from the instance variables self.bandpassDict,
        self.phiArray, and self.waveLenStep

        @param [in] sedobj is an Sed object

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

            magArray = sedobj.manyMagCalc(self.phiArray, self.waveLenStep)

            for i in range(self.nBandpasses):
                magList.append(magArray[i])
        else:
            for i in range(self.nBandpasses):
                magList.append(None)

        return magList

    def normalizeAtmosphereSED(self, obs_metadata):
        fNorm = self.atmosphereSED.calcFluxNorm(obs_metadata.skyBrightness, self.hardwareBandpassDict[obs_metadata.bandpass])
        self.atmosphereSED.multiplyFluxNorm(fNorm)

    def calculatePhotometricUncertainty(self, magnitudes, obs_metadata=None, sig2sys=None):

        if obs_metadata is None:
            raise RuntimeError("Need to pass an ObservationMetaData into calculatePhotometricUncertainty")

        if len(magnitudes) != self.nBandpasses:
            raise RuntimeError("Passed %d magnitudes to " %len(magnitudes) + \
                                " PhotometryBase.calculatePhotometricUncertainty; " + \
                                "needed %d " % self.nBandpasses)

        if self.gammaDict is None or len(self.gammaDict) != self.nBandpasses:
            self.gammaDict = OrderedDict()
            for i in range(self.nBandpasses):
                b = self.bandpassDict.keys()[i]
                self.gammaDict[b] = self.bandpassDict[b].calcGamma(obs_metadata.m5(b),
                                                                   expTime=PhotometricDefaults.exptime,
                                                                   nexp=PhotometricDefaults.nexp,
                                                                   gain=PhotometricDefaults.gain,
                                                                   effarea=PhotometricDefaults.effarea)

        sigma = numpy.zeros(self.nBandpasses, numpy.float64)

        for i in range(self.nBandpasses):
            gg = self.gammaDict.values()[i]
            xx = numpy.power(10.0,0.4*(magnitudes[i] - obs_metadata.m5(self.bandpassDict.keys()[i])))
            sigmaSq = (0.04-gg)*xx+gg*xx*xx
            if sig2sys is None:
                noiseOverSignal = numpy.sqrt(sigmaSq)
            else:
                noiseOverSignal = numpy.sqrt(sigmaSq + sig2sys)

            #see www.ucolick.org/~bolte/AY257/s_n.pdf section 3.1
            sigma[i] = 2.5*numpy.log10(1.0+noiseOverSignal)

        return sigma

    def calculatePhotometricUncertaintyFromColumnNames(self, magnitudeColumnNames):
        magnitudes = []
        for cc in magnitudeColumnNames:
            magnitudes.append(self.column_by_name(cc));

        rows = len(magnitudes[0])
        cols = len(magnitudeColumnNames)

        output = numpy.zeros((cols, rows), numpy.float64)
        for i in range(rows):
            magDummy = []

            for j in range(cols):
                magDummy.append(magnitudes[j][i])
            sigDummy = self.calculatePhotometricUncertainty(magDummy, obs_metadata=self.obs_metadata, sig2sys=self.sig2sys)

            for j in range(cols):
                output[j][i] = sigDummy[j]

        return output

class PhotometryGalaxies(PhotometryBase):
    """
    This mixin provides the code necessary for calculating the component magnitudes associated with
    galaxies.  It assumes that we want LSST filters.
    """

    def calculate_component_magnitudes(self,objectNames, componentNames, \
                                       magNorm = None, internalAv = None, redshift = None,
                                       cosmologicalDistanceModulus = None, specFileMap=None):

        """
        Calculate the magnitudes for different components (disk, bulge, agn, etc) of galaxies.
        This method is designed to be used such that you feed it all of the disk Seds from your data
        base and it returns the associated magnitudes.  Then you feed it all of the bulge Seds, etc.

        @param [in] objectNames is the name of the galaxies (the whole galaxies)

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

        @param [out] componentMags is a dict of lists such that
        magnitude["objectname"][i] will return the magnitude in the ith
        for the associated component Sed

        """


        componentMags = {}

        if componentNames != [] and componentNames is not None:
            componentSed = self.loadSeds(componentNames, magNorm = magNorm, specFileMap=specFileMap)
            self.applyAvAndRedshift(componentSed, internalAv = internalAv, redshift = redshift)

            for i in range(len(objectNames)):
                subList = self.manyMagCalc_list(componentSed[i])

                if isinstance(cosmologicalDistanceModulus, numpy.ndarray):
                    for j in range(len(subList)):
                        subList[j] += cosmologicalDistanceModulus[i]

                componentMags[objectNames[i]] = subList

        else:
            subList=[]
            for i in range(self.nBandpasses):
                subList.append(None)
            for i in range(len(objectNames)):
                componentMags[objectNames[i]]=subList

        return componentMags

    def sum_magnitudes(self, disk = None, bulge = None, agn = None):
        """
        Sum the component magnitudes of a galaxy and return the answer

        @param [in] disk is the disk magnitude

        @param [in] bulge is the bulge magnitude

        @param [in] agn is the agn magnitude

        @param [out] outMag is the total magnitude of the galaxy
        """

        mm_o = 22.

        nn=0.0
        if disk is not None and (not numpy.isnan(disk)):
            nn+=numpy.power(10, (disk - mm_o)/-2.5)

        if bulge is not None and (not numpy.isnan(bulge)):
            nn+=numpy.power(10, (bulge - mm_o)/-2.5)

        if agn is not None and (not numpy.isnan(agn)):
            nn+=numpy.power(10, (agn - mm_o)/-2.5)

        if nn>0.0:
            outMag = -2.5*numpy.log10(nn) + mm_o
        else:
            outMag = None

        return outMag

    def calculate_magnitudes(self, idNames, diskNames=None, diskMagNorm=None, diskAv=None,
                             bulgeNames=None, bulgeMagNorm=None, bulgeAv=None,
                             agnNames=None, agnMagNorm=None,
                             redshift=None, cosmologicalDistanceModulus=None, specFileMap=None):
        """
        Take the array of bandpasses in self.bandpassDict and the array of galaxy
        names idNames ane return a dict of dicts of lists of magnitudes

        the first level key is galid (the name of the galaxy)

        the second level key is "total", "bulge", "disk", or "agn"

        this yields a list of magnitudes corresponding to the bandPasses in self.bandpassDict

        We need to index the galaxies by some unique identifier, such as galid
        because it is possible for galaxies to have the same sed filenames but
        different normalizations

        @param [in] idNames is a list of names uniquely identifying the objects whose magnitudes
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

        @param [out] masterDict is a dict of magnitudes such that
        masterDict['AAA']['BBB'][i] is the magnitude in the ith bandPass of component BBB of galaxy AAA


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

            if len(diskNames) != len(idNames):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d diskNames'
                                   % (len(diskNames), len(idNames)))
            if len(diskNames) != len(diskAv) or len(diskNames) != len(diskMagNorm) or len(diskMagNorm) != len(diskAv):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d diskNames, %d diskAvs, and %d diskMagNorms'
                                   % (len(diskNames), len(diskAv), len(diskMagNorm)))

        if bulgeNames is not None:
            if bulgeAv is None:
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need bulgeAv')

            if bulgeMagNorm is None:
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need bulgeMagNorm')

            if len(bulgeNames) != len(idNames):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d bulgeNames'
                                   % (len(bulgeNames), len(idNames)))
            if len(bulgeNames) != len(bulgeAv) or len(bulgeNames) != len(bulgeMagNorm) or len(bulgeMagNorm) != len(bulgeAv):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d bulgeNames, %d bulgeAvs, and %d bulgeMagNorms'
                                   % (len(bulgeNames), len(bulgeAv), len(bulgeMagNorm)))

        if agnNames is not None:
            if agnMagNorm is None:
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need agnMagNorm')

            if len(agnNames) != len(idNames):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d agnNames'
                                   % (len(agnNames), len(idNames)))
            if len(agnNames) != len(agnMagNorm):
                raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d agnNames and %d agnMagNorms'
                                   % (len(agnNames), len(agnMagNorm)))

        if redshift is None:
            raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes need redshift')

        if len(idNames) != len(redshift):
            raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d redshifts'
                               % (len(idNames), len(redshift)))


        if cosmologicalDistanceModulus is not None and len(idNames) != len(cosmologicalDistanceModulus):
            raise RuntimeError('In PhotometryGalaxies.calculate_magnitudes have %d galaxies and %d cosmologicalDistanceModuli'
                               % (len(idNames), len(cosmologicalDistanceModulus)))

        diskMags = self.calculate_component_magnitudes(idNames,diskNames,magNorm = diskMagNorm, \
                        internalAv = diskAv, redshift = redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulus,
                        specFileMap=specFileMap)

        bulgeMags = self.calculate_component_magnitudes(idNames,bulgeNames,magNorm = bulgeMagNorm, \
                        internalAv = bulgeAv, redshift = redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulus,
                        specFileMap=specFileMap)

        agnMags = self.calculate_component_magnitudes(idNames,agnNames,magNorm = agnMagNorm, \
                        redshift = redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulus,
                        specFileMap=specFileMap)

        total_mags = []
        masterDict = {}

        for i in range(len(idNames)):
            total_mags=[]

            for j in range(self.nBandpasses):
                total_mags.append(self.sum_magnitudes(disk = diskMags[idNames[i]][j],
                                bulge = bulgeMags[idNames[i]][j], agn = agnMags[idNames[i]][j]))

            subDict={}
            subDict["total"] = total_mags
            subDict["bulge"] = bulgeMags[idNames[i]]
            subDict["disk"] = diskMags[idNames[i]]
            subDict["agn"] = agnMags[idNames[i]]

            masterDict[idNames[i]] = subDict


        return masterDict


    def meta_magnitudes_getter(self, idNames):
        """
        This method will return the magnitudes for galaxies in the bandpasses stored in self.bandpassDict

        @param [in] idNames is a list of object IDs

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

        magDict=self.calculate_magnitudes(idNames,
                                          diskNames=diskNames, diskMagNorm=diskmn, diskAv=diskAv,
                                          bulgeNames=bulgeNames, bulgeMagNorm=bulgemn, bulgeAv=bulgeAv,
                                          agnNames=agnNames, agnMagNorm=agnmn,
                                          redshift=redshift, cosmologicalDistanceModulus=cosmologicalDistanceModulus,
                                          specFileMap=self.specFileMap)

        firstRowTotal = []
        firstRowDisk = []
        firstRowBulge = []
        firstRowAgn = []

        failure = None

        outputTotal = None
        outputBulge = None
        outputDisk = None
        outputAgn = None

        for i in range(self.nBandpasses):
            rowTotal = []
            rowDisk = []
            rowBulge = []
            rowAgn = []

            for name in idNames:
                rowTotal.append(magDict[name]["total"][i])

                if magDict[name]["bulge"]:
                    rowBulge.append(magDict[name]["bulge"][i])
                else:
                    rowBulge.append(failure)

                if magDict[name]["disk"]:
                    rowDisk.append(magDict[name]["disk"][i])
                else:
                    rowDisk.append(failure)

                if magDict[name]["agn"]:
                    rowAgn.append(magDict[name]["agn"][i])
                else:
                    rowAgn.append(failure)

            if outputTotal is None:
                outputTotal = numpy.array(rowTotal)
                outputBulge = numpy.array(rowBulge)
                outputDisk = numpy.array(rowDisk)
                outputAgn = numpy.array(rowAgn)
            else:
                outputTotal = numpy.vstack([outputTotal,rowTotal])
                outputBulge = numpy.vstack([outputBulge,rowBulge])
                outputDisk = numpy.vstack([outputDisk,rowDisk])
                outputAgn = numpy.vstack([outputAgn,rowAgn])



        outputTotal = numpy.vstack([outputTotal,outputBulge])
        outputTotal = numpy.vstack([outputTotal,outputDisk])
        outputTotal = numpy.vstack([outputTotal,outputAgn])

        return outputTotal




    @compound('sigma_uRecalc','sigma_gRecalc','sigma_rRecalc',
              'sigma_iRecalc','sigma_zRecalc','sigma_yRecalc',
              'sigma_uBulge','sigma_gBulge','sigma_rBulge',
              'sigma_iBulge','sigma_zBulge','sigma_yBulge',
              'sigma_uDisk','sigma_gDisk','sigma_rDisk',
              'sigma_iDisk','sigma_zDisk','sigma_yDisk',
              'sigma_uAgn','sigma_gAgn','sigma_rAgn',
              'sigma_iAgn','sigma_zAgn','sigma_yAgn')
    def get_photometric_uncertainties(self):
        """
        Getter for photometric uncertainties associated with galaxies
        """

        columnNames = ['uRecalc', 'gRecalc', 'rRecalc', 'iRecalc', 'zRecalc', 'yRecalc']
        output = self.calculatePhotometricUncertaintyFromColumnNames(columnNames)

        columnNames = ['uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge']
        output = numpy.vstack([output, self.calculatePhotometricUncertaintyFromColumnNames(columnNames)])

        columnNames = ['uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk']
        output = numpy.vstack([output, self.calculatePhotometricUncertaintyFromColumnNames(columnNames)])

        columnNames = ['uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn']
        output = numpy.vstack([output, self.calculatePhotometricUncertaintyFromColumnNames(columnNames)])

        return output

    @compound('uRecalc', 'gRecalc', 'rRecalc', 'iRecalc', 'zRecalc', 'yRecalc',
              'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
              'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
              'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn')
    def get_all_mags(self):
        """
        Getter for LSST galaxy magnitudes

        """
        idNames = self.column_by_name('galid')

        """
        Here is where we need some code to load a list of bandPass objects
        into self.bandpassDict so that the bandPasses are available to the
        mixin.  Ideally, we would only do this once for the whole catalog
        """
        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandPassesFromFiles()

        return self.meta_magnitudes_getter(idNames)



class PhotometryStars(PhotometryBase):
    """
    This mixin provides the infrastructure for doing photometry on stars

    It assumes that we want LSST filters.
    """

    def calculate_magnitudes(self, idNames, magNorm, sedNames, specFileMap=None):
        """
        Take the bandpasses in bandpassDict and the array of
        star names idNames and return a dict of lists of magnitudes

        The first level key will be the name of the star (idName)

        This will give you a list of magnitudes corresponding to self.bandpassDict

        As with galaxies, it is important that we identify stars by a unique
        identifier, rather than their sedFilename, because different stars
        can have identical SEDs but different magnitudes.

        @param [in] idNames is a list of names uniquely identifying the objects being considered

        @param [in] magNorm is a list of magnitude normalizations

        @param [in] sedNames is a list of sed file names

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

        if len(idNames) != len(magNorm) or len(idNames) != len(sedNames) or len(sedNames) != len(magNorm):
            raise RuntimeError('In PhotometryStars.calculate_magnitudes, had %d idnames, %d magNorms, and %d sedNames '
                                % (len(idNames), len(magNorm), len(sedNames)))

        sedList = self.loadSeds(sedNames, magNorm=magNorm, specFileMap=specFileMap)

        magDict = {}
        for (name,sed) in zip(idNames,sedList):
            subList = self.manyMagCalc_list(sed)
            magDict[name] = subList

        return magDict


    def meta_magnitudes_getter(self, idNames):
        """
        This method does most of the work for stellar magnitude getters

        @param [in] idNames is a list of object names

        @param [out] output is a 2d numpy array in which the rows are the bandpasses
        from bandpassDict and the columns are the objects from idNames

        """

        magNorm = self.column_by_name('magNorm')
        sedNames = self.column_by_name('sedFilename')
        magDict = self.calculate_magnitudes(idNames, magNorm=magNorm, sedNames=sedNames)
        output = None

        for i in range(self.nBandpasses):
            row = []
            for name in idNames:
                row.append(magDict[name][i])

            if output is None:
                output = numpy.array(row)
            else:
                output=numpy.vstack([output,row])

        return output

    @compound('sigma_lsst_u','sigma_lsst_g','sigma_lsst_r','sigma_lsst_i',
              'sigma_lsst_z','sigma_lsst_y')
    def get_photometric_uncertainties(self):
        """
        Getter for photometric uncertainties associated with stellar
        magnitudes
        """

        magnitudeColumnNames = ['lsst_u', 'lsst_g', 'lsst_r', 'lsst_i',
                                'lsst_z', 'lsst_y']

        return self.calculatePhotometricUncertaintyFromColumnNames(magnitudeColumnNames)


    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y')
    def get_magnitudes(self):
        """
        getter for LSST stellar magnitudes

        """
        idNames = self.column_by_name('id')

        """
        Here is where we need some code to load a list of bandPass objects
        into self.bandpassDict so that the bandPasses are available to the
        mixin.  Ideally, we would only do this once for the whole catalog
        """
        if self.bandpassDict is None or self.phiArray is None:
            self.loadTotalBandPassesFromFiles()

        return self.meta_magnitudes_getter(idNames)

