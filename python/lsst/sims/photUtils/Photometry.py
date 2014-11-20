"""
photUtils -


ljones@astro.washington.edu  (and ajc@astro.washington.edu)

and now (2014 March 28): scott.f.daniel@gmail.com

Collection of utilities to aid usage of Sed and Bandpass with dictionaries.

"""

import os
import numpy
from lsst.sims.photUtils import Sed
from lsst.sims.photUtils import Bandpass
from lsst.sims.catalogs.measures.instance import compound

__all__ = ["PhotometryBase", "PhotometryGalaxies", "PhotometryStars"]

class PhotometryBase(object):
    """
    This mixin provides the basic infrastructure for photometry.
    It can read in SEDs and bandpasses, apply extinction and redshift, and, given
    an SED object it can calculate magnitudes.

    In order to avoid duplication of work, the bandPasses, wavelength array, and phi array
    are stored as instance variables once they are read in by self.loadBandPassesFromFiles()

    To initiailize a different set of bandPasses, call self.loadBandPassesFromFiles() with a different
    set of arguments.

    Once self.loadBandPassesFromFiles() as been called, self.loadSeds() can be used to return an array
    of SED objects.  These objects can be passed to self.manyMagCalc_list() which will calculate
    the magnitudes of the the SEDs, integrated over the loaded bandPasses, and return them as a
    dict keeyed to the array of bandpass keys stored in self.bandPassKey
    """

    bandPassList = None #bandpasses loaded in this particular catalog
    phiArray = None #the response curves for the bandpasses
    waveLenStep = None

    def setupPhiArray_dict(self):
        """
        Generate 2-dimensional numpy array for Phi values associated with the bandpasses in
        self.bandPasses

        The results from this calculation will be stored in the instance variables
        self.phiArray and self.waveLenStep for future use by self.manyMagCalc_list()
        """

        sedobj = Sed()
        self.phiArray, self.waveLenStep = sedobj.setupPhiArray(self.bandPassList)

    def loadBandPassesFromFiles(self,bandPassNames, bandPassDir = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline'),
    bandPassRoot = 'total_'):
        """
        This will take the list of band passes named by bandPassNames and use them to set up
        self.bandPassList (which is being cached so that
        it does not have to be loaded again unless we change which bandpasses we want)

        bandPassRoot contains the first part of the bandpass file name, i.e., it is assumed
        that the bandPasses are stored in files of the type

        $LSST_THROUGHPUTS_DEFAULT/bandPassRoot_bandPassList[i].dat

        if we want to load bandpasses for a telescope other than LSST, we would do so
        by altering bandPassRoot (currently no infrastructure exists for altering the directory
        in which bandpass files are stored)
        """

        self.bandPassList = []

        for w in bandPassNames:
            bandPassDummy = Bandpass()
            bandPassDummy.readThroughput(os.path.join(bandPassDir,"%s.dat" % (bandPassRoot + w)))
            self.bandPassList.append(bandPassDummy)

        self.phiArray = None
        self.waveLenStep = None

    # Handy routines for handling Sed/Bandpass routines with sets of dictionaries.
    def loadSeds(self, sedList, magNorm=15.0, resample_same=False):
        """
        Takes the list of filename sedList and returns an array of SED objects.

        This code will load identical SEDs twice because it is possible for
        (astronomical) objects to have the same SEDs but different magNorms

        @param [in] sedList is a list of file names containing Seds

        @param [in] magNorm is the magnitude normalization

        @param [in] resample_same governs whether or not to resample the Seds
        so that they are all on the same wavelength grid

        @param [out] sedOut is a list of Sed objects

        """

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
                sed.readSED_flambda(os.path.join(dataDir, self.specFileMap[sedName]))

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
            if sedList[i].wavelen != None:
                if internalAv != None:
                    #setupCCMab only depends on the wavelen array
                    #because this is supposed to be the same for every
                    #SED object in sedList, it is only called once for
                    #each invocation of applyAvAndRedshift
                    if wavelen_sampled == [] or (sedList[i].wavelen!=wavelen_sampled).any():
                        a_int, b_int = sedList[i].setupCCMab()
                        wavelen_sampled=sedList[i].wavelen

                    sedList[i].addCCMDust(a_int, b_int, A_v=internalAv[i])
                if redshift != None:
                    sedList[i].redshiftSED(redshift[i], dimming=True)
                    sedList[i].name = sedList[i].name + '_Z' + '%.2f' %(redshift[i])
                    sedList[i].resampleSED(wavelen_match=self.bandPassList[0].wavelen)

    def manyMagCalc_list(self, sedobj):
        """
        Return a list of magnitudes for a single Sed object.

        Bandpass information is taken from the instance variables self.bandPassList,
        self.phiArray, and self.waveLenStep

        @param [in] sedobj is an Sed object

        @param [out] magList is a list of magnitudes in the bandpasses stored in self.bandPassList
        """
        # Set up the SED for using manyMagCalc - note that this CHANGES sedobj
        # Have to check that the wavelength range for sedobj matches bandpass - this is why the dictionary is passed in.

        magList = []
        if sedobj.wavelen != None:
            sedobj.resampleSED(wavelen_match=self.bandPassList[0].wavelen)

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
            i = 0
            for f in self.bandPassList:
                magList.append(magArray[i])
                i = i + 1
        else:
            for f in self.bandPassList:
                magList.append(None)

        return magList

    def calculatePhotometricUncertaintyFromColumn(self, nameTag, columnNames):
        """
        This method reads in a dict of column names and passes out
        the associated photometric uncertainties.  The output will be
        a dict of lists.

        @param [in] nameTag is the name of the column used to identify each object

        @param [in] columnNames is a dict associating filter names with column names,
        e.g. columnName['u'] = 'lsst_u' if the u magnitude is stored in the column
        'lsst_u'

        @param [out] outputDict is a dict of lists such that outputDict['u'] is a list
        of the u band photometric uncertainties for all of the objects queried

        """

        inputDict={}

        idNames = self.column_by_name(nameTag)

        magnitudes = {}

        for filterName in columnNames:
            magnitudes[filterName] = self.column_by_name(columnNames[filterName])

        outputDict = self.calculatePhotometricUncertainty(magnitudes)

        return outputDict

    def calculatePhotometricUncertainty(self, magnitudes):
        """
        This method is based on equations 3.1, 3.2 and Table 3.2
        of the LSST Science Book (version 2.0)

        @param [in] magnitudes will be a dict of lists such that
        magnitudes['A'] will be a list of all the magnitudes in filter A

        @param [out] sigOut is a dict of lists such that sigOut['A'] is
        a list of the photometric uncertainties in filter A
        """
        sigma2Sys = 0.003*0.003 #also taken from the Science Book
                         #see the paragraph between equations 3.1 and 3.2

        gamma = {}

        gamma['u'] = 0.037
        gamma['g'] = 0.038
        gamma['r'] = 0.039
        gamma['i'] = 0.039
        gamma['z'] = 0.040
        gamma['y'] = 0.040

        sigOut={}

        for filterName in magnitudes:

            subList = []

            for i in range(len(magnitudes[filterName])):
                mm = magnitudes[filterName][i]

                if mm != None:

                    xx=10**(0.4*(mm - self.obs_metadata.m5(filterName)))
                    ss = (0.04 - gamma[filterName])*xx + \
                         gamma[filterName]*xx*xx

                    sigmaSquared = ss + sigma2Sys

                    subList.append(numpy.sqrt(sigmaSquared))

                else:
                    subList.append(None)

            sigOut[filterName] = subList

        return sigOut

class PhotometryGalaxies(PhotometryBase):
    """
    This mixin provides the code necessary for calculating the component magnitudes associated with
    galaxies.  It assumes that we want LSST filters.
    """

    def calculate_component_magnitudes(self,objectNames, componentNames, \
                                       magNorm = 15.0, internalAv = None, redshift = None):

        """
        Calculate the magnitudes for different components (disk, bulge, agn, etc) of galaxies.
        This method is designed to be used such that you feed it all of the disk Seds from your data
        base and it returns the associated magnitudes.  Then you feed it all of the bulge Seds, etc.

        @param [in] objectNames is the name of the galaxies (the whole galaxies)

        @param [in] componentNames gives the name of the SED filenames

        @param [in] magNorm is the normalizing magnitude

        @param [in] internalAv is the internal Av extinction

        @param [in] redshift is pretty self-explanatory

        @param [out] componentMags is a dict of lists such that
        magnitude["objectname"][i] will return the magnitude in the ith
        for the associated component Sed

        """

        componentMags = {}

        if componentNames != []:
            componentSed = self.loadSeds(componentNames, magNorm = magNorm)
            self.applyAvAndRedshift(componentSed, internalAv = internalAv, redshift = redshift)

            for i in range(len(objectNames)):
                subList = self.manyMagCalc_list(componentSed[i])
                componentMags[objectNames[i]] = subList

        else:
            subList=[]
            for b in self.bandPassList:
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

    def calculate_magnitudes(self, idNames):
        """
        Take the array of bandpasses in self.bandPassList and the array of galaxy
        names idNames ane return a dict of dicts of lists of magnitudes

        the first level key is galid (the name of the galaxy)

        the second level key is "total", "bulge", "disk", or "agn"

        this yields a list of magnitudes corresponding to the bandPasses in self.bandPassList

        We need to index the galaxies by some unique identifier, such as galid
        because it is possible for galaxies to have the same sed filenames but
        different normalizations

        @param [in] idNames is a list of names uniquely identifying the objects whose magnitudes
        are being calculated


        @param [out] masterDict is a dict of magnitudes such that
        masterDict['AAA']['BBB'][i] is the magnitude in the ith bandPass of component BBB of galaxy AAA


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

        diskMags = self.calculate_component_magnitudes(idNames,diskNames,magNorm = diskmn, \
                        internalAv = diskAv, redshift = redshift)

        bulgeMags = self.calculate_component_magnitudes(idNames,bulgeNames,magNorm = bulgemn, \
                        internalAv = bulgeAv, redshift = redshift)

        agnMags = self.calculate_component_magnitudes(idNames,agnNames,magNorm = agnmn, \
                        redshift = redshift)

        total_mags = []
        masterDict = {}

        for i in range(len(idNames)):
            total_mags=[]
            j=0
            for ff in self.bandPassList:
                total_mags.append(self.sum_magnitudes(disk = diskMags[idNames[i]][j],
                                bulge = bulgeMags[idNames[i]][j], agn = agnMags[idNames[i]][j]))

                j += 1

            subDict={}
            subDict["total"] = total_mags
            subDict["bulge"] = bulgeMags[idNames[i]]
            subDict["disk"] = diskMags[idNames[i]]
            subDict["agn"] = agnMags[idNames[i]]

            masterDict[idNames[i]] = subDict


        return masterDict


    def meta_magnitudes_getter(self, idNames):
        """
        This method will return the magnitudes for galaxies in the bandpasses stored in self.bandPassList

        @param [in] idNames is a list of object IDs

        """

        magDict=self.calculate_magnitudes(idNames)

        firstRowTotal = []
        firstRowDisk = []
        firstRowBulge = []
        firstRowAgn = []

        failure = None

        outputTotal = None
        outputBulge = None
        outputDisk = None
        outputAgn = None

        for i in range(len(self.bandPassList)):
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

        columnNames = {}
        columnNames['u'] = 'uRecalc'
        columnNames['g'] = 'gRecalc'
        columnNames['r'] = 'rRecalc'
        columnNames['i'] = 'iRecalc'
        columnNames['z'] = 'zRecalc'
        columnNames['y'] = 'yRecalc'

        totalDict = self.calculatePhotometricUncertaintyFromColumn('galid',columnNames)

        columnNames = {}
        columnNames['u'] = 'uDisk'
        columnNames['g'] = 'gDisk'
        columnNames['r'] = 'rDisk'
        columnNames['i'] = 'iDisk'
        columnNames['z'] = 'zDisk'
        columnNames['y'] = 'yDisk'

        diskDict = self.calculatePhotometricUncertaintyFromColumn('galid',columnNames)

        columnNames = {}
        columnNames['u'] = 'uBulge'
        columnNames['g'] = 'gBulge'
        columnNames['r'] = 'rBulge'
        columnNames['i'] = 'iBulge'
        columnNames['z'] = 'zBulge'
        columnNames['y'] = 'yBulge'

        bulgeDict = self.calculatePhotometricUncertaintyFromColumn('galid',columnNames)

        columnNames = {}
        columnNames['u'] = 'uAgn'
        columnNames['g'] = 'gAgn'
        columnNames['r'] = 'rAgn'
        columnNames['i'] = 'iAgn'
        columnNames['z'] = 'zAgn'
        columnNames['y'] = 'yAgn'

        agnDict = self.calculatePhotometricUncertaintyFromColumn('galid',columnNames)

        return numpy.array([totalDict['u'],totalDict['g'],totalDict['r'],
                            totalDict['i'],totalDict['z'],totalDict['y'],
                            bulgeDict['u'],bulgeDict['g'],bulgeDict['r'],
                            bulgeDict['i'],bulgeDict['z'],bulgeDict['y'],
                            diskDict['u'],diskDict['g'],diskDict['r'],
                            diskDict['i'],diskDict['z'],diskDict['y'],
                            agnDict['u'],agnDict['g'],agnDict['r'],
                            agnDict['i'],agnDict['z'],agnDict['y']])

    @compound('uRecalc', 'gRecalc', 'rRecalc', 'iRecalc', 'zRecalc', 'yRecalc',
              'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
              'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
              'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn')
    def get_all_mags(self):
        """
        Getter for LSST galaxy magnitudes

        """
        idNames = self.column_by_name('galid')
        bandPassNames = ['u','g','r','i','z','y']

        """
        Here is where we need some code to load a list of bandPass objects
        into self.bandPassList and then call self.setupPhiArray_dict()
        so that the bandPasses are available to the mixin.  Ideally, we
        would only do this once for the whole catalog
        """
        if self.bandPassList is None or self.phiArray is None:
            self.loadBandPassesFromFiles(bandPassNames)
            self.setupPhiArray_dict()

        return self.meta_magnitudes_getter(idNames)



class PhotometryStars(PhotometryBase):
    """
    This mixin provides the infrastructure for doing photometry on stars

    It assumes that we want LSST filters.
    """

    def calculate_magnitudes(self, sedList):
        """
        Take the array of bandpass keys bandPassList and the array of
        star names idNames and return a dict of lists of magnitudes

        The first level key will be the name of the star (idName)

        This will give you a list of magnitudes corresponding to self.bandPassList

        As with galaxies, it is important that we identify stars by a unique
        identifier, rather than their sedFilename, because different stars
        can have identical SEDs but different magnitudes.

        @param [in] sedList is a list of Sed objects corresponding to the objects
        for which you want to do photometry

        @param [out] magList will be a list of lists corresponding to the magnitudes
        of the objects you want output
        """

        magList = []
        for sed in sedList:
            subList = self.manyMagCalc_list(sed)
            magList.append(subList)

        return magList


    def meta_magnitudes_getter(self, idNames):
        """
        This method does most of the work for stellar magnitude getters

        @param [in] idNames is a list of object names

        @param [out] output is a 2d numpy array in which the rows are the bandpasses
        from bandPassList and the columns are the objects from idNames

        """

        sedNames = self.column_by_name('sedFilename')
        magNorm = self.column_by_name('magNorm')
        sedList = self.loadSeds(sedNames, magNorm=magNorm)

        magList = self.calculate_magnitudes(sedList)
        output = None

        for i in range(len(self.bandPassList)):
            row = []
            for iObject in range(len(idNames)):
                row.append(magList[iObject][i])

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

        columnNames = {}
        columnNames['u'] = 'lsst_u'
        columnNames['g'] = 'lsst_g'
        columnNames['r'] = 'lsst_r'
        columnNames['i'] = 'lsst_i'
        columnNames['z'] = 'lsst_z'
        columnNames['y'] = 'lsst_y'

        outputDict = self.calculatePhotometricUncertaintyFromColumn('id',columnNames)


        return numpy.array([outputDict['u'],outputDict['g'],outputDict['r'],
                            outputDict['i'],outputDict['z'],outputDict['y']])


    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y')
    def get_magnitudes(self):
        """
        getter for LSST stellar magnitudes

        """
        idNames = self.column_by_name('id')
        bandPassNames = ['u','g','r','i','z','y']

        """
        Here is where we need some code to load a list of bandPass objects
        into self.bandPassList and then call self.setupPhiArray_dict()
        so that the bandPasses are available to the mixin.  Ideally, we
        would only do this once for the whole catalog
        """
        if self.bandPassList is None or self.phiArray is None:
            self.loadBandPassesFromFiles(bandPassNames)
            self.setupPhiArray_dict()

        return self.meta_magnitudes_getter(idNames)

