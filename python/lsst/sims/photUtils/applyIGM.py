import os
import numpy as np
import warnings
from .Sed import Sed
from lsst.utils import getPackageDir

__all__ = ["ApplyIGM"]

class ApplyIGM(object):

    """
    This class applies IGM to SED objects using lookup tables. If users want to enter their
    own lookup tables they can do that by specifying self.meanLookups and self.varLookups
    which are dictionaries containing redshift as the keys with (wavelength, transmission)
    arrays as the values.
    """

    IGMisInitialized = False
    tablesPresent = False

    def initializeIGM(self, zMin = 1.5, zMax = 2.9, zDelta = 0.1, minWavelen = 300):

        """
        Initialize an applyIGM object with the desired redshift grid.
        If lookup tables are not evenly spaced in redshift then input manually
        desired zRange array.

        @param [in] zMin is the minimum redshift.

        @param [in] zMax is the maximum redshift.

        @param [in] zDelta is the redshift spacing.

        @param [in] minWavelen is the minimum wavelength in the lookup tables
        """
        self.zMin = zMin
        self.zMax = zMax
        self.zDelta = zDelta
        self.minWavelen = minWavelen
        #Don't have max wavelength since transmission goes to 1.0 at longest wavelengths
        self.zRange = np.arange(zMin, zMax + (zDelta/2.), zDelta)

        if self.tablesPresent == False:
            table_dir = getPackageDir('sims_photUtils')
            table_dir = os.path.join(table_dir,
                                     'python/lsst/sims/photUtils/igm_tables')
            self.loadTables(table_dir)

        self.IGMisInitialized = True


    def loadTables(self, filesDir, varianceTbl = True):

        """
        Read in and store in dictionary the IGM Lookup Tables that contain IGM transmission
        for a given redshift and must be formatted in two columns:
        (wavelength (nm), IGM Transmission %) or for variance
        (wavelength (nm), IGM Transmission % Variance). Variance tables are not required and
        can be turned off as a requirement. Names in directory formatted as
        'MeanLookupTable_zSourceX.X.tbl' or 'VarLookupTable_zSourceX.X.tbl' where X.X is the redshift
        of the given lookup table.

        @param [in] filesDir is the location of the directory where lookup table are stored

        @param [in] varianceTbl is a boolean that is True if variance tables are present in dir
        for loading.
        """

        self.meanLookups = {}
        self.varLookups = {}

        for zValue in self.zRange:
            self.meanLookups[str(zValue)] = np.genfromtxt(str(filesDir + '/MeanLookupTable_zSource' +
                                                              str(zValue) + '.tbl.gz'))
            if varianceTbl == True:
                try:
                    self.varLookups[str(zValue)] = np.genfromtxt(str(filesDir + '/VarLookupTable_zSource' +
                                                                     str(zValue) + '.tbl.gz'))
                except IOError:
                    raise IOError("Cannot find variance tables.")

        self.tablesPresent = True

    def applyIGM(self, redshift, sedobj):

        """
        Apply IGM extinction to already redshifted sed with redshift
        between zMin and zMax defined by range of lookup tables

        @param [in] redshift is the redshift of the incoming SED object

        @param [in] sedobj is the SED object to which IGM extinction will be applied. This object
        will be modified as a result of this.
        """

        if self.IGMisInitialized == False:
            self.initializeIGM()

        #First make sure redshift is in range of lookup tables.
        if (redshift < self.zMin) or (redshift > self.zMax):
            warnings.warn(str("IGM Lookup tables only applicable for " + str(self.zMin) + " < z < " + str(self.zMax) + ". No action taken"))
            return

        #Now read in closest two lookup tables for given redshift
        lowerSed = Sed()
        upperSed = Sed()
        for lower, upper in zip(self.zRange[:-1], self.zRange[1:]):
            if lower <= redshift <= upper:
                lowerSed.setSED(self.meanLookups[str(lower)][:,0],
                                flambda = self.meanLookups[str(lower)][:,1])
                upperSed.setSED(self.meanLookups[str(upper)][:,0],
                                flambda = self.meanLookups[str(upper)][:,1])
                break

        #Redshift lookup tables to redshift of source, i.e. if source redshift is 1.78 shift lookup
        #table for 1.7 and lookup table for 1.8 to up and down to 1.78, respectively
        zLowerShift = ((1.0 + redshift)/(1.0 + lower)) - 1.0
        zUpperShift = ((1.0 + redshift)/(1.0 + upper)) - 1.0
        lowerSed.redshiftSED(zLowerShift)
        upperSed.redshiftSED(zUpperShift)

        #Resample lower and upper transmission data onto same wavelength grid.
        minWavelen = 300. #All lookup tables are usable above 300nm
        maxWavelen = np.amin([lowerSed.wavelen[-1], upperSed.wavelen[-1]]) - 0.01
        lowerSed.resampleSED(wavelen_min = minWavelen, wavelen_max = maxWavelen, wavelen_step = 0.01)
        upperSed.resampleSED(wavelen_match = lowerSed.wavelen)

        #Now insert this into a transmission array of 1.0 beyond the limits of current application
        #So that we can get an sed back that extends to the longest wavelengths of the incoming sed
        finalWavelen = np.arange(300., sedobj.wavelen[-1]+0.01, 0.01)
        finalFlambdaExtended = np.ones(len(finalWavelen))

        #Weighted Average of Transmission from each lookup table to get final transmission
        #table at desired redshift
        dzGrid = self.zDelta #Step in redshift between transmission lookup table files
        finalSed = Sed()
        finalFlambda = (lowerSed.flambda*(1.0 - ((redshift - lower)/dzGrid)) +
                        upperSed.flambda*(1.0 - ((upper - redshift)/dzGrid)))
        finalFlambdaExtended[0:len(finalFlambda)] = finalFlambda
        finalSed.setSED(wavelen = finalWavelen, flambda = finalFlambdaExtended)

        #Resample incoming sed to new grid so that we don't get warnings from multiplySED
        #about matching wavelength grids
        sedobj.resampleSED(wavelen_match=finalSed.wavelen)

        #Now multiply transmission curve by input SED to get final result and make it the new flambda
        #data in the original sed which also is now on a new grid starting at 300 nm
        test = sedobj.multiplySED(finalSed)
        sedobj.flambda = test.flambda
