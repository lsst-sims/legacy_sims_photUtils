import os
import numpy as np
import warnings
from lsst.sims.photUtils.Sed import Sed

class applyIGM(object):

    def __init__(self):

        self.zRange = np.arange(1.5, 3.0, 0.1)
        pathToTable = str(os.environ['SIMS_PHOTUTILS_DIR'] + '/python/lsst/sims/photUtils/' + 
                          'IGMLookupTables/')
        self.meanLookups = {}
        self.varLookups = {}
        for zValue in self.zRange:
            self.meanLookups[str(zValue)] = np.genfromtxt(str(pathToTable + 'MeanLookupTable_zSource' + 
                                                              str(zValue) + '.tbl'))
            self.varLookups[str(zValue)] = np.genfromtxt(str(pathToTable + 'VarLookupTable_zSource' + 
                                                             str(zValue) + '.tbl'))

    def applyIGM(self, redshift, sedobj):

        """Apply IGM extinction to already redshifted sed with redshift  
        between z=1.5-2.9 using transmission lookup tables provided by Alex Abate"""

        #First make sure redshift is in range of lookup tables.            
        if (redshift < 1.5) or (redshift > 2.9):
            warnings.warn("IGM Lookup tables only applicable for 1.5 < z < 2.9. No action taken")
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
        dzGrid = 0.1 #Step in redshift between transmission lookup table files 
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
