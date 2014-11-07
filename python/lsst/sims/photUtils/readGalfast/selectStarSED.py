import os
import gzip
import pyfits
import numpy as np
import re

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.photUtils import Photometry as phot
from lsst.sims.catalogs.measures.instance.fileMaps import SpecMap

__all__ = ["selectStarSED"]

class selectStarSED():

    """
    This class provides a way to match incoming star colors to those of the approriate SED.

    In addition, provides ways to load in SED colors.

    Can easily be used for bandpasses other than sdss by redefining self.filterList, 
    self.bandpassDict, self.phiArray, self.wavelenstep before using other methods after initalizing.
    """

    def __init__(self):
        self.sEDDir = os.environ['SIMS_SED_LIBRARY_DIR']
        #Use SpecMap to pull the directory locations
        specMap = SpecMap()
        specMapDict = {}
        specFileStart = ['kp', 'burrows', 'bergeron'] #The beginning of filenames of different SED types
        specFileTypes = ['kurucz', 'mlt', 'wd']
        for specStart, specKey in zip(specFileStart, specFileTypes):
            for key, val in sorted(specMap.subdir_map.iteritems()):
                if re.match(key, specStart):
                    specMapDict[specKey] = str(val)
        
        self.kuruczDir = str(self.sEDDir + '/' + specMapDict['kurucz'] + '/')
        self.mltDir = str(self.sEDDir + '/' + specMapDict['mlt'] + '/')
        self.wdDir = str(self.sEDDir + '/' + specMapDict['wd'] + '/')

        #Load Bandpasses for SDSS colors to match to galfast output.
        #If somebody wants to use this with bandpasses other than sdss all they have to do is redefine
        #self.filterList, self.bandpassDict, self.phiArray self.wavelenstep before using other methods
        starPhot = phot()
        self.filterList = ('u', 'g', 'r', 'i', 'z')
        self.bandpassDict = starPhot.loadBandpasses(filterlist=self.filterList, 
                                                    dataDir = os.getenv("SDSS_THROUGHPUTS"), 
                                                    filterroot='sdss_')
        self.phiArray, self.wavelenstep = starPhot.setupPhiArray_dict(self.bandpassDict, 
                                                                      self.filterList)

        #Load Bandpasses for LSST colors to get colors from matched SEDs
        lsstPhot = phot()
        self.lsstFilterList = ('u', 'g', 'r', 'i', 'z', 'y')
        self.lsstBandpassDict = lsstPhot.loadBandpasses()
        self.lsstPhiArray, self.lsstWavelenstep = lsstPhot.setupPhiArray_dict(self.lsstBandpassDict, 
                                                                              self.lsstFilterList)

    def loadKuruczSEDs(self, subset = None):
        """
        By default will load all seds in kurucz directory, otherwise a subset can be loaded if a 
        list of filenames of only seds wanted in the folder is provided. Will skip over extraneous
        files in sed folder.

        @param [in] subset is the list of the subset of files wanted if one doesn't want all files
        in the kurucz directory.

        @param [out] kDict is a dictionary with all the necessary kurucz information to match SEDs

        """

        files = []
        #Load Kurucz
        kDict = {}
        
        if subset is None:
            for fileName in os.listdir(self.kuruczDir):
                files.append(fileName)
        else:
            for fileName in subset:
                files.append(fileName)
        
        kFiles = []; klogZ = []; kTemp = []; klogg = []
        kumg = []; kgmr = []; krmi = []; kimz = []; krMag = {}
        lsstumg = {}; lsstgmr = {}; lsstrmi = {}; lsstimz = {}; lsstzmy = {}; lsstrMag = {}
        numFiles = len(files)
        numOn = 0

        for fileName in files:
            if numOn % 100 == 0:
                print 'Loading %i of %i: Kurucz SEDs' % (numOn, numFiles)
            fileSED = Sed()

            try:
                fileSED.readSED_flambda(str(self.kuruczDir + fileName))
                
                #Read in Kurucz SED properties from filename
                logZTimesTen, temp, gravity, fineTemp = [x.split(".")[0] for x in fileName.split("_")]

                if logZTimesTen[1] == 'm':
                    fileSED.logZ = -1.0 * float(logZTimesTen[2:]) * 0.1
                else:
                    fileSED.logZ = float(logZTimesTen[2:]) * 0.1

                fileSED.logg = float(gravity[1:]) * 0.1
                fileSED.temp = float(fineTemp)
                
                sEDPhotometry = phot()
                sEDMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, self.phiArray, 
                                                            self.wavelenstep, self.bandpassDict, 
                                                            self.filterList)
                lsstSedMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, self.lsstPhiArray,
                                                                self.lsstWavelenstep, 
                                                                self.lsstBandpassDict,
                                                                self.lsstFilterList)
                
                kFiles.append(fileName)
                kumg.append(sEDMagDict['u']-sEDMagDict['g']) #Calculate colors for SED matching
                kgmr.append(sEDMagDict['g']-sEDMagDict['r'])
                krmi.append(sEDMagDict['r']-sEDMagDict['i'])
                kimz.append(sEDMagDict['i']-sEDMagDict['z'])
                lsstumg[fileName] = (lsstSedMagDict['u']-lsstSedMagDict['g']) #lsst for catalog entry
                lsstgmr[fileName] = (lsstSedMagDict['g']-lsstSedMagDict['r'])
                lsstrmi[fileName] = (lsstSedMagDict['r']-lsstSedMagDict['i'])
                lsstimz[fileName] = (lsstSedMagDict['i']-lsstSedMagDict['z'])
                lsstzmy[fileName] = (lsstSedMagDict['z']-lsstSedMagDict['y'])
                krMag[fileName] = sEDMagDict['r']
                lsstrMag[fileName] = lsstSedMagDict['r']
                klogZ.append(fileSED.logZ)
                klogg.append(fileSED.logg)
                kTemp.append(fileSED.temp)
                
                numOn += 1
                
            except:
                continue

        kDict['sEDName'] = kFiles
        kDict['umg'] = kumg
        kDict['gmr'] = kgmr
        kDict['rmi'] = krmi
        kDict['imz'] = kimz
        kDict['lsstumg'] = lsstumg
        kDict['lsstgmr'] = lsstgmr
        kDict['lsstrmi'] = lsstrmi
        kDict['lsstimz'] = lsstimz
        kDict['lsstzmy'] = lsstzmy
        kDict['rMags'] = krMag
        kDict['lsstrMags'] = lsstrMag
        kDict['logZ'] = klogZ
        kDict['logg'] = klogg
        kDict['temp'] = kTemp

        return kDict

    def loadmltSEDs(self, subset = None):

        """
        By default will load all seds in mlt directory, otherwise a subset can be loaded if a 
        list of filenames of only seds wanted in the folder is provided. Will skip over extraneous
        files in sed folder.

        @param [in] subset is the list of the subset of files wanted if one doesn't want all files
        in the mlt directory.

        @param [out] mltDict is a dictionary with all the necessary mlt information to match SEDs

        """

        files = []
        mltDict = {}
        
        if subset is None:
            for fileName in os.listdir(self.mltDir):
                files.append(fileName)
        else:
            for fileName in subset:
                files.append(fileName)
        
        mltFiles = [] 
        mltumg = []; mltgmr = []; mltrmi = []; mltimz = []; mrMag = {}
        lsstumg = {}; lsstgmr = {}; lsstrmi = {}; lsstimz = {}; lsstzmy = {}; lsstrMag = {}
        numFiles = len(files)
        numOn = 0

        for fileName in files:
            if numOn % 100 == 0:
                print 'Loading %i of %i: MLT SEDs' % (numOn, numFiles)
            
            fileSED = Sed()
            try:
                fileSED.readSED_flambda(str(self.mltDir + fileName))
                
                sEDPhotometry = phot()
                sEDMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, self.phiArray, self.wavelenstep,
                                                            self.bandpassDict, self.filterList)
                lsstSedMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, self.lsstPhiArray,
                                                                self.lsstWavelenstep,
                                                                self.lsstBandpassDict, self.lsstFilterList)

                mltFiles.append(fileName)
                mltumg.append(sEDMagDict['u']-sEDMagDict['g'])
                mltgmr.append(sEDMagDict['g']-sEDMagDict['r'])
                mltrmi.append(sEDMagDict['r']-sEDMagDict['i'])
                mltimz.append(sEDMagDict['i']-sEDMagDict['z'])
                lsstumg[fileName] = (lsstSedMagDict['u']-lsstSedMagDict['g']) #lsst for catalog entry
                lsstgmr[fileName] = (lsstSedMagDict['g']-lsstSedMagDict['r'])
                lsstrmi[fileName] = (lsstSedMagDict['r']-lsstSedMagDict['i'])
                lsstimz[fileName] = (lsstSedMagDict['i']-lsstSedMagDict['z'])
                lsstzmy[fileName] = (lsstSedMagDict['z']-lsstSedMagDict['y'])
                mrMag[fileName] = sEDMagDict['r']
                lsstrMag[fileName] = lsstSedMagDict['r']
                
                numOn += 1
            
            except:
                continue

        mltDict['sEDName'] = mltFiles
        mltDict['umg'] = mltumg
        mltDict['gmr'] = mltgmr
        mltDict['rmi'] = mltrmi
        mltDict['imz'] = mltimz
        mltDict['lsstumg'] = lsstumg
        mltDict['lsstgmr'] = lsstgmr
        mltDict['lsstrmi'] = lsstrmi
        mltDict['lsstimz'] = lsstimz
        mltDict['lsstzmy'] = lsstzmy
        mltDict['rMags'] = mrMag
        mltDict['lsstrMags'] = lsstrMag
        
        return mltDict


    def loadwdSEDs(self, subset = None):

        """
        By default will load all seds in kurucz directory, otherwise a subset can be loaded if a 
        list of filenames of only seds wanted in the folder is provided. Will skip over extraneous
        files in sed folder.

        @param [in] subset is the list of the subset of files wanted if one doesn't want all files
        in the kurucz directory.

        @param [out] wdAllDict is a dictionary with all the necessary wd information to match SEDs.
        Notice that there is an extra layer of the dictionary due to the separation of H and HE WD SEDs.

        """
        files = []
        wdAllDict = {}; wdDict = {}; wdHEDict = {}
        
        if subset is None:
            for fileName in os.listdir(self.wdDir):
                files.append(fileName)
        else:
            for fileName in subset:
                files.append(fileName)
        
        wdFiles = []; wdHEFiles = [] 
        wdumg = []; wdHEumg = []; lsstumg = {}; lsstHEumg = {}
        wdgmr = []; wdHEgmr = []; lsstgmr = {}; lsstHEgmr = {}
        wdrmi = []; wdHErmi = []; lsstrmi = {}; lsstHErmi = {}
        wdimz = []; wdHEimz = []; lsstimz = {}; lsstHEimz = {}
        lsstzmy = {}; lsstHEzmy = {}
        wdrMag = {}; wdHErMag = {}; wdlsstrMag = {}; wdHElsstrMag = {}
        numFiles = len(files)
        numOn = 0

        for fileName in files:
            if numOn % 100 == 0:
                print 'Loading %i of %i: WD SEDs' % (numOn, numFiles)
            
            fileSED = Sed()
            try:
                fileSED.readSED_flambda(str(self.wdDir + fileName))
                
                sEDPhotometry = phot()
                sEDMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, self.phiArray, self.wavelenstep,
                                                        self.bandpassDict, self.filterList)
                lsstSedMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, self.lsstPhiArray,
                                                                self.lsstWavelenstep,
                                                                self.lsstBandpassDict, self.lsstFilterList)

                if fileName.split("_")[1] == 'He':    
                    wdHEFiles.append(fileName)
                    wdHEumg.append(sEDMagDict['u']-sEDMagDict['g'])
                    wdHEgmr.append(sEDMagDict['g']-sEDMagDict['r'])
                    wdHErmi.append(sEDMagDict['r']-sEDMagDict['i'])
                    wdHEimz.append(sEDMagDict['i']-sEDMagDict['z'])
                    lsstHEumg[fileName] = (lsstSedMagDict['u']-lsstSedMagDict['g']) 
                    lsstHEgmr[fileName] = (lsstSedMagDict['g']-lsstSedMagDict['r'])
                    lsstHErmi[fileName] = (lsstSedMagDict['r']-lsstSedMagDict['i'])
                    lsstHEimz[fileName] = (lsstSedMagDict['i']-lsstSedMagDict['z'])
                    lsstHEzmy[fileName] = (lsstSedMagDict['z']-lsstSedMagDict['y'])
                    wdHErMag[fileName] = sEDMagDict['r']
                    wdHElsstrMag[fileName] = lsstSedMagDict['r']
                else:
                    wdFiles.append(fileName)
                    wdumg.append(sEDMagDict['u']-sEDMagDict['g'])
                    wdgmr.append(sEDMagDict['g']-sEDMagDict['r'])
                    wdrmi.append(sEDMagDict['r']-sEDMagDict['i'])
                    wdimz.append(sEDMagDict['i']-sEDMagDict['z'])
                    lsstumg[fileName] = (lsstSedMagDict['u']-lsstSedMagDict['g'])
                    lsstgmr[fileName] = (lsstSedMagDict['g']-lsstSedMagDict['r'])
                    lsstrmi[fileName] = (lsstSedMagDict['r']-lsstSedMagDict['i'])
                    lsstimz[fileName] = (lsstSedMagDict['i']-lsstSedMagDict['z'])
                    lsstzmy[fileName] = (lsstSedMagDict['z']-lsstSedMagDict['y'])
                    wdrMag[fileName] = sEDMagDict['r']
                    wdlsstrMag[fileName] = lsstSedMagDict['r']
                    
                numOn += 1

            except:
                continue

        wdDict['sEDName'] = wdFiles; wdHEDict['sEDName'] = wdHEFiles
        wdDict['umg'] = wdumg; wdHEDict['umg'] = wdHEumg
        wdDict['gmr'] = wdgmr; wdHEDict['gmr'] = wdHEgmr
        wdDict['rmi'] = wdrmi; wdHEDict['rmi'] = wdHErmi
        wdDict['imz'] = wdimz; wdHEDict['imz'] = wdHEimz
        wdDict['lsstumg'] = lsstumg; wdHEDict['lsstumg'] = lsstHEumg
        wdDict['lsstgmr'] = lsstgmr; wdHEDict['lsstgmr'] = lsstHEgmr
        wdDict['lsstrmi'] = lsstrmi; wdHEDict['lsstrmi'] = lsstHErmi
        wdDict['lsstimz'] = lsstimz; wdHEDict['lsstimz'] = lsstHEimz
        wdDict['lsstzmy'] = lsstzmy; wdHEDict['lsstzmy'] = lsstHEzmy
        wdDict['rMags'] = wdrMag; wdHEDict['rMags'] = wdHErMag
        wdDict['lsstrMags'] = wdlsstrMag; wdHEDict['lsstrMags'] = wdHElsstrMag

        wdAllDict['H'] = wdDict
        wdAllDict['HE'] = wdHEDict

        return wdAllDict

    def deReddenGalfast(self, am, magU, magG, magR, magI, magZ, 
                        coeffs=np.array([1.8551, 1.4455, 1.0, 0.7431, 0.5527])):
        """
        
        This will take Galfast magnitudes and remove the reddening effects before matching the SED.

        NOTE: Make sure coeffs match those used in galfast photometry.conf file to add reddening.

        @param [in] am is the extinction parameter for the object from galfast output

        @param [in] magU is the U-band magnitude from galfast output

        @param [in] magG is the G-band magnitude from galfast output

        @param [in] magR is the R-band magnitude from galfast output

        @param [in] magI is the I-band magnitude from galfast output

        @param [in] magZ is the Z-band magnitude form galfast output

        @param [in] coeffs is the set of coefficients from galfast's photometry.conf file that scale
        the extinction in each band, usually calibrated to 1.0 in R-band

        @param [out] galumg is the U-G color after dereddening

        @param [out] galgmr is the G-R color after dereddening

        @param [out] galrmi is the R-I color after dereddening

        @param [out] galimz is the I-Z color after dereddening

        """

        uExt = am * coeffs[0]
        gExt = am * coeffs[1]
        rExt = am * coeffs[2]
        iExt = am * coeffs[3]
        zExt = am * coeffs[4]
        magUCorr = magU - uExt
        magGCorr = magG - gExt
        magRCorr = magR - rExt
        magICorr = magI - iExt
        magZCorr = magZ - zExt
        galumg = magUCorr - magGCorr
        galgmr = magGCorr - magRCorr
        galrmi = magRCorr - magICorr
        galimz = magICorr - magZCorr
        
        return galumg, galgmr, galrmi, galimz

    def findSED(self, sEDDict, magU, magG, magR, magI, magZ, am, comp, reddening = True, 
                coeffs=np.array([1.8551, 1.4455, 1.0, 0.7431, 0.5527])):

        """
        This will find the closest match to an input SED. sEDDict must have 'kurucz', 'mlt', 'wdH', or 'wdHE'
        depending on which types you are inputting and will determine which it is based upon its 
        galfast component value

        @param [in] sEDDict is a dictionary from the load(type) routines that are a part of this class

        @param [in] magU is the U-band magnitude from galfast output

        @param [in] magG is the G-band magnitude from galfast output

        @param [in] magR is the R-band magnitude from galfast output

        @param [in] magI is the I-band magnitude from galfast output

        @param [in] magZ is the Z-band magnitude form galfast output

        @param [in] am is the extinction parameter for the object from galfast output

        @param [in] comp is the component number from the galfast output. Used to determine the type
        of SED to use. See galfast documentation for more info.

        @param [in] reddening indicates whether there is extinction and reddening included
        in the galfast output

        @param [in] coeffs is the set of coefficients from galfast's photometry.conf file that scale
        the extinction in each band, usually calibrated to 1.0 in R-band        

        @param [out] sEDName is the name of the SED file that most closely matches the input mags
        accounting for type of star
        
        """
        
        if reddening == True:
            umg, gmr, rmi, imz = self.deReddenGalfast(am, magU, magG, magR, magI, magZ, coeffs)
        else:
            umg = magU - magG
            gmr = magG - magR
            rmi = magR - magI
            imz = magI - magZ

        if 10 <= comp < 20:
            #This is a Galfast outputted WD
            if 10 <= comp < 15:
                wdDict = sEDDict['wdH']
            else:
                wdDict = sEDDict['wdHE']

            wdumg = np.array(wdDict['umg'])
            wdgmr = np.array(wdDict['gmr'])
            wdrmi = np.array(wdDict['rmi'])
            wdimz = np.array(wdDict['imz'])
            
            sEDName = wdDict['sEDName']
            
            distance = np.power((wdumg - umg),2) + np.power((wdgmr - gmr),2) +\
                np.power((wdrmi - rmi),2) + np.power((wdimz - imz),2)
            
        else:
            #For stars that are not WDs

            if rmi > 0.59:
            #Use mlt SEDs

                mltDict = sEDDict['mlt']
                
                mltumg = np.array(mltDict['umg'])
                mltgmr = np.array(mltDict['gmr'])
                mltrmi = np.array(mltDict['rmi'])
                mltimz = np.array(mltDict['imz'])
                
                sEDName = mltDict['sEDName']

                #u,g mags unreliable for cool stars
                distance = np.power((mltrmi - rmi), 2) + np.power((mltimz - imz), 2)

            else:
            #Use Kurucz otherwise

                kDict = sEDDict['kurucz']
                
                kumg = np.array(kDict['umg'])
                kgmr = np.array(kDict['gmr'])
                krmi = np.array(kDict['rmi'])
                kimz = np.array(kDict['imz'])
                
                sEDName = kDict['sEDName']
                
                distance = np.power((kumg - umg),2) + np.power((kgmr - gmr),2) +\
                    np.power((krmi - rmi),2) + np.power((kimz - imz),2)
            
        return sEDName[np.argmin(distance)]
