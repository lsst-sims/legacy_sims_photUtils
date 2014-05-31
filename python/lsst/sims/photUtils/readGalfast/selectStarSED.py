import os
import gzip
import pyfits
import numpy as np

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.photUtils import Photometry as phot

class selectStarSED():
    def __init__(self):
        self.sEDDir = os.environ['SIMS_SED_LIBRARY_DIR']
        self.kuruczDir = str(self.sEDDir + '/starSED/kurucz/')
        self.mltDir = str(self.sEDDir + '/starSED/mlt/')
        self.wdDir = str(self.sEDDir + '/starSED/wDs/')


        #Load Bandpasses for SDSS colors to match to galfast output
        sdssPhot = phot()
        self.sdssFilterList = ('u', 'g', 'r', 'i', 'z')
        self.sdssBandpassDict = sdssPhot.loadBandpasses(filterlist=self.sdssFilterList, dataDir = os.getenv("SDSS_THROUGHPUTS"), filterroot='sdss_')
        self.sdssPhiArray, self.sdssWavelenstep = sdssPhot.setupPhiArray_dict(self.sdssBandpassDict, self.sdssFilterList)

        #Load Bandpasses for LSST colors to get colors from matched SEDs
        lsstPhot = phot()
        self.lsstFilterList = ('u', 'g', 'r', 'i', 'z', 'y')
        self.lsstBandpassDict = lsstPhot.loadBandpasses()
        self.lsstPhiArray, self.lsstWavelenstep = lsstPhot.setupPhiArray_dict(self.lsstBandpassDict, bandpassKeys = self.lsstFilterList)

    def loadKuruczSEDs(self, subset = None):
        #By default will load all seds in kurucz directory, otherwise a subset can be loaded if a list of filenames of only seds wanted in the folder is provided.

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
        kmagu = []; kumg = []; kgmr = []; krmi = []; kimz = []
        numFiles = len(files)
        numOn = 0

        for file in files:
            if numOn % 100 == 0:
                print 'Loading %i of %i: Kurucz SEDs' % (numOn, numFiles)
            fileSED = Sed()
            fileSED.readSED_flambda(str(self.kuruczDir + file))

            #Read in Kurucz SED properties from filename
            logZTimesTen, temp, gravity, fineTemp = [x.split(".")[0] for x in file.split("_")]

            if logZTimesTen[1] == 'm':
                fileSED.logZ = -1.0 * float(logZTimesTen[2:]) * 0.1
            else:
                fileSED.logZ = float(logZTimesTen[2:]) * 0.1

            fileSED.logg = float(gravity[1:]) * 0.1
            fileSED.temp = float(fineTemp)

            sEDPhotometry = phot()
            sEDMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, self.sdssPhiArray, self.sdssWavelenstep, self.sdssBandpassDict, self.sdssFilterList)

            kFiles.append(file)
            kumg.append(sEDMagDict['u']-sEDMagDict['g']) #Calculate colors for SED matching
            kgmr.append(sEDMagDict['g']-sEDMagDict['r'])
            krmi.append(sEDMagDict['r']-sEDMagDict['i'])
            kimz.append(sEDMagDict['i']-sEDMagDict['z'])
            klogZ.append(fileSED.logZ)
            klogg.append(fileSED.logg)
            kTemp.append(fileSED.temp)
            
            numOn += 1

        kDict['sEDName'] = kFiles
        kDict['umg'] = kumg
        kDict['gmr'] = kgmr
        kDict['rmi'] = krmi
        kDict['imz'] = kimz
        kDict['logZ'] = klogZ
        kDict['logg'] = klogg
        kDict['temp'] = kTemp

        return kDict

    def loadmltSEDs(self, subset = None):

        #Load mlt SEDs. Follow same method as above.
        files = []
        mltDict = {}
        
        if subset is None:
            for fileName in os.listdir(self.mltDir):
                files.append(fileName)
        else:
            for fileName in subset:
                files.append(fileName)
        
        mltFiles = [] 
        mltumg = []; mltgmr = []; mltrmi = []; mltimz = []
        numFiles = len(files)
        numOn = 0

        for file in files:
            if numOn % 100 == 0:
                print 'Loading %i of %i: MLT SEDs' % (numOn, numFiles)
            
            fileSED = Sed()
            fileSED.readSED_flambda(str(self.mltDir + file))

            sEDPhotometry = phot()
            sEDMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, self.sdssPhiArray, self.sdssWavelenstep, self.sdssBandpassDict, self.sdssFilterList)

            mltFiles.append(file)
            mltumg.append(sEDMagDict['u']-sEDMagDict['g'])
            mltgmr.append(sEDMagDict['g']-sEDMagDict['r'])
            mltrmi.append(sEDMagDict['r']-sEDMagDict['i'])
            mltimz.append(sEDMagDict['i']-sEDMagDict['z'])
            
            numOn += 1

        mltDict['sEDName'] = mltFiles
        mltDict['umg'] = mltumg
        mltDict['gmr'] = mltgmr
        mltDict['rmi'] = mltrmi
        mltDict['imz'] = mltimz
        
        return mltDict


    def loadwdSEDs(self, subset = None):

        #Load WD SEDs. Once again following pattern above.
        files = []
        wdAllDict = {}; wdDict = {}; wdHEDict = {}
        
        if subset is None:
            for fileName in os.listdir(self.wdDir):
                files.append(fileName)
        else:
            for fileName in subset:
                files.append(fileName)
        
        wdFiles = []; wdHEFiles = [] 
        wdumg = []; wdHEumg = []
        wdgmr = []; wdHEgmr = []
        wdrmi = []; wdHErmi = [] 
        wdimz = []; wdHEimz = []
        numFiles = len(files)
        numOn = 0

        for file in files:
            if numOn % 100 == 0:
                print 'Loading %i of %i: WD SEDs' % (numOn, numFiles)
            
            fileSED = Sed()
            fileSED.readSED_flambda(str(self.wdDir + file))

            sEDPhotometry = phot()
            sEDMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, self.sdssPhiArray, self.sdssWavelenstep, self.sdssBandpassDict, self.sdssFilterList)

            if file.split("_")[1] == 'He':

                wdHEFiles.append(file)
                wdHEumg.append(sEDMagDict['u']-sEDMagDict['g'])
                wdHEgmr.append(sEDMagDict['g']-sEDMagDict['r'])
                wdHErmi.append(sEDMagDict['r']-sEDMagDict['i'])
                wdHEimz.append(sEDMagDict['i']-sEDMagDict['z'])
            else:
                wdFiles.append(file)
                wdumg.append(sEDMagDict['u']-sEDMagDict['g'])
                wdgmr.append(sEDMagDict['g']-sEDMagDict['r'])
                wdrmi.append(sEDMagDict['r']-sEDMagDict['i'])
                wdimz.append(sEDMagDict['i']-sEDMagDict['z'])

            numOn += 1

        wdDict['sEDName'] = wdFiles; wdHEDict['sEDName'] = wdHEFiles
        wdDict['umg'] = wdumg; wdHEDict['umg'] = wdHEumg
        wdDict['gmr'] = wdgmr; wdHEDict['gmr'] = wdHEgmr
        wdDict['rmi'] = wdrmi; wdHEDict['rmi'] = wdHErmi
        wdDict['imz'] = wdimz; wdHEDict['imz'] = wdHEimz

        wdAllDict['H'] = wdDict
        wdAllDict['HE'] = wdHEDict

        return wdAllDict

    def deReddenGalfast(self, am, magU, magG, magR, magI, magZ, coeffs=np.array([1.8551, 1.4455, 1.0, 0.7431, 0.5527])):
        #Make sure coeffs match those used in galfast photometry.conf file
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

    def findSED(self, sEDDict, magU, magG, magR, magI, magZ, am, comp, reddening = True, coeffs=np.array([1.8551, 1.4455, 1.0, 0.7431, 0.5527])):

        #This will find the closest match to an input SED. sEDDict must have 'kurucz', 'mlt', 'wdH', or 'wdHE' 
        #depending on which types you are inputting and will determine which it is based upon it's galfast component value
        
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

            #Zeljko told me 0.6 and Rob used 0.62 was the cutoff but the lowest r-i in the mlt models is ~0.59 
            #so I went with that otherwise that particular sed (m0.0Full) would never be used and that didn't 
            #seem to make sense. Also I checked and all kurucz currently in sedDir are all r-i<0.59.
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
