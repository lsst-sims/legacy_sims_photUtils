import os
import numpy as np

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.photUtils import Photometry as phot

class selectStarSED():
    def __init__(self):
        self.sEDDir = os.environ['SIMS_SED_LIBRARY_DIR']

    def setupSDSS(self):
        #Load Bandpasses for SDSS colors to match to galfast output
        sdssPhot = phot()
        filterlist = ('u', 'g', 'r', 'i', 'z')
        sdssBandpassDict = sdssPhot.loadBandpasses(filterlist=filterlist, dataDir = os.getenv("SDSS_THROUGHPUTS"), filterroot='sdss_')
        sdssPhiArray, wavelenstep = sdssPhot.setupPhiArray_dict(sdssBandpassDict, filterlist)
        return filterlist, sdssBandpassDict, sdssPhiArray, wavelenstep

    def setupLSST(self):
        #Load Bandpasses for LSST colors to get colors from matched SEDs
        lsstPhot = phot()
        lsstFilterList = ('u', 'g', 'r', 'i', 'z', 'y')
        lsstBandpassDict = lsstPhot.loadBandpasses()
        lsstPhiArray, lsstWavelenstep = lsstPhot.setupPhiArray_dict(lsstBandpassDict, bandpassKeys = lsstFilterList)
        return lsstFilterList, lsstBandpassDict, lsstPhiArray, lsstWavelenstep

    def loadKuruczSEDs(self):
        files = []
        #Load Kurucz
        kDict = {}
        kuruczDir = str(self.sEDDir + '/starSED/kurucz/')
        
        for fileName in os.listdir(kuruczDir):
            files.append(fileName)
        
        filterlist, sdssBandpassDict, sdssPhiArray, wavelenstep = self.setupSDSS()

        kFiles = []; klogZ = []; kTemp = []; klogg = []
        kumg = []; kgmr = []; krmi = []; kimz = []
        numFiles = len(files)
        numOn = 0

        for file in files[0:201]:
            if numOn % 100 == 0:
                print 'Loading %i of %i: Kurucz SEDs' % (numOn, numFiles)
            fileSED = Sed()
            fileSED.readSED_flambda(str(kuruczDir + file))

            #Read in Kurucz SED properties from filename
            logZTimesTen, temp, gravity, fineTemp = [x.split(".")[0] for x in file.split("_")]

            if logZTimesTen[1] == 'm':
                fileSED.logZ = -1.0 * float(logZTimesTen[2:]) * 0.1
            else:
                fileSED.logZ = float(logZTimesTen[2]) * 0.1

            fileSED.logg = float(gravity[1:]) * 0.1
            fileSED.temp = fineTemp

            sEDPhotometry = phot()
            sEDMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, sdssPhiArray, wavelenstep, sdssBandpassDict, filterlist)

            kFiles.append(file)
            kumg.append(sEDMagDict['u']-sEDMagDict['g'])
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

    def loadmltSEDs(self):

        #Load WD SEDs
        files = []
        mltDict = {}
        mltDir = str(self.sEDDir + '/starSED/mlt/')
        
        for fileName in os.listdir(mltDir):
            files.append(fileName)
        
        filterlist, sdssBandpassDict, sdssPhiArray, wavelenstep = self.setupSDSS()

        mltFiles = [] 
        mltumg = []; mltgmr = []; mltrmi = []; mltimz = []
        numFiles = len(files)
        numOn = 0

        for file in files:
            if numOn % 100 == 0:
                print 'Loading %i of %i: MLT SEDs' % (numOn, numFiles)
            
            fileSED = Sed()
            fileSED.readSED_flambda(str(mltDir + file))

            sEDPhotometry = phot()
            sEDMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, sdssPhiArray, wavelenstep, sdssBandpassDict, filterlist)

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


    def loadwdSEDs(self):

        #Load WD SEDs
        files = []
        wdAllDict = {}; wdDict = {}; wdHEDict = {}
        wdDir = str(self.sEDDir + '/starSED/wDs/')
        
        for fileName in os.listdir(wdDir):
            files.append(fileName)
        
        filterlist, sdssBandpassDict, sdssPhiArray, wavelenstep = self.setupSDSS()

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
            fileSED.readSED_flambda(str(wdDir + file))

            sEDPhotometry = phot()
            sEDMagDict = sEDPhotometry.manyMagCalc_dict(fileSED, sdssPhiArray, wavelenstep, sdssBandpassDict, filterlist)

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

    def parseGalfast(self, headerLine):
        galfastDict = {}
        header = headerLine.split(' ')
        colNo = 0
        for title in header:
            if title == 'lb[2]':
                galfastDict['l'] = colNo
                colNo += 1
                galfastDict['b'] = colNo
                colNo += 1
            elif title == 'radec[2]':
                galfastDict['ra'] = colNo
                colNo += 1
                galfastDict['dec'] = colNo
                colNo += 1
            elif title == 'XYZ[3]':
                galfastDict['X'] = colNo
                colNo += 1
                galfastDict['Y'] = colNo
                colNo += 1
                galfastDict['Z'] = colNo
                colNo += 1
            elif title == 'DM':
                galfastDict['DM'] = colNo
                colNo += 1
            elif title == 'absSDSSr{alias=M1;alias=absmag;band=SDSSr;}':
                galfastDict['absSDSSr'] = colNo
                colNo += 1
            elif title == 'comp':
                galfastDict['comp'] = colNo
                colNo += 1
            elif title == 'FeH':
                galfastDict['FeH'] = colNo
                colNo += 1
            elif title == 'vcyl[3]':
                galfastDict['Vr'] = colNo
                colNo += 1
                galfastDict['Vphi'] = colNo
                colNo += 1
                galfastDict['Vz'] = colNo
                colNo += 1
            elif title == 'pmlb[3]':
                galfastDict['pml'] = colNo
                colNo += 1
                galfastDict['pmb'] = colNo
                colNo += 1
                galfastDict['vRadlb'] = colNo
                colNo += 1
            elif title == 'pmradec[3]':
                galfastDict['pmra'] = colNo
                colNo += 1
                galfastDict['pmdec'] = colNo
                colNo += 1
                galfastDict['vRad'] = colNo
                colNo += 1
            elif title == 'Am':
                galfastDict['Am'] = colNo
                colNo += 1
            elif title == 'AmInf':
                galfastDict['AmInf'] = colNo
                colNo += 1
            elif title.startswith('SDSSugriz['):
                bandString = title.split('=')[2]
                bandString1 = bandString.split(',')
                for band in bandString1:
                    band = band.rstrip(';}')
                    bandName = band.split(':')[1]
                    galfastDict[bandName] = colNo
                    colNo += 1
            elif title == 'SDSSugrizPhotoFlags{class=flags;}':
                galfastDict['SDSSPhotoFlags'] = colNo
                colNo += 1
            elif title == '#': pass
            elif title.isspace(): pass
            elif len(title) < 1: pass
            else:
                raise RuntimeError, '*** Unknown field: %s' % (title)
        return galfastDict

    def deReddenGalfast(self, am, magU, magG, magR, magI, magZ, coeffs=np.array([1.8551, 1.4455, 1.0, 0.7431, 0.5527])):
        #Make sure coeffs match those used in galfast photometry.conf file
        uExt = am * coeffs[0]
        gExt = am * coeffs[1]
        rExt = am
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

    def findSED(self, sEDDict, magU, magG, magR, magI, magZ, absMagR, comp):
        
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
            
            distance = np.power((mltumg - (magU - magG)),2) + np.power((mltgmr - (magG - magR)),2) +\
                np.power((mltrmi - (magR - magI)),2) + np.power((mltimz - (magI - magZ)),2)
            
        else:
            #For stars that are not WDs
            if (magR - magI) > 0.6:
            #Use mlt SEDs

                mltDict = sEDDict['mlt']
                
                mltumg = np.array(mltDict['umg'])
                mltgmr = np.array(mltDict['gmr'])
                mltrmi = np.array(mltDict['rmi'])
                mltimz = np.array(mltDict['imz'])
                
                sEDName = mltDict['sEDName']

                distance = np.power((mltumg - (magU - magG)),2) + np.power((mltgmr - (magG - magR)),2) +\
                    np.power((mltrmi - (magR - magI)),2) + np.power((mltimz - (magI - magZ)),2)

            else:
            #Use Kurucz otherwise

                kDict = sEDDict['kurucz']
                
                kumg = np.array(kDict['umg'])
                kgmr = np.array(kDict['gmr'])
                krmi = np.array(kDict['rmi'])
                kimz = np.array(kDict['imz'])
                
                sEDName = kDict['sEDName']
                
                distance = np.power((kumg - (magU - magG)),2) + np.power((kgmr - (magG - magR)),2) +\
                    np.power((krmi - (magR - magI)),2) + np.power((kimz - (magI - magZ)),2)
            
            #Equation B2 from Ivezic et al. 2008, only good for T>4000K, if we want this info need\
                #to find another equation for mlt
        #logT = 3.882 - (magG - magR)*(0.316 - (magG - magR)*(0.0488 + (magG - magR)*0.0283))
        #tEff = np.power(10, logT)

        return sEDName[np.argmin(distance)]

    def findLSSTMags(self, sEDName, absMagR):
        sEDObj = Sed()
        if sEDName.startswith('k'):
            #Kurucz SED
            sEDObj.readSED_flambda(self.sEDDir + '/starSED/kurucz/' + sEDName)
        else:
            #mlt SED
            sEDObj.readSED_flambda(self.sEDDir + '/starSED/mlt/' + sEDName)
        filterlist, sdssBandpassDict, sdssPhiArray, wavelenstep = self.setupSDSS()
        sED_fluxnorm = sEDObj.calcFluxNorm(absMagR, sdssBandpassDict['r'])
        sEDObj.multiplyFluxNorm(sED_fluxnorm)
        
        lsstFilterList, lsstBandpassDict, lsstPhiArray, lsstWavelenstep = self.setupLSST()
        
        sEDPhot = phot()
        lsstMagDict = sEDPhot.manyMagCalc_dict(sEDObj, lsstPhiArray, lsstWavelenstep, lsstBandpassDict, bandpassKeys = ('u', 'g', 'r', 'i', 'z', 'y'))
        
            #For Testing only
        sdssMagDict = sEDPhot.manyMagCalc_dict(sEDObj, sdssPhiArray, wavelenstep, sdssBandpassDict, bandpassKeys = ('u', 'g', 'r', 'i', 'z'))
        
        return lsstMagDict, sdssMagDict #Test with SDSSr which should be same as galfast output

    def loadGalfast(self, filename):
        #Only adding support for .txt file at the moment, will add .fits support later

        sEDDict = {}

        kDict = self.loadKuruczSEDs()
        mltDict = self.loadmltSEDs()
        wdDict = self.loadwdSEDs()

        sEDDict['kurucz'] = kDict
        sEDDict['mlt'] = mltDict
        sEDDict['wdH'] = wdDict['H']
        sEDDict['wdHE'] = wdDict['HE']

        if filename.endswith('.txt'):
            galfastIn = open(filename, 'r')
            lineNum = 0
            for line in galfastIn:
                if lineNum == 0:
                    galfastDict = self.parseGalfast(line)
                    lineNum += 1
                if line[0] == '#': continue
                lineNum += 1
                oID = float(lineNum)
                lineData = line.split()
                DM = float(lineData[galfastDict['DM']])
                absSDSSr = float(lineData[galfastDict['absSDSSr']])
                comp = float(lineData[galfastDict['comp']])
                feH = float(lineData[galfastDict['FeH']])
                am = float(lineData[galfastDict['Am']])
                sDSSu = float(lineData[galfastDict['SDSSu']])
                sDSSg = float(lineData[galfastDict['SDSSg']])
                sDSSr = float(lineData[galfastDict['SDSSr']])
                sDSSi = float(lineData[galfastDict['SDSSi']])
                sDSSz = float(lineData[galfastDict['SDSSz']])
                print sDSSu, sDSSg, sDSSr, sDSSi, sDSSz
                sEDName = self.findSED(sEDDict, sDSSu, sDSSg, sDSSr, sDSSi, sDSSz, absSDSSr, comp)
                lsstMagDict, sdssMagDict = self.findLSSTMags(sEDName, absSDSSr)
                print sdssMagDict
                print float(sdssMagDict['r'] + DM + am), float(sdssMagDict['i'] + DM + (am*0.7431))
        else:
            raise RuntimeError, '*** Can only import .txt Galfast Catalogs at the moment'
