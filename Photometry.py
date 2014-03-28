"""
photUtils - 


ljones@astro.washington.edu  (and ajc@astro.washington.edu)

Collection of utilities to aid usage of Sed and Bandpass with dictionaries.

"""

import os
import numpy
import palpy as pal
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
from lsst.sims.catalogs.measures.instance import compound

class PhotometryBase(object):
    
    bandPasses = {}
    bandPassKey = []   
    phiArray = None
    waveLenStep = None
    subDir = None
    
    def setupPhiArray_dict(self,bandpassDict, bandpassKeys):
        """ Generate 2-dimensional numpy array for Phi values in the bandpassDict.

        You must pass in bandpassKeys so that the ORDER of the phiArray and the order of the magnitudes returned by
        manyMagCalc can be preserved. You only have to do this ONCE and can then reuse phiArray many times for many
        manyMagCalculations."""
        # Make a list of the bandpassDict for phiArray - in the ORDER of the bandpassKeys
        bplist = []
        for f in bandpassKeys:
            bplist.append(bandpassDict[f])
        sedobj = Sed()
        self.phiArray, self.waveLenStep = sedobj.setupPhiArray(bplist)

    def loadBandPasses(self,bandPassList,bandPassRoot="total_"):
        """
        This will take the list of band passes in bandPassList and use them to set up
        self.phiArray and self.waveLenStep (which is being cached so that it does not have
        to be loaded again unless we change which bandpasses we want)
        """
        if self.bandPassKey != bandPassList:
            self.bandPassKey=[]
            self.bandPasses={}
            path = os.getenv('LSST_THROUGHPUTS_DEFAULT')
            for i in range(len(bandPassList)):
                self.bandPassKey.append(bandPassList[i])
            
            for w in self.bandPassKey:    
                self.bandPasses[w] = Bandpass()
                self.bandPasses[w].readThroughput(os.path.join(path,"%s%s.dat"%bandPassRoot%w))
        
            self.setupPhiArray_dict(self.bandPasses,self.bandPassKey)
            
    # Handy routines for handling Sed/Bandpass routines with sets of dictionaries.
    def loadSeds(self,sedList, magNorm=15.0, resample_same=False):
        """
        Generate list of SEDs required for generating magnitudes
        """    
        
        dataDir=os.getenv('SED_DATA')+self.subDir
        
        imsimband = Bandpass()
        imsimband.imsimBandpass()
        sedOut=[]
        firstsed = True
        for i in range(len(sedList)):
            sedName = sedList[i]
            if sedName == None:
                continue
            elif sedName in sedDict:
                continue
            else:            
                sed = Sed()
                sed.readSED_flambda(os.path.join(dataDir, sedName))
                if resample_same:
                    if firstsed:
                        wavelen_same = sed.wavelen
                        firstsed = False
                    else:
                        if sed.needResample(wavelen_same):
                            sed.resampleSED(wavelen_same)
                
                fNorm = sed.calcFluxNorm(magNorm[i], imsimband)
                sed.multiplyFluxNorm(fNorm)

                sedOut.append(sed)
                
        return sedOut
    
    def applyAvAndRedshift(self,sedList,internalAv=None,redshift=None):
        
        for i in range(len(sedList)):
            if internalAv:
                a_int, b_int = sedList[i].setupCCMab()
                sedList[i].addCCMDust(a_int, b_int, A_v=internalAv[i])
            if redshift:
                sedList[i].redshiftSED(redshift[i], dimming=True)
                sedList[i].resampleSED(wavelen_match=self.bandPasses[self.bandPassKey[0]].waveln)

    def manyMagCalc_dict(self,sedobj):
        """Return a dictionary of magnitudes for a single Sed object.

        You must pass the sed itself, phiArray and wavelenstep for the manyMagCalc itself, but
        you must also pass the bandpassDictionary and keys so that the sedobj can be resampled onto
        the correct wavelength range for the bandpass (i.e. maybe you redshifted sedobj??) and so that the
        order of the magnitude calculation / dictionary assignment is preserved from the phiArray setup previously.
        Note that THIS WILL change sedobj by resampling it onto the required wavelength range. """
        # Set up the SED for using manyMagCalc - note that this CHANGES sedobj
        # Have to check that the wavelength range for sedobj matches bandpass - this is why the dictionary is passed in.
        if sedobj.needResample(wavelen_match=self.bandPasses[self.bandPassKey[0]].wavelen):
            sedobj.resampleSED(wavelen_match=self.bandPasses[self.bandPassKey[0]].wavelen)
        sedobj.flambdaTofnu()
        magArray = sedobj.manyMagCalc(self.phiArray, self.waveLenStep)
        magDict = {}
        i = 0
        for f in bandpassKeys:
            magDict[f] = magArray[i]
            i = i + 1
        return magDict

######################################
class PhotometryGalaxies(PhotometryBase):
    
    def calculate_magnitudes(self,bandPassList,idNames):
        """
        will return a dict of magntiudes which is indexed by
        galid
        (total, bulge, disk, agn)
        bandPass
        """
        self.loadBandPasses(bandPassList)
  
        #idNames=self.column_by_name('galid')
        
        diskNames=self.column_by_name('sedFilenameDisk')
        bulgeNames=self.column_by_name('sedFilenameBulge')
        agnNames=self.column_by_name('sedFilenameAgn')

        diskmn = self.column_by_name('magNormDisk')
        bulgemn = self.column_by_name('magNormBulge')
        agnmn = self.column_by_name('magNormAgn')
        
        bulgeAv = self.column_by_name('internalAvBulge')
        diskAv = self.column_by_name('internalAvDisk')

        redshift = self.column_by_name('redshift')
        
        diskMags=[]
        bulgeMags=[]
        agnMags=[]
            
        if diskNames:
            self.subDir="/galaxySED/"
            diskSed = self.loadSeds(diskNames,magNorm = diskmn)
            self.applyAvAndRedshift(diskSed,internalAv = bulgeAv, redshift = redshift)
        
            for dd in diskSed:
                subDict=self.manyMagCalc_dict(dd)
                diskMags.append(subDict)
        
        else:
            subDict={}
            for i in range(len(idNames)):
                diskMags.append(subDict)
         
        if bulgeNames:
            self.subDir="/galaxySED/"
            bulgeSed = self.loadSeds(bulgeNames,magNorm = bulgemn)
            self.applyAvAndRedshift(bulgeSed,internalAv = diskAv, redshift = redshift)
            
            for bb in bulgeSed:
                subDict=self.manyMagCalc_dict(bb)
                bulgeMags.append(subDict)
        else:
            subDict={}
            for i in range(len(idNames)):
                bulgeMags.append(subDict)
        
        if agnNames: 
            self.subDir="/agnSED/"    
            agnSed = self.loadSeds(agnNames,magNorm = agnmn)
            self.applyAvAndRedshift(agnSed,redshift = redshift)
            
            for aa in agnSed:
                subDict=self.manyMagCalc_dict(aa)
                agnMags.append(subDict)
        
        else:
            subDict={}
            for i in range(len(idNames)):
                agnMags.append(subDict)
        
        total_mags = {}
        masterList = {}

        m_o = 22.
        
        for i in len(idNames):
            total_mags={}
            for ff in bandPassList:
                nn=0.0
                if diskMags[i]:
                    nn+=numpy.power(10, (diskMags[i][ff] - m_o)/-2.5)
                
                if bulgeMags[i]:
                    nn+=numpy.power(10, (bulgeMags[i][ff] - m_o)/-2.5)
            
                if agnMags[i]:
                    nn+=numpy.power(10, (agnMags[i][ff] - m_o)/-2.5)
                
                total_mags[ff] = -2.5*log10(nn) + m_o
            
            subDict={}
            subDict["total"] = total_mags
            subDict["bulge"] = bulgeMags[i]
            subDict["disk"] = diskMags[i]
            subDict["agn"] = agnMags[i]
            
            masterDict[idNames[i]] = subDict


        return masterDict
     
    @compound('uRecalc', 'gRecalc', 'rRecalc', 'iRecalc', 'zRecalc', 'yRecalc',
              'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
              'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
              'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn')
    def get_allMags(self):
        bandPassList=['u','g','r','i','z','y']
        idNames=self.column_by_name('galid')
        magDict=self.calculate_magnitudes(bandPassList,idNames)
        
        utotal=numpy.zeros(len(idNames),dtype=float)
        gtotal=numpy.zeros(len(idNames),dtype=float)
        rtotal=numpy.zeros(len(idNames),dtype=float)
        itotal=numpy.zeros(len(idNames),dtype=float)
        ztotal=numpy.zeros(len(idNames),dtype=float)
        ytotal=numpy.zeros(len(idNames),dtype=float)
        
        ubulge=numpy.zeros(len(idNames),dtype=float)
        gbulge=numpy.zeros(len(idNames),dtype=float)
        rbulge=numpy.zeros(len(idNames),dtype=float)
        ibulge=numpy.zeros(len(idNames),dtype=float)
        zbulge=numpy.zeros(len(idNames),dtype=float)
        ybulge=numpy.zeros(len(idNames),dtype=float)
        
        udisk=numpy.zeros(len(idNames),dtype=float)
        gdisk=numpy.zeros(len(idNames),dtype=float)
        rdisk=numpy.zeros(len(idNames),dtype=float)
        idisk=numpy.zeros(len(idNames),dtype=float)
        zdisk=numpy.zeros(len(idNames),dtype=float)
        ydisk=numpy.zeros(len(idNames),dtype=float)
        
        uagn=numpy.zeros(len(idNames),dtype=float)
        gagn=numpy.zeros(len(idNames),dtype=float)
        ragn=numpy.zeros(len(idNames),dtype=float)
        iagn=numpy.zeros(len(idNames),dtype=float)
        zagn=numpy.zeros(len(idNames),dtype=float)
        yagn=numpy.zeros(len(idNames),dtype=float)
        
        i=0
        failure=-999.0
        for name in idNames:
            utotal[i]=magDict[name]["total"]["u"]
            gtotal[i]=magDict[name]["total"]["g"]
            rtotal[i]=magDict[name]["total"]["r"]
            itotal[i]=magDict[name]["total"]["i"]
            ztotal[i]=magDict[name]["total"]["z"]
            ytotal[i]=magDict[name]["total"]["y"]
           
            if magDict[name]["bulge"]:
                ubulge[i]=magDict[name]["bulge"]["u"]
                gbulge[i]=magDict[name]["bulge"]["g"]
                rbulge[i]=magDict[name]["bulge"]["r"]
                ibulge[i]=magDict[name]["bulge"]["i"]
                zbulge[i]=magDict[name]["bulge"]["z"]
                ybulge[i]=magDict[name]["bulge"]["y"]
            else:
                ubulge[i]=failure
                gbulge[i]=failure
                rbulge[i]=failure
                ibulge[i]=failure
                zbulge[i]=failure
                ybulge[i]=failure
           
            if magDict[name]["disk"]:
                udisk[i]=magDict[name]["disk"]["u"]
                gdisk[i]=magDict[name]["disk"]["g"]
                rdisk[i]=magDict[name]["disk"]["r"]
                idisk[i]=magDict[name]["disk"]["i"]
                zdisk[i]=magDict[name]["disk"]["z"]
                ydisk[i]=magDict[name]["disk"]["y"]
            else:
                udisk[i]=failure
                gdisk[i]=failure
                rdisk[i]=failure
                idisk[i]=failure
                zdisk[i]=failure
                ydisk[i]=failure
           
            if magDict[name]["agn"]:
                uagn[i]=magDict[name]["agn"]["u"]
                gagn[i]=magDict[name]["agn"]["g"]
                ragn[i]=magDict[name]["agn"]["r"]
                iagn[i]=magDict[name]["agn"]["i"]
                zagn[i]=magDict[name]["agn"]["z"]
                yagn[i]=magDict[name]["agn"]["y"]
            else:
                uagn[i]=failure
                gagn[i]=failure
                ragn[i]=failure
                iagn[i]=failure
                zagn[i]=faiure
                yagn[i]=failure
        
        return numpy.array([utotal,gtotal,rtotal,itotal,ztotal,ytotal,\
            ubulge,gbulge,rbulge,ibulge,zbulge,ybulge,\
            udisk,gdisk,rdisk,idisk,zdisk,ydisk,\
            uagn,gagn,ragn,iagn,zagn,yagn])



        
    def calculate_star_mags(self, bandPassList):
        """
        return a dictionary of magnitudes indexed first by sedName and then bandpass
        """
        sedNames = self.column_by_name('sedFilename')
        magDict=self.calculate_magnitudes(bandPassList,sedNames)

        return magDict
