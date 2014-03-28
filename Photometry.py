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
                self.bandPasses[w].readThroughput(os.path.join(path,"%s.dat" % (bandPassRoot + w)))
        
            self.setupPhiArray_dict(self.bandPasses,self.bandPassKey)
            
    # Handy routines for handling Sed/Bandpass routines with sets of dictionaries.
    def loadSeds(self,sedList, magNorm=15.0, resample_same=False):
        """
        Generate list of SEDs required for generating magnitudes
        """    
        
        if self.subDir:
            dataDir=os.getenv('SED_DATA')+self.subDir
        else:
            dataDir=os.getenv('SED_DATA')
        
        imsimband = Bandpass()
        imsimband.imsimBandpass()
        sedOut=[]
        firstsed = True
        for i in range(len(sedList)):
            sedName = sedList[i]
            if sedName == "None":
                #assign an empty Sed (one with wavelen==None)
                sed = Sed()
                #continue
            #removed the code below because we now want to load all SEDs
            #since objects can have identical SEDs, but different magNorms
            #elif sedName in sedDict:
                #continue
            else:          
                #print "opening SED ",sedName  
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
        
        if len(sedOut) != len(sedList):
            print "WARNING ",len(sedOut)," ",len(sedList)        
        return sedOut
    
    def applyAvAndRedshift(self,sedList,internalAv=None,redshift=None):
        
        for i in range(len(sedList)):
            if sedList[i].wavelen != None:
                if internalAv != None:
                    a_int, b_int = sedList[i].setupCCMab()
                    sedList[i].addCCMDust(a_int, b_int, A_v=internalAv[i])
                if redshift != None:
                    sedList[i].redshiftSED(redshift[i], dimming=True)
                    sedList[i].resampleSED(wavelen_match=self.bandPasses[self.bandPassKey[0]].wavelen)

    def manyMagCalc_dict(self,sedobj):
        """Return a dictionary of magnitudes for a single Sed object.

        You must pass the sed itself, phiArray and wavelenstep for the manyMagCalc itself, but
        you must also pass the bandpassDictionary and keys so that the sedobj can be resampled onto
        the correct wavelength range for the bandpass (i.e. maybe you redshifted sedobj??) and so that the
        order of the magnitude calculation / dictionary assignment is preserved from the phiArray setup previously.
        Note that THIS WILL change sedobj by resampling it onto the required wavelength range. """
        # Set up the SED for using manyMagCalc - note that this CHANGES sedobj
        # Have to check that the wavelength range for sedobj matches bandpass - this is why the dictionary is passed in.
        
        magDict={}
        if sedobj.wavelen != None:
            if sedobj.needResample(wavelen_match=self.bandPasses[self.bandPassKey[0]].wavelen):
                sedobj.resampleSED(wavelen_match=self.bandPasses[self.bandPassKey[0]].wavelen)
            sedobj.flambdaTofnu()
            magArray = sedobj.manyMagCalc(self.phiArray, self.waveLenStep)
            i = 0
            for f in self.bandPassKey:
                magDict[f] = magArray[i]
                i = i + 1
        else:
            for f in self.bandPassKey:
                magDict[f] = None
                  
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
        
        diskMags={}
        bulgeMags={}
        agnMags={}
            
        if diskNames != []:
            self.subDir="/galaxySED/"
            diskSed = self.loadSeds(diskNames,magNorm = diskmn)
            self.applyAvAndRedshift(diskSed,internalAv = bulgeAv, redshift = redshift)
            
            print "lenidNames ",len(idNames)," lenDiskSed ",len(diskSed)," lenDisknames ",len(diskNames)
            
            for i in range(len(idNames)):
                subDict=self.manyMagCalc_dict(diskSed[i])
                diskMags[idNames[i]]=subDict
        
        else:
            subDict={}
            for b in bandPassList:
                subDict[b]=None
            for i in range(len(idNames)):
                distMags[idNames[i]]=subDict
         
        if bulgeNames != []:
            self.subDir="/galaxySED/"
            bulgeSed = self.loadSeds(bulgeNames,magNorm = bulgemn)
            self.applyAvAndRedshift(bulgeSed,internalAv = diskAv, redshift = redshift)
            
            for i in range(len(idNames)):
                subDict=self.manyMagCalc_dict(bulgeSed[i])
                bulgeMags[idNames[i]]=subDict
        else:
            subDict={}
            for b in bandPassList:
                subDict[b]=None
            for i in range(len(idNames)):
                bulgeMags[idNames[i]]=subDict
        
        if agnNames != []: 
            self.subDir="/agnSED/"    
            agnSed = self.loadSeds(agnNames,magNorm = agnmn)
            self.applyAvAndRedshift(agnSed,redshift = redshift)
            
            for i in range(len(idNames)):
                subDict=self.manyMagCalc_dict(agnSed[i])
                agnMags[idNames[i]]=subDict
        
        else:
            subDict={}
            for b in bandPassList:
                subDict[b]=None
            for i in range(len(idNames)):
                agnMags[idNames[i]]=subDict
        
        total_mags = {}
        masterDict = {}

        m_o = 22.
        
       # print len(diskMags),len(bulgeMags),len(agnMags),len(idNames)
       # print diskMags
       # print bulgeMags
       # print agnMags
        for i in range(len(idNames)):
            total_mags={}
            for ff in bandPassList:
                nn=0.0
                if diskMags[idNames[i]][ff]!=None:
                    nn+=numpy.power(10, (diskMags[idNames[i]][ff] - m_o)/-2.5)
                
                if bulgeMags[idNames[i]][ff]!=None:
                    nn+=numpy.power(10, (bulgeMags[idNames[i]][ff] - m_o)/-2.5)
            
                if agnMags[idNames[i]][ff]!=None:
                    nn+=numpy.power(10, (agnMags[idNames[i]][ff] - m_o)/-2.5)
                
                if nn>0.0:
                    total_mags[ff] = -2.5*numpy.log10(nn) + m_o
                else:
                    total_mags[ff] = None
                
            subDict={}
            subDict["total"] = total_mags
            subDict["bulge"] = bulgeMags[idNames[i]]
            subDict["disk"] = diskMags[idNames[i]]
            subDict["agn"] = agnMags[idNames[i]]
            
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
        for i in range(len(idNames)):
            name=idNames[i]
            
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


class PhotometryStars(PhotometryBase):

    def calculate_magnitudes(self,bandPassList,idNames):
        
        self.subDir="/starSED/kurucz/"
        
        self.loadBandPasses(bandPassList)
        sedNames = self.column_by_name('sedFilename')
        magNorm = self.column_by_name('magNorm')
        sedList = self.loadSeds(sedNames,magNorm = magNorm)
        
        magDict = {}
        for i in range(len(idNames)):
            name = idNames[i]
            subDict = self.manyMagCalc_dict(sedList[i])
            magDict[name] = subDict
        
        return magDict

    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y')
    def get_magnitudes(self):
        idNames = self.column_by_name('id')
        bandPassList = ['u','g','r','i','z','y']
        
        magDict = self.calculate_magnitudes(bandPassList,idNames)
        
        uu = numpy.zeros(len(idNames),dtype=float)
        gg = numpy.zeros(len(idNames),dtype=float)
        rr = numpy.zeros(len(idNames),dtype=float)
        ii = numpy.zeros(len(idNames),dtype=float)
        zz = numpy.zeros(len(idNames),dtype=float)
        yy = numpy.zeros(len(idNames),dtype=float)
        
        for i in range(len(idNames)):
            uu[i] = magDict[idNames[i]]["u"]
            gg[i] = magDict[idNames[i]]["g"]
            rr[i] = magDict[idNames[i]]["r"]
            ii[i] = magDict[idNames[i]]["i"]
            zz[i] = magDict[idNames[i]]["z"]
            yy[i] = magDict[idNames[i]]["y"]
        
        return numpy.array([uu,gg,rr,ii,zz,yy])
      
