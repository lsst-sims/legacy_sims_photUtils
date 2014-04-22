"""
photUtils - 


ljones@astro.washington.edu  (and ajc@astro.washington.edu)

and now (2014 March 28): scott.f.daniel@gmail.com

Collection of utilities to aid usage of Sed and Bandpass with dictionaries.

"""

import os
import numpy
import palpy as pal
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.sims.catalogs.measures.instance import compound

class PhotometryBase(object):
    """
    This mixin provides the basic infrastructure for photometry.
    It can read in SEDs and bandpasses, apply extinction and redshift, and, given
    an SED object it can calculate magnitudes.
    
    In order to avoid duplication of work, the bandPasses, wavelength array, and phi array
    are stored as instance variables once they are read in by self.loadBandPasses()
    
    To initiailize a different set of bandPasses, call self.loadBandPasses() with a different
    set of arguments.
    
    Once self.loadBandPasses() as been called, self.loadSeds() can be used to return an array
    of SED objects.  These objects can be passed to self.manyMagCalc_dict() which will calculate
    the magnitudes of the the SEDs, integrated over the loaded bandPasses, and return them as a 
    dict keeyed to the array of bandpass keys stored in self.bandPassKey
    """
    
    bandPasses = {}
    bandPassKey = []   
    phiArray = None
    waveLenStep = None
        
    def setupPhiArray_dict(self):
        """ 
        Generate 2-dimensional numpy array for Phi values associated with the bandpasses in
        self.bandPasses

        self.bandpassKey is used so that the ORDER of the phiArray and the order of the magnitudes returned by
        manyMagCalc can be preserved. 
        
        The results from this calculation will be stored in the instance variables
        self.phiArray and self.waveLenStep for future use by self.manyMagCalc_dict()
        """
        # Make a list of the bandpassDict for phiArray - in the ORDER of the bandpassKeys
        bplist = []
        for f in self.bandPassKey:
            bplist.append(self.bandPasses[f])
        sedobj = Sed()
        self.phiArray, self.waveLenStep = sedobj.setupPhiArray(bplist)

    def loadBandPasses(self,bandPassList,bandPassRoot="total_"):
        """
        This will take the list of band passes in bandPassList and use them to set up
        self.bandPasses, self.phiArray and self.waveLenStep (which are being cached so that 
        they do not have to be loaded again unless we change which bandpasses we want)
        
        bandPassRoot contains the first part of the bandpass file name, i.e., it is assumed
        that the bandPasses are stored in files of the type
        
        $LSST_THROUGHPUTS_DEFAULT/bandPassRoot_bandPassKey.dat
        
        if we want to load bandpasses for a telescope other than LSST, we would do so
        by altering bandPassRoot (currently no infrastructure exists for altering the directory
        in which bandpass files are stored)
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
        
            self.setupPhiArray_dict()
            
    # Handy routines for handling Sed/Bandpass routines with sets of dictionaries.
    def loadSeds(self,sedList, magNorm=15.0, resample_same=False):
        """
        Takes the list of filename sedList and returns an array of SED objects.
        
        This code will load identical SEDs twice because it is possible for
        (astronomical) objects to have the same SEDs but different magNorms
        """    
        
        dataDir=os.getenv('SED_DATA')
        
        #initialize a delta function bandpass for use in applying magNorm
        imsimband = Bandpass()
        imsimband.imsimBandpass()
        
        
        sedOut=[]
        firstsed = True
        for i in range(len(sedList)):
            sedName = sedList[i]
            if sedName == "None":
                #assign an empty Sed (one with wavelen==None)
                sed = Sed()
            else:          
                sed = Sed()
                sed.readSED_flambda(os.path.join(dataDir, self.specFileMap[sedName]))
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
        """
        Take the array of SED objects sedList and apply the arrays of extinction and redshift
        (internalAV and redshift)
        """
        for i in range(len(sedList)):
            if sedList[i].wavelen != None:
                if internalAv != None:
                    a_int, b_int = sedList[i].setupCCMab()
                    sedList[i].addCCMDust(a_int, b_int, A_v=internalAv[i])
                if redshift != None:
                    sedList[i].redshiftSED(redshift[i], dimming=True)
                    sedList[i].resampleSED(wavelen_match=self.bandPasses[self.bandPassKey[0]].wavelen)

    def manyMagCalc_dict(self,sedobj):
        """
        Return a dictionary of magnitudes for a single Sed object.
        
        Bandpass information is taken from the instance variables self.bandPasses, self.bandPassKey,
        self.phiArray, and self.waveLenStep
        
        Returns a dictionary of magnitudes keyed on self.bandPassKey
        """
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

    def calculatePhotometricUncertaintyFromColumn(self,nameTag,columnNames):
        """
        This method reads in a dict of column names and passes out
        the associated photometric uncertainties.  The output will be
        a dict of lists.
        
        The dict is such that  
        
        columnNames[filterName] gives the name of the column corresponding the 
        filter denoted by filterName
   
        nameTag indicates what column is used for object names
        """
        
        inputDict={}
        
        idNames = self.column_by_name(nameTag)
        
        magnitudes = {}
        
        for filterName in columnNames:
            magnitudes[filterName] = self.column_by_name(columnNames[filterName])
        
        i=0
        for name in idNames:
            subDict={}
            
            for filterName in columnNames:
                subDict[filterName] = magnitudes[filterName][i]
            
            inputDict[name] = subDict
            i = i+1
        
        outputDict=self.calculatePhotometricUncertainty(inputDict)
        
        finalDict = {}
        for filterName in columnNames:
            subList = []
            for name in idNames:
                subList.append(outputDict[name][filterName])
            
            finalDict[filterName] = subList
    
        return finalDict
        
    def calculatePhotometricUncertainty(self,magDict):
        """
        This method is based on equations 3.1, 3.2 and Table 3.2
        of the LSST Science Book (version 2.0)
        
        magDict will be two-level dict of magnitudes, i.e.
        
        magDict['name']['filter'] will be the magnitude of object 'name'
        in the appropriate filter
        
        This method will return a similar dict of photometric uncertainties
        """
        sigma2Sys = 0.003*0.003 #also taken from the Science Book
                         #see the paragraph between equations 3.1 and 3.2
        
        gamma = {}
        m5 = {}
        
        gamma['u'] = 0.037
        gamma['g'] = 0.038
        gamma['r'] = 0.039
        gamma['i'] = 0.039
        gamma['z'] = 0.040
        gamma['y'] = 0.040
        
        m5['u'] = 23.9
        m5['g'] = 25.0
        m5['r'] = 24.7
        m5['i'] = 24.0
        m5['z'] = 23.3
        m5['y'] = 22.1
        
        sigOut={}
        
        for name in magDict:
            filterDict = magDict[name]
            
            subDict = {}
            
            for filterName in filterDict:
                mm = filterDict[filterName]

                xx=10**(0.4*(mm - m5[filterName]))
                ss = (0.04 - gamma[filterName])*xx + \
                     gamma[filterName]*xx*xx
                
                sigmaSquared = ss + sigma2Sys
                
                subDict[filterName] = numpy.sqrt(sigmaSquared)
            
            sigOut[name] = subDict
        
        return sigOut

class PhotometryGalaxies(PhotometryBase):
    """
    This mixin provides the code necessary for calculating the component magnitudes associated with
    galaxies.  It assumes that we want LSST filters.
    """
    
    def calculate_component_magnitudes(self,objectNames, componentNames, bandPassList, \
                                       magNorm = 15.0, internalAv = None, redshift = None):
        
        """
        Calculate the magnitudes for different components (disk, bulge, agn, etc) of galaxies.
        
        @param [in] objectNames is the name of the galaxies (the whole galaxies)
        
        @param [in] componentNames gives the name of the SED files for the component in question
        
        @param [in] bandPassList lists the bandpasses for which we want magnitudes (this will come
        from calculate_magnitudes()
        
        @param [in] magNorm is the normalizing magnitude
        
        @param [in] internalAv is the internal Av extinction
        
        @param [in] redshift is pretty self-explanatory
        
        This will return a dict of dicts such that
        
        magnitude["objectname"]["filter label"] will return the magnitude in that filter
        for the component being calculated
        
        """
        
        componentMags = {}
        
        if componentNames != []:
            componentSed = self.loadSeds(componentNames,magNorm = magNorm)
            self.applyAvAndRedshift(componentSed,internalAv = internalAv, redshift = redshift)
            
            for i in range(len(objectNames)):
                subDict=self.manyMagCalc_dict(componentSed[i])
                componentMags[objectNames[i]]=subDict
        
        else:
            subDict={}
            for b in bandPassList:
                subDict[b]=None
            for i in range(len(objectNames)):
                componentMags[objectNames[i]]=subDict
    
        return componentMags
    
    def sum_magnitudes(self,disk = None, bulge = None, agn = None):
        mm_o = 22.
        
        nn=0.0
        if disk is not None:
            nn+=numpy.power(10, (disk - mm_o)/-2.5)
                
        if bulge is not None:
            nn+=numpy.power(10, (bulge - mm_o)/-2.5)
            
        if agn is not None:
            nn+=numpy.power(10, (agn - mm_o)/-2.5)
                
        if nn>0.0:
            outMag = -2.5*numpy.log10(nn) + mm_o
        else:
            outMag = None
        
        return outMag
    
    def calculate_magnitudes(self,bandPassList,idNames):
        """
        Take the array of bandpass keys bandPassList and the array of galaxy
        names idNames ane return a dict of dicts of dicts of magnitudes
        
        the first level key is galid (the name of the galaxy)
        
        the second level key is "total", "bulge", "disk", or "agn"
        
        the third level key is bandPassList
        
        We need to index the galaxies by some unique identifier, such as galid
        because it is possible for galaxies to have the same sed filenames but 
        different normalizations
        
        """
        self.loadBandPasses(bandPassList)
        
        diskNames=self.column_by_name('sedFilenameDisk')
        bulgeNames=self.column_by_name('sedFilenameBulge')
        agnNames=self.column_by_name('sedFilenameAgn')

        diskmn = self.column_by_name('magNormDisk')
        bulgemn = self.column_by_name('magNormBulge')
        agnmn = self.column_by_name('magNormAgn')
        
        bulgeAv = self.column_by_name('internalAvBulge')
        diskAv = self.column_by_name('internalAvDisk')

        redshift = self.column_by_name('redshift')
         
        diskMags = self.calculate_component_magnitudes(idNames,diskNames,bandPassList,magNorm = diskmn, \
                        internalAv = diskAv, redshift = redshift)
                        
        bulgeMags = self.calculate_component_magnitudes(idNames,bulgeNames,bandPassList,magNorm = bulgemn, \
                        internalAv = bulgeAv, redshift = redshift)
                        
        agnMags = self.calculate_component_magnitudes(idNames,agnNames,bandPassList,magNorm = agnmn, \
                        redshift = redshift)
        
        total_mags = {}
        masterDict = {}

        for i in range(len(idNames)):
            total_mags={}
            for ff in bandPassList:
                total_mags[ff]=self.sum_magnitudes(disk = diskMags[idNames[i]][ff],
                                bulge = bulgeMags[idNames[i]][ff], agn = agnMags[idNames[i]][ff])
                
                
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
    
    @compound('sigma_uRecalc','sigma_gRecalc','sigma_rRecalc',
              'sigma_iRecalc','sigma_zRecalc','sigma_yRecalc')
    def get_photometric_uncertainties(self):
        
        columnNames = {}
        columnNames['u'] = 'uRecalc'
        columnNames['g'] = 'gRecalc'
        columnNames['r'] = 'rRecalc'
        columnNames['i'] = 'iRecalc'
        columnNames['z'] = 'zRecalc'
        columnNames['y'] = 'yRecalc'
        
        outputDict = self.calculatePhotometricUncertaintyFromColumn('galid',columnNames)
        
        return numpy.array([outputDict['u'],outputDict['g'],outputDict['r'],
                            outputDict['i'],outputDict['z'],outputDict['y']])
        
        
        

class PhotometryStars(PhotometryBase):
    """
    This mixin provides the infrastructure for doing photometry on stars
    
    It assumes that we want LSST filters.
    """
                         
    def calculate_magnitudes(self,bandPassList,idNames):
        """
        Take the array of bandpass keys bandPassList and the array of
        star names idNames and return a dict of dicts of magnitudes
        
        The first level key will be the name of the star (idName)
        
        The second level key will be the name of the filter (bandPassList)
        
        As with galaxies, it is important that we identify stars by a unique
        identifier, rather than their sedFilename, because different stars
        can have identical SEDs but different magnitudes.
        """

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
      
    @compound('sigma_lsst_u','sigma_lsst_g','sigma_lsst_r','sigma_lsst_i',
              'sigma_lsst_z','sigma_lsst_y')
    def get_photometric_uncertainties(self):
        idNames = self.column_by_name('id')
        
        uu = self.column_by_name('lsst_u')
        gg = self.column_by_name('lsst_g')
        rr = self.column_by_name('lsst_r')
        ii = self.column_by_name('lsst_i')
        zz = self.column_by_name('lsst_z')
        yy = self.column_by_name('lsst_y')
        
        inputDict={}
        i = 0
        for name in idNames:
            subDict={}
            subDict['u'] = uu[i]
            subDict['g'] = gg[i]
            subDict['r'] = rr[i]
            subDict['i'] = ii[i]
            subDict['z'] = zz[i]
            subDict['y'] = yy[i]
          
            inputDict[name] = subDict
          
            i += 1
        
        outputDict = self.calculatePhotometricUncertainty(inputDict)
        
        uuOut = []
        ggOut = []
        rrOut = []
        iiOut = []
        zzOut = []
        yyOut = []
        
        for name in idNames:
            uuOut.append(outputDict[name]['u'])
            ggOut.append(outputDict[name]['g'])
            rrOut.append(outputDict[name]['r'])
            iiOut.append(outputDict[name]['i'])
            zzOut.append(outputDict[name]['z'])
            yyOut.append(outputDict[name]['y'])    

        return numpy.array([uuOut,ggOut,rrOut,iiOut,zzOut,yyOut])
