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
import lsst.sims.catalogs.measures.photometry.EBV as EBV
from lsst.sims.catalogs.measures.instance import compound

class Photometry(object):
    
    #these variables will tell the mixin where to get the dust maps
    ebvDataDir=os.environ.get("CAT_SHARE_DATA")
    ebvMapNorthName="data/Dust/SFD_dust_4096_ngp.fits"
    ebvMapSouthName="data/Dust/SFD_dust_4096_sgp.fits"
    ebvMapNorth=None
    ebvMapSouth=None
    
    initializedPhotometry=False
    masterBandpassDict=None
    masterFilterList=None
    masterPhiArray=None
    masterWavelenStep=None
    masterSedDirectory=None
    
    def initializePhotometry(self,filterList=None,filterDir=None,filterRoot=None,sedDir=None):
        self.initializedPhotometry=True
        
        if filterList == None:
            filterList=('u','g','r','i','z','y')
        
        self.setFilters(filterList=filterList,filterDir=filterDir,filterRoot=filterRoot)
        
        if sedDir == None:
            sedDir = os.environ.get("CAT_SHARE_DATA")+"data/starSED/kurucz/"
        
        self.setSedDir(sedDir)
        
    
    def setFilters(self,filterList,filterDir,filterRoot):
        self.masterFilterList=filterList
        if filterRoot == None:
            self.masterBandpassDict=self.loadBandpasses(filterlist=self.masterFilterList,dataDir=filterDir)
        else:
            self.masterBandpassDict=self.loadBandpasses(filterlist=self.masterFilterList,dataDir=filterDir,filterroot=filterRoot)
        self.masterPhiArray, self.masterWavelenStep = self.setupPhiArray_dict(self.masterBandpassDict,self.masterFilterList)
    
    def setSedDir(self,sedDir):
        self.masterSedDirectory=sedDir
        
    
    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y')
    def get_LSSTmagnitudes(self):
        if self.initializedPhotometry == False:
            self.initializePhotometry()
        
        """    
        if self.masterBandpassDict == None:
            print "cannot get magnitudes; BandpassDict is None"
        
        if self.masterFilterList == None:
            print "cannot get magnitudes; FilterList is None"
        
        if self.masterPhiArray == None:
            print "cannot get magnitudes; PhiArray is None"
        
        if self.masterWavelenStep == None:
            print "cannot get magnitudes; WavelenStep is None"
        
        if self.masterSedDirectory == None:
            print "cannot get magnitudes; Sed Dir is None"
        """
        
        sedname=self.column_by_name('sedFilename')
        print "sedname ",sedname,len(sedname),self.masterSedDirectory
        sedObj=self.loadSeds(sedname,self.masterSedDirectory)
        print "type sedObj ",type(sedObj)
        
        uu=numpy.zeros(len(sedname),dtype=float)
        gg=numpy.zeros(len(sedname),dtype=float)
        rr=numpy.zeros(len(sedname),dtype=float)
        ii=numpy.zeros(len(sedname),dtype=float)
        zz=numpy.zeros(len(sedname),dtype=float)
        yy=numpy.zeros(len(sedname),dtype=float)
        for i in range(len(sedname)):
            magDict = self.manyMagCalc_dict(sedObj[sedname[i]],self.masterPhiArray,self.masterWavelenStep,self.masterBandpassDict,self.masterFilterList)
            uu[i]=magDict['u']
            gg[i]=magDict['g']
            rr[i]=magDict['r']
            ii[i]=magDict['i']
            zz[i]=magDict['z']
            yy[i]=magDict['y']
            
        #print magDict

        return numpy.array([uu,gg,rr,ii,zz,yy])
     
        
    def get_ug_color(self):
        u = self.column_by_name('lsst_u')
        g = self.column_by_name('lsst_g')
        return u - g
    
    def get_gr_color(self):
        g = self.column_by_name('lsst_g')
        r = self.column_by_name('lsst_r')
        return g - r    
    
    #the set_xxxx routines below will allow the user to point elsewhere for the dust maps
    def set_ebvMapNorth(self,word):
        self.ebvMapNorthName=word
    
    def set_ebvMapSouth(self,word):
        self.ebvMapSouthName=word
    
    #these routines will load the dust maps for the galactic north and south hemispheres
    def load_ebvMapNorth(self):
        self.ebvMapNorth=EBV.EbvMap()
        self.ebvMapNorth.readMapFits(os.path.join(self.ebvDataDir,self.ebvMapNorthName))
    
    def load_ebvMapSouth(self):
        self.ebvMapSouth=EBV.EbvMap()
        self.ebvMapSouth.readMapFits(os.path.join(self.ebvDataDir,self.ebvMapSouthName))
    
    #and finally, here is the getter
    #it relies ont he calculateEbv routine defined in EBV.py
    def get_EBV(self):
        if self.ebvMapNorth==None:
            self.load_ebvMapNorth()
        
        if self.ebvMapSouth==None:
            self.load_ebvMapSouth()
        
        glon=self.column_by_name("glon")
        glat=self.column_by_name("glat")
        
        EBV_out=numpy.array(EBV.calculateEbv(glon,glat,self.ebvMapNorth,self.ebvMapSouth,interp=True))
        #print EBV_out
        return EBV_out
    
    @compound('glon','glat')
    def get_galactic_coords(self):
        ra=self.column_by_name('raJ2000')
        dec=self.column_by_name('decJ2000')
        
        glon=numpy.zeros(len(ra))
        glat=numpy.zeros(len(ra))
        for i in range(len(ra)):
            gg=pal.eqgal(ra[i],dec[i])
            glon[i]=gg[0]
            glat[i]=gg[1]
        
        return numpy.array([glon,glat])
    
    
    def get_galacticRv(self):
        av = self.column_by_name('galacticAv')
        ee = self.column_by_name('EBV')
        return av/ee
    
    # Handy routines for handling Sed/Bandpass routines with sets of dictionaries.
    def loadSeds(self,sedList, dataDir = "./", resample_same=False):
        """Generate dictionary of SEDs required for generating magnitudes

        Given a dataDir and a list of seds return a dictionary with sedName and sed as key, value
        """    
        sedDict={}
        firstsed = True
        for sedName in sedList:
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
                sedDict[sedName] = sed
        return sedDict

    def loadBandpasses(self,filterlist=('u', 'g', 'r', 'i', 'z', 'y'), dataDir=None, filterroot='total_'):
        """ Generate dictionary of bandpasses for the LSST nominal throughputs

        Given a list of filter keys (like u,g,r,i,z,y), return a dictionary of the total bandpasses.
        dataDir is the directory where these bandpasses are stored; leave blank to use environment
         variable 'LSST_THROUGHPUTS_DEFAULT' which is set if throughputs is setup using eups.
        This routine uses the 'total' bandpass values by default, but can be changed (such as to 'filter') using
         the filterroot option (filename = filterroot + filterkey + '.dat'). 
        """
        bandpassDict = {}
        if dataDir == None:
            dataDir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
            if dataDir == None:
                raise Exception("dataDir not given and unable to access environment variable 'LSST_THROUGHPUTS_DEFAULT'")
        for f in filterlist:
            bandpassDict[f] = Bandpass()
            bandpassDict[f].readThroughput(os.path.join(dataDir, filterroot + f + ".dat"))
        
        return bandpassDict

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
        phiArray, wavelenstep = sedobj.setupPhiArray(bplist)
        return phiArray, wavelenstep

    def manyMagCalc_dict(self,sedobj, phiArray, wavelenstep, bandpassDict, bandpassKeys):
        """Return a dictionary of magnitudes for a single Sed object.

        You must pass the sed itself, phiArray and wavelenstep for the manyMagCalc itself, but
        you must also pass the bandpassDictionary and keys so that the sedobj can be resampled onto
        the correct wavelength range for the bandpass (i.e. maybe you redshifted sedobj??) and so that the
        order of the magnitude calculation / dictionary assignment is preserved from the phiArray setup previously.
        Note that THIS WILL change sedobj by resampling it onto the required wavelength range. """
        # Set up the SED for using manyMagCalc - note that this CHANGES sedobj
        # Have to check that the wavelength range for sedobj matches bandpass - this is why the dictionary is passed in.
        if sedobj.needResample(wavelen_match=bandpassDict[bandpassKeys[0]].wavelen):
            sedobj.resampleSED(wavelen_match=bandpassDict[bandpassKeys[0]].wavelen)
        sedobj.flambdaTofnu()
        magArray = sedobj.manyMagCalc(phiArray, wavelenstep)
        magDict = {}
        i = 0
        for f in bandpassKeys:
            magDict[f] = magArray[i]
            i = i + 1
        return magDict


