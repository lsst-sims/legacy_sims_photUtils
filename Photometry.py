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

    def loadBandPasses(self,bandPassList):
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
                self.bandPasses[w].readThroughput(os.path.join(path,"total_%s.dat"%w))
        
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
    
        def calculate_magnitudes(self,bandPassList):
            """
            will return a dict of magntiudes which is indexed by
            galid
            (total, bulge, disk, agn)
            bandPass
            """
            self.loadBandPasses(bandPassList)
  
            idNames=self.column_by_name('galid')
        
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
                
        
            if bulgeNames:
                self.subDir="/galaxySED/"
                bulgeSed = self.loadSeds(bulgeNames,magNorm = bulgemn)
                self.applyAvAndRedshift(bulgeSed,internalAv = diskAv, redshift = redshift)
            
                for bb in bulgeSed:
                    subDict=self.manyMagCalc_dict(bb)
                    bulgeMags.append(subDict)
        
            if agnNames: 
                self.subDir="/agnSED/"    
                agnSed = self.loadSeds(agnNames,magNorm = agnmn)
                self.applyAvAndRedshift(agnSed,redshift = redshift)
            
                for aa in agnSed:
                    subDict=self.manyMagCalc_dict(aa)
                    agnMags.append(subDict)
       
            total_mags = {}
            masterList={}

            m_o = 22.
        
            for i in len(idNames):
                total_mags={}
                for ff in bandPassList:
                    nn=0.0
                    if diskMags:
                        nn+=numpy.power(10, (diskMags[i][ff] - m_o)/-2.5)
                
                    if bulgeMags:
                        nn+=numpy.power(10, (bulgeMags[i][ff] - m_o)/-2.5)
                
                    if agnMags:
                        nn+=numpy.power(10, (agnMags[i][ff] - m_o)/-2.5)
                
                    total_mags[ff] = -2.5*log10(nn) + m_o
              
                masterDict[idNames[i]]["total"] = total_mags
                masterDict[idNames[i]]["bulge"] = bulgeMags[bb]
                masterDict[idNames[i]]["disk"] = diskMags[dd]
                masterDict[idNames[i]]["agn"] = agnMags[aa]

            return masterDict
     
     #spock put the galaxy getter here
        
    def calculate_star_mags(self, bandPassList):
        """
        return a dictionary of magnitudes indexed first by sedName and then bandpass
        """
        sedNames = self.column_by_name('sedFilename')
        magDict=self.calculate_magnitudes(bandPassList,sedNames)

        return magDict
    

    
    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y',
    'lsst_bulge_u','lsst_bulge_g','lsst_bulge_r','lsst_bulge_i',
    'lsst_bulge_z','lsst_bulge_y','lsst_disk_u','lsst_disk_g','lsst_disk_r',
    'lsst_disk_i','lsst_disk_z','lsst_disk_y','lsst_agn_u','lsst_agn_g','lsst_agn_r',
    'lsst_agn_i','lsst_agn_z','lsst_agn_y')
    def get_magnitudes(self):
        bandPassList=['u','g','r','i','z','y']
        sedNames=self.column_by_name('sedFilename')
        magDict=self.calculate_magnitudes(bandPassList,sedNames)
        
        uu=numpy.zeros(len(sedNames),dtype=float)
        gg=numpy.zeros(len(sedNames),dtype=float)
        rr=numpy.zeros(len(sedNames),dtype=float)
        ii=numpy.zeros(len(sedNames),dtype=float)
        zz=numpy.zeros(len(sedNames),dtype=float)
        yy=numpy.zeros(len(sedNames),dtype=float)
        buu=numpy.zeros(len(sedNames),dtype=float)
        bgg=numpy.zeros(len(sedNames),dtype=float)
        brr=numpy.zeros(len(sedNames),dtype=float)
        bii=numpy.zeros(len(sedNames),dtype=float)
        bzz=numpy.zeros(len(sedNames),dtype=float)
        byy=numpy.zeros(len(sedNames),dtype=float)
        duu=numpy.zeros(len(sedNames),dtype=float)
        dgg=numpy.zeros(len(sedNames),dtype=float)
        drr=numpy.zeros(len(sedNames),dtype=float)
        dii=numpy.zeros(len(sedNames),dtype=float)
        dzz=numpy.zeros(len(sedNames),dtype=float)
        dyy=numpy.zeros(len(sedNames),dtype=float)
        auu=numpy.zeros(len(sedNames),dtype=float)
        agg=numpy.zeros(len(sedNames),dtype=float)
        arr=numpy.zeros(len(sedNames),dtype=float)
        aii=numpy.zeros(len(sedNames),dtype=float)
        azz=numpy.zeros(len(sedNames),dtype=float)
        ayy=numpy.zeros(len(sedNames),dtype=float)
        
        print "sedNames ",sedNames
        for i in range(len(sedNames)):
           name=sedNames[i]
           uu[i]=magDict[name]['u']
           gg[i]=magDict[name]['g']
           rr[i]=magDict[name]['r']
           ii[i]=magDict[name]['i']
           zz[i]=magDict[name]['z']
           yy[i]=magDict[name]['y']
           
           buu[i]=0.0
           bgg[i]=0.0
           brr[i]=0.0
           bii[i]=0.0
           bzz[i]=0.0
           byy[i]=0.0
           
           duu[i]=0.0
           dgg[i]=0.0
           drr[i]=0.0
           dii[i]=0.0
           dzz[i]=0.0
           dyy[i]=0.0
           
           auu[i]=0.0
           agg[i]=0.0
           arr[i]=0.0
           aii[i]=0.0
           azz[i]=0.0
           ayy[i]=0.0
       
       
        return numpy.array([uu,gg,rr,ii,zz,yy,buu,bgg,brr,bii,bzz,byy,duu,dgg,drr,dii,dzz,dyy,auu,agg,arr,aii,azz,ayy])
