import os
import gzip
import pyfits
import eups

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.readGalfast.selectStarSED import selectStarSED
from lsst.sims.photUtils.Photometry import PhotometryBase as phot

__all__ = ["readGalfast"]

class readGalfast():

    def parseGalfast(self, headerLine):

        """
        Use galfast header line to organize input

        @param [in] headerLine is the first line from a galfast catalog output file

        @param [out] galfastDict is a dictionary relating parameter name to input column

        """

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
            elif len(title) < 1: pass
            elif title.isspace(): pass
            else:
                raise RuntimeError, '*** Unknown field: %s' % (title)
        return galfastDict

    def convDMtoKpc(self, DM): 
        """
        Change from distance modulus to distance in kiloparsecs

        @param [in] DM is the distance modulus

        @param [out] distanceKpc is the distance in kiloparsecs

        """

        distancePc = 10**((0.2*DM) + 1)
        distanceKpc = distancePc / 1000.
        return distanceKpc

    def loadGalfast(self, filenameList, outFileList, sEDPath = None, kuruczPath = None, 
                    mltPath = None, wdPath = None, kuruczSubset = None, 
                    mltSubset = None, wdSubset = None):
        """
        This is customized for the outputs we currently need for the purposes of consistent output
        It will read in a galfast output file and output desired values for database input into a file

        @param [in] filenameList is a list of the galfast output files that will be loaded and processed

        @param [in] outFileList is a list of the names of the output files that will be created

        @param [in] kuruczPath is a place to specify a different path to kurucz SED files than the
        files in the LSST sims_sed_library. If set to None it will default to the LSST library. 
        Will probably be most useful for those who want to use loadGalfast without downloading the
        entire LSST sims_sed_library which contains much more than just the star SEDs.

        @param [in] mltPath is the same as kuruczPath except that it specifies a directory for the 
        mlt SEDs

        @param [in] wdPath is the same as the previous two except that it specifies a path to an
        alternate white dwarf SED directory.

        @param [in] kuruczSubset is a list which provides a subset of the kurucz files within the
        kurucz folder that one wants to use

        @param [in] mltSubset is a list which provides a subset of the mlt files within the
        mlt folder that one wants to use

        @param [in] wdSubset is a list which provides a subset of the wd files within the
        wd folder that one wants to use

        """

        for filename in filenameList:
            #Make sure input file exists and is readable format before doing anything else
            if os.path.isfile(filename) == False:
                raise RuntimeError, '*** File does not exist'
            
            #Process various possible galfast outputs
            if filename.endswith(('.txt', '.gz', '.fits')):
                continue
            else:
                raise RuntimeError, str('*** Unsupported File Format in file: ' + str(filename))

        #If all files exist and are in proper formats then load seds
        
        selectStarSED0 = selectStarSED(sEDDir = sEDPath, kuruczDir = kuruczPath, 
                                       mltDir = mltPath, wdDir = wdPath)

        if kuruczSubset is None:
            kuruczList = selectStarSED0.loadKuruczSEDs()
        else:
            kuruczList = selectStarSED0.loadKuruczSEDs(subset = kuruczSubset)

        #Only need one dictionary since none of the names overlap
        positionDict = {}
        for kuruczSED, kNum in zip(kuruczList, range(0, len(kuruczList))):
            positionDict[kuruczSED.name] = kNum
            
        if mltSubset is None:
            mltList = selectStarSED0.loadmltSEDs()
        else:
            mltList = selectStarSED0.loadmltSEDs(subset = mltSubset)

        for mltSED, mltNum in zip(mltList, range(0, len(mltList))):
            positionDict[mltSED.name] = mltNum

        if wdSubset is None:
            wdListH, wdListHE = selectStarSED0.loadwdSEDs()
        else:
            wdListH, wdListHE = selectStarSED0.loadwdSEDs(subset = wdSubset)

        for hSED, hNum in zip(wdListH, range(0, len(wdListH))):
            positionDict[hSED.name] = hNum

        for heSED, heNum in zip(wdListHE, range(0, len(wdListHE))):
            positionDict[heSED.name] = heNum

        #For adding/subtracting extinction when calculating colors
        sdssExtCoeffs = [1.8551, 1.4455, 1.0, 0.7431, 0.5527]
        lsstExtCoeffs = [1.8140, 1.4166, 0.9947, 0.7370, 0.5790, 0.4761]

        sdssPhot = phot()
        sdssPhot.loadBandPassesFromFiles(['u','g','r','i','z'], 
                                         bandPassDir = os.path.join(eups.productDir('throughputs'),
                                                                    'sdss'),
                                         bandPassRoot = 'sdss_')
        sdssPhot.setupPhiArray_dict()

        #Calculate colors and add them to the SED objects
        kuruczColors = selectStarSED0.calcBasicColors(kuruczList, sdssPhot)
        mltColors = selectStarSED0.calcBasicColors(mltList, sdssPhot)
        hColors = selectStarSED0.calcBasicColors(wdListH, sdssPhot)
        heColors = selectStarSED0.calcBasicColors(wdListHE, sdssPhot)

        listDict = {'kurucz':kuruczList, 'mlt':mltList, 'H':wdListH, 'HE':wdListHE}
        colorDict = {'kurucz':kuruczColors, 'mlt':mltColors, 'H':hColors, 'HE':heColors}

        for filename, outFile in zip(filenameList, outFileList):

            if filename.endswith('.txt'):
                galfastIn = open(filename, 'r')
                inFits = False
            elif filename.endswith('.gz'):
                galfastIn = gzip.open(filename, 'r')
                inFits = False
            elif filename.endswith('fits'):
                hdulist = pyfits.open(filename)
                galfastIn = hdulist[1].data
                inFits = True

            fOut = open(outFile, 'w')
            fOut.write('#oID, ra, dec, gall, galb, coordX, coordY, coordZ, sEDName, magNorm, ' +\
                       'LSSTugrizy, SDSSugriz, absSDSSr, pmRA, pmDec, vRad, pml, pmb, vRadlb, ' +\
                       'vR, vPhi, vZ, FeH, pop, distKpc, ebv, ebvInf\n')

            lineNum = 0

            for line in galfastIn:
                if inFits:
                    oID = int(lineNum)
                    starData = galfastIn[lineNum]
                    gall, galb = starData.field('lb')
                    ra, dec = starData.field('radec')
                    coordX, coordY, coordZ = starData.field('XYZ')
                    DM = starData.field('DM')
                    absSDSSr = starData.field('absSDSSr')
                    pop = starData.field('comp')
                    FeH = starData.field('FeH')
                    vR, vPhi, vZ = starData.field('vcyl')
                    pml, pmb, vRadlb = starData.field('pmlb')
                    pmRA, pmDec, vRad = starData.field('pmradec')
                    am = starData.field('Am')
                    amInf = starData.field('AmInf')
                    sDSSu, sDSSg, sDSSr, sDSSi, sDSSz = starData.field('SDSSugriz')
                    sdssPhotoFlags = starData.field('SDSSugrizPhotoFlags')
                else:
                    if lineNum == 0:
                        galfastDict = self.parseGalfast(line)
                        lineNum += 1
                    if line[0] == '#': continue
                    oID = int(lineNum)
                    lineData = line.split()
                    gall = float(lineData[galfastDict['l']])
                    galb = float(lineData[galfastDict['b']])
                    ra = float(lineData[galfastDict['ra']])
                    dec = float(lineData[galfastDict['dec']])
                    coordX = float(lineData[galfastDict['X']])
                    coordY = float(lineData[galfastDict['Y']])
                    coordZ = float(lineData[galfastDict['Z']])
                    DM = float(lineData[galfastDict['DM']])
                    absSDSSr = float(lineData[galfastDict['absSDSSr']])
                    pop = float(lineData[galfastDict['comp']])
                    FeH = float(lineData[galfastDict['FeH']])
                    vR = float(lineData[galfastDict['Vr']])
                    vPhi = float(lineData[galfastDict['Vphi']])
                    vZ = float(lineData[galfastDict['Vz']])
                    pml = float(lineData[galfastDict['pml']])
                    pmb = float(lineData[galfastDict['pmb']])
                    vRadlb = float(lineData[galfastDict['vRadlb']])
                    pmRA = float(lineData[galfastDict['pmra']])
                    pmDec = float(lineData[galfastDict['pmdec']])
                    vRad = float(lineData[galfastDict['vRad']])
                    am = float(lineData[galfastDict['Am']])
                    amInf = float(lineData[galfastDict['AmInf']])
                    sDSSu = float(lineData[galfastDict['SDSSu']])
                    sDSSg = float(lineData[galfastDict['SDSSg']])
                    sDSSr = float(lineData[galfastDict['SDSSr']])
                    sDSSi = float(lineData[galfastDict['SDSSi']])
                    sDSSz = float(lineData[galfastDict['SDSSz']])
                    sDSSPhotoFlags = float(lineData[galfastDict['SDSSPhotoFlags']])
                    
                #End of input, now onto processing and output
                sdssMags = selectStarSED0.deReddenMags(am, [sDSSu, sDSSg, sDSSr, sDSSi, sDSSz], 
                                                       sdssExtCoeffs)[0]

                if 10 <= pop < 15:
                    sedType = 'H'
                elif 15 <= pop < 20:
                    sedType = 'HE'
                else:
                    #Use mlt is sdss r-i > 0.59
                    if (sdssMags[2] - sdssMags[3]) > 0.59:
                        sedType = 'mlt'
                    else:
                        sedType = 'kurucz'

                sEDName, magNorm = selectStarSED0.findSED(listDict[sedType], [sdssMags], ra, dec,
                                                          reddening = False, colors = colorDict[sedType])
                testSED = Sed()
                testSED.setSED(listDict[sedType][positionDict[sEDName[0]]].wavelen, 
                               flambda = listDict[sedType][positionDict[sEDName[0]]].flambda)

                #Load Bandpasses for LSST colors to get colors from matched SEDs
                lsstPhot = phot()
                lsstFilterList = ('u', 'g', 'r', 'i', 'z', 'y')
                lsstPhot.loadBandPassesFromFiles(lsstFilterList)
                lsstPhot.setupPhiArray_dict()

                imSimBand = Bandpass()
                imSimBand.imsimBandpass()
                fluxNorm = testSED.calcFluxNorm(magNorm[0], imSimBand)
                testSED.multiplyFluxNorm(fluxNorm)
                
                lsstMagsUnred = lsstPhot.manyMagCalc_list(testSED)
                #If the extinction value is negative then it will add the reddening back in
                lsstMags = selectStarSED0.deReddenMags((-1.0*am), lsstMagsUnred, lsstExtCoeffs)[0]

                distKpc = self.convDMtoKpc(DM)
                ebv = am / 2.285 #From Schlafly and Finkbeiner for sdssr
                ebvInf = amInf / 2.285
                outFmt = '%i,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%s,%3.7f,' +\
                         '%3.7f,%3.7f,%3.7f,' +\
                         '%3.7f,%3.7f,%3.7f,' +\
                         '%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,' +\
                         '%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,' +\
                         '%3.7f,%i,%3.7f,%3.7f,%3.7f\n'
                outDat = (oID, ra, dec, gall, galb, coordX, coordY, coordZ, sEDName[0], magNorm[0],
                          lsstMags[0], lsstMags[1], lsstMags[2], 
                          lsstMags[3], lsstMags[4], lsstMags[5],
                          sDSSu, sDSSg, sDSSr, sDSSi, sDSSz, absSDSSr,
                          pmRA, pmDec, vRad, pml, pmb, vRadlb, vR, vPhi, vZ,
                          FeH, pop, distKpc, ebv, ebvInf)
                fOut.write(outFmt % outDat)
                lineNum += 1
                if lineNum % 10000 == 0:
                    print str(str(lineNum) + ' done')
