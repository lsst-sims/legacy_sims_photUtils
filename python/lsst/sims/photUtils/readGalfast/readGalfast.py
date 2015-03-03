import os
import gzip
import pyfits
import eups
import itertools
import numpy as np

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
                    mltSubset = None, wdSubset = None, magNormAcc = 1, chunkSize = 10000):
        """
        This is customized for the outputs we currently need for the purposes of consistent output
        It will read in a galfast output file and output desired values for database input into a file

        @param [in] filenameList is a list of the galfast output files that will be loaded and processed.
        Can process fits, gzipped, or txt output from galfast.

        @param [in] outFileList is a list of the names of the output files that will be created. If gzipped
        output is desired simply write the filenames with .gz at the end.

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
        
        @param [in] magNormAcc is the number of decimal places within the magNorm result will be accurate.
        
        @param [in] chunkSize is the size of chunks of lines to be read from the catalog at one time.
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
        
        #Load Bandpasses for LSST colors to get colors from matched SEDs        
        lsstPhot = phot()
        lsstFilterList = ('u', 'g', 'r', 'i', 'z', 'y')
        lsstPhot.loadBandPassesFromFiles(lsstFilterList)
        lsstPhot.setupPhiArray_dict()
        imSimBand = Bandpass()
        imSimBand.imsimBandpass()

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
                gzFile = False
                num_lines = sum(1 for line in open(filename))
            elif filename.endswith('.gz'):
                galfastIn = gzip.open(filename, 'r')
                inFits = False
                gzFile = True
                num_lines = sum(1 for line in gzip.open(filename))
            elif filename.endswith('fits'):
                hdulist = pyfits.open(filename)
                galfastIn = hdulist[1].data
                num_lines = len(galfastIn)
                gzFile = False
                inFits = True

            if outFile.endswith('.txt'):
                fOut = open(outFile, 'w')
            elif outFile.endswith('.gz'):
                fOut = gzip.open(outFile, 'w')
            fOut.write('#oID, ra, dec, gall, galb, coordX, coordY, coordZ, sEDName, magNorm, ' +\
                       'LSSTugrizy, SDSSugriz, absSDSSr, pmRA, pmDec, vRad, pml, pmb, vRadlb, ' +\
                       'vR, vPhi, vZ, FeH, pop, distKpc, ebv, ebvInf\n')
            header_length = 0
            numChunks = 0
            if inFits == False:            
                galfastDict = self.parseGalfast(galfastIn.readline())
                header_length += 1
                header_status = True
                while header_status == True:
                    newLine = galfastIn.readline()
                    if newLine[0] != '#':
                        header_status = False
                    else:
                        header_length += 1
            print 'Total objects = %i' % (num_lines - header_length)
            numChunks = ((num_lines-header_length)/chunkSize) + 1

            for chunk in range(0,numChunks):
                if chunk == numChunks-1:
                    lastChunkSize = (num_lines - header_length) % chunkSize
                    readSize = lastChunkSize
                else:
                    readSize = chunkSize
                oID = np.arange(readSize*chunk, readSize*(chunk+1))
                if inFits:
                    starData = galfastIn[readSize*chunk:(readSize*chunk + readSize)]
                    sDSS = starData.field('SDSSugriz')
                    gall, galb = np.transpose(starData.field('lb'))
                    ra, dec = np.transpose(starData.field('radec'))
                    coordX, coordY, coordZ = np.transpose(starData.field('XYZ'))
                    DM = starData.field('DM')
                    absSDSSr = starData.field('absSDSSr')
                    pop = starData.field('comp')
                    FeH = starData.field('FeH')
                    vR, vPhi, vZ = np.transpose(starData.field('vcyl'))
                    pml, pmb, vRadlb = np.transpose(starData.field('pmlb'))
                    pmRA, pmDec, vRad = np.transpose(starData.field('pmradec'))
                    am = starData.field('Am')
                    amInf = starData.field('AmInf')
                    sdssPhotoFlags = starData.field('SDSSugrizPhotoFlags')
                else:
                    if gzFile == False:
                        with open(filename) as t_in:
                            starData = np.loadtxt(itertools.islice(t_in,((readSize*chunk)+header_length),
                                                                   ((readSize*(chunk+1))+header_length)))
                    else:
                        with gzip.open(filename) as t_in:
                            starData = np.loadtxt(itertools.islice(t_in,((readSize*chunk)+header_length),
                                                                   ((readSize*(chunk+1))+header_length)))
                    starData = np.transpose(starData)
                    gall = starData[galfastDict['l']]
                    galb = starData[galfastDict['b']]
                    ra = starData[galfastDict['ra']]
                    dec = starData[galfastDict['dec']]
                    coordX = starData[galfastDict['X']]
                    coordY = starData[galfastDict['Y']]
                    coordZ = starData[galfastDict['Z']]
                    DM = starData[galfastDict['DM']]
                    absSDSSr = starData[galfastDict['absSDSSr']]
                    pop = starData[galfastDict['comp']]
                    FeH = starData[galfastDict['FeH']]
                    vR = starData[galfastDict['Vr']]
                    vPhi = starData[galfastDict['Vphi']]
                    vZ = starData[galfastDict['Vz']]
                    pml = starData[galfastDict['pml']]
                    pmb = starData[galfastDict['pmb']]
                    vRadlb = starData[galfastDict['vRadlb']]
                    pmRA = starData[galfastDict['pmra']]
                    pmDec = starData[galfastDict['pmdec']]
                    vRad = starData[galfastDict['vRad']]
                    am = starData[galfastDict['Am']]
                    amInf = starData[galfastDict['AmInf']]
                    sDSS = np.transpose(starData[galfastDict['SDSSu']:galfastDict['SDSSz']+1])
                    sDSSPhotoFlags = starData[galfastDict['SDSSPhotoFlags']]
                    
                #End of input, now onto processing and output
                sDSSunred = selectStarSED0.deReddenMags(am, sDSS, sdssExtCoeffs)
                if readSize == 1:
                    ra = np.array([ra])
                    dec = np.array([dec])
                mIn = np.where(((pop < 10) | (pop >= 20)) & (sDSSunred[:,2] - sDSSunred[:,3] > 0.59))
                kIn = np.where(((pop < 10) | (pop >= 20)) & (sDSSunred[:,2] - sDSSunred[:,3] <= 0.59))
                hIn = np.where((pop >= 10) & (pop < 15))
                heIn = np.where((pop >= 15) & (pop < 20))

                sEDNameK, magNormK = selectStarSED0.findSED(listDict['kurucz'], 
                                                            sDSSunred[kIn], ra[kIn], dec[kIn],
                                                            reddening = False, 
                                                            magNormAcc = magNormAcc, 
                                                            colors = colorDict['kurucz'])
                sEDNameM, magNormM = selectStarSED0.findSED(listDict['mlt'], 
                                                            sDSSunred[mIn], ra[mIn], dec[mIn],
                                                            reddening = False, 
                                                            magNormAcc = magNormAcc, 
                                                            colors = colorDict['mlt'])
                sEDNameH, magNormH = selectStarSED0.findSED(listDict['H'], 
                                                            sDSSunred[hIn], ra[hIn], dec[hIn],
                                                            reddening = False, 
                                                            magNormAcc = magNormAcc, 
                                                            colors = colorDict['H'])
                sEDNameHE, magNormHE = selectStarSED0.findSED(listDict['HE'], 
                                                              sDSSunred[heIn], ra[heIn], dec[heIn],
                                                              reddening = False, 
                                                              magNormAcc = magNormAcc, 
                                                              colors = colorDict['HE'])
                chunkNames = np.empty(readSize, dtype = 'S32')
                chunkTypes = np.empty(readSize, dtype = 'S8')
                chunkMagNorms = np.zeros(readSize)
                chunkNames[kIn] = sEDNameK
                chunkTypes[kIn] = 'kurucz'
                chunkMagNorms[kIn] = magNormK
                chunkNames[mIn] = sEDNameM
                chunkTypes[mIn] = 'mlt'
                chunkMagNorms[mIn] = magNormM
                chunkNames[hIn] = sEDNameH
                chunkTypes[hIn] = 'H'
                chunkMagNorms[hIn] = magNormH
                chunkNames[heIn] = sEDNameHE
                chunkTypes[heIn] = 'HE'
                chunkMagNorms[heIn] = magNormHE
                lsstMagsUnred = []
                for sedName, sedType, magNorm in zip(chunkNames, chunkTypes, chunkMagNorms):
                    testSED = Sed()
                    testSED.setSED(listDict[sedType][positionDict[sedName]].wavelen, 
                                   flambda = listDict[sedType][positionDict[sedName]].flambda)
                    fluxNorm = testSED.calcFluxNorm(magNorm, imSimBand)
                    testSED.multiplyFluxNorm(fluxNorm)
                    lsstMagsUnred.append(lsstPhot.manyMagCalc_list(testSED))
                #If the extinction value is negative then it will add the reddening back in
                lsstMags = selectStarSED0.deReddenMags((-1.0*am), lsstMagsUnred, 
                                                       lsstExtCoeffs)
                distKpc = self.convDMtoKpc(DM)
                ebv = am / 2.285 #From Schlafly and Finkbeiner for sdssr
                ebvInf = amInf / 2.285
                for line in range(0, readSize):
                    outFmt = '%i,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,' +\
                             '%3.7f,%3.7f,%s,%3.7f,' +\
                             '%3.7f,%3.7f,%3.7f,' +\
                             '%3.7f,%3.7f,%3.7f,' +\
                             '%3.7f,%3.7f,%3.7f,%3.7f,' +\
                             '%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,' +\
                             '%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,' +\
                             '%3.7f,%i,%3.7f,%3.7f,%3.7f\n'
                    if readSize == 1:
                        if inFits == True:
                            sDSS = sDSS[0]
                        outDat = (oID, ra[line], dec[line], gall, galb, coordX, 
                                  coordY, coordZ, chunkNames, chunkMagNorms,
                                  lsstMags[line][0], lsstMags[line][1], lsstMags[line][2], 
                                  lsstMags[line][3], lsstMags[line][4], lsstMags[line][5], 
                                  sDSS[0], sDSS[1], sDSS[2], sDSS[3], 
                                  sDSS[4], absSDSSr, pmRA, pmDec, vRad, 
                                  pml, pmb, vRadlb, vR, vPhi, vZ,
                                  FeH, pop, distKpc, ebv, ebvInf)
                    else:
                        outDat = (oID[line], ra[line], dec[line], gall[line], galb[line], coordX[line], 
                                  coordY[line], coordZ[line], chunkNames[line], chunkMagNorms[line],
                                  lsstMags[line][0], lsstMags[line][1], lsstMags[line][2], 
                                  lsstMags[line][3], lsstMags[line][4], lsstMags[line][5], 
                                  sDSS[line][0], sDSS[line][1], sDSS[line][2], sDSS[line][3], 
                                  sDSS[line][4], absSDSSr[line], pmRA[line], pmDec[line], vRad[line], 
                                  pml[line], pmb[line], vRadlb[line], vR[line], vPhi[line], vZ[line],
                                  FeH[line], pop[line], distKpc[line], ebv[line], ebvInf[line])
                    fOut.write(outFmt % outDat)
                print 'Chunk Num Done = %i out of %i' % (chunk+1, numChunks)