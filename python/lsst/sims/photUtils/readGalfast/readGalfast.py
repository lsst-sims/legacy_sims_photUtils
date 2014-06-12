import os
import gzip
import pyfits
import numpy as np

from lsst.sims.photUtils import Sed
from lsst.sims.photUtils import Bandpass
from lsst.sims.photUtils.photUtils import Photometry as phot
from lsst.sims.photUtils.readGalfast.selectStarSED import selectStarSED

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

    def findLSSTMags(self, sEDName, absMagR, DM, am, reddening = True, 
                     lsstExtCoeffs = [1.8140, 1.4166, 0.9947, 0.7370, 0.5790, 0.4761]):

        """
        Takes input SED and object's galfast output properties and generates appropriate LSST mags

        @param [in] sEDname is the matched SED file to the galfast object

        @param [in] absMagR is the galfast output's absolute SDSS R magnitude

        @param [in] DM is the distance modulus from the galfast output

        @param [in] am is the extinction from the galfast output

        @param [in] reddening indicates whether extinction and reddening were included in the
        galfast output

        @param [in] lsstExtCoeffs are the appropriately scaled reddening coefficients compared to
        SDSS R from schlafly and finkbeiner 2011, table 2

        @param [out] lsstFinalMagDict is a dictionary with the filter names and lsst magnitudes

        @param [out] sED_fluxnorm is the flux normalization for the object and the given absMagR
        """

        sEDObj = Sed()
        setupSEDParams = selectStarSED()
        if sEDName.startswith('k'):
            #Kurucz SED
            sEDObj.readSED_flambda(setupSEDParams.kuruczDir + sEDName)
        elif sEDName.startswith('bergeron'):
            #WD SED
            sEDObj.readSED_flambda(setupSEDParams.wdDir + sEDName)
        else:
            #mlt SED
            sEDObj.readSED_flambda(setupSEDParams.mltDir + sEDName)

        sED_fluxnorm = sEDObj.calcFluxNorm(absMagR, setupSEDParams.sdssBandpassDict['r'])
        sEDObj.multiplyFluxNorm(sED_fluxnorm)
        
        sEDPhot = phot()
        lsstMagDict = sEDPhot.manyMagCalc_dict(sEDObj, setupSEDParams.lsstPhiArray, 
                                               setupSEDParams.lsstWavelenstep, 
                                               setupSEDParams.lsstBandpassDict, 
                                               setupSEDParams.lsstFilterList)
        
        #Add DM and reddening the way mario does
        #Find proportional reddening coeffecients to sdssr reddening 
        #from schlafly and finkbeiner 2011, table 2

        lsstFinalMagDict = {}
        for filterName, ext in zip(setupSEDParams.lsstFilterList, lsstExtCoeffs):
            mag = lsstMagDict[filterName]
            if reddening == True:
                lsstFinalMagDict[filterName] = mag + DM + (am * ext)
            else:
                lsstFinalMagDict[filterName] = mag + DM

        return lsstFinalMagDict, sED_fluxnorm 

    def convDMtoKpc(self, DM): 
        """
        Change from distance modulus to distance in kiloparsecs

        @param [in] DM is the distance modulus

        @param [out] distanceKpc is the distance in kiloparsecs

        """

        distancePc = 10**((0.2*DM) + 1)
        distanceKpc = distancePc / 1000.
        return distanceKpc

    def loadGalfast(self, filenameList, outFileList):
        """
        This is customized for the outputs we currently need for the purposes of consistent output
        It will read in a galfast output file and output desired values for database input into a file

        @param [in] filenameList is a list of the galfast output files that will be loaded and processed

        @param [in] outFileList is a list of the names of the output files that will be created
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
        sEDDict = {}
        
        selectStarSED0 = selectStarSED()

        sEDDict['kurucz'] = selectStarSED0.loadKuruczSEDs()
        sEDDict['mlt'] = selectStarSED0.loadmltSEDs()
        wdDict = selectStarSED0.loadwdSEDs()
        sEDDict['wdH'] = wdDict['H']
        sEDDict['wdHE'] = wdDict['HE']

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
            fOut.write('#ra, dec, gall, galb, coordX, coordY, coordZ, sEDName, fluxNorm, ' +\
                       'LSSTugrizy, SDSSugriz, pmRA, pmDec, vRad, pml, pmb, vRadlb, ' +\
                       'vR, vPhi, vZ, FeH, pop, distKpc, ebv, ebvInf\n')

            lineNum = 0
            for line in galfastIn:
                if inFits:
                    oID = float(lineNum)
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
                    oID = float(lineNum)
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
                sEDName = selectStarSED0.findSED(sEDDict, sDSSu, sDSSg, sDSSr, sDSSi, sDSSz, am, pop)
                lsstMagDict, fluxNorm = self.findLSSTMags(sEDName, absSDSSr, DM, am)
                distKpc = self.convDMtoKpc(DM)
                ebv = am / 2.285 #From Schlafly and Finkbeiner for sdssr
                ebvInf = amInf / 2.285
                outFmt = '%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%s,%3.7e,' +\
                         '%3.7f,%3.7f,%3.7f,' +\
                         '%3.7f,%3.7f,%3.7f,' +\
                         '%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,' +\
                         '%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f,%3.7f' +\
                         '%3.7f,%i,%3.7f,%3.7f,%3.7f\n'
                outDat = (ra, dec, gall, galb, coordX, coordY, coordZ, sEDName, fluxNorm,
                          lsstMagDict['u'], lsstMagDict['g'], lsstMagDict['r'], 
                          lsstMagDict['i'], lsstMagDict['z'], lsstMagDict['y'],
                          sDSSu, sDSSg, sDSSr, sDSSi, sDSSz, absSDSSr,
                          pmRA, pmDec, vRad, pml, pmb, vRadlb, vR, vPhi, vZ,
                          FeH, pop, distKpc, ebv, ebvInf)
                fOut.write(outFmt % outDat)
                lineNum += 1
