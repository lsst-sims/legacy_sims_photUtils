import unittest
import os
import numpy as np
import gzip
import pyfits
from lsst.sims.photUtils.readGalfast.selectStarSED import selectStarSED
from lsst.sims.photUtils.readGalfast.readGalfast import readGalfast
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.photUtils import Photometry as phot

class TestSelectStarSED(unittest.TestCase):

    def testLoadKurucz(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for"""
        #Test Matching to Kurucz SEDs
        testSEDs = selectStarSED().loadKuruczSEDs()

        #Read in a list of the SEDs in the kurucz sims sed directory
        testKuruczList = os.listdir(os.environ['SIMS_SED_LIBRARY_DIR'] + '/starSED/kurucz/')

        #First make sure that all SEDs are correctly accounted for if no subset provided
        self.assertEqual(testSEDs['sEDName'], testKuruczList)

        #Then test to make sure no values are null in any of the lists within the dictionary
        for key in testSEDs.keys():
            for item in testSEDs[key]:
                self.assertIsNotNone(item)

        #Test same condition if subset is provided
        testSubsetList = ['km01_7000.fits_g40_7140.gz', 'kp01_7000.fits_g40_7240.gz']
        testSEDsSubset = selectStarSED().loadKuruczSEDs(subset = testSubsetList)

        #Next make sure that correct subset loads if subset is provided
        self.assertEqual(testSEDsSubset['sEDName'], testSubsetList)

        #Finally make sure that values are being entered into dictionaries correctly
        self.assertEqual(testSEDsSubset['logZ'], [-0.1, 0.1]) #Test both pos. and neg. get in right
        self.assertEqual(testSEDsSubset['logg'], [4.0, 4.0])
        self.assertEqual(testSEDsSubset['temp'], [7140, 7240])

        #Now test colors are being calculated correctly by loading a flat sed into it
        testSED = Sed()
        testSEDsColors = selectStarSED()

        testSED.setFlatSED()
        testSED.multiplyFluxNorm(testSED.calcFluxNorm(10, testSEDsColors.sdssBandpassDict['r']))
        #Give kurucz like filename just so it works
        testName = 'km99_9999.fits_g99_9999'
        testSED.writeSED(testName)

        testSEDsColors.kuruczDir = str(os.getcwd() + '/')
        testColors = testSEDsColors.loadKuruczSEDs(subset = [testName])
        self.assertAlmostEqual(testColors['umg'][0], 0.0)
        self.assertAlmostEqual(testColors['gmr'][0], 0.0)
        self.assertAlmostEqual(testColors['rmi'][0], 0.0)
        self.assertAlmostEqual(testColors['imz'][0], 0.0)

    def testLoadMLT(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for"""
        #Test Matching to mlt SEDs
        testSEDs = selectStarSED().loadmltSEDs()

        #Read in a list of the SEDs in the kurucz sims sed directory
        testMLTList = os.listdir(os.environ['SIMS_SED_LIBRARY_DIR'] + '/starSED/mlt/')

        #First make sure that all SEDs are correctly accounted for if no subset provided
        self.assertItemsEqual(testSEDs['sEDName'], testMLTList)

        #Then test to make sure no values are null in any of the lists within the dictionary
        for key in testSEDs.keys():
            for item in testSEDs[key]:
                self.assertIsNotNone(item)

        #Test same condition if subset is provided
        testSubsetList = ['m0.0Full.dat.gz']
        testSEDsSubset = selectStarSED().loadmltSEDs(subset = testSubsetList)

        #Next make sure that correct subset loads if subset is provided
        self.assertItemsEqual(testSEDsSubset['sEDName'], testSubsetList)

        #Now test colors are being calculated correctly by loading a flat sed into it
        testSED = Sed()
        testSEDsColors = selectStarSED()

        testSED.setFlatSED()
        testSED.multiplyFluxNorm(testSED.calcFluxNorm(10, testSEDsColors.sdssBandpassDict['r']))
        #Give mlt like filename just so it works
        testName = 'm99.99Full.dat'
        testSED.writeSED(testName)

        testSEDsColors.mltDir = str(os.getcwd() + '/')
        testColors = testSEDsColors.loadmltSEDs(subset = [testName])
        self.assertAlmostEqual(testColors['umg'][0], 0.0)
        self.assertAlmostEqual(testColors['gmr'][0], 0.0)
        self.assertAlmostEqual(testColors['rmi'][0], 0.0)
        self.assertAlmostEqual(testColors['imz'][0], 0.0)

    def testLoadWD(self):
        """Test SED loading algorithm by making sure SEDs are all accounted for and
        values for each property have been calculated."""
        #Test Matching to WD SEDs
        testSEDs = selectStarSED().loadwdSEDs()

        #Add extra step because WD SEDs are separated into helium and hydrogen
        testSEDNamesLists = []
        testSEDNamesLists.append(testSEDs['H']['sEDName'])
        testSEDNamesLists.append(testSEDs['HE']['sEDName'])
        testSEDNames = [name for nameList in testSEDNamesLists for name in nameList]

        #Read in a list of the SEDs in the kurucz sims sed directory
        testWDList = os.listdir(os.environ['SIMS_SED_LIBRARY_DIR'] + '/starSED/wDs/')

        #First make sure that all SEDs are correctly accounted for if no subset provided
        self.assertItemsEqual(testSEDNames, testWDList)

        #Then test to make sure no values are null in any of the lists within the dictionary
        for wdType in ['H', 'HE']:
            for key in testSEDs[wdType].keys():
                for item in testSEDs[wdType][key]:
                    self.assertIsNotNone(item)

        #Test same condition if subset is provided
        testSubsetList = ['bergeron_10000_75.dat_10100.gz', 'bergeron_He_9000_80.dat_9400.gz']

        testSEDsSubset = selectStarSED().loadwdSEDs(subset = testSubsetList)

        #Add extra step because WD SEDs are separated into helium and hydrogen
        testSEDNamesSubsetLists = []
        testSEDNamesSubsetLists.append(testSEDsSubset['H']['sEDName'])
        testSEDNamesSubsetLists.append(testSEDsSubset['HE']['sEDName'])
        testSEDNamesSubset = [name for nameList in testSEDNamesSubsetLists for name in nameList]

        #Next make sure that correct subset loads if subset is provided
        self.assertItemsEqual(testSEDNamesSubset, testSubsetList)

        #Make sure that the names get separated into correct wd type
        self.assertItemsEqual(testSEDsSubset['H']['sEDName'], [testSubsetList[0]])
        self.assertItemsEqual(testSEDsSubset['HE']['sEDName'], [testSubsetList[1]])

        #Now test colors are being calculated correctly by loading a flat sed into it
        testSED = Sed()
        testSEDsColors = selectStarSED()

        testSED.setFlatSED()
        testSED.multiplyFluxNorm(testSED.calcFluxNorm(10, testSEDsColors.sdssBandpassDict['r']))
        #Give WD like filename just so it works
        testName = 'bergeron_9999_99.dat_9999'
        testNameHe = 'bergeron_He_9999_99.dat_9999'
        testSED.writeSED(testName)
        testSED.writeSED(testNameHe)

        testSEDsColors.wdDir = str(os.getcwd() + '/')
        testColors = testSEDsColors.loadwdSEDs(subset = [testName, testNameHe])
        self.assertAlmostEqual(testColors['H']['umg'][0], 0.0)
        self.assertAlmostEqual(testColors['H']['gmr'][0], 0.0)
        self.assertAlmostEqual(testColors['H']['rmi'][0], 0.0)
        self.assertAlmostEqual(testColors['H']['imz'][0], 0.0)
        self.assertAlmostEqual(testColors['HE']['umg'][0], 0.0)
        self.assertAlmostEqual(testColors['HE']['gmr'][0], 0.0)
        self.assertAlmostEqual(testColors['HE']['rmi'][0], 0.0)
        self.assertAlmostEqual(testColors['HE']['imz'][0], 0.0)

    def testDeReddenGalfast(self):
        
        """Test that consistent numbers come out of deReddening procedure"""

        am = 0.5
        coeffs = np.ones(5)
        mags = np.arange(2,-3,-1)

        testDeRedColors = selectStarSED().deReddenGalfast(am, mags[0], mags[1], mags[2], mags[3], mags[4], coeffs)
        
        #Test Output
        expectedDeRed = np.ones(4)
        np.testing.assert_equal(testDeRedColors, expectedDeRed)

    def testFindSED(self):
        
        """Pull one SED from each type and make sure that it gets matched to."""

        testMatching = selectStarSED()

        testSubsetList = []
        testSubsetType = []
        testSubsetComp = []
        testOutputName = []
        testOutputNameReddened = []
        #Populate Lists
        for fileName in os.listdir(testMatching.kuruczDir):
            testSubsetList.append(fileName)
            testSubsetType.append('kurucz')
            testSubsetComp.append(1)
        for fileName in os.listdir(testMatching.mltDir):
            testSubsetList.append(fileName)
            testSubsetType.append('mlt')
            testSubsetComp.append(1)
        for fileName in os.listdir(testMatching.wdDir):
            testSubsetList.append(fileName)
            testSubsetType.append('wDs')
            if 'He' in fileName:
                testSubsetComp.append(15)
            else:
                testSubsetComp.append(10)
        
        testMatchingDict = {}
        testMatchingDict['kurucz'] = testMatching.loadKuruczSEDs()
        testMatchingDict['mlt'] = testMatching.loadmltSEDs()
        testMatchingWDs = testMatching.loadwdSEDs()
        testMatchingDict['wdH'] = testMatchingWDs['H']
        testMatchingDict['wdHE'] = testMatchingWDs['HE']

        for testSEDNum in range(0, len(testSubsetList)):

            if testSEDNum % 100 == 0:
                print 'Calculating test colors for SED %i of %i' % (testSEDNum, len(testSubsetList))

            getSEDColors = Sed()
            getSEDColors.readSED_flambda(str(testMatching.sEDDir + '/starSED/' + testSubsetType[testSEDNum] + '/' + testSubsetList[testSEDNum]))
            getSEDColors.multiplyFluxNorm(getSEDColors.calcFluxNorm(10, testMatching.sdssBandpassDict['r']))
            testSEDPhotometry = phot()
            testMagDict = testSEDPhotometry.manyMagCalc_dict(getSEDColors, testMatching.sdssPhiArray, testMatching.sdssWavelenstep, testMatching.sdssBandpassDict, testMatching.sdssFilterList)

            #First test without reddening
            testOutputName.append(testMatching.findSED(testMatchingDict, testMagDict['u'], testMagDict['g'], testMagDict['r'], testMagDict['i'], testMagDict['z'], 0, testSubsetComp[testSEDNum], reddening = False))
            #Next test with reddening and custom coeffs
            am = 0.5
            reddenCoeffs = np.array([1.2, 1.1, 1.0, 0.9, 0.8])
            testReddening = am * reddenCoeffs
            testReddenedMagDict = {}
            for filter, coeffNum in zip(testMatching.sdssFilterList, range(0, len(testReddening))):
                testReddenedMagDict[filter] = testMagDict[filter] + testReddening[coeffNum]
            testOutputNameReddened.append(testMatching.findSED(testMatchingDict, testReddenedMagDict['u'], testReddenedMagDict['g'], testReddenedMagDict['r'], testReddenedMagDict['i'], testReddenedMagDict['z'], am, testSubsetComp[testSEDNum], True, reddenCoeffs))
        self.assertEqual(testOutputName, testSubsetList)
        self.assertEqual(testOutputNameReddened, testSubsetList)

class TestReadGalfast(unittest.TestCase):

    def testParseGalfast(self):

        """Test Read-in of Galfast Header"""
        
        testRG = readGalfast()
        #First test that exception is raised when an invalid header label is used
        testInvalidHeader = 'lb[2] header'
        self.assertRaises(RuntimeError, testRG.parseGalfast, testInvalidHeader)

        #Next test that '#' is ignored
        testSymbolHeader = '# lb[2]'
        testSymbolDict = testRG.parseGalfast(testSymbolHeader)
        self.assertEqual(testSymbolDict, {'l':0, 'b':1})

        #Test that new line marker at end of line is ignored
        testNewLineHeader = '# lb[2] \n '
        testNewLineDict = testRG.parseGalfast(testNewLineHeader)
        self.assertEqual(testNewLineDict, {'l':0, 'b':1})

        #Next test that extra spaces are ignored
        testSpaceHeader = 'lb[2]     radec[2]'
        testSpaceDict = testRG.parseGalfast(testSpaceHeader)
        self.assertEqual(testSpaceDict, {'l':0, 'b':1, 'ra':2, 'dec':3})

        #Now test that every header value gets put in the right place using different order than usually come out
        testFullHeader = '# lb[2] XYZ[3] radec[2] absSDSSr{alias=M1;alias=absmag;band=SDSSr;} DM comp FeH vcyl[3] pmlb[3] pmradec[3] Am AmInf SDSSugriz[5]{class=magnitude;fieldNames=0:SDSSu,1:SDSSg,2:SDSSr,3:SDSSi,4:SDSSz;} SDSSugrizPhotoFlags{class=flags;}'
        actualFullHeaderDict = {'l':0, 'b':1, 'X':2, 'Y':3, 'Z':4, 'ra':5, 'dec':6, 'absSDSSr':7, 'DM':8, 'comp':9, 'FeH':10, 'Vr':11, 'Vphi':12, 'Vz':13, 'pml':14, 'pmb':15, 'vRadlb':16, 'pmra':17, 'pmdec':18, 'vRad':19, 'Am':20, 'AmInf':21, 'SDSSu':22, 'SDSSg':23, 'SDSSr':24, 'SDSSi':25, 'SDSSz':26, 'SDSSPhotoFlags':27}
        testFullHeaderDict = testRG.parseGalfast(testFullHeader)
        self.assertEqual(testFullHeaderDict, actualFullHeaderDict)

    def testFindLSSTMags(self):
        
        """Make sure that the Magnitudes you get from spectra are correct"""
        
        testRG = readGalfast()
        sEDParams = selectStarSED()
        testPhot = phot()
        absLSSTr = 10.
        DM = 10.
        am = 0.5
        lsstExtCoords = [1.25, 1.15, 1.05, 0.95, 0.85, 0.75]
        
        #Test a kurucz, an mlt, and a wd
        testSpectra = ['bergeron_10000_75.dat_10100.gz', 'm0.0Full.dat.gz', 'km01_7000.fits_g40_7140.gz']
        testTypes = ['wDs', 'mlt', 'kurucz']
        for testSpectrum, testType in zip(testSpectra, testTypes):
            sEDObj = Sed()
            sEDObj.readSED_flambda(sEDParams.sEDDir + '/starSED/' + testType + '/' + testSpectrum)
            #Calculate the SDSSr if the SED's LSSTr is set to absLSSTr above
            #Then Make sure with this input you get the same LSSTr out of findLSSTMags
            lsstFluxNorm = sEDObj.calcFluxNorm(absLSSTr, sEDParams.lsstBandpassDict['r'])
            sEDObj.multiplyFluxNorm(lsstFluxNorm)
            absSDSSr = sEDObj.calcMag(sEDParams.sdssBandpassDict['r'])
            testMags, testFluxNorm = testRG.findLSSTMags(testSpectrum, absSDSSr, DM, am, reddening=False)
            self.assertAlmostEqual(testMags['r'] - DM, absLSSTr)
            self.assertAlmostEqual(testFluxNorm, lsstFluxNorm)
            testMagsReddened, testFluxNormReddened = testRG.findLSSTMags(testSpectrum, absSDSSr, DM, am, True, lsstExtCoords)
            self.assertAlmostEqual(testMagsReddened['r'] - DM - (am * lsstExtCoords[2]), absLSSTr)
            self.assertAlmostEqual(testFluxNormReddened, lsstFluxNorm)
            

    def testConvDMtoKpc(self):

        """Make sure Distance Modulus get correctly converted to distance in kpc"""

        testRG = readGalfast()
        dm = 10.
        actualKpcDistance = 1.
        testKpcDistance = testRG.convDMtoKpc(dm)
        self.assertEqual(testKpcDistance, actualKpcDistance)
        
    def testLoadGalfast(self):

        """Make sure all desired input file types are correctly processed and return correct output files"""
        
        testRG = readGalfast()
        #First test that it makes sure file exists
        self.assertRaises(RuntimeError, testRG.loadGalfast, 'notarealfile.txt', 'noOutput.txt')

        #Next test that if an unknown file format is entered it exits
        self.assertRaises(RuntimeError, testRG.loadGalfast, 'notarealfile.dat', 'noOutput.txt')
        
        #Write example files and then load in and make sure example output files are created
        #First .txt
        exampleIn = open('example.txt', 'w')
        inHeader = '# lb[2] XYZ[3] radec[2] absSDSSr{alias=M1;alias=absmag;band=SDSSr;} DM comp FeH vcyl[3] pmlb[3] pmradec[3] Am AmInf SDSSugriz[5]{class=magnitude;fieldNames=0:SDSSu,1:SDSSg,2:SDSSr,3:SDSSi,4:SDSSz;} SDSSugrizPhotoFlags{class=flags;} \n'
        testComment = '# Comment\n'
        inData = '   1.79371816  -89.02816704   11.92064832  -27.62775082       7.15       0.22    -421.87   8.126   4.366   0 -0.095    13.7  -183.4    -6.2   -20.58   -12.60    13.02    21.34   -11.26    13.02  0.037  0.037  14.350  12.949  12.529  12.381  12.358 0'
        exampleIn.write(inHeader)
        exampleIn.write(testComment)
        exampleIn.write(inData)
        exampleIn.close()
        testRG.loadGalfast('example.txt', 'exampleOutput.txt')
        self.assertTrue(os.path.isfile('exampleOutput.txt'))
        #Then gzipped
        exampleGzipIn = gzip.open('gzipExample.txt.gz', 'w')
        exampleGzipIn.write(inHeader)
        exampleGzipIn.write(testComment)
        exampleGzipIn.write(inData)
        exampleGzipIn.close()
        testRG.loadGalfast('gzipExample.txt.gz', 'exampleOutputGzip.txt')
        self.assertTrue(os.path.isfile('exampleOutputGzip.txt'))
        #Finally a fits file
        columnNames = ['lb', 'XYZ', 'radec', 'absSDSSr', 'DM', 'comp', 'FeH', 'vcyl', 'pmlb', 'pmradec', 'Am', 'AmInf', 'SDSSugriz', 'SDSSugrizPhotoFlags']
        columnArrays = [[[1.79371816, -89.02816704]],  [[7.15, 0.22, -421.87]], [[11.92064832, -27.62775082]], [[8.126]], [[4.366]], [[0]], [[-0.095]], [[13.7,  -183.4, -6.2]], [[-20.58, -12.60, 13.02]], [[21.34, -11.26, 13.02]],[[0.037]], [[0.037]],  [[14.350, 12.949, 12.529, 12.381, 12.358]], [[0]]]
        columnFormats = ['2E', '3E', '2E', 'E', 'E', 'E', 'E', '3E', '3E', '3E', 'E', 'E', '5E', 'E']
        cols = pyfits.ColDefs([pyfits.Column(name = columnNames[0], format = columnFormats[0], array = columnArrays[0])])
        for colName, colArray, colFormat in zip(columnNames[1:], columnArrays[1:], columnFormats[1:]):
            cols.add_col(pyfits.Column(name = colName, format = colFormat, array = colArray))
        exampleTable = pyfits.new_table(cols)
        exampleTable.writeto('exampleFits.fits')
        testRG.loadGalfast('exampleFits.fits', 'exampleOutputFits.txt')
        self.assertTrue(os.path.isfile('exampleOutputFits.txt'))        

if __name__ == "__main__":

    suite = unittest.TestLoader().loadTestsFromTestCase(TestSelectStarSED)
    unittest.TextTestRunner(verbosity=2).run(suite)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestReadGalfast)
    unittest.TextTestRunner(verbosity=2).run(suite)

