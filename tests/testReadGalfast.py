from builtins import zip
from builtins import str
import unittest
import os
import gzip
import re
from astropy.io import fits
import lsst.utils
import lsst.utils.tests
from lsst.sims.photUtils.readGalfast.readGalfast import readGalfast
from lsst.utils import getPackageDir


def setup_module(module):
    lsst.utils.tests.init()


class TestReadGalfast(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # Set up Test Spectra Directory
        cls.testSpecDir = os.path.join(getPackageDir('sims_photUtils'),
                                       'tests/cartoonSedTestData/starSed/')
        cls.testKDir = str(cls.testSpecDir + '/kurucz/')
        cls.testMLTDir = str(cls.testSpecDir + '/mlt/')
        cls.testWDDir = str(cls.testSpecDir + '/wDs/')

    @classmethod
    def tearDownClass(cls):
        del cls.testSpecDir
        del cls.testKDir
        del cls.testMLTDir
        del cls.testWDDir

    def testParseGalfast(self):

        """Test Read-in of Galfast Header"""

        testRG = readGalfast()
        # First test that exception is raised when an invalid header label is used
        testInvalidHeader = 'lb[2] header'
        self.assertRaises(RuntimeError, testRG.parseGalfast, testInvalidHeader)

        # Next test that '#' is ignored
        testSymbolHeader = '# lb[2]'
        testSymbolDict = testRG.parseGalfast(testSymbolHeader)
        self.assertEqual(testSymbolDict, {'l': 0, 'b': 1})

        # Test that new line marker at end of line is ignored
        testNewLineHeader = '# lb[2] \n '
        testNewLineDict = testRG.parseGalfast(testNewLineHeader)
        self.assertEqual(testNewLineDict, {'l': 0, 'b': 1})

        # Next test that extra spaces are ignored
        testSpaceHeader = 'lb[2]     radec[2]'
        testSpaceDict = testRG.parseGalfast(testSpaceHeader)
        self.assertEqual(testSpaceDict, {'l': 0, 'b': 1, 'ra': 2, 'dec': 3})

        # Test that every header value gets put in the right place using different order than usually come out
        testFullHeader = '# lb[2] XYZ[3] radec[2] absSDSSr{alias=M1;alias=absmag;band=SDSSr;} ' +\
                         'DM comp FeH vcyl[3] pmlb[3] pmradec[3] Am AmInf SDSSugriz[5]{class=magnitude;' +\
                         'fieldNames=0:SDSSu,1:SDSSg,2:SDSSr,3:SDSSi,4:SDSSz;} ' +\
                         'SDSSugrizPhotoFlags{class=flags;}'
        actualFullHeaderDict = {'l': 0, 'b': 1, 'X': 2, 'Y': 3, 'Z': 4,
                                'ra': 5, 'dec': 6, 'absSDSSr': 7, 'DM': 8,
                                'comp': 9, 'FeH': 10, 'Vr': 11, 'Vphi': 12,
                                'Vz': 13, 'pml': 14, 'pmb': 15, 'vRadlb': 16,
                                'pmra': 17, 'pmdec': 18, 'vRad': 19, 'Am': 20,
                                'AmInf': 21, 'SDSSu': 22, 'SDSSg': 23,
                                'SDSSr': 24, 'SDSSi': 25, 'SDSSz': 26,
                                'SDSSPhotoFlags': 27}
        testFullHeaderDict = testRG.parseGalfast(testFullHeader)
        self.assertEqual(testFullHeaderDict, actualFullHeaderDict)

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
        # First test that it makes sure file exists
        self.assertRaises(RuntimeError, testRG.loadGalfast, ['notarealfile.txt'], ['noOutput.txt'])

        # Next test that if an unknown file format is entered it exits
        self.assertRaises(RuntimeError, testRG.loadGalfast, ['notarealfile.dat'], ['noOutput.txt'])

        # Write example files and then load in and make sure example output files are created
        with lsst.utils.tests.getTempFilePath('.in.txt') as inTxtName:
            with lsst.utils.tests.getTempFilePath('.in.txt.gz') as inGzipName:
                with lsst.utils.tests.getTempFilePath('.in.fits') as inFitsName:

                    # First write .txt
                    with open(inTxtName, 'w') as exampleIn:
                        inHeader = '# lb[2] radec[2] XYZ[3] absSDSSr{alias=M1;alias=absmag;band=SDSSr;} DM comp FeH ' +\
                                   'vcyl[3] pmlb[3] pmradec[3] Am AmInf SDSSugriz[5]{class=magnitude;fieldNames=0:SDSSu,' +\
                                   '1:SDSSg,2:SDSSr,3:SDSSi,4:SDSSz;} SDSSugrizPhotoFlags{class=flags;} \n'
                        testComment = '# Comment\n'
                        inData = '   1.79371816  -89.02816704   11.92064832  -27.62775082       7.15       0.22   ' +\
                                 '-421.87   8.126   4.366   0 -0.095    13.7  -183.4    -6.2   -20.58   -12.60    ' +\
                                 '13.02    21.34   -11.26    13.02  0.037  0.037  14.350  12.949  12.529  12.381  12.358 0\n'
                        exampleIn.write(inHeader)
                        exampleIn.write(testComment)
                        exampleIn.write(inData)

                    # Then gzipped. Also testing multiple lines in catalog.
                    with gzip.open(inGzipName, 'wt') as exampleGzipIn:
                        exampleGzipIn.write(inHeader)
                        exampleGzipIn.write(testComment)
                        exampleGzipIn.write(inData)
                        exampleGzipIn.write(inData)

                    # Finally a fits file, but first make sure to remove pre-existing file
                    columnNames = ['lb', 'XYZ', 'radec', 'absSDSSr', 'DM', 'comp', 'FeH', 'vcyl', 'pmlb', 'pmradec',
                                   'Am', 'AmInf', 'SDSSugriz', 'SDSSugrizPhotoFlags']
                    columnArrays = [[[1.79371816, -89.02816704]], [[7.15, 0.22, -421.87]],
                                    [[11.92064832, -27.62775082]], [[8.126]], [[4.366]], [[0]], [[-0.095]],
                                    [[13.7, -183.4, -6.2]], [[-20.58, -12.60, 13.02]], [[21.34, -11.26, 13.02]],
                                    [[0.037]], [[0.037]], [[14.350, 12.949, 12.529, 12.381, 12.358]], [[0]]]
                    columnFormats = ['2E', '3E', '2E', 'E', 'E', 'E', 'E', '3E', '3E', '3E', 'E', 'E', '5E', 'E']
                    cols = fits.ColDefs([fits.Column(name = columnNames[0], format = columnFormats[0],
                                                     array = columnArrays[0])])
                    for colName, colArray, colFormat in zip(columnNames[1:], columnArrays[1:], columnFormats[1:]):
                        cols.add_col(fits.Column(name = colName, format = colFormat, array = colArray))
                    exampleTable = fits.BinTableHDU.from_columns(cols)

                    exampleTable.writeto(inFitsName)

                    with lsst.utils.tests.getTempFilePath('.fits.txt') as outFitsName:
                        with lsst.utils.tests.getTempFilePath('.txt') as outTxtName:
                            with lsst.utils.tests.getTempFilePath('.gz') as outGzipName:
                                testRG.loadGalfast([inTxtName, inGzipName, inFitsName],
                                                   [outTxtName, outGzipName, outFitsName],
                                                   kuruczPath = self.testKDir,
                                                   mltPath = self.testMLTDir,
                                                   wdPath = self.testWDDir)
                                self.assertTrue(os.path.isfile(outTxtName),
                                                msg='file .txt output file was not created')
                                self.assertTrue(os.path.isfile(outGzipName),
                                                msg='file gzip output file was not created')
                                self.assertTrue(os.path.isfile(outFitsName),
                                                msg='file fit.txt output file was not created')


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
