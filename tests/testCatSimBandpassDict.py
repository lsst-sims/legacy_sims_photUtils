from __future__ import with_statement
import unittest
import os
import copy
import numpy
import lsst.utils.tests as utilsTests
from lsst.utils import getPackageDir

from lsst.sims.photUtils import Bandpass, Sed, CatSimBandpassDict, \
                                CatSimSedList

class BandpassDictTest(unittest.TestCase):

    def setUp(self):
        numpy.random.seed(32)
        self.bandpassPossibilities = ['u', 'g', 'r', 'i', 'z', 'y']
        self.bandpassDir = os.path.join(getPackageDir('throughputs'), 'baseline')
        self.sedDir = os.path.join(getPackageDir('sims_sed_library'))
        self.sedDir = os.path.join(self.sedDir, 'galaxySED')
        self.sedPossibilities = os.listdir(self.sedDir)


    def getListOfSedNames(self, nNames):
        return [self.sedPossibilities[ii].replace('.gz','') \
                for ii in \
                numpy.random.random_integers(0, len(self.sedPossibilities)-1, nNames)]


    def getListOfBandpasses(self, nBp):
        """
        Generate a list of nBp bandpass names and bandpasses

        Intentionally do so a nonsense order so that we can test
        that order is preserved in the CatSimBandpassDict
        """
        dexList = numpy.random.random_integers(0, len(self.bandpassPossibilities)-1, nBp)
        bandpassNameList = []
        bandpassList = []
        for dex in dexList:
            name = self.bandpassPossibilities[dex]
            bp = Bandpass()
            bp.readThroughput(os.path.join(self.bandpassDir,'total_%s.dat' % name))
            while name in bandpassNameList:
                name += '0'
            bandpassNameList.append(name)
            bandpassList.append(bp)

        return bandpassNameList, bandpassList


    def testInitialization(self):
        """
        Test that all of the member variables of CatSimBandpassDict are set
        to the correct value upon construction.
        """

        for nBp in range(3,10,1):
            nameList, bpList = self.getListOfBandpasses(nBp)
            testDict = CatSimBandpassDict(bpList, nameList)

            self.assertEqual(testDict.nBandpasses, nBp)
            self.assertEqual(len(testDict), nBp)

            for controlName, testName in zip(nameList, testDict):
                self.assertEqual(controlName, testName)

            for controlName, testName in zip(nameList, testDict.keys()):
                self.assertEqual(controlName, testName)

            for name, bp in zip(nameList, bpList):
                numpy.testing.assert_array_equal(bp.wavelen, testDict[name].wavelen)
                numpy.testing.assert_array_equal(bp.sb, testDict[name].sb)

            for bpControl, bpTest in zip(bpList, testDict.values()):
                numpy.testing.assert_array_equal(bpControl.wavelen, bpTest.wavelen)
                numpy.testing.assert_array_equal(bpControl.sb, bpTest.sb)


    def testWavelenMatch(self):
        """
        Test that when you load bandpasses sampled over different
        wavelength grids, they all get sampled to the same wavelength
        grid.
        """
        dwavList = numpy.arange(5.0,25.0,5.0)
        bpList = []
        bpNameList = []
        for ix, dwav in enumerate(dwavList):
            name = 'bp_%d' % ix
            wavelen = numpy.arange(10.0, 1500.0, dwav)
            sb = numpy.exp(-0.5*(numpy.power((wavelen-100.0*ix)/100.0,2)))
            bp = Bandpass(wavelen=wavelen, sb=sb)
            bpList.append(bp)
            bpNameList.append(name)

        # First make sure that we have created distinct wavelength grids
        for ix in range(len(bpList)):
            for iy in range(ix+1,len(bpList)):
                self.assertTrue(len(bpList[ix].wavelen)!=len(bpList[iy].wavelen))

        testDict = CatSimBandpassDict(bpList, bpNameList)

        # Now make sure that the wavelength grids in the dict were resampled, but that
        # the original wavelength grids were not changed
        for ix in range(len(bpList)):
            numpy.testing.assert_array_equal(testDict.values()[ix].wavelen, testDict.wavelenMatch)
            if ix!=0:
                self.assertTrue(len(testDict.wavelenMatch)!=len(bpList[ix].wavelen))


    def testPhiArray(self):
        """
        Test that the phi array is correctly calculated by CatSimBandpassDict
        upon construction.
        """

        for nBp in range(3, 10, 1):
            nameList, bpList  = self.getListOfBandpasses(nBp)
            testDict = CatSimBandpassDict(bpList, nameList)
            dummySed = Sed()
            controlPhi, controlWavelenStep = dummySed.setupPhiArray(bpList)
            numpy.testing.assert_array_equal(controlPhi, testDict.phiArray)
            self.assertAlmostEqual(controlWavelenStep, testDict.wavelenStep, 10)


    def testExceptions(self):
        """
        Test that the correct exceptions are thrown by CatSimBandpassDict
        """

        nameList, bpList = self.getListOfBandpasses(4)
        dummyNameList = copy.deepcopy(nameList)
        dummyNameList[1] = dummyNameList[0]

        with self.assertRaises(RuntimeError) as context:
            testDict = CatSimBandpassDict(bpList, dummyNameList)

        self.assertTrue('occurs twice' in context.exception.message)


        testDict = CatSimBandpassDict(bpList, nameList)

        with self.assertRaises(RuntimeError) as context:
            testDict.nBandpasses = 6
        self.assertTrue('You should not be setting nBandpasses' in context.exception.message)

        with self.assertRaises(RuntimeError) as context:
            testDict.phiArray = None
        self.assertTrue('You should not be setting phiArray' in context.exception.message)

        with self.assertRaises(RuntimeError) as context:
            testDict.wavelenStep = 0.9
        self.assertTrue('You should not be setting wavelenStep' in context.exception.message)

        with self.assertRaises(RuntimeError) as context:
            testDict.wavelenMatch = numpy.arange(10.0,100.0,1.0)
        self.assertTrue('You should not be setting wavelenMatch' in context.exception.message)


    def testCalcMagListFromSed(self):
        """
        Test that calcMagListFromSed calculates the correct magnitude
        """

        wavelen = numpy.arange(10.0,2000.0,1.0)
        flux = (wavelen*2.0-5.0)*1.0e-6
        spectrum = Sed(wavelen=wavelen, flambda=flux)

        for nBp in range(3, 10, 1):

            nameList, bpList = self.getListOfBandpasses(nBp)
            testDict = CatSimBandpassDict(bpList, nameList)
            self.assertFalse(len(testDict.values()[0].wavelen)==len(spectrum.wavelen))

            magList = testDict.calcMagListFromSed(spectrum)
            for ix, (name, bp, magTest) in enumerate(zip(nameList, bpList, magList)):
                magControl = spectrum.calcMag(bp)
                self.assertAlmostEqual(magTest, magControl, 5)


    def testCalcMagListFromSedList(self):
        """
        Test that calcMagListFromSedList calculates the correct magnitude
        """

        nBandpasses = 7
        bpNameList, bpList = self.getListOfBandpasses(nBandpasses)
        testBpDict = CatSimBandpassDict(bpList, bpNameList)

        nSed = 20
        sedNameList = self.getListOfSedNames(nSed)
        magNormList = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList = numpy.random.random_sample(nSed)*5.0
        galacticAvList = numpy.random.random_sample(nSed)*0.3 + 0.1

        # first, test on an SedList without a wavelenMatch
        testSedList = CatSimSedList(sedNameList, magNormList,
                                    internalAvList=internalAvList,
                                    redshiftList=redshiftList,
                                    galacticAvList=galacticAvList)

        magList = testBpDict.calcMagListFromSedList(testSedList)
        self.assertEqual(magList.shape[0], nSed)
        self.assertEqual(magList.shape[1], nBandpasses)

        imsimBand = Bandpass()
        imsimBand.imsimBandpass()

        for ix, sedObj in enumerate(testSedList):
            dummySed = Sed(wavelen=copy.deepcopy(sedObj.wavelen),
                           flambda=copy.deepcopy(sedObj.flambda))

            for iy, bp in enumerate(testBpDict):
                mag = dummySed.calcMag(testBpDict[bp])
                self.assertAlmostEqual(mag, magList[ix][iy], 3)

        # now use wavelenMatch
        testSedList = CatSimSedList(sedNameList, magNormList,
                                    internalAvList=internalAvList,
                                    redshiftList=redshiftList,
                                    galacticAvList=galacticAvList,
                                    wavelenMatch=testBpDict.wavelenMatch)

        magList = testBpDict.calcMagListFromSedList(testSedList)
        self.assertEqual(magList.shape[0], nSed)
        self.assertEqual(magList.shape[1], nBandpasses)

        imsimBand = Bandpass()
        imsimBand.imsimBandpass()

        for ix, sedObj in enumerate(testSedList):
            dummySed = Sed(wavelen=copy.deepcopy(sedObj.wavelen),
                           flambda=copy.deepcopy(sedObj.flambda))

            for iy, bp in enumerate(testBpDict):
                mag = dummySed.calcMag(testBpDict[bp])
                self.assertAlmostEqual(mag, magList[ix][iy], 3)



def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(BandpassDictTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
