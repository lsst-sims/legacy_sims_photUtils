from __future__ import with_statement
import unittest
import os
import copy
import numpy
import lsst.utils.tests as utilsTests
from lsst.utils import getPackageDir

from lsst.sims.photUtils import Bandpass, Sed, SedList

class SedListTest(unittest.TestCase):

    def setUp(self):
        numpy.random.seed(18233)
        self.sedDir = os.path.join(getPackageDir('sims_sed_library'))
        self.sedDir = os.path.join(self.sedDir, 'galaxySED')
        self.sedPossibilities = os.listdir(self.sedDir)

    def getListOfSedNames(self, nNames):
        return [self.sedPossibilities[ii].replace('.gz','') \
                for ii in \
                numpy.random.random_integers(0, len(self.sedPossibilities)-1, nNames)]


    def testExceptions(self):
        nSed = 10
        sedNameList = self.getListOfSedNames(nSed)
        magNormList = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList = numpy.random.random_sample(nSed)*5.0
        galacticAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        wavelen_match = numpy.arange(300.0, 1500.0, 10.0)
        testList = SedList(sedNameList, magNormList, internalAvList=internalAvList,
                                 redshiftList=redshiftList, galacticAvList=galacticAvList,
                                 wavelenMatch=wavelen_match)


        with self.assertRaises(AttributeError) as context:
            testList.wavelenMatch = numpy.arange(10.0, 1000.0, 1000.0)

        with self.assertRaises(AttributeError) as context:
            testList.cosmologicalDimming = False

        with self.assertRaises(AttributeError) as context:
            testList.redshiftList = [1.8]

        with self.assertRaises(AttributeError) as context:
            testList.internalAvList = [2.5]

        with self.assertRaises(AttributeError) as context:
            testList.galacticAvList = [1.9]

        testList = SedList(sedNameList, magNormList)

        with self.assertRaises(RuntimeError) as context:
            testList.loadSedsFromList(sedNameList, magNormList, internalAvList=internalAvList)

        self.assertTrue('does not contain internalAvList' in context.exception.message)

        with self.assertRaises(RuntimeError) as context:
            testList.loadSedsFromList(sedNameList, magNormList, galacticAvList=galacticAvList)

        self.assertTrue('does not contain galacticAvList' in context.exception.message)

        with self.assertRaises(RuntimeError) as context:
            testList.loadSedsFromList(sedNameList, magNormList, redshiftList=redshiftList)

        self.assertTrue('does not contain redshiftList' in context.exception.message)



    def testSetUp(self):

        ############## Try just reading in an normalizing some SEDs
        nSed = 10
        sedNameList = self.getListOfSedNames(nSed)
        magNormList = numpy.random.random_sample(nSed)*5.0 + 15.0
        testList = SedList(sedNameList, magNormList)
        self.assertEqual(len(testList), nSed)
        self.assertTrue(testList.internalAvList is None)
        self.assertTrue(testList.galacticAvList is None)
        self.assertTrue(testList.redshiftList is None)
        self.assertTrue(testList.wavelenMatch is None)
        self.assertTrue(testList.cosmologicalDimming is True)

        imsimBand = Bandpass()
        imsimBand.imsimBandpass()

        for name, norm, sedTest in zip(sedNameList, magNormList, testList):
            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))
            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)

        ################# now add an internalAv
        sedNameList = self.getListOfSedNames(nSed)
        magNormList = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        testList = SedList(sedNameList, magNormList, internalAvList=internalAvList)
        self.assertTrue(testList.galacticAvList is None)
        self.assertTrue(testList.redshiftList is None)
        self.assertTrue(testList.wavelenMatch is None)
        self.assertTrue(testList.cosmologicalDimming is True)
        for avControl, avTest in zip(internalAvList, testList.internalAvList):
            self.assertAlmostEqual(avControl, avTest, 10)

        for name, norm, av, sedTest in zip(sedNameList, magNormList, internalAvList, testList):
            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))
            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=av)

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


        ################ now add redshift
        sedNameList = self.getListOfSedNames(nSed)
        magNormList = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList = numpy.random.random_sample(nSed)*5.0
        testList = SedList(sedNameList, magNormList, internalAvList=internalAvList,
                                 redshiftList=redshiftList)
        self.assertTrue(testList.galacticAvList is None)
        self.assertTrue(testList.wavelenMatch is None)
        self.assertTrue(testList.cosmologicalDimming is True)
        for avControl, avTest in zip(internalAvList, testList.internalAvList):
            self.assertAlmostEqual(avControl, avTest, 10)

        for zControl, zTest in zip(redshiftList, testList.redshiftList):
            self.assertAlmostEqual(zControl, zTest, 10)

        for name, norm, av, zz, sedTest in \
        zip(sedNameList, magNormList, internalAvList, redshiftList, testList):
            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))
            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=av)

            sedControl.redshiftSED(zz, dimming=True)

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


        ################# without cosmological dimming
        sedNameList = self.getListOfSedNames(nSed)
        magNormList = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList = numpy.random.random_sample(nSed)*5.0
        testList = SedList(sedNameList, magNormList, internalAvList=internalAvList,
                                 redshiftList=redshiftList, cosmologicalDimming=False)
        self.assertTrue(testList.galacticAvList is None)
        self.assertTrue(testList.wavelenMatch is None)
        self.assertTrue(testList.cosmologicalDimming is False)
        for avControl, avTest in zip(internalAvList, testList.internalAvList):
            self.assertAlmostEqual(avControl, avTest, 10)

        for zControl, zTest in zip(redshiftList, testList.redshiftList):
            self.assertAlmostEqual(zControl, zTest, 10)

        for name, norm, av, zz, sedTest in \
        zip(sedNameList, magNormList, internalAvList, redshiftList, testList):
            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))
            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=av)

            sedControl.redshiftSED(zz, dimming=False)

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


        ################ now add galacticAv
        sedNameList = self.getListOfSedNames(nSed)
        magNormList = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList = numpy.random.random_sample(nSed)*5.0
        galacticAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        testList = SedList(sedNameList, magNormList, internalAvList=internalAvList,
                                 redshiftList=redshiftList, galacticAvList=galacticAvList)
        self.assertTrue(testList.wavelenMatch is None)
        self.assertTrue(testList.cosmologicalDimming is True)
        for avControl, avTest in zip(internalAvList, testList.internalAvList):
            self.assertAlmostEqual(avControl, avTest, 10)

        for zControl, zTest in zip(redshiftList, testList.redshiftList):
            self.assertAlmostEqual(zControl, zTest, 10)

        for avControl, avTest in zip(galacticAvList, testList.galacticAvList):
            self.assertAlmostEqual(avControl, avTest, 10)

        for name, norm, av, zz, gav, sedTest in \
        zip(sedNameList, magNormList, internalAvList, \
            redshiftList, galacticAvList, testList):

            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))
            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=av)

            sedControl.redshiftSED(zz, dimming=True)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


        ################ now use a wavelen_match
        sedNameList = self.getListOfSedNames(nSed)
        magNormList = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList = numpy.random.random_sample(nSed)*5.0
        galacticAvList = numpy.random.random_sample(nSed)*0.3 + 0.1
        wavelen_match = numpy.arange(300.0, 1500.0, 10.0)
        testList = SedList(sedNameList, magNormList, internalAvList=internalAvList,
                                 redshiftList=redshiftList, galacticAvList=galacticAvList,
                                 wavelenMatch=wavelen_match)

        self.assertTrue(testList.cosmologicalDimming is True)
        for avControl, avTest in zip(internalAvList, testList.internalAvList):
            self.assertAlmostEqual(avControl, avTest, 10)

        for zControl, zTest in zip(redshiftList, testList.redshiftList):
            self.assertAlmostEqual(zControl, zTest, 10)

        for avControl, avTest in zip(galacticAvList, testList.galacticAvList):
            self.assertAlmostEqual(avControl, avTest, 10)

        numpy.testing.assert_array_equal(wavelen_match, testList.wavelenMatch)

        for name, norm, av, zz, gav, sedTest in \
        zip(sedNameList, magNormList, internalAvList, \
            redshiftList, galacticAvList, testList):

            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))

            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=av)

            sedControl.redshiftSED(zz, dimming=True)
            sedControl.resampleSED(wavelen_match=wavelen_match)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


    def testAddingToList(self):
        """
        Test that we can add Seds to an already instantiated SedList
        """
        imsimBand = Bandpass()
        imsimBand.imsimBandpass()
        nSed = 10
        sedNameList_0 = self.getListOfSedNames(nSed)
        magNormList_0 = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList_0 = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList_0 = numpy.random.random_sample(nSed)*5.0
        galacticAvList_0 = numpy.random.random_sample(nSed)*0.3 + 0.1
        wavelen_match = numpy.arange(300.0, 1500.0, 10.0)
        testList = SedList(sedNameList_0, magNormList_0, internalAvList=internalAvList_0, \
                                 redshiftList=redshiftList_0, galacticAvList=galacticAvList_0,
                                 wavelenMatch=wavelen_match)


        # experiment with adding different combinations of physical parameter lists
        # as None and not None
        for addIav in [True, False]:
            for addRedshift in [True, False]:
                for addGav in [True, False]:

                    testList = SedList(sedNameList_0, magNormList_0, internalAvList=internalAvList_0, \
                                             redshiftList=redshiftList_0, galacticAvList=galacticAvList_0,
                                             wavelenMatch=wavelen_match)

                    sedNameList_1 = self.getListOfSedNames(nSed)
                    magNormList_1 = numpy.random.random_sample(nSed)*5.0 + 15.0

                    if addIav:
                        internalAvList_1 = numpy.random.random_sample(nSed)*0.3 + 0.1
                    else:
                        internalAvList_1 = None

                    if addRedshift:
                        redshiftList_1 = numpy.random.random_sample(nSed)*5.0
                    else:
                        redshiftList_1 = None

                    if addGav:
                        galacticAvList_1 = numpy.random.random_sample(nSed)*0.3 + 0.1
                    else:
                        galacticAvList_1 = None


                    testList.loadSedsFromList(sedNameList_1, magNormList_1,
                                              internalAvList=internalAvList_1,
                                              galacticAvList=galacticAvList_1,
                                              redshiftList=redshiftList_1)

                    self.assertEqual(len(testList), 2*nSed)
                    numpy.testing.assert_array_equal(wavelen_match, testList.wavelenMatch)

                    for ix in range(len(sedNameList_0)):
                        self.assertAlmostEqual(internalAvList_0[ix], testList.internalAvList[ix], 10)
                        self.assertAlmostEqual(galacticAvList_0[ix], testList.galacticAvList[ix], 10)
                        self.assertAlmostEqual(redshiftList_0[ix], testList.redshiftList[ix], 10)


                    for ix in range(len(sedNameList_1)):
                        if addIav:
                            self.assertAlmostEqual(internalAvList_1[ix], testList.internalAvList[ix+nSed], 10)
                        else:
                            self.assertTrue(testList.internalAvList[ix+nSed] is None)

                        if addGav:
                            self.assertAlmostEqual(galacticAvList_1[ix], testList.galacticAvList[ix+nSed], 10)
                        else:
                            self.assertTrue(testList.galacticAvList[ix+nSed] is None)

                        if addRedshift:
                            self.assertAlmostEqual(redshiftList_1[ix], testList.redshiftList[ix+nSed], 10)
                        else:
                            self.assertTrue(testList.redshiftList[ix+nSed] is None)

                    for ix, (name, norm, iav, gav, zz) in \
                      enumerate(zip(sedNameList_0, magNormList_0, internalAvList_0, \
                                  galacticAvList_0, redshiftList_0)):

                        sedControl = Sed()
                        sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))

                        fnorm = sedControl.calcFluxNorm(norm, imsimBand)
                        sedControl.multiplyFluxNorm(fnorm)

                        a_coeff, b_coeff = sedControl.setupCCMab()
                        sedControl.addCCMDust(a_coeff, b_coeff, A_v=iav)

                        sedControl.redshiftSED(zz, dimming=True)
                        sedControl.resampleSED(wavelen_match=wavelen_match)

                        a_coeff, b_coeff = sedControl.setupCCMab()
                        sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

                        sedTest = testList[ix]

                        numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
                        numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
                        numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


                    if not addIav:
                        internalAvList_1 = [None] * nSed

                    if not addRedshift:
                        redshiftList_1 = [None] * nSed

                    if not addGav:
                        galacticAvList_1 = [None] * nSed

                    for ix, (name, norm, iav, gav, zz) in \
                        enumerate(zip(sedNameList_1, magNormList_1, internalAvList_1, \
                                      galacticAvList_1, redshiftList_1)):

                        sedControl = Sed()
                        sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))

                        fnorm = sedControl.calcFluxNorm(norm, imsimBand)
                        sedControl.multiplyFluxNorm(fnorm)

                        if addIav:
                            a_coeff, b_coeff = sedControl.setupCCMab()
                            sedControl.addCCMDust(a_coeff, b_coeff, A_v=iav)

                        if addRedshift:
                            sedControl.redshiftSED(zz, dimming=True)

                        sedControl.resampleSED(wavelen_match=wavelen_match)

                        if addGav:
                            a_coeff, b_coeff = sedControl.setupCCMab()
                            sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

                        sedTest = testList[ix+nSed]

                        numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
                        numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
                        numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


    def testAddingNonesToList(self):
        """
        Test what happens if you add SEDs to an SedList that have None for
        one or more of the physical parameters (i.e. galacticAv, internalAv, or redshift)
        """
        imsimBand = Bandpass()
        imsimBand.imsimBandpass()
        nSed = 10
        sedNameList_0 = self.getListOfSedNames(nSed)
        magNormList_0 = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList_0 = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList_0 = numpy.random.random_sample(nSed)*5.0
        galacticAvList_0 = numpy.random.random_sample(nSed)*0.3 + 0.1
        wavelen_match = numpy.arange(300.0, 1500.0, 10.0)
        testList = SedList(sedNameList_0, magNormList_0, internalAvList=internalAvList_0, \
                                 redshiftList=redshiftList_0, galacticAvList=galacticAvList_0,
                                 wavelenMatch=wavelen_match)


        sedNameList_1 = self.getListOfSedNames(nSed)
        magNormList_1 = list(numpy.random.random_sample(nSed)*5.0 + 15.0)
        internalAvList_1 = list(numpy.random.random_sample(nSed)*0.3 + 0.1)
        redshiftList_1 = list(numpy.random.random_sample(nSed)*5.0)
        galacticAvList_1 = list(numpy.random.random_sample(nSed)*0.3 + 0.1)

        internalAvList_1[0] = None
        redshiftList_1[1] = None
        galacticAvList_1[2] = None

        internalAvList_1[3] = None
        redshiftList_1[3] = None

        internalAvList_1[4] = None
        galacticAvList_1[4] = None

        redshiftList_1[5] = None
        galacticAvList_1[5] = None

        internalAvList_1[6] = None
        redshiftList_1[6] = None
        galacticAvList_1[6] = None

        testList.loadSedsFromList(sedNameList_1, magNormList_1,
                                  internalAvList=internalAvList_1,
                                  galacticAvList=galacticAvList_1,
                                  redshiftList=redshiftList_1)

        self.assertEqual(len(testList), 2*nSed)
        numpy.testing.assert_array_equal(wavelen_match, testList.wavelenMatch)

        for ix in range(len(sedNameList_0)):
            self.assertAlmostEqual(internalAvList_0[ix], testList.internalAvList[ix], 10)
            self.assertAlmostEqual(galacticAvList_0[ix], testList.galacticAvList[ix], 10)
            self.assertAlmostEqual(redshiftList_0[ix], testList.redshiftList[ix], 10)


        for ix in range(len(sedNameList_1)):
            self.assertAlmostEqual(internalAvList_1[ix], testList.internalAvList[ix+nSed], 10)
            self.assertAlmostEqual(galacticAvList_1[ix], testList.galacticAvList[ix+nSed], 10)
            self.assertAlmostEqual(redshiftList_1[ix], testList.redshiftList[ix+nSed], 10)

        for ix, (name, norm, iav, gav, zz) in \
        enumerate(zip(sedNameList_0, magNormList_0, internalAvList_0, \
                      galacticAvList_0, redshiftList_0)):

            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))

            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=iav)

            sedControl.redshiftSED(zz, dimming=True)
            sedControl.resampleSED(wavelen_match=wavelen_match)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

            sedTest = testList[ix]

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


        for ix, (name, norm, iav, gav, zz) in \
        enumerate(zip(sedNameList_1, magNormList_1, internalAvList_1, \
                      galacticAvList_1, redshiftList_1)):

            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))

            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            if iav is not None:
                a_coeff, b_coeff = sedControl.setupCCMab()
                sedControl.addCCMDust(a_coeff, b_coeff, A_v=iav)

            if zz is not None:
                sedControl.redshiftSED(zz, dimming=True)

            sedControl.resampleSED(wavelen_match=wavelen_match)

            if gav is not None:
                a_coeff, b_coeff = sedControl.setupCCMab()
                sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

            sedTest = testList[ix+nSed]

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


    def testAlternateNormalizingBandpass(self):
        """
        A reiteration of testAddingToList, but testing with a non-imsimBandpass
        normalizing bandpass
        """
        normalizingBand = Bandpass()
        normalizingBand.readThroughput(os.path.join(getPackageDir('throughputs'),'baseline','total_r.dat'))
        nSed = 10
        sedNameList_0 = self.getListOfSedNames(nSed)
        magNormList_0 = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList_0 = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList_0 = numpy.random.random_sample(nSed)*5.0
        galacticAvList_0 = numpy.random.random_sample(nSed)*0.3 + 0.1
        wavelen_match = numpy.arange(300.0, 1500.0, 10.0)
        testList = SedList(sedNameList_0, magNormList_0,
                           normalizingBandpass=normalizingBand,
                           internalAvList=internalAvList_0,
                           redshiftList=redshiftList_0, galacticAvList=galacticAvList_0,
                           wavelenMatch=wavelen_match)


        sedNameList_1 = self.getListOfSedNames(nSed)
        magNormList_1 = numpy.random.random_sample(nSed)*5.0 + 15.0

        internalAvList_1 = numpy.random.random_sample(nSed)*0.3 + 0.1

        redshiftList_1 = numpy.random.random_sample(nSed)*5.0

        galacticAvList_1 = numpy.random.random_sample(nSed)*0.3 + 0.1


        testList.loadSedsFromList(sedNameList_1, magNormList_1,
                                  internalAvList=internalAvList_1,
                                  galacticAvList=galacticAvList_1,
                                  redshiftList=redshiftList_1)

        self.assertEqual(len(testList), 2*nSed)
        numpy.testing.assert_array_equal(wavelen_match, testList.wavelenMatch)

        for ix in range(len(sedNameList_0)):
            self.assertAlmostEqual(internalAvList_0[ix], testList.internalAvList[ix], 10)
            self.assertAlmostEqual(galacticAvList_0[ix], testList.galacticAvList[ix], 10)
            self.assertAlmostEqual(redshiftList_0[ix], testList.redshiftList[ix], 10)


        for ix in range(len(sedNameList_1)):
            self.assertAlmostEqual(internalAvList_1[ix], testList.internalAvList[ix+nSed], 10)
            self.assertAlmostEqual(galacticAvList_1[ix], testList.galacticAvList[ix+nSed], 10)
            self.assertAlmostEqual(redshiftList_1[ix], testList.redshiftList[ix+nSed], 10)

        for ix, (name, norm, iav, gav, zz) in \
          enumerate(zip(sedNameList_0, magNormList_0, internalAvList_0, \
                     galacticAvList_0, redshiftList_0)):

            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))

            fnorm = sedControl.calcFluxNorm(norm, normalizingBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=iav)

            sedControl.redshiftSED(zz, dimming=True)
            sedControl.resampleSED(wavelen_match=wavelen_match)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

            sedTest = testList[ix]

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)

        for ix, (name, norm, iav, gav, zz) in \
            enumerate(zip(sedNameList_1, magNormList_1, internalAvList_1, \
                          galacticAvList_1, redshiftList_1)):

            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))

            fnorm = sedControl.calcFluxNorm(norm, normalizingBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=iav)

            sedControl.redshiftSED(zz, dimming=True)

            sedControl.resampleSED(wavelen_match=wavelen_match)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

            sedTest = testList[ix+nSed]

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)


    def testFlush(self):
        imsimBand = Bandpass()
        imsimBand.imsimBandpass()
        nSed = 10
        sedNameList_0 = self.getListOfSedNames(nSed)
        magNormList_0 = numpy.random.random_sample(nSed)*5.0 + 15.0
        internalAvList_0 = numpy.random.random_sample(nSed)*0.3 + 0.1
        redshiftList_0 = numpy.random.random_sample(nSed)*5.0
        galacticAvList_0 = numpy.random.random_sample(nSed)*0.3 + 0.1
        wavelen_match = numpy.arange(300.0, 1500.0, 10.0)
        testList = SedList(sedNameList_0, magNormList_0, internalAvList=internalAvList_0, \
                                 redshiftList=redshiftList_0, galacticAvList=galacticAvList_0,
                                 wavelenMatch=wavelen_match)

        self.assertEqual(len(testList), nSed)
        numpy.testing.assert_array_equal(wavelen_match, testList.wavelenMatch)

        for ix in range(len(sedNameList_0)):
            self.assertAlmostEqual(internalAvList_0[ix], testList.internalAvList[ix], 10)
            self.assertAlmostEqual(galacticAvList_0[ix], testList.galacticAvList[ix], 10)
            self.assertAlmostEqual(redshiftList_0[ix], testList.redshiftList[ix], 10)

        for ix, (name, norm, iav, gav, zz) in \
        enumerate(zip(sedNameList_0, magNormList_0, internalAvList_0, \
                      galacticAvList_0, redshiftList_0)):

            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))

            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=iav)

            sedControl.redshiftSED(zz, dimming=True)
            sedControl.resampleSED(wavelen_match=wavelen_match)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

            sedTest = testList[ix]

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)



        testList.flush()

        sedNameList_1 = self.getListOfSedNames(nSed/2)
        magNormList_1 = numpy.random.random_sample(nSed/2)*5.0 + 15.0
        internalAvList_1 = numpy.random.random_sample(nSed/2)*0.3 + 0.1
        redshiftList_1 = numpy.random.random_sample(nSed/2)*5.0
        galacticAvList_1 = numpy.random.random_sample(nSed/2)*0.3 + 0.1

        testList.loadSedsFromList(sedNameList_1, magNormList_1,
                                  internalAvList=internalAvList_1,
                                  galacticAvList=galacticAvList_1,
                                  redshiftList=redshiftList_1)

        self.assertEqual(len(testList), nSed/2)
        self.assertEqual(len(testList.redshiftList), nSed/2)
        self.assertEqual(len(testList.internalAvList), nSed/2)
        self.assertEqual(len(testList.galacticAvList), nSed/2)
        numpy.testing.assert_array_equal(wavelen_match, testList.wavelenMatch)


        for ix in range(len(sedNameList_1)):
            self.assertAlmostEqual(internalAvList_1[ix], testList.internalAvList[ix], 10)
            self.assertAlmostEqual(galacticAvList_1[ix], testList.galacticAvList[ix], 10)
            self.assertAlmostEqual(redshiftList_1[ix], testList.redshiftList[ix], 10)

        for ix, (name, norm, iav, gav, zz) in \
        enumerate(zip(sedNameList_1, magNormList_1, internalAvList_1, \
                      galacticAvList_1, redshiftList_1)):

            sedControl = Sed()
            sedControl.readSED_flambda(os.path.join(self.sedDir, name+'.gz'))

            fnorm = sedControl.calcFluxNorm(norm, imsimBand)
            sedControl.multiplyFluxNorm(fnorm)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=iav)

            sedControl.redshiftSED(zz, dimming=True)
            sedControl.resampleSED(wavelen_match=wavelen_match)

            a_coeff, b_coeff = sedControl.setupCCMab()
            sedControl.addCCMDust(a_coeff, b_coeff, A_v=gav)

            sedTest = testList[ix]

            numpy.testing.assert_array_equal(sedControl.wavelen, sedTest.wavelen)
            numpy.testing.assert_array_equal(sedControl.flambda, sedTest.flambda)
            numpy.testing.assert_array_equal(sedControl.fnu, sedTest.fnu)




def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(SedListTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
