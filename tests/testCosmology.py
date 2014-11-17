import os
import unittest
import lsst.utils.tests as utilsTests
import numpy
import scipy

from lsst.sims.photUtils import CosmologyObject, CosmologyWrapper, PhotometryGalaxies
from lsst.sims.catalogs.measures.instance import InstanceCatalog

from lsst.sims.catalogs.generation.utils import myTestGals, makeGalTestDB
from lsst.sims.photUtils.utils import testGalaxies, comovingDistanceIntegrand, \
                                      cosmologicalOmega

class cosmologicalGalaxyCatalog(testGalaxies, CosmologyWrapper):
    column_outputs = ['uRecalc', 'gRecalc', 'rRecalc', 'iRecalc', 'zRecalc', 'yRecalc',
                      'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
                      'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
                      'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn',
                      'redshift', 'distanceModulus']


class CosmologyUnitTest(unittest.TestCase):

    def setUp(self):
        self.speedOfLight = 2.9979e5 #in km/sec

    def tearDown(self):
        del self.speedOfLight

    def testFlatLCDM(self):
        """
        Test the evolution of H and Omega_i as a function of redshift for
        flat Lambda CDM models
        """
        H0 = 50.0
        for Om0 in numpy.arange(start=0.1, stop=0.91, step=0.4):
            universe = CosmologyObject(H0=H0, Om0=Om0)

            Og0 = universe.OmegaPhotons(redshift=0.0)
            Onu0 = universe.OmegaNeutrinos(redshift=0.0)

            self.assertTrue(Og0 < 1.0e-4)
            self.assertTrue(Onu0 < 1.0e-4)
            self.assertAlmostEqual(universe.OmegaMatter(redshift=0.0), Om0, 10)
            self.assertAlmostEqual(1.0 - Om0 - universe.OmegaDarkEnergy(redshift=0.0), Og0+Onu0, 6)
            self.assertAlmostEqual(universe.H(redshift=0.0),H0,10)
            self.assertEqual(universe.OmegaCurvature(),0.0)

            Om0 = universe.OmegaMatter(redshift=0.0)
            Ode0 = universe.OmegaDarkEnergy(redshift=0.0)

            for zz in numpy.arange(start=0.0, stop=4.1, step=2.0):
               aa = (1.0+zz)

               Hcontrol, OmControl, OdeControl, OgControl, OnuControl, \
                   OkControl, = cosmologicalOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0)

               self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
               self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
               self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
               self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
               self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

            del universe

    def testFlatW0Wa(self):
        """
        Test the evolution of H and Omega_i as a function of redshift for
        flat models with w = w0 + wa * z / (1 + z)
        """

        H0 = 96.0
        for Om0 in numpy.arange(start=0.1, stop=0.95, step=0.4):
            for w0 in numpy.arange(start=-1.1, stop=-0.89, step=0.2):
                for wa in numpy.arange(start=-0.1, stop=0.11, step=0.2):

                    universe = CosmologyObject(H0=H0, Om0=Om0, w0=w0, wa=wa)

                    Og0 = universe.OmegaPhotons(redshift=0.0)
                    Onu0 = universe.OmegaNeutrinos(redshift=0.0)

                    self.assertTrue(Og0 < 1.0e-4)
                    self.assertTrue(Onu0 < 1.0e-4)
                    self.assertAlmostEqual(universe.OmegaMatter(redshift=0.0), Om0, 10)
                    self.assertAlmostEqual(1.0 - Om0 - universe.OmegaDarkEnergy(redshift=0.0), Og0+Onu0, 6)
                    self.assertAlmostEqual(universe.H(redshift=0.0),H0,10)
                    self.assertEqual(universe.OmegaCurvature(),0.0)

                    Om0 = universe.OmegaMatter(redshift=0.0)
                    Ode0 = universe.OmegaDarkEnergy(redshift=0.0)

                    for zz in numpy.arange(start=0.0, stop=4.1, step=2.0):

                       wControl = w0 + wa*(1.0 - 1.0/(1.0+zz))
                       self.assertAlmostEqual(wControl, universe.w(redshift=zz), 6)

                       Hcontrol, OmControl, OdeControl, OgControl, OnuControl, \
                           OkControl = cosmologicalOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0,
                                                    w0=w0, wa=wa)

                       self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
                       self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
                       self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
                       self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
                       self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

                    del universe

    def testFlatW0(self):
        """
        Test the evolution of H and Omega_i as a function of redshift for flat
        models with constant w
        """

        H0 = 96.0
        for Om0 in numpy.arange(start=0.1, stop=0.95, step=0.4):
            for w0 in numpy.arange(start=-1.5, stop=-0.49, step=1.0):

                universe = CosmologyObject(H0=H0, Om0=Om0, w0=w0)

                Og0 = universe.OmegaPhotons(redshift=0.0)
                Onu0 = universe.OmegaNeutrinos(redshift=0.0)

                self.assertTrue(Og0 < 1.0e-4)
                self.assertTrue(Onu0 < 1.0e-4)
                self.assertAlmostEqual(universe.OmegaMatter(redshift=0.0), Om0, 10)
                self.assertAlmostEqual(1.0 - Om0 - universe.OmegaDarkEnergy(redshift=0.0), Og0+Onu0, 6)
                self.assertAlmostEqual(universe.H(redshift=0.0),H0,10)
                self.assertEqual(universe.OmegaCurvature(),0.0)

                Om0 = universe.OmegaMatter(redshift=0.0)
                Ode0 = universe.OmegaDarkEnergy(redshift=0.0)

                for zz in numpy.arange(start=0.0, stop=4.1, step=2.0):

                   self.assertAlmostEqual(w0, universe.w(redshift=zz), 6)

                   Hcontrol, OmControl, OdeControl, OgControl, OnuControl, \
                       OkControl = cosmologicalOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0,
                                                w0=w0, wa=0.0)

                   self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
                   self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
                   self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
                   self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
                   self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

                del universe

    def testNonFlatLCDM(self):
        """
        Test the evolution of H and Omega_i as a function of redshift for non-flat
        Lambda CDM models
        """
        w0 = -1.0
        wa = 0.0
        H0 = 77.0

        for Om0 in numpy.arange(start=0.15, stop=0.96, step=0.4):
            for Ode0 in numpy.arange(start=1.0-Om0-0.1, stop=1.0-Om0+0.11, step=0.2):

                universe = CosmologyObject(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                Og0 = universe.OmegaPhotons(redshift=0.0)
                Onu0 = universe.OmegaNeutrinos(redshift=0.0)

                self.assertTrue(Og0 < 1.0e-4)
                self.assertTrue(Onu0 < 1.0e-4)
                self.assertAlmostEqual(universe.OmegaMatter(redshift=0.0), Om0, 10)
                self.assertAlmostEqual(universe.OmegaDarkEnergy(redshift=0.0), Ode0, 10)
                self.assertAlmostEqual(1.0 - Ode0 - Om0 - universe.OmegaCurvature(redshift=0.0), Og0+Onu0, 6)
                self.assertAlmostEqual(universe.H(redshift=0.0),H0,10)

                Om0 = universe.OmegaMatter(redshift=0.0)
                Ode0 = universe.OmegaDarkEnergy(redshift=0.0)
                Ok0 = universe.OmegaCurvature(redshift=0.0)

                for zz in numpy.arange(start=0.0, stop=4.0, step=2.0):
                    Hcontrol, OmControl, OdeControl, OgControl, OnuControl, \
                        OkControl = cosmologicalOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0,
                                                 Ode0=Ode0)

                    self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
                    self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
                    self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
                    self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
                    self.assertAlmostEqual(OkControl, universe.OmegaCurvature(redshift=zz), 6)
                    self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

                del universe

    def testNonFlatW0Wa(self):
        """
        Test the evolution of H and Omega_i as a function of redshift for non-flat
        models with w = w0 + wa * z / (1+z)
        """

        H0 = 60.0

        for Om0 in numpy.arange(start=0.15, stop=0.76, step=0.3):
            for Ode0 in numpy.arange(1.0-Om0-0.1, stop = 1.0-Om0+0.11, step=0.2):
                for w0 in numpy.arange(start=-1.1, stop = -0.89, step=0.2):
                    for wa in numpy.arange(start=-0.1, stop=0.15, step=0.2):

                        universe = CosmologyObject(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                        Og0 = universe.OmegaPhotons(redshift=0.0)
                        Onu0 = universe.OmegaNeutrinos(redshift=0.0)

                        self.assertTrue(Og0 < 1.0e-4)
                        self.assertTrue(Onu0 < 1.0e-4)
                        self.assertAlmostEqual(universe.OmegaMatter(redshift=0.0), Om0, 10)
                        self.assertAlmostEqual(Ode0, universe.OmegaDarkEnergy(redshift=0.0), 10)
                        self.assertAlmostEqual(1.0 - Om0 -Ode0 - universe.OmegaCurvature(redshift=0.0),
                                                   Og0+Onu0, 10)
                        self.assertAlmostEqual(universe.H(redshift=0.0),H0,10)

                        Om0 = universe.OmegaMatter(redshift=0.0)
                        Ode0 = universe.OmegaDarkEnergy(redshift=0.0)

                        for zz in numpy.arange(start=0.0, stop=4.0, step=2.0):

                           wControl = w0 + wa*(1.0 - 1.0/(1.0+zz))
                           self.assertAlmostEqual(wControl, universe.w(redshift=zz), 6)

                           Hcontrol, OmControl, OdeControl, OgControl, OnuControl, \
                               OkControl = cosmologicalOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0,
                                                        w0=w0, wa=wa, Ode0=Ode0)

                           self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
                           self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
                           self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
                           self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
                           self.assertAlmostEqual(OkControl, universe.OmegaCurvature(redshift=zz), 6)
                           self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

                        del universe

    def testNonFlatW0(self):
        """
        Test the evolution of H and Omega_i as a function of redshift for non-flat
        models with constant w
        """

        H0 = 60.0

        for Om0 in numpy.arange(start=0.15, stop=0.76, step=0.3):
            for Ode0 in numpy.arange(1.0-Om0-0.1, stop = 1.0-Om0+0.11, step=0.2):
                for w0 in numpy.arange(start=-1.1, stop = -0.89, step=0.2):

                    universe = CosmologyObject(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0)

                    Og0 = universe.OmegaPhotons(redshift=0.0)
                    Onu0 = universe.OmegaNeutrinos(redshift=0.0)

                    self.assertTrue(Og0 < 1.0e-4)
                    self.assertTrue(Onu0 < 1.0e-4)
                    self.assertAlmostEqual(universe.OmegaMatter(redshift=0.0), Om0, 10)
                    self.assertAlmostEqual(Ode0, universe.OmegaDarkEnergy(redshift=0.0), 10)
                    self.assertAlmostEqual(1.0 - Om0 -Ode0 - universe.OmegaCurvature(redshift=0.0),
                                               Og0+Onu0, 10)
                    self.assertAlmostEqual(universe.H(redshift=0.0),H0,10)

                    Om0 = universe.OmegaMatter(redshift=0.0)
                    Ode0 = universe.OmegaDarkEnergy(redshift=0.0)

                    for zz in numpy.arange(start=0.0, stop=4.0, step=2.0):

                       self.assertAlmostEqual(w0, universe.w(redshift=zz), 6)

                       Hcontrol, OmControl, OdeControl, OgControl, OnuControl, \
                           OkControl = cosmologicalOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0,
                                                    w0=w0, wa=0.0, Ode0=Ode0)

                       self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
                       self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
                       self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
                       self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
                       self.assertAlmostEqual(OkControl, universe.OmegaCurvature(redshift=zz), 6)
                       self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

                    del universe

    def testComovingDistance(self):
        """
        Test comoving distance calculation

        Note: this is comoving distance defined as X in the FRW metric

        ds^2 = -c^2 dt^2 + a^2 dX^2 + sin^2(X) dOmega^2

        where spatial curvature is accounted for in the sin function
        """


        H0 = 73.0
        for Om0 in numpy.arange(start=0.15, stop=0.56, step=0.2):
            for Ode0 in numpy.arange(start=1.0-Om0-0.1, stop=1.0-Om0+0.11, step=0.2):
                for w0 in numpy.arange(start=-1.1, stop=-0.85, step=0.2):
                    for wa in numpy.arange(start=-0.1, stop=0.115, step=0.02):

                        universe = CosmologyObject(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)
                        Og0 = universe.OmegaPhotons()
                        Onu0 = universe.OmegaNeutrinos()

                        for zz in numpy.arange(start=0.1, stop=4.2, step=2.0):
                            comovingControl = universe.comovingDistance(redshift=zz)
                            comovingTest = \
                                self.speedOfLight*scipy.integrate.quad(comovingDistanceIntegrand, 0.0, zz,
                                                                   args=(H0, Om0, Ode0, Og0, Onu0, w0, wa))[0]
                            self.assertAlmostEqual(comovingControl/comovingTest,1.0,4)

    def testLuminosityDistance(self):
        """
        Test the calculation of the luminosity distance
        """


        H0 = 73.0

        for Om0 in numpy.arange(start=0.15, stop=0.56, step=0.2):
            for Ode0 in numpy.arange(start=1.0-Om0-0.1, stop=1.0-Om0+0.11, step=0.2):
                for w0 in numpy.arange(start=-1.1, stop=-0.85, step=0.2):
                    for wa in numpy.arange(start=-0.1, stop=0.11, step=0.2):

                        universe = CosmologyObject(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                        sqrtkCurvature = numpy.sqrt(numpy.abs(universe.OmegaCurvature()))*universe.H()/self.speedOfLight
                        Og0 = universe.OmegaPhotons()
                        Onu0 = universe.OmegaNeutrinos()

                        for zz in numpy.arange(start=0.1, stop=4.2, step=2.0):
                            luminosityControl = universe.luminosityDistance(redshift=zz)
                            comovingDistance = self.speedOfLight*scipy.integrate.quad(comovingDistanceIntegrand, 0.0, zz,
                                                                                 args=(H0, Om0, Ode0, Og0, Onu0, w0, wa))[0]

                            if universe.OmegaCurvature()<0.0:
                                nn =sqrtkCurvature*comovingDistance
                                nn = numpy.sin(nn)
                                luminosityTest = (1.0+zz)*nn/sqrtkCurvature
                            elif universe.OmegaCurvature()>0.0:
                                nn = sqrtkCurvature*comovingDistance
                                nn = numpy.sinh(nn)
                                luminosityTest = (1.0+zz)*nn/sqrtkCurvature
                            else:
                                luminosityTest = (1.0+zz)*comovingDistance
                            self.assertAlmostEqual(luminosityControl/luminosityTest,1.0,4)

    def testAngularDiameterDistance(self):
        """
        Test the calculation of the angular diameter distance
        """

        H0 = 56.0
        universe=CosmologyObject()
        for Om0 in numpy.arange(start=0.15, stop=0.56, step=0.2):
            for Ode0 in numpy.arange(start=1.0-Om0-0.1, stop=1.0-Om0+0.11, step=0.2):
                for w0 in numpy.arange(start=-1.1, stop=-0.85, step=0.2):
                    for wa in numpy.arange(start=-0.1, stop=0.11, step=0.2):

                        universe = CosmologyObject(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                        sqrtkCurvature = numpy.sqrt(numpy.abs(universe.OmegaCurvature()))*universe.H()/self.speedOfLight
                        Og0 = universe.OmegaPhotons()
                        Onu0 = universe.OmegaNeutrinos()

                        for zz in numpy.arange(start=0.1, stop=4.2, step=2.0):
                            angularControl = universe.angularDiameterDistance(redshift=zz)
                            comovingDistance = self.speedOfLight*scipy.integrate.quad(comovingDistanceIntegrand, 0.0, zz,
                                                                                 args=(H0, Om0, Ode0, Og0, Onu0, w0, wa))[0]

                            if universe.OmegaCurvature()<0.0:
                                nn =sqrtkCurvature*comovingDistance
                                nn = numpy.sin(nn)
                                angularTest = nn/sqrtkCurvature
                            elif universe.OmegaCurvature()>0.0:
                                nn = sqrtkCurvature*comovingDistance
                                nn = numpy.sinh(nn)
                                angularTest = nn/sqrtkCurvature
                            else:
                                angularTest = comovingDistance
                            angularTest /= (1.0+zz)
                            self.assertAlmostEqual(angularControl/angularTest,1.0,4)


    def testDistanceModulus(self):
        """
        Test the calculation of the distance modulus out to a certain redshift
        """
        H0 = 73.0

        universe=CosmologyObject()
        for Om0 in numpy.arange(start=0.15, stop=0.56, step=0.2):
            for Ode0 in numpy.arange(start=1.0-Om0-0.1, stop=1.0-Om0+0.11, step=0.2):
                for w0 in numpy.arange(start=-1.1, stop=-0.85, step=0.2):
                    for wa in numpy.arange(start=-0.1, stop=0.11, step=0.2):

                        universe = CosmologyObject(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                        sqrtkCurvature = numpy.sqrt(numpy.abs(universe.OmegaCurvature()))*universe.H()/self.speedOfLight
                        Og0 = universe.OmegaPhotons()
                        Onu0 = universe.OmegaNeutrinos()

                        for zz in numpy.arange(start=0.1, stop=4.2, step=2.0):
                            modulusControl = universe.distanceModulus(redshift=zz)
                            comovingDistance = self.speedOfLight*scipy.integrate.quad(comovingDistanceIntegrand, 0.0, zz,
                                                                                 args=(H0, Om0, Ode0, Og0, Onu0, w0, wa))[0]

                            if universe.OmegaCurvature()<0.0:
                                nn =sqrtkCurvature*comovingDistance
                                nn = numpy.sin(nn)
                                luminosityDistance = (1.0+zz)*nn/sqrtkCurvature
                            elif universe.OmegaCurvature()>0.0:
                                nn = sqrtkCurvature*comovingDistance
                                nn = numpy.sinh(nn)
                                luminosityDistance = (1.0+zz)*nn/sqrtkCurvature
                            else:
                                luminosityDistance = (1.0+zz)*comovingDistance

                            modulusTest = 5.0*numpy.log10(luminosityDistance) + 25.0
                            self.assertAlmostEqual(modulusControl/modulusTest,1.0,4)

    def testDistanceModulusAtZero(self):
        """
        Test to make sure that the distance modulus is set to zero if the distance modulus method
        returns a negative number
        """
        universe = CosmologyObject()
        ztest = [0.0, 1.0, 2.0, 0.0, 3.0]
        mm = universe.distanceModulus(redshift=ztest)
        self.assertEqual(mm[0], 0.0)
        self.assertEqual(mm[3], 0.0)
        self.assertEqual(mm[1], 5.0*numpy.log10(universe.luminosityDistance(ztest[1])) + 25.0)
        self.assertEqual(mm[2], 5.0*numpy.log10(universe.luminosityDistance(ztest[2])) + 25.0)
        self.assertEqual(mm[4], 5.0*numpy.log10(universe.luminosityDistance(ztest[4])) + 25.0)

    def testGetCurrent(self):
        """
        Test to make sure that getCurrent returns the activeCosmology
        """

        for Om0 in numpy.arange(start=0.2, stop=0.5, step=0.29):
            for Ode0 in numpy.arange(start=1.0-Om0-0.2, stop=1.0-Om0+0.2, step=0.39):
                for w0 in numpy.arange(start=-1.2, stop=-0.7, step=0.49):
                    for wa in numpy.arange(start=-0.2, stop=0.2, step=0.39):
                        universe = CosmologyObject(Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)
                        testUniverse = universe.getCurrent()

                        for zz in numpy.arange(start=1.0, stop=2.1, step=1.0):
                            self.assertEqual(universe.OmegaMatter(redshift=zz),
                                             testUniverse.Om(zz))
                            self.assertEqual(universe.OmegaDarkEnergy(redshift=zz),
                                             testUniverse.Ode(zz))
                            self.assertEqual(universe.OmegaPhotons(redshift=zz),
                                             testUniverse.Ogamma(zz))
                            self.assertEqual(universe.OmegaNeutrinos(redshift=zz),
                                             testUniverse.Onu(zz))
                            self.assertEqual(universe.OmegaCurvature(redshift=zz),
                                             testUniverse.Ok(zz))

class CosmologyMixinUnitTest(unittest.TestCase):
    """
    This class will test to make sure that our example CosmologyMixin
    (defined in lsst/sims/photUtils/examples/CosmologyMixin.py)
    can produce a catalog
    """

    @classmethod
    def setUpClass(cls):
        cls.dbName = 'cosmologyTestDB.db'
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)
        makeGalTestDB(size=100, seedVal=1, filename=cls.dbName)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        del cls.dbName

    def setUp(self):
        self.catName = 'cosmologyCatalog.txt'
        if os.path.exists(self.catName):
            os.unlink(self.catName)

    def tearDown(self):
        pass
        #if os.path.exists(self.catName):
        #    os.unlink(self.catName)
        #del self.catName

    def testCosmologyCatalog(self):
        """
        Does a catalog get written?
        """
        address = 'sqlite:///' + self.dbName
        dbObj = myTestGals(address=address)
        cat = cosmologicalGalaxyCatalog(dbObj)
        cat.write_catalog(self.catName)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(CosmologyUnitTest)
    suites += unittest.makeSuite(CosmologyMixinUnitTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
