import os
import unittest
import lsst.utils.tests as utilsTests
import numpy
import scipy

from lsst.sims.photUtils import CosmologyWrapper, PhotometryGalaxies
from lsst.sims.photUtils.examples import ExampleCosmologyMixin
from lsst.sims.catalogs.measures.instance import InstanceCatalog

from lsst.sims.catalogs.generation.utils import myTestGals, makeGalTestDB
from lsst.sims.photUtils.utils import testGalaxies

class cosmologicalGalaxyCatalog(testGalaxies, ExampleCosmologyMixin):
    column_outputs = ['uAbs', 'gAbs', 'rAbs', 'iAbs', 'zAbs', 'yAbs',
                      'uBulgeAbs', 'gBulgeAbs', 'rBulgeAbs', 'iBulgeAbs', 'zBulgeAbs', 'yBulgeAbs',
                      'uDiskAbs', 'gDiskAbs', 'rDiskAbs', 'iDiskAbs', 'zDiskAbs', 'yDiskAbs',
                      'uAgnAbs', 'gAgnAbs', 'rAgnAbs', 'iAgnAbs', 'zAgnAbs', 'yAgnAbs',
                      'redshift', 'distanceModulus']

def controlOmega(redshift, H0, Om0, Ode0 = None, Og0=0.0, Onu0=0.0, w0=-1.0, wa=0.0):
    if Ode0 is None:
        Ode0 = 1.0 - Om0 - Og0 - Onu0

    Ok0 = 1.0 - Om0 - Ode0 - Og0 - Onu0

    aa = 1.0/(1.0+redshift)
    Omz = Om0 * numpy.power(1.0+redshift, 3)
    Ogz = Og0 * numpy.power(1.0+redshift, 4)
    Onuz = Onu0 * numpy.power(1.0+redshift, 4)
    Okz = Ok0 * numpy.power(1.0+redshift, 2)
    Odez = Ode0 * numpy.exp(-3.0*(numpy.log(aa)*(w0 + wa +1.0) - wa*(aa - 1.0)))

    Ototal = Omz + Ogz + Onuz + Odez + Okz


    return Omz/Ototal, Ogz/Ototal, Onuz/Ototal, Odez/Ototal, Okz/Ototal, H0*numpy.sqrt(Ototal)

class CosmologyUnitTest(unittest.TestCase):

    def setUp(self):
        self.speedOfLight = 2.9979e5

    def tearDown(self):
        del self.speedOfLight

    #@unittest.skip("fornow")
    def testFlatLCDM(self):
        H0 = 50.0
        for Om0 in numpy.arange(start=0.1, stop=0.91, step=0.4):
            universe = CosmologyWrapper()
            universe.initializeCosmology(H0=H0, Om0=Om0)

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

            ztest = numpy.arange(start=0.0, stop=4.0, step=0.5)
            for zz in ztest:
               aa = (1.0+zz)

               OmControl, OgControl, OnuControl, \
                   OdeControl, OkControl, Hcontrol = controlOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0)

               self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
               self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
               self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
               self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
               self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

            del universe

    #@unittest.skip("fornow")
    def testFlatW0Wa(self):

        H0 = 96.0
        for Om0 in numpy.arange(start=0.1, stop=0.95, step=0.4):
            for w0 in numpy.arange(start=-1.1, stop=-0.89, step=0.2):
                for wa in numpy.arange(start=-0.1, stop=0.11, step=0.2):

                    universe = CosmologyWrapper()
                    universe.initializeCosmology(H0=H0, Om0=Om0, w0=w0, wa=wa)

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

                    ztest = numpy.arange(start=0.0, stop=4.0, step=0.5)
                    for zz in ztest:

                       wControl = w0 + wa*(1.0 - 1.0/(1.0+zz))
                       self.assertAlmostEqual(wControl, universe.w(redshift=zz), 6)

                       OmControl, OgControl, OnuControl, \
                       OdeControl, OkControl, Hcontrol = controlOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0,
                                                                          w0=w0, wa=wa)

                       self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
                       self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
                       self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
                       self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
                       self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

                    del universe

    #@unittest.skip("fornow")
    def testNonFlatLCDM(self):
        w0 = -1.0
        wa = 0.0
        H0 = 77.0

        for Om0 in numpy.arange(start=0.15, stop=0.96, step=0.4):
            for Ode0 in numpy.arange(start=1.0-Om0-0.1, stop=1.0-Om0+0.11, step=0.2):

                universe = CosmologyWrapper()
                universe.initializeCosmology(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

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

                ztest = numpy.arange(start=0.0, stop=4.0, step=0.5)
                for zz in ztest:
                    OmControl, OgControl, OnuControl, \
                    OdeControl, OkControl, Hcontrol = controlOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0,
                                                                          Ode0=Ode0)

                    self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
                    self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
                    self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
                    self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
                    self.assertAlmostEqual(OkControl, universe.OmegaCurvature(redshift=zz), 6)
                    self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

                del universe

    #@unittest.skip("fornow")
    def testNonFlatW0Wa(self):

        H0 = 60.0

        for Om0 in numpy.arange(start=0.15, stop=0.76, step=0.3):
            for Ode0 in numpy.arange(1.0-Om0-0.1, stop = 1.0-Om0+0.11, step=0.2):
                for w0 in numpy.arange(start=-1.1, stop = -0.89, step=0.1):
                    for wa in numpy.arange(start=-0.1, stop=0.15, step=0.1):

                        universe = CosmologyWrapper()
                        universe.initializeCosmology(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

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

                        ztest = numpy.arange(start=0.0, stop=4.0, step=1.0)
                        for zz in ztest:

                           wControl = w0 + wa*(1.0 - 1.0/(1.0+zz))
                           self.assertAlmostEqual(wControl, universe.w(redshift=zz), 6)

                           OmControl, OgControl, OnuControl, \
                           OdeControl, OkControl, Hcontrol = controlOmega(zz, H0, Om0, Og0=Og0, Onu0=Onu0,
                                                                          w0=w0, wa=wa, Ode0=Ode0)

                           self.assertAlmostEqual(OmControl, universe.OmegaMatter(redshift=zz), 6)
                           self.assertAlmostEqual(OdeControl, universe.OmegaDarkEnergy(redshift=zz), 6)
                           self.assertAlmostEqual(OgControl, universe.OmegaPhotons(redshift=zz), 6)
                           self.assertAlmostEqual(OnuControl, universe.OmegaNeutrinos(redshift=zz), 6)
                           self.assertAlmostEqual(OkControl, universe.OmegaCurvature(redshift=zz), 6)
                           self.assertAlmostEqual(Hcontrol, universe.H(redshift=zz), 6)

                        del universe

    #@unittest.skip("fornow")
    def testComovingDistance(self):

        universe = CosmologyWrapper()
        H0 = 73.0
        for Om0 in numpy.arange(start=0.1, stop=0.55, step=0.2):
            for Ode0 in numpy.arange(start=1.0-Om0-0.05, stop=1.0-Om0+0.06, step=0.05):
                for w0 in numpy.arange(start=-1.1, stop=-0.85, step=0.1):
                    for wa in numpy.arange(start=-0.1, stop=0.115, step=0.05):

                        universe.initializeCosmology(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                        ztest = numpy.arange(start=0.1, stop=2.0, step=0.3)
                        for zz in ztest:
                            comovingControl = universe.comovingDistance(redshift=zz)
                            comovingTest = self.speedOfLight*scipy.integrate.quad(lambda z: 1.0/universe.H(z), 0.0, zz)[0]
                            self.assertAlmostEqual(comovingControl/comovingTest,1.0,4)

    #@unittest.skip("fornow")
    def testLuminosityDistance(self):

        H0 = 73.0

        universe=CosmologyWrapper()
        for Om0 in numpy.arange(start=0.1, stop=0.55, step=0.2):
            for Ode0 in numpy.arange(start=1.0-Om0-0.05, stop=1.0-Om0+0.06, step=0.05):
                for w0 in numpy.arange(start=-1.1, stop=-0.85, step=0.1):
                    for wa in numpy.arange(start=-0.1, stop=0.115, step=0.05):

                        universe.initializeCosmology(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                        ztest = numpy.arange(start=0.1, stop=2.0, step=0.3)

                        sqrtkCurvature = numpy.sqrt(numpy.abs(universe.OmegaCurvature()))*universe.H()/self.speedOfLight

                        for zz in ztest:
                            luminosityControl = universe.luminosityDistance(redshift=zz)
                            comovingDistance = self.speedOfLight*scipy.integrate.quad(lambda z: 1.0/universe.H(z), 0.0, zz)[0]

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

    #@unittest.skip("fornow")
    def testAngularDiameterDistance(self):

        H0 = 56.0
        universe=CosmologyWrapper()
        for Om0 in numpy.arange(start=0.1, stop=0.55, step=0.2):
            for Ode0 in numpy.arange(start=1.0-Om0-0.05, stop=1.0-Om0+0.06, step=0.05):
                for w0 in numpy.arange(start=-1.1, stop=-0.85, step=0.1):
                    for wa in numpy.arange(start=-0.1, stop=0.115, step=0.05):

                        universe.initializeCosmology(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                        ztest = numpy.arange(start=0.1, stop=2.0, step=0.3)

                        sqrtkCurvature = numpy.sqrt(numpy.abs(universe.OmegaCurvature()))*universe.H()/self.speedOfLight

                        for zz in ztest:
                            angularControl = universe.angularDiameterDistance(redshift=zz)
                            comovingDistance = self.speedOfLight*scipy.integrate.quad(lambda z: 1.0/universe.H(z), 0.0, zz)[0]

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

        H0 = 73.0

        universe=CosmologyWrapper()
        for Om0 in numpy.arange(start=0.1, stop=0.55, step=0.2):
            for Ode0 in numpy.arange(start=1.0-Om0-0.05, stop=1.0-Om0+0.06, step=0.05):
                for w0 in numpy.arange(start=-1.1, stop=-0.85, step=0.1):
                    for wa in numpy.arange(start=-0.1, stop=0.115, step=0.05):

                        universe.initializeCosmology(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                        ztest = numpy.arange(start=0.1, stop=2.0, step=0.3)

                        sqrtkCurvature = numpy.sqrt(numpy.abs(universe.OmegaCurvature()))*universe.H()/self.speedOfLight

                        for zz in ztest:
                            modulusControl = universe.distanceModulus(redshift=zz)
                            comovingDistance = self.speedOfLight*scipy.integrate.quad(lambda z: 1.0/universe.H(z), 0.0, zz)[0]

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

class CosmologyMixinUnitTest(unittest.TestCase):

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
        if os.path.exists(self.catName):
            os.unlink(self.catName)
        del self.catName

    def testCosmologyCatalog(self):
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
