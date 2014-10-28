import unittest
import lsst.utils.tests as utilsTests
import numpy

from lsst.sims.photUtils import CosmologyWrapper

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

    def testFlatLCDM(self):
        w0=-1.0
        wa=0.0

        matter = numpy.arange(start=0.1, stop=0.5, step=0.1)
        hubble = numpy.arange(start=50.0, stop=90.0, step=5.0)

        for Om0 in matter:
            for H0 in hubble:
                Ode0 = 1.0 - Om0

                universe = CosmologyWrapper()
                universe.Initialize(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                Og0 = universe.OmegaPhotons(redshift=0.0)
                Onu0 = universe.OmegaNeutrinos(redshift=0.0)

                self.assertTrue(Og0 < 1.0e-4)
                self.assertTrue(Onu0 < 1.0e-4)
                self.assertAlmostEqual(universe.OmegaMatter(redshift=0.0), Om0, 10)
                self.assertAlmostEqual(Ode0 - universe.OmegaDarkEnergy(redshift=0.0), Og0+Onu0, 6)
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

    def testFlatW0Wa(self):

        matter = numpy.arange(start=0.1, stop=0.7, step=0.2)
        ww0 = numpy.arange(start=-1.0, stop = -0.5, step=0.1)
        wwa = numpy.arange(start = -0.3, stop = 0.3, step=0.1)
        hubble = numpy.arange(start=50.0, stop=90.0, step=10.0)

        for Om0 in matter:
            for H0 in hubble:
                for w0 in ww0:
                    for wa in wwa:
                        Ode0 = 1.0 - Om0

                        universe = CosmologyWrapper()
                        universe.Initialize(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

                        Og0 = universe.OmegaPhotons(redshift=0.0)
                        Onu0 = universe.OmegaNeutrinos(redshift=0.0)

                        self.assertTrue(Og0 < 1.0e-4)
                        self.assertTrue(Onu0 < 1.0e-4)
                        self.assertAlmostEqual(universe.OmegaMatter(redshift=0.0), Om0, 10)
                        self.assertAlmostEqual(Ode0 - universe.OmegaDarkEnergy(redshift=0.0), Og0+Onu0, 6)
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

    def testNonFlatLCDM(self):
        w0=-1.0
        wa=0.0

        matter = numpy.arange(start=0.1, stop=1.0, step=0.3)
        darkEnergy = numpy.arange(start=0.1, stop=1.0, step=0.3)
        hubble = numpy.arange(start=50.0, stop=90.0, step=10.0)

        for Om0 in matter:
            for H0 in hubble:
                for Ode0 in darkEnergy:

                    universe = CosmologyWrapper()
                    universe.Initialize(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

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

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(CosmologyUnitTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
