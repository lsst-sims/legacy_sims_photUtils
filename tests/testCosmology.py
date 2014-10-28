import unittest
import lsst.utils.tests as utilsTests
import numpy

from lsst.sims.photUtils import CosmologyWrapper

class CosmologyUnitTest(unittest.TestCase):

    def testParams(self):
        universe = CosmologyWrapper()
        Om0=0.24
        Ode0=0.76
        w0=-1.0
        wa=0.0
        H0=86.0
        universe.Initialize(H0=H0, Om0=Om0, Ode0=Ode0, w0=w0, wa=wa)

        Omg = universe.OmegaPhotons(redshift=0.0)
        Onu = universe.OmegaNeutrinos(redshift=0.0)

        self.assertTrue(Omg < 1.0e-4)
        self.assertTrue(Onu < 1.0e-4)
        self.assertTrue(numpy.abs(universe.OmegaMatter(redshift=0.0)-Om0) < Omg+Onu)
        self.assertTrue(numpy.abs(universe.OmegaDarkEnergy(redshift=0.0)-Ode0) < Omg+Onu)
        self.assertEqual(universe.H(redshift=0.0),H0)
        self.assertEqual(universe.OmegaCurvature(),0.0)



def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(CosmologyUnitTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
