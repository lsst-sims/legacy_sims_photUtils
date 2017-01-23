import unittest
import numpy as np
import lsst.utils.tests

from lsst.sims.photUtils import read_close_Kuruz


def setup_module(module):
    lsst.utils.tests.init()


class readKuruzTest(unittest.TestCase):

    def testRead(self):
        sed1, paramDict = read_close_Kuruz(6000., 0., 4.4)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()

