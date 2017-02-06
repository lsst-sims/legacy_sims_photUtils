# This unittest is not a part of sims_sed_library because introducing code into
# sims_sed_library would force users to install a new copy every time any
# upstream lsst utility code changed.  This is the lowest level package that
# depends on sims_sed_library, so this is where we are putting the unit test.

import unittest
import os

import lsst.utils.tests
from lsst.utils import getPackageDir


def setup_module(module):
    lsst.utils.tests.init()


class SedLibraryContents(unittest.TestCase):
    """
    This TestCase will verify that the contents of sims_sed_library were
    correctly loaded.
    """

    longMessage = True

    def verify_dir(self, dir_name, n_files, min_size=10):
        """
        verify the contents of a sims_sed_library sub directory

        dir_name is the name of the sub-directory under sims_sed_library

        n_files is the number of files meant to be in that directory

        min_size is the size (in kb) that we demand all files be greater than
        """
        msg = 'failed on %s ' % dir_name
        kb = 1024

        target_dir = os.path.join(getPackageDir('sims_sed_library'),
                                  dir_name)

        list_of_files = os.listdir(target_dir)
        self.assertEqual(len(list_of_files), n_files, msg=msg)

        for file_name in list_of_files:
            full_name = os.path.join(target_dir, file_name)
            msg = 'failed on %s' % full_name
            self.assertGreater(os.path.getsize(full_name), min_size*kb,
                               msg=msg)

    def test_directories(self):
        self.verify_dir('starSED/kurucz', 4885)
        self.verify_dir('starSED/mlt', 869)
        self.verify_dir('starSED/phoSimMLT', 869)
        self.verify_dir('starSED/wDs', 1333)
        self.verify_dir('galaxySED', 959)
        self.verify_dir('agnSED', 1)
        self.verify_dir('igm', 30)
        self.verify_dir('cepheid_lc', 5, min_size=5.0)
        self.verify_dir('flatSED', 1)
        self.verify_dir('eb_lc', 1842)
        self.verify_dir('mflare', 50)
        self.verify_dir('microlens/bh_binary_source', 71, min_size=0.001)
        self.verify_dir('rrly_lc/RRab', 758)
        self.verify_dir('rrly_lc/RRc', 208)
        self.verify_dir('ssmSED', 26, min_size=9)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
