import numpy as np

import unittest as ut
import espressomd.utils


class UtilsTest(ut.TestCase):

    def test_is_valid_type(self):
        self.assertFalse(espressomd.utils.is_valid_type(np.inf, float))
        self.assertFalse(espressomd.utils.is_valid_type(np.nan, float))
        self.assertFalse(espressomd.utils.is_valid_type(None, float))
        self.assertFalse(espressomd.utils.is_valid_type(float('nan'), float))
        self.assertFalse(espressomd.utils.is_valid_type(float('inf'), float))
        self.assertFalse(espressomd.utils.is_valid_type(np.inf, int))
        self.assertFalse(espressomd.utils.is_valid_type(np.nan, int))
        self.assertFalse(espressomd.utils.is_valid_type(None, int))

if __name__ == "__main__":
    ut.main()
