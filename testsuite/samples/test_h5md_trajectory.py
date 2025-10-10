#
# Copyright (C) 2022 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import unittest as ut
import numpy as np
import importlib_wrapper

sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/h5md_trajectory.py")


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_data_integrity(self):
        n_frames = 4
        self.assertEqual(sample.pos_folded.shape[0], n_frames)
        self.assertEqual(sample.xyz_folded.shape[0], n_frames)
        self.assertEqual(sample.pos_unfolded.shape[0], n_frames)
        self.assertEqual(sample.xyz_unfolded.shape[0], n_frames)
        np.testing.assert_allclose(
            sample.pos_folded,
            sample.xyz_folded,
            atol=1e-12)
        np.testing.assert_allclose(
            sample.pos_unfolded,
            sample.xyz_unfolded,
            atol=1e-12)


if __name__ == "__main__":
    ut.main()
