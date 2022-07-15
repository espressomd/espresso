#
# Copyright (C) 2019-2022 The ESPResSo project
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
import espressomd


class Features(ut.TestCase):

    def test_has_features(self):
        for feature in espressomd.code_info.features():
            self.assertTrue(espressomd.has_features(feature))

        for feature in espressomd.code_info.all_features(
        ) - set(espressomd.code_info.features()):
            self.assertFalse(espressomd.has_features(feature))

        with self.assertRaises(RuntimeError) as _:
            espressomd.has_features("NotAFeature")


if __name__ == '__main__':
    ut.main()
