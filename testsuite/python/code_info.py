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
import espressomd
import espressomd.code_info
import espressomd.version


class Test(ut.TestCase):

    def test_code_info(self):
        # check CMake build type
        build_types = {
            "Debug", "Release", "RelWithDebInfo", "MinSizeRel", "Coverage",
            "RelWithAssert"}
        self.assertIn(espressomd.code_info.build_type(), build_types)

        # check features
        features = espressomd.code_info.features()
        all_features = espressomd.code_info.all_features()
        self.assertTrue(set(features).issubset(all_features))

        # check arrays are sorted
        scafacos_methods = espressomd.code_info.scafacos_methods()
        self.assertEqual(features, sorted(features))
        self.assertEqual(all_features, sorted(all_features))
        self.assertEqual(scafacos_methods, sorted(scafacos_methods))

    def test_version(self):
        version_full = espressomd.version.version()
        version_major_minor = (espressomd.version.major(),
                               espressomd.version.minor())
        self.assertTrue(all(x >= 0 for x in version_full))
        self.assertIn(len(version_full), (2, 3))
        self.assertEqual(version_full[:2], version_major_minor)
        self.assertEqual(".".join(map(str, espressomd.version.version())),
                         espressomd.version.friendly())

    def test_git_info(self):
        git_states = {"CLEAN", "DIRTY"}
        commit_charset = set("abcdef0123456789")
        self.assertIn(espressomd.version.git_state(), git_states)
        self.assertIsInstance(espressomd.version.git_branch(), str)
        self.assertIsInstance(espressomd.version.git_commit(), str)
        git_commit = espressomd.version.git_commit()
        self.assertLessEqual(len(git_commit), 40)
        self.assertTrue(set(git_commit).issubset(commit_charset))


if __name__ == "__main__":
    ut.main()
