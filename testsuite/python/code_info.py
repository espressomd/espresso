#
# Copyright (C) 2022-2026 The ESPResSo project
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

import os
import unittest as ut
import packaging.specifiers
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

        self.assertIsNone(
            espressomd.code_info._CodeInfo().call_method("unknown"))

    def test_version(self):
        version_full = espressomd.version.version()
        version_major_minor = (espressomd.version.major(),
                               espressomd.version.minor())
        self.assertTrue(all(x >= 0 for x in version_full))
        self.assertIn(len(version_full), (2, 3))
        self.assertEqual(version_full[:2], version_major_minor)
        self.assertEqual(espressomd.version.friendly(), espressomd.__version__)
        self.assertEqual(".".join(map(str, espressomd.version.version())),
                         espressomd.version.friendly())
        self.assertIsNone(espressomd.version._Version().call_method("unknown"))

    def test_git_info(self):
        git_states = {"CLEAN", "DIRTY"}
        commit_charset = set("abcdef0123456789")
        self.assertIsInstance(espressomd.version.git_branch(), str)
        self.assertIsInstance(espressomd.version.git_commit(), str)
        self.assertIsInstance(espressomd.version.git_state(), str)
        git_commit = espressomd.version.git_commit()
        self.assertLessEqual(len(git_commit), 40)
        self.assertTrue(set(git_commit).issubset(commit_charset))
        if git_commit:
            self.assertIn(espressomd.version.git_state(), git_states)

    def test_build_info(self):
        def check_specs(version, specs):
            return packaging.specifiers.SpecifierSet(specs).contains(version)
        toolchain = espressomd.code_info.toolchain()
        self.assertIsNotNone(toolchain["CMAKE_C_COMPILER_ID"])
        self.assertIsNotNone(toolchain["CMAKE_C_COMPILER_VERSION"])
        self.assertIsNotNone(toolchain["CMAKE_CXX_COMPILER_ID"])
        self.assertIsNotNone(toolchain["CMAKE_CXX_COMPILER_VERSION"])
        if espressomd.has_features(["CUDA"]):
            self.assertIsNotNone(toolchain["CMAKE_CUDA_COMPILER_ID"])
            self.assertIsNotNone(toolchain["CMAKE_CUDA_COMPILER_VERSION"])
        else:
            self.assertIsNone(toolchain["CMAKE_CUDA_COMPILER_ID"])
            self.assertIsNone(toolchain["CMAKE_CUDA_COMPILER_VERSION"])
            self.assertIsNone(toolchain["CMAKE_CUDA_HOST_COMPILER_ID"])
            self.assertIsNone(toolchain["CMAKE_CUDA_HOST_COMPILER_VERSION"])
        self.assertIsNotNone(toolchain["ESPRESSO_MPIEXEC_VENDOR"])
        self.assertIsNotNone(toolchain["ESPRESSO_MPIEXEC_VERSION"])
        self.assertIsNotNone(toolchain["OpenMP_VERSION"])
        self.assertIsNotNone(toolchain["OpenMP_CXX_VERSION"])
        self.assertTrue(check_specs(toolchain["OpenMP_CXX_VERSION"], ">=4.5"))
        icp_ci_infrastructure = "https://gitlab.icp.uni-stuttgart.de/espressomd/espresso"
        if os.environ.get("CI_PROJECT_URL") == icp_ci_infrastructure:
            specs = {
                "GNU": ">=12.2.0",
                "Clang": ">=18.1.0",
                "NVHPC": ">=25.5",
                "NVIDIA": ">=12.0",
            }
            self.assertIn(
                toolchain["CMAKE_C_COMPILER_ID"], ["GNU", "Clang", "NVHPC"])
            self.assertIn(
                toolchain["CMAKE_CXX_COMPILER_ID"], ["GNU", "Clang", "NVHPC"])
            if espressomd.has_features(["CUDA"]):
                self.assertIn(
                    toolchain["CMAKE_CUDA_COMPILER_ID"], ["Clang", "NVIDIA"])
                self.assertIn(
                    toolchain["CMAKE_CUDA_HOST_COMPILER_ID"],
                    ["GNU", "NVHPC", None])
            self.assertTrue(check_specs(
                toolchain["CMAKE_C_COMPILER_VERSION"],
                specs[toolchain["CMAKE_C_COMPILER_ID"]]))
            self.assertTrue(check_specs(
                toolchain["CMAKE_CXX_COMPILER_VERSION"],
                specs[toolchain["CMAKE_CXX_COMPILER_ID"]]))
            if espressomd.has_features(["CUDA"]):
                self.assertTrue(check_specs(
                    toolchain["CMAKE_CUDA_COMPILER_VERSION"],
                    specs[toolchain["CMAKE_CUDA_COMPILER_ID"]]))
            specs = {"OpenMPI": ">=4.1.6", "MPICH": ">=4.2"}
            self.assertIn(
                toolchain["ESPRESSO_MPIEXEC_VENDOR"], ["OpenMPI", "MPICH"])
            self.assertTrue(check_specs(
                toolchain["ESPRESSO_MPIEXEC_VERSION"],
                specs[toolchain["ESPRESSO_MPIEXEC_VENDOR"]]))


if __name__ == "__main__":
    ut.main()
