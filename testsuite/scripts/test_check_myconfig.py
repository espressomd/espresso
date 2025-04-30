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

import sys
import pathlib
import tempfile
import importlib
import unittest as ut

sys.path.insert(0, "@CMAKE_SOURCE_DIR@/src/config")
module = importlib.import_module("check_myconfig")


class Test(ut.TestCase):

    def test_exceptions(self):
        tmp_directory = tempfile.TemporaryDirectory()
        path_root = pathlib.Path(tmp_directory.name)
        path_myconfig = path_root / "myconfig.hpp"
        path_cmake_config = path_root / "cmake_config.hpp"
        path_cmakedefine = path_root / "cmake_config.cmakein"
        with open("@CMAKE_SOURCE_DIR@/cmake/espresso_cmake_config.cmakein") as fp:
            default_cmakedefine = fp.read()

        def run(myconfig="", cmake_config="", cmakedefine=default_cmakedefine):
            path_myconfig.write_text(myconfig)
            path_cmake_config.write_text(cmake_config)
            path_cmakedefine.write_text(cmakedefine)
            module.check_myconfig(
                "@CMAKE_CXX_COMPILER@",
                "@CMAKE_SOURCE_DIR@/src/config/features.def",
                str(path_cmakedefine),
                str(path_myconfig),
                str(path_cmake_config))

        # no error for empty configurations
        run()
        # no error for multi-line comments
        run(cmakedefine=default_cmakedefine + "/*\n#cmake" + "define ERROR\n*/")

        # check all errors
        with self.assertRaisesRegex(RuntimeError, "unknown feature 'ROTATIONN', did you mean 'ROTATION'"):
            run(myconfig="#define ROTATIONN")
        with self.assertRaisesRegex(RuntimeError, "unknown feature 'UNKNOWN_FEATURE'"):
            run(myconfig="#define UNKNOWN_FEATURE")
        with self.assertRaisesRegex(RuntimeError, "external feature 'FFTW' cannot be defined in myconfig"):
            run(myconfig="#define FFTW")
        with self.assertRaisesRegex(RuntimeError, "cmake_config.hpp` returned non-zero exit code 1, output:"):
            run(cmake_config="#define")
        with self.assertRaisesRegex(RuntimeError, "myconfig.hpp` returned non-zero exit code 1, output:"):
            run(myconfig="#define")
        with self.assertRaisesRegex(RuntimeError, "external feature 'UNKNOWN' is missing from '.*features.def'"):
            run(cmakedefine="#cmake" + "define ESPRESSO_BUILD_WITH_UNKNOWN")
        with self.assertRaisesRegex(RuntimeError, "cmakedefine 'FFTW' is missing from '.*cmake_config.cmakein'"):
            run(cmakedefine="")

        tmp_directory.cleanup()


if __name__ == "__main__":
    ut.main()
