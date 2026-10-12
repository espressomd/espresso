#
# Copyright (C) 2026 The ESPResSo project
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

import pathlib
import tempfile
import unittest

import vendor_pystencils_runtime as vendor


def write(path, content=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


class VendorPystencilsRuntimeTest(unittest.TestCase):

    def setUp(self):
        tmp = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.include_dir = tmp / "include"
        runtime = self.include_dir / "pystencils_runtime"
        write(runtime / "generic_cpu.hpp", '#include "./random.hpp"\n')
        write(runtime / "random.hpp",
              '#include "./bits/philox_rand.h"\n#include <cmath>\n')
        write(runtime / "bits" / "philox_rand.h",
              '#include "./myintrin.h"\n#include "opencl_stdint.h"\n')
        write(runtime / "bits" / "myintrin.h")
        write(runtime / "unused.hpp")
        self.kernels = tmp / "generated_kernels"
        write(self.kernels / "Kernel.cpp",
              '#include "Kernel.h"\n'
              '#include "pystencils_runtime/generic_cpu.hpp"\n')
        write(self.kernels / "Kernel.h", '#include "core/DataTypes.h"\n')
        self.dest = tmp / "third_party"

    def test_needed_headers_follow_includes_transitively(self):
        roots = vendor.runtime_includes([self.kernels])
        self.assertEqual(roots, {"pystencils_runtime/generic_cpu.hpp"})
        needed = vendor.needed_headers(self.include_dir, roots)
        self.assertEqual(needed, [
            "pystencils_runtime/bits/myintrin.h",
            "pystencils_runtime/bits/philox_rand.h",
            "pystencils_runtime/generic_cpu.hpp",
            "pystencils_runtime/random.hpp",
        ])

    def test_vendor_replaces_destination(self):
        write(self.dest / "pystencils_runtime" / "stale.hpp")
        copied = vendor.vendor(self.include_dir, [self.kernels], self.dest)
        self.assertEqual(len(copied), 4)
        found = sorted(str(p.relative_to(self.dest))
                       for p in self.dest.rglob("*") if p.is_file())
        self.assertEqual(found, copied)

    def test_missing_root_is_an_error(self):
        write(self.kernels / "Other.cu",
              '#include "pystencils_runtime/cuda.cuh"\n')
        with self.assertRaisesRegex(FileNotFoundError, "cuda.cuh"):
            vendor.vendor(self.include_dir, [self.kernels], self.dest)


if __name__ == "__main__":
    unittest.main()
