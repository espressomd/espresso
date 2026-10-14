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

"""
The CPU and the CUDA pack infos are linked into the same library when
ESPResSo is built with CUDA. The kernels of both targets have the same
kinds of parameters (pointers, sizes and strides), so a kernel of a CPU
file and a kernel of a CUDA file must never have the same qualified name
and parameter types: the linker would bind one to the other (the
``__global__`` kernels have host stubs) without a diagnostic. This is the
check of ``test_accessor_symbols`` for the pack infos.
"""

import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

from tests.test_accessor_symbols import definitions

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "generate_lb_kernels.py"
STEMS = ("PackInfoPdf", "PackInfoVec")
PRECISIONS = ("Double", "Single")
# distinct sets of populations of D3Q19: 6 faces and 12 edges, and
# one set for all directions of a vector field; each has 3 kernels
KERNEL_COUNTS = {"PackInfoPdf": 3 * 18, "PackInfoVec": 3 * 1}


class PackInfoSymbolsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.directory = pathlib.Path(tempfile.mkdtemp())
        subprocess.run(
            [sys.executable, str(SCRIPT), "--kernels", "pack_info_pdf",
             "pack_info_vec", "--target", "cpu", "cuda",
             "--output-dir", str(cls.directory)],
            capture_output=True, text=True, check=True)
        cls.cpu, cls.cuda = {}, {}
        for stem in STEMS:
            for precision in PRECISIONS:
                name = f"{stem}{precision}Precision"
                cls.cpu[name] = cls.read(name, ".h", ".cpp")
                cls.cuda[name] = cls.read(f"{name}CUDA", ".h", ".cu")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.directory)

    @classmethod
    def read(cls, name, header, source):
        texts = [(cls.directory / f"{name}{extension}").read_text()
                 for extension in (header, source)]
        return definitions("\n".join(texts))

    def test_functions_are_found(self):
        for stem in STEMS:
            for precision in PRECISIONS:
                name = f"{stem}{precision}Precision"
                for target, functions in (("CPU", self.cpu[name]),
                                          ("CUDA", self.cuda[name])):
                    kernels = {n for n, _ in functions if "_kernels::" in n}
                    self.assertEqual(len(kernels), KERNEL_COUNTS[stem],
                                     f"{name} {target}")

    def test_kernels_are_in_the_namespace_of_their_class(self):
        for name, functions in self.cpu.items():
            for qualified, _ in functions:
                if "_kernels::" in qualified:
                    self.assertTrue(qualified.startswith(
                        f"walberla::pystencils::{name}_kernels::"),
                        qualified)
        for name, functions in self.cuda.items():
            for qualified, _ in functions:
                if "_kernels::" in qualified:
                    self.assertTrue(qualified.startswith(
                        f"walberla::pystencils::{name}CUDA_kernels::"),
                        qualified)

    def test_cpu_and_cuda_functions_are_disjoint(self):
        cpu = set().union(*self.cpu.values())
        cuda = set().union(*self.cuda.values())
        self.assertEqual(
            cpu & cuda, set(),
            "functions with the same qualified name and parameter types "
            "on both targets")

    def test_kernel_names_are_unique_among_all_pack_infos(self):
        # the kernels of two classes of the same target are also never
        # defined twice
        names = [qualified for functions in
                 (*self.cpu.values(), *self.cuda.values())
                 for qualified, _ in functions if "_kernels::" in qualified]
        self.assertEqual(len(names), len(set(names)))


if __name__ == "__main__":
    unittest.main()
