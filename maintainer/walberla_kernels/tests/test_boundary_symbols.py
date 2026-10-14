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
The CPU and the CUDA boundary classes are linked into the same library
when ESPResSo is built with CUDA. The link kernels of both targets have the
same kinds of parameters (pointers, sizes and strides), so a kernel of a CPU
file and a kernel of a CUDA file must never have the same qualified name and
parameter types: the linker would bind one to the other (the ``__global__``
kernels have host stubs) without a diagnostic. This is the check of
``test_pack_info_symbols`` for the boundaries. The hand-written base class of
the boundaries is included by both targets: it must only define templates,
whose instances for the host and the device are different types.
"""

import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

from tests.test_accessor_symbols import definitions

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "generate_lb_kernels.py"
BASE_HEADER = (pathlib.Path(__file__).resolve().parents[3] / "src" /
               "walberla_bridge" / "src" / "utils" /
               "IndexVectorBoundary.hpp")
PRECISIONS = ("Double", "Single")


class BoundarySymbolsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.directory = pathlib.Path(tempfile.mkdtemp())
        subprocess.run(
            [sys.executable, str(SCRIPT), "--kernels", "dynamic_ubb",
             "--target", "cpu", "cuda", "--output-dir", str(cls.directory)],
            capture_output=True, text=True, check=True)
        cls.cpu, cls.cuda, cls.sources = {}, {}, {}
        for precision in PRECISIONS:
            name = f"DynamicUBB{precision}Precision"
            cls.cpu[name] = cls.read(name, ".cpp")
            cls.cuda[f"{name}CUDA"] = cls.read(f"{name}CUDA", ".cu")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.directory)

    @classmethod
    def read(cls, name, source):
        texts = [(cls.directory / f"{name}{extension}").read_text()
                 for extension in (".h", source)]
        cls.sources[name] = texts[1]
        return definitions("\n".join(texts))

    def test_each_class_defines_one_kernel_in_its_namespace(self):
        for functions in (self.cpu, self.cuda):
            for name, defined in functions.items():
                kernels = {n for n, _ in defined if "_kernels::" in n}
                self.assertEqual(
                    kernels, {f"walberla::lbm::{name}_kernels::linkKernel"})

    def test_each_class_defines_its_own_method(self):
        for name, source in self.sources.items():
            self.assertEqual(
                len(re.findall(rf"\b{name}::linkKernel\(", source)), 1, name)

    def test_cpu_and_cuda_functions_are_disjoint(self):
        cpu = set().union(*self.cpu.values())
        cuda = set().union(*self.cuda.values())
        self.assertEqual(
            cpu & cuda, set(),
            "functions with the same qualified name and parameter types "
            "on both targets")

    def test_kernel_names_are_unique_among_all_classes(self):
        names = [qualified for functions in
                 (*self.cpu.values(), *self.cuda.values())
                 for qualified, _ in functions if "_kernels::" in qualified]
        self.assertEqual(len(names), 4)
        self.assertEqual(len(names), len(set(names)))

    def test_base_header_defines_no_function(self):
        # only templates: no function can be defined by the CPU and by the
        # CUDA translation units with the same name and signature
        text = BASE_HEADER.read_text()
        self.assertEqual(definitions(text), set())
        for keyword in ("__global__", "__device__", "__host__"):
            self.assertNotIn(keyword, text)


if __name__ == "__main__":
    unittest.main()
