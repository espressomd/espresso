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
The CPU and the CUDA flux sweeps of the electrokinetics are linked into the
same library when ESPResSo is built with CUDA. The kernels of both targets
have the same kinds of parameters (pointers, sizes and strides), so a kernel
of a CPU file and a kernel of a CUDA file must never have the same qualified
name and parameter types: the linker would bind one to the other (the
``__global__`` kernels have host stubs) without a diagnostic (see
``test_boundary_symbols``). The kernels of a class are in the namespace
``<ClassName>_kernels`` of ``walberla::pystencils``.
"""

import re
import unittest

from tests import ek_generated
from tests.test_accessor_symbols import definitions

PRECISIONS = (("Double", "double"), ("Single", "single"))


class SweepSymbolsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cpu, cls.cuda, cls.names = {}, {}, {}
        for suffix, precision in PRECISIONS:
            directory = ek_generated.generated(precision)
            for kernel, stem in ek_generated.CLASS_NAMES.items():
                if kernel == "fixed_flux":
                    continue
                name = f"{stem}{suffix}Precision"
                if not (directory / f"{name}.cpp").exists():
                    continue
                cls.cpu[name] = cls.read(directory, name, ".cpp")
                cls.cuda[f"{name}CUDA"] = cls.read(
                    directory, f"{name}CUDA", ".cu")

    @classmethod
    def read(cls, directory, name, source):
        return definitions("\n".join(
            (directory / f"{name}{extension}").read_text()
            for extension in (".h", source)))

    def test_classes(self):
        # all the sweeps in double precision, the ones of
        # `ek_generated.SINGLE_PRECISION` in single precision
        self.assertEqual(len(self.cpu), 7 + 2)
        self.assertEqual(len(self.cuda), 7 + 2)

    def test_each_class_defines_one_kernel_in_its_namespace(self):
        for functions in (self.cpu, self.cuda):
            for name, defined in functions.items():
                kernels = {n for n, _ in defined if "_kernels::" in n}
                self.assertEqual(
                    kernels,
                    {f"walberla::pystencils::{name}_kernels::{name}"}, name)

    def test_cpu_and_cuda_functions_are_disjoint(self):
        cpu = set().union(*self.cpu.values())
        cuda = set().union(*self.cuda.values())
        self.assertEqual(
            cpu & cuda, set(),
            "functions with the same qualified name and parameter types "
            "on both targets")

    def test_kernel_names_are_unique(self):
        names = [qualified for functions in
                 (*self.cpu.values(), *self.cuda.values())
                 for qualified, _ in functions if "_kernels::" in qualified]
        self.assertEqual(len(names), len(self.cpu) + len(self.cuda))
        self.assertEqual(len(names), len(set(names)))

    def test_each_class_defines_its_own_methods(self):
        for suffix, precision in PRECISIONS:
            directory = ek_generated.generated(precision)
            for name in self.cpu:
                if not name.endswith(f"{suffix}Precision"):
                    continue
                for target, extension in (("", ".cpp"), ("CUDA", ".cu")):
                    source = (directory
                              / f"{name}{target}{extension}").read_text()
                    for method in (r"operator\(\)", "runOnCellInterval"):
                        self.assertEqual(len(re.findall(
                            rf"\b{name}{target}::{method}\(", source)), 1,
                            (name, target, method))
