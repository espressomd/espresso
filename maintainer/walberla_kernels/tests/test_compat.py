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

import pystencils as ps
import sympy as sp

import compat
import sfg_driver


def toy_assignments(ctx):
    a, b = ps.fields(f"a, b: {ctx.data_type.c_string()}[3D]", layout="fzyx")
    alpha, omega = sp.symbols("alpha omega")
    return [ps.Assignment(a.center, alpha * b.center + omega)]


class LegacyAccessorSuffixTest(unittest.TestCase):

    def test_suffix(self):
        self.assertEqual(compat.legacy_accessor_suffix("kT"), "Kt")
        self.assertEqual(compat.legacy_accessor_suffix("rho_0"), "Rho_0")
        self.assertEqual(compat.legacy_accessor_suffix("time_step"),
                         "Time_step")


class LegacySweepTest(unittest.TestCase):

    def generate(self, target, **kwargs):
        tmp = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        ctx = sfg_driver.CodegenContext("double", target)
        with sfg_driver.KernelFileGenerator("Toy", ctx, tmp, "test.py") as sfg:
            compat.legacy_sweep(sfg, ctx, "Toy", toy_assignments(ctx),
                                **kwargs)
        return (tmp / "Toy.h").read_text()

    def test_cpu_interface(self):
        header = self.generate("cpu")
        expected = [
            "class Toy {",
            "Toy (walberla::BlockDataID aId, walberla::BlockDataID bId, "
            "double alpha, double omega)",
            "void run(walberla::IBlock * block)",
            "impl_(block);",
            "void runOnCellInterval(const walberla::shared_ptr< "
            "walberla::StructuredBlockStorage >& blocks, "
            "const walberla::CellInterval& globalCellInterval, "
            "walberla::cell_idx_t ghostLayers, walberla::IBlock * block)",
            "impl_.runOnCellInterval(block, ci);",
            "void operator()(walberla::IBlock * block)",
            "static std::function< void(walberla::IBlock *) > getSweep(",
            "getSweepOnCellInterval(",
            "void configure(",
            "double getAlpha() const",
            "void setAlpha(double value)",
            "double getOmega() const",
            "void setOmega(double value)",
            "internal::Toy impl_;",
        ]
        for snippet in expected:
            self.assertIn(snippet, header)
        self.assertNotIn("gpuStream_t", header)

    def test_cuda_interface(self):
        header = self.generate("cuda")
        expected = [
            "void run(walberla::IBlock * block, gpuStream_t stream = nullptr)",
            "impl_(block, walberla::v8::sweep::exectag::GPU{stream});",
            "static std::function< void(walberla::IBlock *, gpuStream_t) > "
            "getSweepOnCellInterval(",
            "std::function< void(walberla::IBlock *) > "
            "getSweep(gpuStream_t stream = nullptr)",
        ]
        for snippet in expected:
            self.assertIn(snippet, header)

    def test_explicit_constructor_order(self):
        header = self.generate(
            "cpu", ctor_order=["bId", "aId", "omega", "alpha"])
        self.assertIn("Toy (walberla::BlockDataID bId, "
                      "walberla::BlockDataID aId, double omega, double alpha)",
                      header)

    def test_wrong_constructor_order(self):
        with self.assertRaisesRegex(ValueError, "constructor order"):
            self.generate("cpu", ctor_order=["aId", "bId", "alpha"])


if __name__ == "__main__":
    unittest.main()
