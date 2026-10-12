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
from unittest import mock

import pystencils as ps
import sympy as sp
from pystencils import DynamicType
from pystencils.sympyextensions import tcast

import compat
import sfg_driver
import sweepgen
import sweepgen.symbolic


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
            "impl_.runOnCellInterval(block, ci, "
            "walberla::v8::sweep::exectag::GPU{stream});",
            "(walberla::IBlock *b, gpuStream_t stream = nullptr)",
        ]
        for snippet in expected:
            self.assertIn(snippet, header)

    def test_explicit_constructor_order(self):
        header = self.generate(
            "cpu", ctor_order=["bId", "aId", "omega", "alpha"])
        self.assertIn("Toy (walberla::BlockDataID bId, "
                      "walberla::BlockDataID aId, double omega, double alpha)",
                      header)

    def test_scalars_sorted_alphabetically(self):
        a, b = ps.fields("a, b: double[3D]", layout="fzyx")
        zeta, beta = sp.symbols("zeta beta")
        assignments = [ps.Assignment(a.center, zeta * b.center + beta)]
        tmp = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        ctx = sfg_driver.CodegenContext("double", "cpu")
        with sfg_driver.KernelFileGenerator("Toy", ctx, tmp, "test.py") as sfg:
            compat.legacy_sweep(sfg, ctx, "Toy", assignments)
        header = (tmp / "Toy.h").read_text()
        self.assertIn("double beta, double zeta)", header)

    def test_blockforest_property_unsupported(self):
        a, b = ps.fields("a, b: double[3D]", layout="fzyx")
        rhs = b.center + sweepgen.symbolic.DomainCoordinates.x_min()
        assignments = [ps.Assignment(a.center, rhs)]
        tmp = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        ctx = sfg_driver.CodegenContext("double", "cpu")
        with self.assertRaisesRegex(
                NotImplementedError, "not supported by the legacy wrapper"):
            with sfg_driver.KernelFileGenerator(
                    "Toy", ctx, tmp, "test.py") as sfg:
                compat.legacy_sweep(sfg, ctx, "Toy", assignments)

    def test_sparse_sweep_unsupported(self):
        tmp = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        ctx = sfg_driver.CodegenContext("double", "cpu")
        with mock.patch.object(sweepgen.Sweep, "sparse",
                               new_callable=mock.PropertyMock,
                               return_value=True):
            with self.assertRaisesRegex(
                    NotImplementedError, "sparse sweeps are not supported"):
                with sfg_driver.KernelFileGenerator(
                        "Toy", ctx, tmp, "test.py") as sfg:
                    compat.legacy_sweep(sfg, ctx, "Toy",
                                        toy_assignments(ctx))

    def test_wrong_constructor_order(self):
        with self.assertRaisesRegex(ValueError, "constructor order"):
            self.generate("cpu", ctor_order=["aId", "bId", "alpha"])


def swap_assignments(ctx, rng_parameters=True):
    """
    Toy stream step ``f_tmp <- f`` with the parameters of a thermalized
    kernel (``kT`` and, optionally, RNG seed, counter and block offset).
    """
    f, f_tmp = ps.fields(f"f, f_tmp: {ctx.data_type.c_string()}[3D]",
                         layout="fzyx")
    rhs = f[1, 0, 0]
    if rng_parameters:
        ints = [ps.TypedSymbol(name, DynamicType.INDEX_TYPE) for name in (
            "block_offset_0", "block_offset_1", "block_offset_2", "seed",
            "time_step")]
        rhs = sp.Symbol("kT") * rhs + tcast(sum(ints), ctx.data_type)
    return [ps.Assignment(f_tmp.center, rhs)], (f, f_tmp)


class SwappingSweepTest(unittest.TestCase):

    def generate(self, target="cpu", **kwargs):
        tmp = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        ctx = sfg_driver.CodegenContext("double", target)
        assignments, swap = swap_assignments(ctx)
        with sfg_driver.KernelFileGenerator("Swap", ctx, tmp, "test.py") as sfg:
            compat.legacy_sweep(sfg, ctx, "Swap", assignments,
                                swap_fields=[swap], **kwargs)
        header = (tmp / "Swap.h").read_text()
        return header[header.rindex("  class Swap {"):]

    def test_shadow_field_and_cell_intervals_are_hidden(self):
        wrapper = self.generate()
        self.assertIn("Swap (walberla::BlockDataID fId, double kT, "
                      "int64_t seed, int64_t time_step)", wrapper)
        self.assertNotIn("f_tmp", wrapper)
        self.assertNotIn("runOnCellInterval", wrapper)
        self.assertIn("static std::function< void(walberla::IBlock *) > "
                      "getSweep(", wrapper)

    def test_block_offset_is_set_by_configure(self):
        wrapper = self.generate()
        self.assertIn(":impl_(fId, walberla::Vector3< int64_t >(0), kT, "
                      "seed, time_step)", wrapper)
        self.assertIn("const walberla::Cell blockMin = "
                      "blocks->getBlockCellBB(*block).min();", wrapper)
        self.assertIn("impl_.block_offset() = walberla::Vector3< int64_t >"
                      "(blockMin[0], blockMin[1], blockMin[2]);", wrapper)
        self.assertIn("int64_t getBlock_offset_2() const", wrapper)
        self.assertIn("void setBlock_offset_0(int64_t value)", wrapper)

    def test_legacy_types(self):
        wrapper = self.generate(legacy_types={
            "seed": "uint32_t", "time_step": "uint32_t",
            "block_offset": "uint32_t"})
        self.assertIn("Swap (walberla::BlockDataID fId, double kT, "
                      "uint32_t seed, uint32_t time_step)", wrapper)
        self.assertIn("return static_cast< uint32_t >(impl_.time_step());",
                      wrapper)
        self.assertIn("void setTime_step(uint32_t value)", wrapper)
        self.assertIn("uint32_t getBlock_offset_1() const", wrapper)

    def test_legacy_type_of_unknown_scalar(self):
        with self.assertRaisesRegex(ValueError, "unknown scalar omega"):
            self.generate(legacy_types={"omega": "uint32_t"})

    def test_athermal_dispatch(self):
        ctx = sfg_driver.CodegenContext("double", "cpu")
        athermal, _ = swap_assignments(ctx, rng_parameters=False)
        wrapper = self.generate(athermal_assignments=athermal)
        self.assertIn("impl_athermal_(fId)", wrapper)
        self.assertIn("if (impl_.kT() > 0) {", wrapper)
        self.assertIn("impl_athermal_(block);", wrapper)
        self.assertIn("internal::SwapAthermal impl_athermal_;", wrapper)

    def test_athermal_dispatch_on_gpu(self):
        ctx = sfg_driver.CodegenContext("double", "cuda")
        athermal, _ = swap_assignments(ctx, rng_parameters=False)
        wrapper = self.generate("cuda", athermal_assignments=athermal)
        self.assertIn("impl_athermal_(block, "
                      "walberla::v8::sweep::exectag::GPU{stream});", wrapper)

    def test_athermal_dispatch_needs_swaps(self):
        tmp = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        ctx = sfg_driver.CodegenContext("double", "cpu")
        with self.assertRaisesRegex(NotImplementedError, "swapping sweep"):
            with sfg_driver.KernelFileGenerator(
                    "Toy", ctx, tmp, "test.py") as sfg:
                compat.legacy_sweep(sfg, ctx, "Toy", toy_assignments(ctx),
                                    athermal_assignments=toy_assignments(ctx))


if __name__ == "__main__":
    unittest.main()
