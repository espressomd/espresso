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

"""
Tests of the generated reaction kernels of the electrokinetics: the
assignments, the structure of the generated classes for 1 to 5 reactants
(bulk and indexed, CPU and CUDA, both precisions), and the symbols of their
kernels, which have to be unique because the CPU and CUDA classes are linked
into the same library.
"""

import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

import pystencils as ps
import sympy as sp

import generate_ek_kernels
import reactions
import sfg_driver
from tests.test_accessor_symbols import definitions

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "generate_ek_kernels.py"
NUMBERS = range(1, 6)
PRECISIONS = (("Double", "double"), ("Single", "float"))


def class_name(kind, num, precision="Double", target=""):
    return f"ReactionKernel{kind}{num}{precision}Precision{target}"


class ReactionAssignmentsTest(unittest.TestCase):

    dtype = ps.create_type("float64")

    def inlined(self, num, indexed):
        collection = reactions.reaction_assignments(self.dtype, num, indexed)
        return collection, collection.new_without_subexpressions()

    def test_reaction_equation(self):
        # rho_i += stoech_i * k * prod_j rho_j ** order_j, with the densities
        # of before the reaction
        for indexed in (False, True):
            for num in NUMBERS:
                _, inlined = self.inlined(num, indexed)
                lhs = [a.lhs for a in inlined.main_assignments]
                symbols = {s.name: s for a in inlined.main_assignments
                           for s in a.rhs.atoms(ps.TypedSymbol)}
                rate = symbols["rate_coefficient"]
                for i in range(num):
                    rate *= sp.Pow(lhs[i], symbols[f"order{i}"])
                for i, assignment in enumerate(inlined.main_assignments):
                    expected = lhs[i] + symbols[f"stoech{i}"] * rate
                    with self.subTest(indexed=indexed, num=num, i=i):
                        self.assertEqual(
                            sp.simplify(assignment.rhs - expected), 0)

    def test_densities_are_read_once_before_they_are_written(self):
        for num in NUMBERS:
            collection = reactions.reaction_assignments(self.dtype, num)
            reads = [a.rhs for a in collection.subexpressions[:num]]
            self.assertEqual(
                [r.field.name for r in reads],
                [f"rho_{i}" for i in range(num)])
            for assignment in collection.main_assignments:
                self.assertEqual(
                    assignment.rhs.atoms(ps.Field.Access), set())

    def test_scalars_have_the_data_type_of_the_kernel(self):
        for dtype in ("float32", "float64"):
            collection = reactions.reaction_assignments(
                ps.create_type(dtype), 3)
            # the scalars are not read from fields, which have typed symbols
            # too
            scalars = {s.name: s.dtype for a in (
                *collection.main_assignments, collection.subexpressions[-1])
                for s in a.rhs.atoms(ps.TypedSymbol)}
            self.assertEqual(sorted(scalars), sorted(
                ["rate_coefficient"] + [f"{n}{i}" for n in ("order", "stoech")
                                        for i in range(3)]))
            self.assertEqual(
                {t.numpy_dtype.name for t in scalars.values()}, {dtype})

    def test_bulk_accesses_the_cell_itself(self):
        collection = reactions.reaction_assignments(self.dtype, 2)
        for assignment in collection.main_assignments:
            self.assertFalse(assignment.lhs.is_absolute_access)
            self.assertEqual(assignment.lhs.offsets, (0, 0, 0))
            self.assertEqual(assignment.lhs.field.field_type,
                             ps.FieldType.GENERIC)

    def test_indexed_accesses_the_cells_of_the_index_list(self):
        collection = reactions.reaction_assignments(self.dtype, 2, True)
        index_list = sfg_driver.index_list_field()
        for assignment in collection.main_assignments:
            self.assertTrue(assignment.lhs.is_absolute_access)
            self.assertEqual(assignment.lhs.field.field_type,
                             ps.FieldType.CUSTOM)
            self.assertEqual(
                {o.field for o in assignment.lhs.offsets}, {index_list})
            self.assertEqual(
                [o.index[0] for o in assignment.lhs.offsets],
                ["x", "y", "z"])

    def test_number_of_reactants(self):
        for num in (0, 6, -1):
            with self.assertRaisesRegex(ValueError, "not implemented"):
                reactions.reaction_assignments(self.dtype, num)


class GeneratedReactionKernelsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.directory = pathlib.Path(tempfile.mkdtemp())
        names = [f"reaction_{kind}_{num}" for kind in ("bulk", "indexed")
                 for num in NUMBERS]
        subprocess.run(
            [sys.executable, str(SCRIPT), "--kernels", *names,
             "--target", "cpu", "cuda", "--output-dir", str(cls.directory)],
            capture_output=True, text=True, check=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.directory)

    def read(self, kind, num, precision="Double", target="", extension="h"):
        name = class_name(kind, num, precision, target)
        return (self.directory / f"{name}.{extension}").read_text()

    def constructor(self, kind, num, precision, target, c_type):
        density_ids = ", ".join(
            f"walberla::BlockDataID rho_{i}Id" for i in range(num))
        scalars = ", ".join(
            [f"{c_type} order{i}" for i in range(num)]
            + [f"{c_type} rate_coefficient"]
            + [f"{c_type} stoech{i}" for i in range(num)])
        memory = "unified" if target == "CUDA" else "host"
        index_list = (
            "walberla::v8::sweep::SparseIndexList< "
            "walberla::v8::sweep::CellIdx, "
            f"walberla::v8::memory::memtag::{memory} > indexList, "
            if kind == "Indexed" else "")
        return (f"{class_name(kind, num, precision, target)} "
                f"({index_list}{density_ids}, {scalars})")

    def test_files(self):
        suffixes = (".h", ".cpp", "CUDA.h", "CUDA.cu")
        expected = sorted(
            f"ReactionKernel{kind}{num}{p}Precision{s}"
            for kind in ("Bulk", "Indexed") for num in NUMBERS
            for p in ("Single", "Double") for s in suffixes)
        self.assertEqual(
            sorted(p.name for p in self.directory.iterdir()), expected)

    def test_constructors(self):
        # the same constructor for all numbers of reactants: the index list
        # (indexed), the density fields, then the scalars in alphabetical
        # order. The scalars of 2 and 3 reactants must not be combined to a
        # Vector2 or Vector3
        for kind in ("Bulk", "Indexed"):
            for num in NUMBERS:
                for precision, c_type in PRECISIONS:
                    for target in ("", "CUDA"):
                        header = self.read(kind, num, precision, target)
                        constructor = self.constructor(
                            kind, num, precision, target, c_type)
                        with self.subTest(kind=kind, num=num,
                                          precision=precision, target=target):
                            self.assertIn(f"    {constructor}\n", header)
                            self.assertNotIn("Vector", header)

    def test_accessors_of_the_scalars(self):
        for num in NUMBERS:
            header = self.read("Bulk", num)
            for i in range(num):
                for name in (f"order{i}", f"stoech{i}"):
                    self.assertIn(f"double& {name}()\n", header)
                    self.assertIn(f"const double& {name}() const", header)
            self.assertIn("double& rate_coefficient()\n", header)
            self.assertIn(f"const walberla::BlockDataID& rho_{num - 1}Id()",
                          header)

    def test_run_methods(self):
        for num in NUMBERS:
            for target in ("", "CUDA"):
                bulk = self.read("Bulk", num, "Single", target)
                indexed = self.read("Indexed", num, "Single", target)
                self.assertIn(
                    "void operator()(walberla::IBlock * block", bulk)
                self.assertIn("runOnCellInterval(", bulk)
                self.assertIn("void operator()(walberla::IBlock * block",
                              indexed)
                self.assertNotIn("runOnCellInterval", indexed)
                self.assertNotIn("SparseIndexList", bulk)

    def test_memory_of_the_index_list(self):
        for num in NUMBERS:
            cpu = self.read("Indexed", num)
            cuda = self.read("Indexed", num, target="CUDA")
            self.assertIn("memtag::host", cpu)
            self.assertNotIn("unified", cpu)
            self.assertIn("memtag::unified", cuda)
            self.assertNotIn("memtag::host", cuda)

    def test_kernels(self):
        for num in NUMBERS:
            for kind in ("Bulk", "Indexed"):
                single = self.read(kind, num, "Single", extension="cpp")
                double = self.read(kind, num, "Double", extension="cpp")
                cuda = self.read(kind, num, "Single", "CUDA", "cu")
                self.assertEqual(single.count("std::pow(local_rho_"), num)
                self.assertEqual(double.count("std::pow(local_rho_"), num)
                self.assertNotRegex(single, r"(?<![:\w])pow\(")
                self.assertIn("#pragma omp", single)
                self.assertIn("__global__", cuda)
                self.assertEqual(cuda.count("powf(local_rho_"), num)
                for i in range(num):
                    self.assertIn(f"local_rho_{i} + rate_factor * stoech{i}",
                                  double)
            # the bulk sweep loops over cells, the indexed sweep over the list
            self.assertNotIn(
                "sparse_idx", self.read("Bulk", num, extension="cpp"))
            self.assertIn(
                "sparse_idx", self.read("Indexed", num, extension="cpp"))

    def test_header_of_the_generator(self):
        header = self.read("Indexed", 3, "Single", "CUDA")
        self.assertIn("Generated by "
                      "maintainer/walberla_kernels/generate_ek_kernels.py",
                      header)
        self.assertIn("namespace walberla::pystencils {", header)


class ReactionRegistryTest(unittest.TestCase):

    def test_registered_kernels(self):
        kernels = generate_ek_kernels.KERNELS
        reaction = {k: v for k, v in kernels.items()
                    if k.startswith("reaction_")}
        self.assertEqual(sorted(reaction), sorted(
            f"reaction_{kind}_{num}" for kind in ("bulk", "indexed")
            for num in NUMBERS))
        for key, kernel in reaction.items():
            _, kind, num = key.split("_")
            self.assertEqual(kernel.stem,
                             f"ReactionKernel{kind.capitalize()}{num}")
            self.assertEqual(kernel.targets, ("cpu", "cuda"))
            self.assertFalse(kernel.header_only)
            self.assertIsInstance(kernel, sfg_driver.Kernel)

    def test_class_names(self):
        ctx = sfg_driver.CodegenContext("single", "cuda")
        kernel = generate_ek_kernels.KERNELS["reaction_indexed_4"]
        self.assertEqual(ctx.class_name(kernel.stem),
                         "ReactionKernelIndexed4SinglePrecisionCUDA")


class ReactionSymbolsTest(unittest.TestCase):
    """
    The CPU and CUDA kernels of the reactions have the same kinds of
    parameters, so a kernel of a CPU file and a kernel of a CUDA file must
    never have the same qualified name and parameter types (see
    ``test_boundary_symbols``). The kernels of a class are in the namespace
    ``<ClassName>_kernels`` of ``walberla::pystencils``.
    """

    @classmethod
    def setUpClass(cls):
        cls.directory = pathlib.Path(tempfile.mkdtemp())
        kernels = ["dirichlet"] + [
            f"reaction_{kind}_{num}" for kind in ("bulk", "indexed")
            for num in NUMBERS]
        subprocess.run(
            [sys.executable, str(SCRIPT), "--kernels", *kernels,
             "--target", "cpu", "cuda", "--output-dir", str(cls.directory)],
            capture_output=True, text=True, check=True)
        cls.cpu, cls.cuda = {}, {}
        for kind in ("Bulk", "Indexed"):
            for num in NUMBERS:
                for precision, _ in PRECISIONS:
                    for target, functions in (("", cls.cpu),
                                              ("CUDA", cls.cuda)):
                        name = class_name(kind, num, precision, target)
                        suffix = ".cu" if target else ".cpp"
                        functions[name] = definitions("\n".join(
                            (cls.directory / f"{name}{s}").read_text()
                            for s in (".h", suffix)))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.directory)

    def test_each_class_defines_one_kernel_in_its_namespace(self):
        for functions in (self.cpu, self.cuda):
            self.assertEqual(len(functions), 20)
            for name, defined in functions.items():
                kernels = {n for n, _ in defined if "_kernels::" in n}
                self.assertEqual(
                    kernels, {f"walberla::pystencils::{name}_kernels::{name}"},
                    name)

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
        self.assertEqual(len(names), 40)
        self.assertEqual(len(names), len(set(names)))

    def test_no_function_is_shared_with_other_kernels(self):
        # the Dirichlet boundary of the EK is in the same library
        for target, functions in ((".cpp", self.cpu), (".cu", self.cuda)):
            reactions_ = set().union(*functions.values())
            for precision, _ in PRECISIONS:
                name = f"Dirichlet{precision}Precision" + (
                    "CUDA" if target == ".cu" else "")
                dirichlet = definitions("\n".join(
                    (self.directory / f"{name}{s}").read_text()
                    for s in (".h", target)))
                self.assertTrue(dirichlet)
                self.assertEqual(reactions_ & dirichlet, set())


if __name__ == "__main__":
    unittest.main()
