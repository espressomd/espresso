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

import functools
import os
import pathlib
import re
import shutil
import tempfile
import unittest

import lbmpy
import lbmpy.boundaries
import pystencils as ps
from pystencils import FieldType
from pystencils.types import PsStructType

import boundaries
import generate_lb_kernels
import pystencils_espresso
import sfg_driver

D3Q19 = lbmpy.LBStencil(lbmpy.Stencil.D3Q19)
D3Q27 = lbmpy.LBStencil(lbmpy.Stencil.D3Q27)
REPO = pathlib.Path(__file__).resolve().parents[3]
BRIDGE_SRC = REPO / "src" / "walberla_bridge" / "src"
BASE_HEADER = BRIDGE_SRC / "utils" / "IndexVectorBoundary.hpp"
GENERATED = BRIDGE_SRC / "lattice_boltzmann" / "generated_kernels"
CPP_TYPES = {"double": "double", "single": "float"}


def generate_uncached(precision, target):
    """Files written for ``DynamicUBB``, as a dict from the file name to the
    text."""
    ctx = sfg_driver.CodegenContext(precision, target)
    name = ctx.class_name("DynamicUBB")
    fields = pystencils_espresso.generate_fields(D3Q19, ctx.data_type)
    method = generate_lb_kernels.create_lb_method(fields)
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        with sfg_driver.KernelFileGenerator(
                name, ctx, tmp, "generate_lb_kernels.py") as sfg:
            boundaries.generate_dynamic_ubb(
                sfg, ctx, name, method, fields["pdfs"])
        return {p.name: p.read_text() for p in sorted(tmp.iterdir())}
    finally:
        shutil.rmtree(tmp)


generate = functools.cache(generate_uncached)


def struct_members(header, name):
    """(type, name) of the members of the struct ``name``, in order."""
    body = re.search(r"struct %s \{(.*?)\n\s*%s ?\(" % (name, name),
                     header, re.S).group(1)
    return re.findall(r"(\w+) (\w+);", body)


class StructTypesTest(unittest.TestCase):

    def test_index_struct(self):
        data = [("vel_0", ps.create_type("float32"))]
        struct = boundaries.index_struct_type("ns::Name", data)
        self.assertIsInstance(struct, PsStructType)
        self.assertEqual(struct.name, "ns::Name")
        self.assertEqual(
            [(m.name, str(m.dtype)) for m in struct.members],
            [("x", "int32"), ("y", "int32"), ("z", "int32"),
             ("dir", "int32"), ("vel_0", "float32")])

    def test_index_coordinates_are_32_bit(self):
        # the kernels index with int32_t, see the module docstring
        self.assertEqual(str(boundaries.INDEX_TYPE), "int32")

    def test_force_struct(self):
        struct = boundaries.force_struct_type(
            "ns::Name", ps.create_type("float64"), 2)
        self.assertEqual(
            [(m.name, str(m.dtype)) for m in struct.members],
            [("F_0", "float64"), ("F_1", "float64")])

    def test_link_fields(self):
        index_field, force_field = boundaries.link_fields(
            "Name", ps.create_type("float64"), 3,
            [("vel_0", ps.create_type("float64"))])
        for field in (index_field, force_field):
            self.assertEqual(field.field_type, FieldType.INDEXED)
            self.assertEqual(field.spatial_dimensions, 1)
        self.assertEqual(index_field.name, "indexVector")
        self.assertEqual(force_field.name, "forceVector")
        self.assertEqual(index_field.dtype.name,
                         "walberla::lbm::NameIndexInfo")
        self.assertEqual(force_field.dtype.name,
                         "walberla::lbm::NameForceStruct")

    def test_offsets_table(self):
        table = boundaries._offsets_table([(0, 0), (1, 0), (-1, 1)])
        self.assertEqual(table, "{0, 1, -1},\n{0, 0, 1}")

    def test_array_extraction(self):
        extraction = boundaries.ArrayExtraction("ptr", "size")
        self.assertEqual(str(extraction._extract_ptr()), "ptr")
        self.assertEqual(str(extraction._extract_size(0)), "size")
        self.assertEqual(str(extraction._extract_stride(0)), "1")
        self.assertIsNone(extraction._extract_size(1))
        self.assertIsNone(extraction._extract_stride(1))


class LinkAssignmentsTest(unittest.TestCase):

    @staticmethod
    def assignments(stencil=D3Q19, boundary=None):
        data_type = ps.create_type("float64")
        fields = pystencils_espresso.generate_fields(stencil, data_type)
        method = lbmpy.creationfunctions.create_lb_method(
            lbm_config=lbmpy.LBMConfig(
                stencil=stencil, method=lbmpy.Method.SRT,
                relaxation_rate=1.0))
        if boundary is None:
            boundary = lbmpy.boundaries.UBB(
                lambda *args: None, dim=3, data_type=data_type,
                calculate_force_on_boundary=True, name="test")
        index_field, force_field = boundaries.link_fields(
            "Test", data_type, 3, boundary.additional_data)
        elements = boundaries.lbm_link_assignments(
            boundary, method, fields["pdfs"], index_field, force_field)
        return elements, fields["pdfs"], index_field, force_field

    @staticmethod
    def accesses(expression, field):
        return {a.field.name + (str(a.index) if a.index else "")
                for a in expression.atoms(ps.Field.Access)
                if a.field.name == field.name}

    def test_force_and_populations_are_written(self):
        elements, pdfs, index_field, force_field = self.assignments()
        written = [e.lhs for e in elements
                   if isinstance(e.lhs, ps.Field.Access)]
        self.assertEqual(
            sorted(str(a.index[0]) for a in written
                   if a.field.name == "forceVector"),
            ["F_0", "F_1", "F_2"])
        # one population is written: the inverse of the link direction
        self.assertEqual(
            len([a for a in written if a.field.name == pdfs.name]), 1)
        self.assertFalse(
            [a for a in written if a.field.name == index_field.name])

    def test_direction_and_velocity_come_from_the_index_field(self):
        elements, _, index_field, _ = self.assignments()
        read = set()
        for element in elements:
            for access in element.rhs.atoms(ps.Field.Access):
                if access.field.name == index_field.name:
                    read.add(str(access.index[0]))
        self.assertEqual(read, {"dir", "vel_0", "vel_1", "vel_2"})
        # the symbol of the direction is defined by the index field
        definitions = [
            e for e in elements if isinstance(e.rhs, ps.Field.Access)
            and e.rhs.field.name == index_field.name]
        self.assertEqual([str(e.rhs.index[0]) for e in definitions], ["dir"])

    def test_other_boundary_conditions(self):
        # nothing is specific to the UBB: no additional data, other stencil
        elements, _, index_field, _ = self.assignments(
            D3Q27, lbmpy.boundaries.NoSlip(calculate_force_on_boundary=True))
        self.assertEqual([m.name for m in index_field.dtype.members],
                         ["x", "y", "z", "dir"])
        self.assertTrue(elements)


class GeneratedBoundaryTest(unittest.TestCase):

    def check_files(self, precision, target):
        ctx = sfg_driver.CodegenContext(precision, target)
        name = ctx.class_name("DynamicUBB")
        files = generate(precision, target)
        impl = "cu" if target == "cuda" else "cpp"
        self.assertEqual(sorted(files), [f"{name}.{impl}", f"{name}.h"])
        return name, files[f"{name}.h"], files[f"{name}.{impl}"]

    def test_index_info_layout(self):
        for precision, cpp in CPP_TYPES.items():
            for target in ("cpu", "cuda"):
                name, header, _ = self.check_files(precision, target)
                self.assertEqual(
                    struct_members(header, f"{name}IndexInfo"),
                    [("int32_t", "x"), ("int32_t", "y"), ("int32_t", "z"),
                     ("int32_t", "dir"), (cpp, "vel_0"), (cpp, "vel_1"),
                     (cpp, "vel_2")])
                self.assertEqual(
                    struct_members(header, f"{name}ForceStruct"),
                    [(cpp, "F_0"), (cpp, "F_1"), (cpp, "F_2")])

    def test_index_info_constructor(self):
        _, header, _ = self.check_files("double", "cpu")
        self.assertRegex(
            header,
            r"DynamicUBBDoublePrecisionIndexInfo ?\(int32_t x_, int32_t y_, "
            r"int32_t z_, int32_t dir_\)\s*:\s*x\(x_\),\s*y\(y_\),\s*z\(z_\),"
            r"\s*dir\(dir_\),\s*vel_0\(double\(0\)\)")
        self.assertIn("walberla::floatIsEqual(vel_2, o.vel_2)", header)
        self.assertIn("dir == o.dir", header)

    def test_class(self):
        for precision, cpp in CPP_TYPES.items():
            for target, memory in (("cpu", "Host"), ("cuda", "Device")):
                name, header, _ = self.check_files(precision, target)
                base = (
                    f"walberla::IndexVectorBoundary< {name}, {name}IndexInfo, "
                    f"{name}ForceStruct, walberla::stencil::D3Q19, "
                    f"walberla::Vector3< {cpp} >, "
                    f"walberla::IndexVectorMemory::{memory} >")
                self.assertIn(f"class {name} : public {base}", header)
                self.assertIn(f"using Base = {base};", header)
                self.assertIn(f"using IndexInfo = {name}IndexInfo;", header)
                self.assertIn(f"using ForceStruct = {name}ForceStruct;",
                              header)
                self.assertIn(f"{name} (const std::shared_ptr< "
                              "walberla::StructuredBlockForest >& blocks, "
                              "walberla::BlockDataID pdfsID_, const "
                              f"std::function< walberla::Vector3< {cpp} >("
                              "const walberla::Cell &, const "
                              "walberla::shared_ptr< "
                              "walberla::StructuredBlockForest > &, "
                              "walberla::IBlock &) >& valueCallback)",
                              header)
                self.assertIn(f':Base(blocks, pdfsID_, valueCallback, '
                              f'"{name}")', header)

    def test_pdf_field_type(self):
        _, cpu, _ = self.check_files("single", "cpu")
        self.assertIn("using PdfField_T = walberla::field::GhostLayerField< "
                      "float, 19 >;", cpu)
        _, cuda, _ = self.check_files("single", "cuda")
        self.assertIn("using PdfField_T = walberla::gpu::GPUField< float >;",
                      cuda)

    def test_assign_value(self):
        _, header, _ = self.check_files("double", "cpu")
        self.assertRegex(
            header,
            r"static void assignValue\(IndexInfo &e, walberla::Vector3< "
            r"double > const &v\) \{\s*e\.vel_0 = v\[0\];\s*e\.vel_1 = v\[1\];"
            r"\s*e\.vel_2 = v\[2\];\s*\}")

    def test_link_kernel_signature(self):
        _, cpu, _ = self.check_files("double", "cpu")
        self.assertIn(
            "void linkKernel(walberla::field::GhostLayerField< double, 19 > "
            "* pdfs, DynamicUBBDoublePrecisionIndexInfo * indexVector_, "
            "DynamicUBBDoublePrecisionForceStruct * forceVector_, "
            "std::size_t size);", cpu)
        _, cuda, _ = self.check_files("double", "cuda")
        self.assertIn(
            "void linkKernel(walberla::gpu::GPUField< double > * pdfs, "
            "DynamicUBBDoublePrecisionCUDAIndexInfo * indexVector_, "
            "DynamicUBBDoublePrecisionCUDAForceStruct * forceVector_, "
            "std::size_t size, cudaStream_t stream);", cuda)

    def test_neighbor_offsets(self):
        # the table of the stencil, independent of the implementation
        _, header, _ = self.check_files("double", "cpu")
        match = re.search(
            r"neighborOffset = \{\{(.*?)\}\};", header, re.S)
        rows = [list(map(int, row.split(","))) for row in
                re.findall(r"\{([-\d, ]+)\}", match.group(1))]
        self.assertEqual(
            rows, [[int(v[i]) for v in D3Q19] for i in range(3)])
        self.assertIn("std::array< std::array< int, 19u >, 3u > "
                      "neighborOffset", header)

    def test_base_header_is_included_by_relative_path(self):
        for target in ("cpu", "cuda"):
            _, header, _ = self.check_files("double", target)
            include = re.search(
                r'#include "([^"]*IndexVectorBoundary[^"]*)"', header).group(1)
            self.assertEqual(
                os.path.normpath(GENERATED / include), str(BASE_HEADER))
        self.assertTrue(BASE_HEADER.is_file())

    def test_base_header_provides_what_the_class_uses(self):
        text = BASE_HEADER.read_text()
        self.assertIn("enum class IndexVectorMemory { Host, Device };", text)
        self.assertIn("class IndexVectorBoundary {", text)
        for member in ("assignValue", "linkKernel", "PdfField_T"):
            self.assertIn(member, text)

    def test_cpu_kernel(self):
        name, _, source = self.check_files("double", "cpu")
        kernels = f"namespace {name}_kernels"
        self.assertIn(kernels, source)
        self.assertIn(f"void linkKernel (walberla::lbm::{name}ForceStruct * "
                      f"RESTRICT const _data_forceVector, walberla::lbm::"
                      f"{name}IndexInfo * RESTRICT const _data_indexVector, ",
                      source)
        self.assertIn("#pragma omp parallel", source)
        self.assertIn("#pragma omp for schedule(static)", source)
        self.assertNotIn("__global__", source)
        # the elements of the vectors are accessed as structs
        self.assertIn("_data_indexVector[sparse_idx * _stride_indexVector_0]"
                      ".vel_0", source)
        self.assertIn("_data_forceVector[sparse_idx * _stride_forceVector_0]"
                      ".F_2 = ", source)
        self.assertIn("pdfs->dataAt(0, 0, 0, 0)", source)
        self.assertIn(f"walberla::lbm::{name}_kernels::linkKernel(", source)

    def test_cuda_kernel(self):
        name, _, source = self.check_files("double", "cuda")
        self.assertIn(f"namespace {name}_kernels", source)
        self.assertIn("__global__ void linkKernel", source)
        self.assertNotIn("#pragma omp", source)
        # one thread for each link, guarded by the number of links
        self.assertIn("const int32_t sparse_idx = (int32_t) blockIdx.x * "
                      "(int32_t) blockDim.x + (int32_t) threadIdx.x;", source)
        self.assertIn("if(sparse_idx < _size_indexVector_0)", source)
        self.assertIn("const dim3 __block_size { 256, 1, 1 };", source)
        self.assertRegex(
            source, r"linkKernel<<< __grid_size, __block_size, 0, stream >>>")

    def test_kernels_use_32_bit_indices(self):
        for target in ("cpu", "cuda"):
            _, _, source = self.check_files("double", target)
            self.assertIn("const int32_t _stride_pdfs_0 {", source)
            self.assertIn("const int32_t _size_indexVector_0 {", source)
            self.assertNotIn("int64_t", source)

    def test_precision(self):
        for target in ("cpu", "cuda"):
            _, _, source = self.check_files("single", target)
            self.assertIn("float * RESTRICT const _data_pdfs", source)
            self.assertNotIn("double", source)

    def test_populations_and_forces_of_the_kernel(self):
        # the boundary conditions is the UBB with the force of the link:
        # f(inverse direction, boundary cell) = f(direction) - 6 w rho c.u
        # and F = (2 f(direction) - 6 w rho c.u) c
        _, _, source = self.check_files("double", "cpu")
        self.assertIn("const double rho = 1.0 + delta_rho;", source)
        self.assertIn("6.0 * _data_indexVector[sparse_idx * "
                      "_stride_indexVector_0].vel_0 * (double) "
                      "neighbour_offset_x[dir]", source)
        self.assertIn("weights[19] = { __c_1_0o3_0, __c_1_0o18_0", source)
        self.assertIn("const double f = 2.0 * _data_pdfs[", source)
        self.assertIn("_data_forceVector[sparse_idx * _stride_forceVector_0]"
                      ".F_0 = f * (double) neighbour_offset_x[dir];", source)
        self.assertIn("f_in_inv_dir_idx[19] = { 0, 2, 1, 4, 3, 6, 5, 10, 9, "
                      "8, 7, 16, 15, 18, 17, 12, 11, 14, 13 }", source)


class GenericBoundaryTest(unittest.TestCase):
    """`generate_link_boundary` is independent of the UBB and of D3Q19."""

    def test_no_slip_d3q27(self):
        ctx = sfg_driver.CodegenContext("double", "cpu")
        data_type = ctx.data_type
        fields = pystencils_espresso.generate_fields(D3Q27, data_type)
        method = lbmpy.creationfunctions.create_lb_method(
            lbm_config=lbmpy.LBMConfig(
                stencil=D3Q27, method=lbmpy.Method.SRT,
                relaxation_rate=1.0))
        boundary = lbmpy.boundaries.NoSlip(calculate_force_on_boundary=True)
        index_field, force_field = boundaries.link_fields(
            "NoSlipTest", data_type, 3, boundary.additional_data)
        assignments = boundaries.lbm_link_assignments(
            boundary, method, fields["pdfs"], index_field, force_field)
        stencil = [tuple(int(c) for c in v) for v in D3Q27]
        tmp = pathlib.Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        with sfg_driver.KernelFileGenerator(
                "NoSlipTest", ctx, tmp, "test") as sfg:
            boundaries.generate_link_boundary(
                sfg, ctx, "NoSlipTest", assignments=assignments,
                index_field=index_field, force_field=force_field,
                pdf_field=fields["pdfs"], stencil=stencil,
                value=boundaries.LinkValue("double", ()))
        header = (tmp / "NoSlipTest.h").read_text()
        source = (tmp / "NoSlipTest.cpp").read_text()
        self.assertEqual(struct_members(header, "NoSlipTestIndexInfo"),
                         [("int32_t", n) for n in ("x", "y", "z", "dir")])
        self.assertIn("walberla::stencil::D3Q27", header)
        self.assertIn("walberla::field::GhostLayerField< double, 27 >",
                      header)
        self.assertIn("std::array< std::array< int, 27u >, 3u >", header)
        self.assertRegex(header, r"static void assignValue\(IndexInfo &e, "
                         r"double const &v\) \{\s*\}")
        self.assertIn("namespace NoSlipTest_kernels", source)
        self.assertIn("#pragma omp parallel", source)


if __name__ == "__main__":
    unittest.main()
