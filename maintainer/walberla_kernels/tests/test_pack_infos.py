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
import warnings

import lbmpy
import numpy as np
import pystencils as ps

import pack_infos
import sfg_driver

D3Q19 = lbmpy.LBStencil(lbmpy.Stencil.D3Q19)
D3Q27 = lbmpy.LBStencil(lbmpy.Stencil.D3Q27)
REPO = pathlib.Path(__file__).resolve().parents[3]
BRIDGE_SRC = REPO / "src" / "walberla_bridge" / "src"
CPU_BASE_HEADER = BRIDGE_SRC / "utils" / "UniformCpuFieldPackInfoBase.hpp"
FACES = {"W", "E", "S", "N", "B", "T"}


def generate_uncached(kind, precision, target):
    """
    Files written for a pack info (``kind`` is ``pdf`` or ``vec``), as a
    dict from the file name to the text.
    """
    stem, generator = {
        "pdf": ("PackInfoPdf", pack_infos.generate_pdf_pack_info),
        "vec": ("PackInfoVec", pack_infos.generate_vector_pack_info),
    }[kind]
    ctx = sfg_driver.CodegenContext(precision, target)
    name = ctx.class_name(stem)
    tmp = pathlib.Path(tempfile.mkdtemp())
    try:
        with sfg_driver.KernelFileGenerator(
                name, ctx, tmp, "generate_lb_kernels.py") as sfg:
            generator(sfg, ctx, name, D3Q19)
        return {p.name: p.read_text() for p in sorted(tmp.iterdir())}
    finally:
        shutil.rmtree(tmp)


generate = functools.cache(generate_uncached)


def elements_per_cell(header):
    """Number of components per direction from ``elementsPerCell()``."""
    return {m.group(1): int(m.group(2)) for m in re.finditer(
        r"case walberla::stencil::Direction::(\w+): \{\s*return (\d+);",
        header)}


class CommunicationSetsTest(unittest.TestCase):

    def test_directions(self):
        directions = pack_infos.COMMUNICATION_DIRECTIONS
        self.assertEqual(len(directions), 26)
        self.assertEqual(len(set(directions)), 26)
        self.assertNotIn((0, 0, 0), directions)
        self.assertEqual(directions[0], (-1, -1, -1))
        self.assertEqual(directions[-1], (1, 1, 1))

    def test_pdf_sets_follow_pull_streaming(self):
        # a population leaves the block through the faces, edges and
        # corners that are contained in its own direction
        for stencil in (D3Q19, D3Q27):
            sets = pack_infos.pdf_streaming_sets(stencil)
            for direction in pack_infos.COMMUNICATION_DIRECTIONS:
                expected = tuple(
                    i for i, e in enumerate(stencil)
                    if all(d in (0, v) for d, v in zip(direction, e)))
                self.assertEqual(sets.get(direction, ()), expected,
                                 (len(stencil), direction))

    def test_pdf_sets_are_physical(self):
        # independent of the implementation: the neighbouring block in
        # direction d has cells at coordinate d along the axes of d. Under
        # pull streaming, a cell x pulls population e from x - e, so the
        # population needs data of this block (coordinate 0 along the axes
        # of d) if x - e is in this block there. Such populations move
        # towards the neighbour, e.d > 0.
        for stencil in (D3Q19, D3Q27):
            sets = pack_infos.pdf_streaming_sets(stencil)
            for d in pack_infos.COMMUNICATION_DIRECTIONS:
                needed = {
                    i for i, e in enumerate(stencil)
                    if all(dc - ec == 0 for dc, ec in zip(d, e) if dc != 0)}
                self.assertEqual(set(sets.get(d, ())), needed, d)
                for i in needed:
                    self.assertGreater(
                        sum(a * b for a, b in zip(stencil[i], d)), 0, d)

    def test_pdf_sets_d3q19(self):
        sets = pack_infos.pdf_streaming_sets(D3Q19)
        self.assertEqual(len(sets), 18)
        for direction, indices in sets.items():
            n_axes = sum(map(abs, direction))
            self.assertIn(n_axes, (1, 2))
            self.assertEqual(len(indices), 5 if n_axes == 1 else 1)
            self.assertEqual(list(indices), sorted(indices))
        self.assertEqual(sum(map(len, sets.values())), 6 * 5 + 12)
        # no corner and no rest population
        self.assertNotIn((1, 1, 1), sets)
        rest = [i for i, e in enumerate(D3Q19) if e == (0, 0, 0)]
        self.assertEqual(len(rest), 1)
        self.assertFalse(any(rest[0] in i for i in sets.values()))
        # the populations of the face in x direction move in +x direction
        for i in sets[(1, 0, 0)]:
            self.assertEqual(D3Q19[i][0], 1)

    def test_pdf_sets_d3q27(self):
        sets = pack_infos.pdf_streaming_sets(D3Q27)
        self.assertEqual(len(sets), 26)
        self.assertEqual(len(sets[(1, 1, 1)]), 1)
        self.assertEqual(sum(map(len, sets.values())), 6 * 9 + 12 * 3 + 8)

    def test_vector_sets(self):
        for stencil, n_directions in ((D3Q19, 18), (D3Q27, 26)):
            sets = pack_infos.vector_sets(stencil, 3)
            self.assertEqual(len(sets), n_directions)
            self.assertEqual(set(sets.values()), {(0, 1, 2)})
            self.assertEqual(set(sets), set(stencil) - {(0, 0, 0)})


class GenerationTest(unittest.TestCase):

    def test_cpu_pdf_files(self):
        files = generate("pdf", "double", "cpu")
        self.assertEqual(sorted(files), ["PackInfoPdfDoublePrecision.cpp",
                                         "PackInfoPdfDoublePrecision.h"])
        header = files["PackInfoPdfDoublePrecision.h"]
        self.assertIn("namespace walberla::pystencils {", header)
        self.assertIn(
            "class PackInfoPdfDoublePrecision : public "
            "walberla::UniformCpuFieldPackInfoBase< "
            "PackInfoPdfDoublePrecision >", header)
        self.assertIn("using Base::Base;", header)
        self.assertIn(
            "using Field_T = walberla::field::GhostLayerField< double, 19 >;",
            header)
        self.assertIn('#include "../../utils/UniformCpuFieldPackInfoBase.hpp"',
                      header)
        self.assertNotIn("UniformGpuFieldPackInfoBase", header)

    def test_cpu_pdf_methods_are_const(self):
        header = generate("pdf", "double", "cpu")[
            "PackInfoPdfDoublePrecision.h"]
        for method in ("doPack", "doUnpack", "doLocalCopy"):
            declaration = re.search(rf"void {method}\([^;]*;", header)
            self.assertIsNotNone(declaration, method)
            self.assertTrue(declaration.group().endswith(") const;"), method)
        self.assertNotIn("cudaStream_t", header)

    def test_pdf_elements_per_cell(self):
        header = generate("pdf", "double", "cpu")[
            "PackInfoPdfDoublePrecision.h"]
        counts = elements_per_cell(header)
        self.assertEqual(len(counts), 18)
        self.assertEqual({d for d, n in counts.items() if n == 5}, FACES)
        self.assertEqual({n for d, n in counts.items() if d not in FACES},
                         {1})

    def test_cpu_pdf_kernels(self):
        source = generate("pdf", "double", "cpu")[
            "PackInfoPdfDoublePrecision.cpp"]
        self.assertIn("namespace PackInfoPdfDoublePrecision_kernels {",
                      source)
        # every PDF direction has a distinct set of populations, so there
        # are 18 kernels per type; the deduplication of kernels is exercised
        # by the vector pack info
        for kernel in ("pack", "unpack", "localCopy"):
            self.assertEqual(len(re.findall(rf"\n    void {kernel}\w+ \(",
                                            source)), 18, kernel)
        self.assertNotIn("omp", source)

    def test_vec_shares_one_kernel_between_directions(self):
        files = generate("vec", "single", "cpu")
        self.assertEqual(sorted(files), ["PackInfoVecSinglePrecision.cpp",
                                         "PackInfoVecSinglePrecision.h"])
        header = files["PackInfoVecSinglePrecision.h"]
        self.assertIn("walberla::field::GhostLayerField< float, 3 >", header)
        counts = elements_per_cell(header)
        self.assertEqual(len(counts), 18)
        self.assertEqual(set(counts.values()), {3})
        source = files["PackInfoVecSinglePrecision.cpp"]
        for kernel in ("pack", "unpack", "localCopy"):
            self.assertEqual(len(re.findall(rf"\n    void {kernel}\w+ \(",
                                            source)), 1, kernel)
        # the single kernel of all directions is called in all cases
        self.assertEqual(source.count("PackInfoVecSinglePrecision_kernels::"
                                      "pack"), 18)

    def test_buffer_is_x_fastest(self):
        # the regular expressions and strings below depend on the format of
        # the pystencils output (``int64_t( ... )``, ``_data_buffer[...]``);
        # they are the only guard of the x-fastest layout of the buffer
        # the axes of the symbolic fields are (z, y, x): the loops run over
        # z, y, x from the outside to the inside, the buffer index grows with
        # the innermost loop (x), and the strides of the waLBerla field are
        # bound in reverse order
        source = generate("vec", "double", "cpu")[
            "PackInfoVecDoublePrecision.cpp"]
        for axis, spatial in ((0, "z"), (2, "x")):
            self.assertIn(f"_size_force_src_{axis} {{ int64_t( "
                          f"srcInterval.{spatial}Size() ) }}", source)
            self.assertIn(f"_stride_force_src_{axis} {{ int64_t( "
                          f"force_src->{spatial}Stride() ) }}", source)
        # the component axis is not reversed
        self.assertIn("_stride_force_src_3 { int64_t( "
                      "force_src->fStride() ) }", source)
        loops = re.findall(r"for\(int64_t (ctr_\d) = ", source)[:3]
        self.assertEqual(loops, ["ctr_0", "ctr_1", "ctr_2"])
        self.assertIn(
            "_data_buffer[((ctr_0 * _size_force_src_1 + ctr_1) * "
            "_size_force_src_2 + ctr_2) * 3LL]", source)

    def test_unpack_binds_the_destination_field(self):
        source = generate("pdf", "double", "cpu")[
            "PackInfoPdfDoublePrecision.cpp"]
        unpack = source[source.index("::doUnpack("):source.index(
            "::doLocalCopy(")]
        self.assertIn("pdfs_dst->dataAt(dstInterval.xMin()", unpack)
        self.assertNotIn("pdfs_src", unpack)

    def test_cuda_pdf(self):
        files = generate("pdf", "double", "cuda")
        self.assertEqual(sorted(files), ["PackInfoPdfDoublePrecisionCUDA.cu",
                                         "PackInfoPdfDoublePrecisionCUDA.h"])
        header = files["PackInfoPdfDoublePrecisionCUDA.h"]
        self.assertIn(
            "class PackInfoPdfDoublePrecisionCUDA : public "
            "walberla::sweepgen::UniformGpuFieldPackInfoBase< "
            "PackInfoPdfDoublePrecisionCUDA >", header)
        self.assertIn(
            '#include "walberla/sweepgen/communication/'
            'UniformGpuFieldPackInfoBase.hpp"', header)
        self.assertIn("using Field_T = walberla::gpu::GPUField< double >;",
                      header)
        self.assertEqual(len(re.findall(r"cudaStream_t stream\)", header)), 3)
        self.assertEqual(elements_per_cell(header),
                         elements_per_cell(generate(
                             "pdf", "double", "cpu")[
                                 "PackInfoPdfDoublePrecision.h"]))
        source = files["PackInfoPdfDoublePrecisionCUDA.cu"]
        self.assertIn(
            "namespace PackInfoPdfDoublePrecisionCUDA_kernels {", source)
        self.assertIn("__global__", source)

    def test_cuda_vec_single(self):
        files = generate("vec", "single", "cuda")
        self.assertEqual(sorted(files), ["PackInfoVecSinglePrecisionCUDA.cu",
                                         "PackInfoVecSinglePrecisionCUDA.h"])
        self.assertIn("using Field_T = walberla::gpu::GPUField< float >;",
                      files["PackInfoVecSinglePrecisionCUDA.h"])

    def test_kernel_namespaces_follow_the_class_names(self):
        for kind in ("pdf", "vec"):
            for precision in ("double", "single"):
                for target in ("cpu", "cuda"):
                    files = generate(kind, precision, target)
                    namespaces = {
                        n for text in files.values() for n in re.findall(
                            r"namespace (\w+_kernels) \{", text)}
                    name = {n.rsplit(".", 1)[0] for n in files}
                    self.assertEqual(len(name), 1)
                    self.assertEqual(namespaces, {f"{name.pop()}_kernels"})

    def test_generation_emits_no_warnings(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            generate_uncached("pdf", "double", "cpu")
            generate_uncached("pdf", "double", "cuda")
        self.assertEqual([str(w.message) for w in caught], [])

    def test_scalar_field(self):
        # fields without index dimension (scalar fields of the EK)
        field = ps.Field.create_generic(
            "phi", 3, np.float64, index_dimensions=0)
        comm_sets = pack_infos.vector_sets(D3Q19, 1)
        for target, suffix in (("cpu", ""), ("cuda", "CUDA")):
            ctx = sfg_driver.CodegenContext("double", target)
            name = ctx.class_name("PackInfoScalar")
            tmp = pathlib.Path(tempfile.mkdtemp())
            try:
                with sfg_driver.KernelFileGenerator(
                        name, ctx, tmp, "generate_lb_kernels.py") as sfg:
                    sfg.generate(
                        pack_infos.FieldPackInfo(name, field, comm_sets))
                files = {p.name: p.read_text() for p in tmp.iterdir()}
            finally:
                shutil.rmtree(tmp)
            header = files[f"{name}.h"]
            source = next(t for n, t in files.items() if n != f"{name}.h")
            field_t = ("walberla::gpu::GPUField< double >" if suffix
                       else "walberla::field::GhostLayerField< double, 1 >")
            self.assertIn(f"using Field_T = {field_t};", header)
            self.assertEqual(set(elements_per_cell(header).values()), {1})
            # no component axis
            self.assertNotIn("fStride", source)
            self.assertNotIn("_stride_phi_src_3", source)
            # one kernel is shared by all directions
            self.assertEqual(len(re.findall(
                r"void pack\w+ ?\(", source)), 1)

    def test_fields_with_two_index_dimensions_are_rejected(self):
        field = ps.Field.create_generic(
            "m", 3, np.float64, index_dimensions=2, index_shape=(3, 3))
        with self.assertRaises(ValueError):
            pack_infos.FieldPackInfo("X", field, {})

    def test_unsupported_target(self):
        with self.assertRaises(AssertionError):
            generate_uncached("pdf", "double", "avx")


class CpuBaseHeaderTest(unittest.TestCase):

    def test_include_of_the_generated_header_resolves(self):
        header = generate("pdf", "double", "cpu")[
            "PackInfoPdfDoublePrecision.h"]
        include = re.search(r'#include "([^"]*UniformCpuFieldPack[^"]*)"',
                            header).group(1)
        generated_kernels = (
            BRIDGE_SRC / "lattice_boltzmann" / "generated_kernels")
        self.assertEqual(
            os.path.normpath(generated_kernels / include),
            str(CPU_BASE_HEADER))
        self.assertTrue(CPU_BASE_HEADER.is_file())


if __name__ == "__main__":
    unittest.main()
