#
# Copyright (C) 2021-2026 The ESPResSo project
# Copyright (C) 2020-2022 The waLBerla project
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

import re
import unittest

import pystencils as ps

import accessor_backend_cpu
import accessor_backend_cuda
import accessor_catalogue as cat
from tests import header_api
from tests.test_accessor_sfg import generate

# The backends know nothing about a lattice model: they work from the
# catalogue. This catalogue of a scalar field (the accessors of the EK
# densities) has no population field and no lattice stencil.
DOUBLE = ps.create_type("float64")
SCALAR_FIELD = cat.FieldArg("scalar", True, "scalar_field")
SCALAR_FIELD_MUT = cat.FieldArg("scalar", False, "scalar_field")


def scalar_catalogue():
    def field(name):
        return ps.Field.create_generic(
            name, 3, DOUBLE, index_dimensions=1, layout="fzyx",
            index_shape=(1,))
    scalar, buf = field("scalar"), field("buf1")

    def symbols(prefix, n):
        return [ps.TypedSymbol(f"{prefix}_{i}", DOUBLE) for i in range(n)]
    value, values = cat.Value(cat.SCALAR, 1, "value"), cat.Value(
        cat.STD_VECTOR, 1, "values")
    out = cat.Value(cat.SCALAR, 1, "value")
    accessors = (
        cat.Accessor("get", cat.Kind.CELL, "get", (SCALAR_FIELD,), out=out,
                     targets={"cpu": cat.Variant(returns="auto"),
                              "cuda": cat.Variant()}),
        cat.Accessor("set", cat.Kind.CELL, "set", (SCALAR_FIELD_MUT,),
                     inp=value),
        cat.Accessor("initialize", cat.Kind.ALL, "broadcast",
                     (SCALAR_FIELD_MUT,), inp=value),
        cat.Accessor("get", cat.Kind.INTERVAL, "get", (SCALAR_FIELD,),
                     out=cat.Value(cat.STD_VECTOR, 1, "out")),
        cat.Accessor("set", cat.Kind.INTERVAL, "set", (SCALAR_FIELD_MUT,),
                     inp=values))
    kernels = (
        cat.Kernel("get", lambda: [ps.Assignment(buf(0), scalar(0))]),
        cat.Kernel("set", lambda: [ps.Assignment(scalar(0), buf(0))]),
        cat.Kernel("broadcast", lambda: [
            ps.Assignment(scalar(0), symbols("value", 1)[0])],
            components=("value", 1)))
    return cat.Catalogue(
        dtype=DOUBLE, namespace="walberla::ek::accessor",
        prelude="Scalar accessors.", includes={},
        fields={"scalar": cat.FieldSpec(scalar, 1)}, dim=3,
        namespaces=(cat.Namespace("Scalar", kernels, accessors),),
        component_symbols=symbols, buffer=lambda n: buf)


def generate_scalar_accessors(backend, target):
    catalogue = scalar_catalogue()

    def body(sfg):
        ctx = type("Ctx", (), {"is_gpu": target == "cuda",
                               "target": target})()
        backend(sfg, ctx, catalogue).generate()
    return generate(body, header_only=target == "cpu", target=target)


def signatures(header, declarations):
    return sorted((f[0], f[1], tuple((p[0], p[1]) for p in f[2]), f[3])
                  for f in header_api.functions(header, declarations))


class ScalarFieldTest(unittest.TestCase):
    """A catalogue of a scalar field generates the accessors on both
    targets."""

    FIELD = "GhostLayerField<double,1>"

    def test_cpu(self):
        files = generate_scalar_accessors(
            accessor_backend_cpu.CpuBackend, "cpu")
        self.assertEqual(sorted(files), ["Foo.h"])
        header = files["Foo.h"]
        f = self.FIELD
        self.assertEqual(signatures(header, False), sorted([
            ("ek::accessor::Scalar", "get",
             ((f"{f} const*", "scalar_field"), ("Cell const&", "cell")),
             "auto"),
            ("ek::accessor::Scalar", "set",
             ((f"{f}*", "scalar_field"), ("double", "value"),
              ("Cell const&", "cell")), "void"),
            ("ek::accessor::Scalar", "initialize",
             ((f"{f}*", "scalar_field"), ("double", "value")), "void"),
            ("ek::accessor::Scalar", "get",
             ((f"{f} const*", "scalar_field"), ("CellInterval const&", "ci")),
             "std::vector<double>"),
            ("ek::accessor::Scalar", "set",
             ((f"{f}*", "scalar_field"), ("std::vector<double> const&",
                                          "values"),
              ("CellInterval const&", "ci")), "void")]))
        self.assertEqual(
            len(re.findall(r"namespace\s+kernels\s*\{", header)), 3)

    def test_cuda(self):
        files = generate_scalar_accessors(
            accessor_backend_cuda.CudaBackend, "cuda")
        self.assertEqual(sorted(files), ["Foo.cu", "Foo.h"])
        self.assertEqual(
            [f[:2] for f in signatures(files["Foo.h"], True)],
            sorted([("ek::accessor::Scalar", n) for n in (
                "get", "set", "initialize", "get", "set")]))

    def test_the_scalar_input_of_a_broadcast_is_the_kernel_parameter(self):
        header = generate_scalar_accessors(
            accessor_backend_cpu.CpuBackend, "cpu")["Foo.h"]
        self.assertIn("kernels::broadcast(", header)
        self.assertNotIn("value[", header)


if __name__ == "__main__":
    unittest.main()
