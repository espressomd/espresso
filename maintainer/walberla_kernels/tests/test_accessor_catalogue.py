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

import unittest

import pystencils as ps

import accessor_catalogue as cat
import accessor_sfg
from tests.test_accessor_sfg import generate

DOUBLE = ps.create_type("float64")
VEC = cat.Value(cat.VECTOR3, 3, "vec")
FIELD = cat.FieldArg("scalar", True, "field")
FIELD_MUT = cat.FieldArg("scalar", False, "field")


def make_catalogue(*namespaces):
    """A catalogue of one scalar field and the given namespaces."""
    field = ps.Field.create_generic("f", 3, DOUBLE, index_dimensions=1,
                                    index_shape=(1,))
    return cat.Catalogue(
        dtype=DOUBLE, namespace="ns", prelude="", includes={},
        fields={"scalar": cat.FieldSpec(field, 1)}, dim=3,
        namespaces=namespaces,
        component_symbols=lambda prefix, n: [], buffer=lambda n: field)


def make_kernel(name, targets=cat.TARGETS):
    return cat.Kernel(name, lambda: [], targets=targets)


class ContainerTest(unittest.TestCase):

    def test_sizes(self):
        self.assertEqual(cat.SCALAR.size, 1)
        self.assertEqual(cat.VECTOR3.size, 3)
        self.assertEqual(cat.MATRIX3.size, 9)
        self.assertEqual(cat.array(19), cat.Container("array", 19))
        self.assertIsNone(cat.STD_VECTOR.size)


class AccessorTest(unittest.TestCase):

    ACCESSOR = cat.Accessor(
        "set", cat.Kind.CELL, "set", (FIELD_MUT,), inp=VEC,
        scalars=("density",), targets={
            "cpu": cat.Variant(),
            "cuda": cat.Variant(names={"field": "vec_field", "vec": "u"})})

    def test_default_names_are_in_the_order_of_the_signature(self):
        self.assertEqual(self.ACCESSOR.default_names(),
                         ("field", "vec", "density", "cell"))
        interval = cat.Accessor(
            "get", cat.Kind.INTERVAL, "get", (FIELD,),
            out=cat.Value(cat.STD_VECTOR, 1, "out"))
        self.assertEqual(interval.default_names(), ("field", "ci"))
        broadcast = cat.Accessor("init", cat.Kind.ALL, "b", (FIELD_MUT,),
                                 inp=VEC)
        self.assertEqual(broadcast.default_names(), ("field", "vec"))

    def test_parameter_names_of_the_targets(self):
        self.assertEqual(self.ACCESSOR.parameter_names("cpu"),
                         ("field", "vec", "density", "cell"))
        self.assertEqual(self.ACCESSOR.parameter_names("cuda"),
                         ("vec_field", "u", "density", "cell"))
        self.assertEqual(self.ACCESSOR.param_name("vec", "cuda"), "u")

    def test_accessors_exist_on_the_targets_of_the_mapping(self):
        self.assertEqual(set(cat.BOTH), {"cpu", "cuda"})
        self.assertEqual(set(cat.CPU_ONLY), {"cpu"})
        self.assertEqual(set(cat.CUDA_ONLY), {"cuda"})
        self.assertEqual(set(cat.Accessor("a", cat.Kind.CELL).targets),
                         {"cpu", "cuda"})


class CatalogueTest(unittest.TestCase):

    @staticmethod
    def namespace(name="N", kernels=None, accessors=None):
        both = cat.Accessor("get", cat.Kind.CELL, "k", (FIELD,))
        cuda_only = cat.Accessor("scatter", cat.Kind.SCATTER, None, (FIELD,),
                                 targets=cat.CUDA_ONLY)
        return cat.Namespace(
            name, kernels=kernels or (make_kernel("k"),),
            accessors=accessors or (both, cuda_only))

    def test_namespaces_of_a_target(self):
        cpu_only = cat.Accessor("point", cat.Kind.POINT,
                                targets=cat.CPU_ONLY)
        catalogue = make_catalogue(
            self.namespace("A"),
            cat.Namespace("OnlyCuda", (), (self.namespace().accessors[1],)),
            cat.Namespace("OnlyCpu", (), (cpu_only,)))
        self.assertEqual([n.name for n in catalogue.for_target("cpu")],
                         ["A", "OnlyCpu"])
        self.assertEqual([n.name for n in catalogue.for_target("cuda")],
                         ["A", "OnlyCuda"])
        self.assertEqual(
            [(n, a.name) for n, a in catalogue.accessors("cuda")],
            [("A", "get"), ("A", "scatter"), ("OnlyCuda", "scatter")])
        cpu_namespace = catalogue.for_target("cpu")[0]
        self.assertEqual([a.name for a in cpu_namespace.accessors], ["get"])

    def test_kernels_of_a_target(self):
        namespace = self.namespace(kernels=(
            make_kernel("a"), make_kernel("b", ("cuda",)),
            make_kernel("c", ("cpu",))))
        self.assertEqual(
            [k.name for k in namespace.for_target("cpu")[0]], ["a", "c"])
        self.assertEqual(
            [k.name for k in namespace.for_target("cuda")[0]], ["a", "b"])

    def test_unknown_target(self):
        with self.assertRaises(AssertionError):
            make_catalogue(self.namespace()).for_target("hip")

    def test_unknown_kernel(self):
        accessor = cat.Accessor("get", cat.Kind.CELL, "missing", (FIELD,))
        with self.assertRaises(AssertionError):
            make_catalogue(self.namespace(accessors=(accessor,)))

    def test_kernel_must_exist_on_the_targets_of_the_accessor(self):
        accessor = cat.Accessor("get", cat.Kind.CELL, "k", (FIELD,))
        namespace = self.namespace(kernels=(make_kernel("k", ("cpu",)),),
                                   accessors=(accessor,))
        with self.assertRaises(AssertionError):
            make_catalogue(namespace)
        accessor = cat.Accessor(
            "get", cat.Kind.CELL, make_kernel("p", ("cpu",)), (FIELD,))
        with self.assertRaises(AssertionError):
            make_catalogue(self.namespace(accessors=(accessor,)))

    def test_unknown_field(self):
        accessor = cat.Accessor(
            "get", cat.Kind.CELL, "k", (cat.FieldArg("x", True, "f"),))
        with self.assertRaises(AssertionError):
            make_catalogue(self.namespace(accessors=(accessor,)))

    def test_rename_of_an_unknown_parameter(self):
        accessor = cat.Accessor(
            "get", cat.Kind.CELL, "k", (FIELD,),
            targets={"cpu": cat.Variant(names={"unknown": "x"})})
        with self.assertRaises(AssertionError):
            make_catalogue(self.namespace(accessors=(accessor,)))

    def test_duplicate_overload(self):
        accessor = cat.Accessor("get", cat.Kind.CELL, "k", (FIELD,))
        with self.assertRaises(AssertionError):
            make_catalogue(self.namespace(accessors=(accessor, accessor)))
        interval = cat.Accessor("get", cat.Kind.INTERVAL, "k", (FIELD,))
        make_catalogue(self.namespace(accessors=(accessor, interval)))

    def test_duplicate_kernel(self):
        with self.assertRaises(AssertionError):
            make_catalogue(self.namespace(
                kernels=(make_kernel("k"), make_kernel("k"))))


class ToyBackend(accessor_sfg.AccessorBackend):
    """Backend that records what it is asked to generate."""

    target = "cpu"
    ps_target = ps.Target.GenericCPU
    events = None

    def create_kernel(self, spec):
        self.events.append(("kernel", spec.name))
        return spec.name

    def emit_cell(self, acc, spec, kernel):
        self.events.append(("cell", acc.name, spec.name, kernel))

    def emit_point(self, acc, spec, kernel):
        self.events.append(("point", acc.name, spec, kernel))


class BackendTest(unittest.TestCase):

    def run_backend(self, namespace, backend=ToyBackend, target="cpu"):
        events = []

        def body(sfg):
            ctx = type("Ctx", (), {"is_gpu": target == "cuda",
                                   "target": target})()
            backend.events = events
            backend(sfg, ctx, make_catalogue(namespace)).generate()
        generate(body, target="cpu")
        return events

    def test_shared_kernels_precede_the_accessors(self):
        both = cat.Accessor("a", cat.Kind.CELL, "k1", (FIELD,))
        other = cat.Accessor("b", cat.Kind.CELL, "k2", (FIELD,))
        namespace = cat.Namespace(
            "N", (make_kernel("k1"), make_kernel("k2")), (both, other))
        self.assertEqual(self.run_backend(namespace), [
            ("kernel", "k1"), ("kernel", "k2"),
            ("cell", "a", "k1", "k1"), ("cell", "b", "k2", "k2")])

    def test_private_kernels_are_created_at_their_accessor(self):
        first = cat.Accessor("a", cat.Kind.CELL, "k1", (FIELD,))
        private = cat.Accessor(
            "b", cat.Kind.CELL, make_kernel("late"), (FIELD,))
        last = cat.Accessor("c", cat.Kind.CELL, "k1", (FIELD,))
        namespace = cat.Namespace("N", (make_kernel("k1"),),
                                  (first, private, last))
        self.assertEqual(self.run_backend(namespace), [
            ("kernel", "k1"), ("cell", "a", "k1", "k1"),
            ("kernel", "late"), ("cell", "b", "late", "late"),
            ("cell", "c", "k1", "k1")])

    def test_accessors_without_kernel(self):
        namespace = cat.Namespace("N", (), (
            cat.Accessor("p", cat.Kind.POINT, point="data"),))
        self.assertEqual(self.run_backend(namespace),
                         [("point", "p", None, None)])

    def test_unsupported_target_kernels_are_skipped(self):
        accessor = cat.Accessor("a", cat.Kind.CELL, "k1", (FIELD,))
        namespace = cat.Namespace(
            "N", (make_kernel("k1"), make_kernel("k2", ("cuda",))),
            (accessor,))
        self.assertEqual(self.run_backend(namespace)[0], ("kernel", "k1"))
        self.assertNotIn(("kernel", "k2"), self.run_backend(namespace))

    def test_unsupported_kind(self):
        namespace = cat.Namespace("N", (make_kernel("k1"),), (
            cat.Accessor("a", cat.Kind.ALL, "k1", (FIELD,)),))
        with self.assertRaisesRegex(ValueError, r"N::a of the kind all"):
            self.run_backend(namespace)

    def test_the_backend_checks_the_target_of_the_context(self):
        with self.assertRaises(ValueError):
            self.run_backend(self.namespace_of_a_cell(), target="cuda")

    @staticmethod
    def namespace_of_a_cell():
        return cat.Namespace("N", (make_kernel("k1"),), (
            cat.Accessor("a", cat.Kind.CELL, "k1", (FIELD,)),))


if __name__ == "__main__":
    unittest.main()
