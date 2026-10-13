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
OUT = cat.Value(cat.SCALAR, 1, "out")


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
        both = cat.Accessor("get", cat.Kind.CELL, "k", (FIELD,), out=OUT)
        cuda_only = cat.Accessor("scatter", cat.Kind.SCATTER, None, (FIELD,),
                                 dims=3, targets=cat.CUDA_ONLY)
        return cat.Namespace(
            name, kernels=kernels or (make_kernel("k"),),
            accessors=accessors or (both, cuda_only))

    def test_namespaces_of_a_target(self):
        cpu_only = cat.Accessor("point", cat.Kind.POINT, point="data",
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
        with self.assertRaises(cat.CatalogueError):
            make_catalogue(self.namespace()).for_target("hip")

    def test_unknown_kernel(self):
        accessor = cat.Accessor("get", cat.Kind.CELL,
                                "missing", (FIELD,), out=OUT)
        with self.assertRaises(cat.CatalogueError):
            make_catalogue(self.namespace(accessors=(accessor,)))

    def test_kernel_must_exist_on_the_targets_of_the_accessor(self):
        accessor = cat.Accessor("get", cat.Kind.CELL, "k", (FIELD,), out=OUT)
        namespace = self.namespace(kernels=(make_kernel("k", ("cpu",)),),
                                   accessors=(accessor,))
        with self.assertRaises(cat.CatalogueError):
            make_catalogue(namespace)
        accessor = cat.Accessor(
            "get", cat.Kind.CELL, make_kernel("p", ("cpu",)), (FIELD,), out=OUT)
        with self.assertRaises(cat.CatalogueError):
            make_catalogue(self.namespace(accessors=(accessor,)))

    def test_unknown_field(self):
        accessor = cat.Accessor(
            "get", cat.Kind.CELL, "k", (cat.FieldArg("x", True, "f"),), out=OUT)
        with self.assertRaises(cat.CatalogueError):
            make_catalogue(self.namespace(accessors=(accessor,)))

    def test_rename_of_an_unknown_parameter(self):
        accessor = cat.Accessor(
            "get", cat.Kind.CELL, "k", (FIELD,),
            out=OUT, targets={"cpu": cat.Variant(names={"unknown": "x"})})
        with self.assertRaises(cat.CatalogueError):
            make_catalogue(self.namespace(accessors=(accessor,)))

    def test_duplicate_overload(self):
        accessor = cat.Accessor("get", cat.Kind.CELL, "k", (FIELD,), out=OUT)
        with self.assertRaises(cat.CatalogueError):
            make_catalogue(self.namespace(accessors=(accessor, accessor)))
        interval = cat.Accessor("get", cat.Kind.INTERVAL, "k", (FIELD,),
                                out=OUT)
        make_catalogue(self.namespace(accessors=(accessor, interval)))

    def test_duplicate_overload_on_a_target(self):
        both = cat.Accessor("get", cat.Kind.CELL, "k", (FIELD,), out=OUT)
        cpu_only = cat.Accessor("get", cat.Kind.CELL, "k", (FIELD,), out=OUT,
                                targets=cat.CPU_ONLY)
        cuda_only = cat.Accessor("get", cat.Kind.CELL, "k", (FIELD,),
                                 out=OUT, targets=cat.CUDA_ONLY)
        for accessors in ((both, cpu_only), (cuda_only, both)):
            with self.assertRaisesRegex(cat.CatalogueError, "defined twice"):
                make_catalogue(self.namespace(accessors=accessors))
        make_catalogue(self.namespace(accessors=(cpu_only, cuda_only)))

    def test_duplicate_kernel(self):
        with self.assertRaises(cat.CatalogueError):
            make_catalogue(self.namespace(
                kernels=(make_kernel("k"), make_kernel("k"))))

    def test_scalars_and_the_domain_cannot_be_renamed(self):
        for name in ("density", "cell"):
            accessor = cat.Accessor(
                "get", cat.Kind.CELL, "k", (FIELD,), out=OUT,
                scalars=("density",),
                targets={"cuda": cat.Variant(names={name: "x"})})
            with self.assertRaisesRegex(cat.CatalogueError, "can't rename"):
                make_catalogue(self.namespace(accessors=(accessor,)))
        renamed = cat.Accessor(
            "get", cat.Kind.CELL, "k", (FIELD,), inp=VEC,
            targets={"cuda": cat.Variant(names={"field": "f", "vec": "v"})})
        make_catalogue(self.namespace(accessors=(renamed,)))

    def check_kind(self, message, kind, *args, kernel="k", **kwargs):
        accessor = cat.Accessor("a", kind, kernel, *args, **kwargs)
        with self.assertRaisesRegex(cat.CatalogueError, message):
            make_catalogue(self.namespace(
                kernels=(make_kernel("k"),
                         cat.ReductionKernel("r", lambda: ([], []))),
                accessors=(accessor,)))

    def test_kind_requirements(self):
        K = cat.Kind
        reduction = cat.Reduction(cat.VECTOR3, "r")
        self.check_kind("requires exactly one of inp and out",
                        K.CELL, (FIELD,))
        self.check_kind("requires exactly one of inp and out",
                        K.INTERVAL, (FIELD,), inp=VEC, out=OUT)
        self.check_kind("requires fields", K.CELL, (), out=OUT)
        self.check_kind("requires inp", K.ALL, (FIELD,))
        self.check_kind("kernel with components",
                        K.ALL, (FIELD,), inp=VEC)
        self.check_kind("requires reduction", K.REDUCE, (FIELD,),
                        kernel="r")
        self.check_kind("requires a ReductionKernel", K.REDUCE, (FIELD,),
                        reduction=reduction)
        self.check_kind("requires a ReductionKernel", K.REDUCE, (FIELD,),
                        kernel="k", reduction=reduction)
        self.check_kind("the other kinds a Kernel", K.CELL, (FIELD,),
                        kernel="r", out=OUT)
        self.check_kind("needs a size", K.REDUCE, (FIELD,), kernel="r",
                        reduction=cat.Reduction(cat.STD_VECTOR, "r"))
        self.check_kind("unknown scale", K.REDUCE, (FIELD,), kernel="r",
                        reduction=cat.Reduction(cat.VECTOR3, "r", "s"))
        self.check_kind("requires point", K.POINT, kernel=None)
        self.check_kind("doesn't take a kernel", K.POINT, point="p")
        self.check_kind("requires dims", K.LIST, (FIELD,))
        self.check_kind("requires dims", K.SCATTER, (FIELD,), kernel=None)

    def test_kind_exclusions(self):
        K = cat.Kind
        self.check_kind("doesn't take reduction", K.CELL, (FIELD,), out=OUT,
                        reduction=cat.Reduction(cat.VECTOR3, "r"))
        self.check_kind("doesn't take dims", K.INTERVAL, (FIELD,), out=OUT,
                        dims=3)
        self.check_kind("doesn't take out", K.ALL, (FIELD,), inp=VEC,
                        out=OUT)
        self.check_kind("doesn't take inp", K.REDUCE, (FIELD,), kernel="r",
                        inp=VEC, reduction=cat.Reduction(cat.VECTOR3, "r"))
        self.check_kind("doesn't take fields", K.POINT, (FIELD,),
                        kernel=None, point="p")
        self.check_kind("doesn't take point", K.GATHER, (FIELD,), dims=3,
                        point="p")

    def test_valid_kinds(self):
        K = cat.Kind
        kernels = (cat.Kernel("b", lambda: [], components=("c", 3)),
                   cat.ReductionKernel("r", lambda: ([], [])))
        make_catalogue(cat.Namespace("N", kernels, (
            cat.Accessor("a", K.ALL, "b", (FIELD_MUT,), inp=VEC),
            cat.Accessor("r", K.REDUCE, "r", (FIELD,),
                         reduction=cat.Reduction(cat.VECTOR3, "r")),
            cat.Accessor("p", K.POINT, point="p"),
            cat.Accessor("s", K.SCATTER, None, (FIELD,), dims=3,
                         targets=cat.CUDA_ONLY))))


class ValidationTest(unittest.TestCase):

    def test_container_kind_is_checked(self):
        with self.assertRaisesRegex(cat.CatalogueError, "container kind"):
            cat.Container("vector", 3)

    def test_value_must_fit_its_container(self):
        cat.Value(cat.VECTOR3, 3, "v")
        cat.Value(cat.STD_VECTOR, 5, "v")
        cat.Value(cat.array(19), 19, "v")
        with self.assertRaisesRegex(cat.CatalogueError, "can't hold 5"):
            cat.Value(cat.VECTOR3, 5, "v")

    def test_return_type_of_a_variant(self):
        cat.Variant(returns="auto")
        with self.assertRaisesRegex(cat.CatalogueError, "return type"):
            cat.Variant(returns="double")

    def test_target_mappings_are_read_only(self):
        for targets in (cat.BOTH, cat.CPU_ONLY, cat.CUDA_ONLY):
            with self.assertRaises(TypeError):
                targets["cpu"] = cat.Variant()
        accessor = cat.Accessor("a", cat.Kind.CELL)
        self.assertEqual(dict(accessor.targets), dict(cat.BOTH))

    def test_catalogue_error_is_a_value_error(self):
        self.assertTrue(issubclass(cat.CatalogueError, ValueError))


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
        both = cat.Accessor("a", cat.Kind.CELL, "k1", (FIELD,), out=OUT)
        other = cat.Accessor("b", cat.Kind.CELL, "k2", (FIELD,), out=OUT)
        namespace = cat.Namespace(
            "N", (make_kernel("k1"), make_kernel("k2")), (both, other))
        self.assertEqual(self.run_backend(namespace), [
            ("kernel", "k1"), ("kernel", "k2"),
            ("cell", "a", "k1", "k1"), ("cell", "b", "k2", "k2")])

    def test_private_kernels_are_created_at_their_accessor(self):
        first = cat.Accessor("a", cat.Kind.CELL, "k1", (FIELD,), out=OUT)
        private = cat.Accessor(
            "b", cat.Kind.CELL, make_kernel("late"), (FIELD,), out=OUT)
        last = cat.Accessor("c", cat.Kind.CELL, "k1", (FIELD,), out=OUT)
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
        accessor = cat.Accessor("a", cat.Kind.CELL, "k1", (FIELD,), out=OUT)
        namespace = cat.Namespace(
            "N", (make_kernel("k1"), make_kernel("k2", ("cuda",))),
            (accessor,))
        self.assertEqual(self.run_backend(namespace)[0], ("kernel", "k1"))
        self.assertNotIn(("kernel", "k2"), self.run_backend(namespace))

    def test_unsupported_kind(self):
        namespace = cat.Namespace("N", (make_kernel("k1"),), (
            cat.Accessor("a", cat.Kind.INTERVAL, "k1", (FIELD,), out=OUT),))
        with self.assertRaisesRegex(ValueError, r"N::a of the kind interval"):
            self.run_backend(namespace)

    def test_the_backend_checks_the_target_of_the_context(self):
        with self.assertRaises(ValueError):
            self.run_backend(self.namespace_of_a_cell(), target="cuda")

    @staticmethod
    def namespace_of_a_cell():
        return cat.Namespace("N", (make_kernel("k1"),), (
            cat.Accessor("a", cat.Kind.CELL, "k1", (FIELD,), out=OUT),))


if __name__ == "__main__":
    unittest.main()
