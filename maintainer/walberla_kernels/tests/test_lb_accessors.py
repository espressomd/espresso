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

import collections
import pathlib
import unittest

import pystencils as ps

import accessor_catalogue as cat
import generate_lb_kernels
import lb_accessors
import pystencils_espresso
import sfg_driver
from tests import header_api

DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"
# the frozen APIs have one function per accessor of the target
N_FUNCTIONS = {"cpu": 30, "cuda": 32}
API_FILES = {"cpu": "field_accessors_api_{}.json",
             "cuda": "field_accessors_cuda_api_{}.json"}
# the kinds whose parameters are all known to the catalogue
EXACT_KINDS = (cat.Kind.CELL, cat.Kind.INTERVAL, cat.Kind.ALL,
               cat.Kind.REDUCE)


def create_catalogue(precision):
    ctx = sfg_driver.CodegenContext(precision, "cpu")
    fields = pystencils_espresso.generate_fields(
        generate_lb_kernels.STENCIL, ctx.data_type)
    method = generate_lb_kernels.create_lb_method(fields)
    return lb_accessors.LbAccessors(ctx.data_type, method).catalogue()


def frozen_api(target, precision):
    """The functions of a target as ``(namespace, name, parameter names,
    return type)``."""
    records = header_api.load_records(
        DATA_DIR / API_FILES[target].format(precision))
    return [(f[0].split("::")[-1], f[1], tuple(p[1] for p in f[2]), f[3])
            for f in records]


class LbCatalogueTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.catalogues = {p: create_catalogue(p) for p in ("double", "single")}

    def test_numbers_of_accessors(self):
        for catalogue in self.catalogues.values():
            for target, n in N_FUNCTIONS.items():
                self.assertEqual(len(catalogue.accessors(target)), n, target)

    def test_accessors_equal_the_frozen_api(self):
        # every accessor of a target is a function of its frozen API and
        # vice versa: the namespaces, the names, the parameter names, and
        # whether the return type is ``auto``
        for precision, catalogue in self.catalogues.items():
            for target in cat.TARGETS:
                where = (precision, target)
                frozen = frozen_api(target, precision)
                self.assertEqual(len(frozen), N_FUNCTIONS[target], where)
                accessors = catalogue.accessors(target)
                exact = sorted(
                    (n, a.name, a.parameter_names(target),
                     a.targets[target].returns == "auto")
                    for n, a in accessors if a.kind in EXACT_KINDS)
                loose = {(n, a.name): a for n, a in accessors
                         if a.kind not in EXACT_KINDS}
                self.assertEqual(
                    exact, sorted(
                        (n, name, names, return_type == "auto")
                        for n, name, names, return_type in frozen
                        if (n, name) not in loose), where)
                self.assertEqual(
                    sorted(loose), sorted(
                        (n, name) for n, name, _, _ in frozen
                        if (n, name) in loose), where)
                # the other accessors have more parameters than the
                # catalogue knows, in the same order
                for n, name, names, _ in frozen:
                    if (n, name) in loose:
                        known = loose[(n, name)].parameter_names(target)
                        self.assertEqual([p for p in names if p in known],
                                         list(known), (where, n, name))

    def test_target_only_accessors(self):
        catalogue = self.catalogues["double"]
        only = {t: {(n, a.name) for n, a in catalogue.accessors(t)}
                for t in cat.TARGETS}
        self.assertEqual(
            only["cpu"] - only["cuda"],
            {("Equilibrium", "set"), ("EquilibriumDistribution", "get")})
        self.assertEqual(
            only["cuda"] - only["cpu"],
            {("Interpolation", "get_rho"), ("Interpolation", "get_vel"),
             ("Interpolation", "add_force"), ("Vector", "set_from_list")})

    def test_namespaces(self):
        catalogue = self.catalogues["double"]
        self.assertEqual(
            [n.name for n in catalogue.for_target("cpu")],
            ["Population", "Vector", "EquilibriumDistribution",
             "Equilibrium", "Density", "Velocity", "Force",
             "MomentumDensity", "PressureTensor"])
        self.assertEqual(
            [n.name for n in catalogue.for_target("cuda")],
            ["Population", "Vector", "Interpolation", "Density", "Velocity",
             "Force", "MomentumDensity", "PressureTensor"])

    def test_kinds(self):
        counts = collections.Counter(
            a.kind for _, a in self.catalogues["double"].accessors("cuda"))
        self.assertEqual(counts, {
            cat.Kind.CELL: 12, cat.Kind.INTERVAL: 11, cat.Kind.ALL: 3,
            cat.Kind.REDUCE: 2, cat.Kind.LIST: 1, cat.Kind.GATHER: 2,
            cat.Kind.SCATTER: 1})

    def test_reductions_are_sums_of_the_values_of_the_cells(self):
        for _, acc in self.catalogues["double"].accessors("cpu"):
            if acc.kind == cat.Kind.REDUCE:
                kernel = [k for ns in self.catalogues["double"].namespaces
                          for k in ns.kernels if k.name == acc.kernel
                          and isinstance(k, cat.ReductionKernel)]
                self.assertEqual(len(kernel), 1, acc.name)
                _, values = kernel[0].build()
                self.assertEqual(len(values), acc.reduction.container.size)

    def test_the_kernels_build_assignments(self):
        for catalogue in self.catalogues.values():
            for ns in catalogue.namespaces:
                for kernel in ns.kernels:
                    built = kernel.build()
                    if isinstance(kernel, cat.ReductionKernel):
                        assignments, values = built
                        self.assertTrue(values, kernel.name)
                    else:
                        assignments = built
                    self.assertTrue(assignments, (ns.name, kernel.name))
                    for a in assignments:
                        self.assertIsInstance(a, ps.Assignment)

    def test_component_symbols(self):
        catalogue = self.catalogues["single"]
        symbols = catalogue.component_symbols("pop", 19)
        self.assertEqual([s.name for s in symbols],
                         [f"pop_{i}" for i in range(19)])
        self.assertTrue(all(s.dtype == catalogue.dtype for s in symbols))

    def test_fields(self):
        catalogue = self.catalogues["double"]
        self.assertEqual(
            {k: f.fsize for k, f in catalogue.fields.items()},
            {"pdf": 19, "vel": 3, "force": 3, "vec": 3, "pos": 3,
             "indices": 3, "values": 3, "out1": 1, "out3": 3})
        self.assertEqual(catalogue.dim, 3)
        for n in (1, 3, 9, 19):
            self.assertEqual(catalogue.buffer(n).index_shape, (n,))

    def test_the_stencil_header_is_a_cpu_include(self):
        self.assertEqual(self.catalogues["double"].includes,
                         {"cpu": ("stencil/D3Q19.h",)})

    def test_working_precision(self):
        self.assertEqual(self.catalogues["double"].dtype,
                         ps.create_type("float64"))
        self.assertEqual(self.catalogues["single"].dtype,
                         ps.create_type("float32"))


if __name__ == "__main__":
    unittest.main()
