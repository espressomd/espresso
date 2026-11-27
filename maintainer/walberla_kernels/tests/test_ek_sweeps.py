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
Tests of the flux sweeps of the electrokinetics: the guard of the faces of
the extended region, the finite-volume stencils, the structure of the
generated classes (all targets and precisions) and the registry.
"""

import itertools
import os
import pathlib
import re
import subprocess
import sys
import unittest

import pystencils as ps
import sympy as sp
from pystencils.defaults import DEFAULTS

import ek_sweeps
import flux_layout
import generate_ek_kernels
import sfg_driver
from tests import ek_generated

DOUBLE = ps.create_type("float64")
KERNELS = pathlib.Path(__file__).resolve().parents[1]


class FaceConditionTest(unittest.TestCase):
    """
    The sweeps run on the interior of a block expanded by one cell. The
    counters of the cells of the interior are 1 to ``shape - 2``.
    """

    shape = (6, 5, 4)

    def computed(self, edge, cell):
        condition = ek_sweeps.face_condition(edge, self.shape)
        values = dict(zip(DEFAULTS.spatial_counters, cell))
        return bool(condition.subs(values))

    def inside(self, cell):
        return all(0 <= c < n for c, n in zip(cell, self.shape))

    def interior(self, cell):
        return all(1 <= c <= n - 2 for c, n in zip(cell, self.shape))

    def test_conditions(self):
        # both cells of the face are in the extended region and, on the axes
        # that the face doesn't cross, the cell is an interior cell
        for edge in flux_layout.STORED_DIRECTIONS:
            for cell in itertools.product(*map(range, self.shape)):
                other = tuple(c + e for c, e in zip(cell, edge))
                expected = self.inside(cell) and self.inside(other) and all(
                    1 <= c <= n - 2 for c, e, n in zip(cell, edge, self.shape)
                    if e == 0)
                self.assertEqual(self.computed(edge, cell), expected,
                                 (edge, cell))

    def test_every_face_of_an_interior_cell_is_computed_once(self):
        # the face between an interior cell and a neighbor is stored by the
        # cell that owns it, in the component of its direction
        for cell in itertools.product(*map(range, self.shape)):
            if not self.interior(cell):
                continue
            for direction in flux_layout.NEIGHBORS:
                offset, component, _ = flux_layout.face_location(direction)
                owner = tuple(c + o for c, o in zip(cell, offset))
                edge = flux_layout.STORED_DIRECTIONS[component]
                self.assertTrue(self.computed(edge, owner),
                                (cell, direction))

    def test_faces_between_ghost_cells_are_not_accessed_outside(self):
        # no face reaches outside of the extended region
        for edge in flux_layout.STORED_DIRECTIONS:
            for cell in itertools.product(*map(range, self.shape)):
                if self.computed(edge, cell):
                    other = tuple(c + e for c, e in zip(cell, edge))
                    self.assertTrue(self.inside(cell) and self.inside(other))

    def test_condition_is_a_range_on_each_axis(self):
        for edge in flux_layout.STORED_DIRECTIONS:
            condition = ek_sweeps.face_condition(edge, self.shape)
            self.assertEqual(len(condition.args), 6)


class DiscreteFluxTest(unittest.TestCase):
    """The stencils of the flux through the faces."""

    # a linear density and a linear potential
    density_slope = (sp.Rational(3, 10), sp.Rational(-1, 5), sp.Rational(1, 2))
    potential_slope = (sp.Rational(1, 10), sp.Rational(1, 4),
                       sp.Rational(-3, 20))
    ext_efield = (sp.Rational(1, 20), sp.Rational(-1, 10), sp.Rational(1, 5))
    diffusion, kT, valency = sp.Rational(37, 100), sp.Rational(13, 10), 2

    def evaluate(self, store):
        """The value of the right-hand side for the linear fields."""
        values = {}
        for access in store.rhs.atoms(ps.Field.Access):
            offset = access.offsets
            slope = (self.density_slope if access.field.name == "density"
                     else self.potential_slope)
            offset_value = 2 + sum(s * o for s, o in zip(slope, offset))
            values[access] = offset_value
        symbols = {s.name: s for s in store.rhs.atoms(ps.TypedSymbol)
                   if not isinstance(s, ps.Field.Access)}
        numbers = {"diffusion": self.diffusion, "kT": self.kT,
                   "valency": self.valency}
        numbers.update({f"ext_efield_{i}": e
                        for i, e in enumerate(self.ext_efield)})
        values.update({symbols[n]: v for n, v in numbers.items()
                       if n in symbols})
        return sp.N(store.rhs.xreplace(values), 30)

    def expected(self, edge, electrostatic):
        a0 = sum(sp.Matrix(e).norm() for e in
                 flux_layout.STORED_DIRECTIONS) / 3
        edge = sp.Matrix(edge)
        gradient = sp.Matrix(self.density_slope)
        flux = -self.diffusion * gradient
        if electrostatic:
            density_at_face = 2 + gradient.dot(edge) / 2
            field = sp.Matrix(self.potential_slope) - sp.Matrix(
                self.ext_efield)
            flux += -self.diffusion / self.kT * density_at_face \
                * self.valency * field
        return sp.N(flux.dot(edge) / edge.norm() / a0, 30)

    def test_the_flux_of_linear_fields_is_exact(self):
        for electrostatic in (False, True):
            stores = ek_sweeps.diffusive_flux_stores(DOUBLE, electrostatic)
            self.assertEqual(len(stores), 13)
            for edge, store in zip(flux_layout.STORED_DIRECTIONS, stores):
                with self.subTest(electrostatic=electrostatic, edge=edge):
                    self.assertAlmostEqual(
                        float(self.evaluate(store)),
                        float(self.expected(edge, electrostatic)),
                        places=13)

    def test_stores_are_to_the_components_of_the_flux_field(self):
        stores = ek_sweeps.diffusive_flux_stores(DOUBLE, False)
        for k, store in enumerate(stores):
            self.assertEqual(store.lhs.index, k)
            self.assertEqual(store.lhs.offsets, (0, 0, 0))

    def test_stencil_of_the_body_diagonals(self):
        # the two-point difference: the diffusive flux has the two cells of
        # the face (and the electrostatic flux the potential in these cells)
        # on the four body diagonals
        for electrostatic, n_accesses in ((False, 2), (True, 4)):
            stores = ek_sweeps.diffusive_flux_stores(DOUBLE, electrostatic)
            for edge, store in zip(flux_layout.STORED_DIRECTIONS, stores):
                if all(edge):
                    accesses = store.rhs.atoms(ps.Field.Access)
                    self.assertEqual(len(accesses), n_accesses, edge)
                    self.assertEqual(
                        {a.offsets for a in accesses},
                        {(0, 0, 0), edge}, edge)

    def test_stencils_do_not_depend_on_the_hash_seed(self):
        # pystencils chooses the stencil of the body diagonals at random;
        # derive only these (the last four stored directions)
        script = (
            "import pystencils as ps, ek_sweeps, flux_layout; "
            "dtype = ps.create_type('float64'); "
            "stores = ek_sweeps.discrete_flux("
            "ek_sweeps.diffusive_flux_expression(dtype, True), "
            "flux_layout.staggered_flux_field(dtype), "
            "flux_layout.STORED_DIRECTIONS[9:]); "
            "print([str(s.rhs) for s in stores])")
        outputs = []
        for seed in ("1", "2"):
            result = subprocess.run(
                [sys.executable, "-c", script], cwd=KERNELS, check=True,
                capture_output=True, text=True,
                env=dict(os.environ, PYTHONHASHSEED=seed))
            outputs.append(result.stdout.strip().splitlines()[-1])
        self.assertEqual(outputs[0], outputs[1])
        self.assertIn("diffusion", outputs[0])


class FluctuationNoiseTest(unittest.TestCase):

    def setUp(self):
        self.density = ek_sweeps.scalar_field("density", DOUBLE)
        self.diffusion = ps.TypedSymbol("diffusion", DOUBLE)
        self.block, self.noise = ek_sweeps.fluctuation_noise(
            DOUBLE, self.density, self.diffusion)

    def test_one_draw_for_two_faces(self):
        # 13 faces need 7 draws of 2 random numbers each
        self.assertEqual(len(self.noise), 13)
        self.assertEqual(sorted(s.name for s in self.block.exports),
                         [f"fluctuation_rng_{i}" for i in range(7)])
        for i, term in enumerate(self.noise):
            numbers = {(str(e.base), e.indices[0])
                       for e in term.atoms(sp.Indexed)}
            self.assertEqual(numbers, {(f"fluctuation_rng_{i // 2}", i % 2)})

    def test_the_amplitude_of_the_noise(self):
        # the uniform random number w has the mean 1 / 2 and the variance
        # 1 / 12, so that (w - 1 / 2) * sqrt(12) has the variance 1: the
        # noise has the variance 2 D rho_face / |e| / A0
        a0 = sum(sp.Matrix(e).norm() for e in
                 flux_layout.STORED_DIRECTIONS) / 3
        one_sigma = sp.Rational(1, 2) + 1 / sp.sqrt(12)
        for edge, term in zip(flux_layout.STORED_DIRECTIONS, self.noise):
            substitutions = {
                self.diffusion: sp.Rational(1, 5),
                self.density.center: sp.Rational(3, 2),
                self.density[tuple(edge)]: sp.Rational(3, 2)}
            amplitude = term.replace(
                lambda e: isinstance(e, sp.Indexed),
                lambda e: one_sigma).xreplace(substitutions)
            expected = sp.sqrt(2 * sp.Rational(1, 5) * sp.Rational(3, 2)
                               / sp.Matrix(edge).norm() / a0)
            self.assertAlmostEqual(float(amplitude), float(expected),
                                   places=12, msg=str(edge))

    def test_the_noise_depends_on_the_density_of_the_face(self):
        for edge, term in zip(flux_layout.STORED_DIRECTIONS, self.noise):
            accesses = term.atoms(ps.Field.Access)
            self.assertEqual(
                {a.offsets for a in accesses}, {(0, 0, 0), tuple(edge)})
            self.assertIn(self.diffusion, term.free_symbols)

    def test_the_counter_is_the_time_step(self):
        self.assertEqual({s.name for s in self.block.free_symbols},
                         {"time_step"})


class GeneratedSweepsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.double = ek_generated.generated("double")
        cls.single = ek_generated.generated("single")

    def read(self, kernel, precision="Double", target="", extension="h"):
        directory = self.double if precision == "Double" else self.single
        name = ek_generated.CLASS_NAMES[kernel]
        return (directory
                / f"{name}{precision}Precision{target}.{extension}"
                ).read_text()

    def constructor(self, text, name):
        match = re.search(rf"\n\s*{name}\s*\((.*?)\)\s*:\s*properties_",
                          text, re.S)
        self.assertIsNotNone(match, name)
        return " ".join(match.group(1).split())

    def test_files(self):
        names = {p.name for p in self.double.iterdir()}
        for kernel in ek_generated.SWEEPS:
            stem = f"{ek_generated.CLASS_NAMES[kernel]}DoublePrecision"
            for suffix in (".h", ".cpp", "CUDA.h", "CUDA.cu"):
                self.assertIn(stem + suffix, names)
        single = {p.name for p in self.single.iterdir()}
        self.assertIn("AdvectiveFluxKernelSinglePrecision.cpp", single)
        self.assertNotIn("ContinuityKernelSinglePrecision.cpp", single)

    def test_constructors(self):
        # the field IDs in alphabetical order, then the scalars in
        # alphabetical order: densityId, fluxId, ...
        ids = "walberla::BlockDataID"
        vec = "walberla::Vector3< int64_t >"
        expected = {
            "continuity": f"{ids} densityId, {ids} fluxId",
            "advective_flux":
                f"{ids} densityId, {ids} fluxId, {ids} velocityId",
            "diffusive_flux": f"{ids} densityId, {ids} fluxId, "
            "double diffusion",
            "diffusive_flux_thermalized":
                f"{ids} densityId, {ids} fluxId, {vec} block_offset, "
                f"double diffusion, {vec} grid_size, int64_t seed, "
                "int64_t time_step",
            "diffusive_flux_electrostatic":
                f"{ids} densityId, {ids} fluxId, {ids} potentialId, "
                "double diffusion, walberla::Vector3< double > ext_efield, "
                "double kT, double valency",
            "diffusive_flux_electrostatic_thermalized":
                f"{ids} densityId, {ids} fluxId, {ids} potentialId, "
                f"{vec} block_offset, double diffusion, "
                f"walberla::Vector3< double > ext_efield, {vec} grid_size, "
                "double kT, int64_t seed, int64_t time_step, double valency",
            "friction_coupling":
                f"{ids} fluxId, {ids} forceId, double diffusion, double kT, "
                "double lb_density",
        }
        for kernel, parameters in expected.items():
            name = f"{ek_generated.CLASS_NAMES[kernel]}DoublePrecision"
            for target in ("", "CUDA"):
                with self.subTest(kernel=kernel, target=target):
                    self.assertEqual(
                        self.constructor(self.read(kernel, target=target),
                                         name + target),
                        parameters)

    def test_single_precision_constructor(self):
        text = self.read("diffusive_flux_electrostatic_thermalized",
                         "Single")
        name = "DiffusiveFluxKernelWithElectrostaticThermalizedSinglePrecision"
        parameters = self.constructor(text, name)
        self.assertIn("float diffusion", parameters)
        self.assertIn("walberla::Vector3< float > ext_efield", parameters)
        self.assertIn("int64_t seed", parameters)
        self.assertNotIn("double", parameters)

    def test_properties_of_the_electrostatic_kernel_are_assignable(self):
        # the bridge updates the scalars of persistent kernel objects; the
        # field IDs are constant
        text = self.read("diffusive_flux_electrostatic_thermalized")
        for member in ("diffusion", "kT", "valency", "ext_efield", "seed",
                       "time_step", "block_offset", "grid_size"):
            self.assertRegex(text, rf"\n\s*[\w:<>, ]+&\s*{member}\(\)\s*\n?"
                             rf"\s*\{{", member)
            self.assertRegex(text, rf"const [\w:<>, ]+&\s*{member}\(\) "
                             "const", member)
        for member in ("densityId", "fluxId", "potentialId"):
            self.assertRegex(text, rf"const walberla::BlockDataID\s*&\s*"
                             rf"{member}\(\) const")
            self.assertNotRegex(text, rf"\n\s*walberla::BlockDataID\s*&\s*"
                                rf"{member}\(\)")

    def test_run_methods(self):
        for kernel in ek_generated.SWEEPS:
            for target in ("", "CUDA"):
                text = self.read(kernel, target=target)
                self.assertRegex(
                    text, r"void operator\(\)\(walberla::IBlock \* block")
                self.assertIn("runOnCellInterval(", text)
                self.assertNotIn("SparseIndexList", text)
        for text in (self.read("advective_flux", "Single", "CUDA"),):
            self.assertIn("exectag::GPU", text)

    def test_only_the_thermalized_kernels_have_an_rng(self):
        for kernel in ek_generated.SWEEPS:
            thermalized = "thermalized" in kernel
            cpu = self.read(kernel, extension="cpp")
            cuda = self.read(kernel, target="CUDA", extension="cu")
            header = self.read(kernel)
            self.assertEqual("philox_float64" in cpu, thermalized, kernel)
            self.assertEqual("curand" in cuda or "philox" in cuda,
                             thermalized, kernel)
            for name in ("seed", "time_step", "block_offset", "grid_size"):
                self.assertEqual(f" {name}_;" in header, thermalized,
                                 (kernel, name))

    def test_electrostatic_properties(self):
        for kernel in ek_generated.SWEEPS:
            electrostatic = "electrostatic" in kernel
            header = self.read(kernel)
            self.assertEqual("potentialId" in header, electrostatic, kernel)
            self.assertEqual("ext_efield" in header, electrostatic, kernel)
            self.assertEqual("valency" in header, electrostatic, kernel)

    def test_rng_coordinates_are_periodic(self):
        # the counters are the global coordinates (offset by the grid size
        # to make them positive) modulo the size of the grid
        cpu = self.read("diffusive_flux_thermalized", extension="cpp")
        self.assertEqual(len(re.findall(
            r"\(block_offset_0 \+ ctr_0 \+ grid_size_0\) % grid_size_0",
            cpu)), 7)
        # 13 faces need 7 draws of 2 numbers
        self.assertEqual(len(re.findall(
            r"philox_float64\s*<\s*int64_t\s*>\s*\(", cpu)), 7)
        single = self.read("diffusive_flux_electrostatic_thermalized",
                           "Single", extension="cpp")
        self.assertEqual(len(re.findall(
            r"philox_float32\s*<\s*int64_t\s*>\s*\(", single)), 4)
        self.assertNotIn("philox_float64", single)

    def test_flux_stores_are_guarded(self):
        # one condition for each of the 13 components, each guarding the
        # store to the component
        for kernel in ("diffusive_flux", "advective_flux"):
            for target, extension in (("", "cpp"), ("CUDA", "cu")):
                text = self.read(kernel, target=target, extension=extension)
                self.assertEqual(
                    len(re.findall(r"ctr_0 >= \d+LL", text)), 13, kernel)
                self.assertEqual(
                    len(re.findall(r"_data_flux\[[^\]]*\] = ", text)), 13,
                    kernel)
                for k in range(1, 13):
                    self.assertIn(f"{k}LL * _stride_flux_3", text)

    def test_the_sweeps_run_on_the_cells_of_the_loop(self):
        for kernel in ("continuity", "friction_coupling"):
            cpu = self.read(kernel, extension="cpp")
            # no condition on the position: the sweep covers the interval
            self.assertNotIn("ctr_0 >=", cpu)
            self.assertIn("#pragma omp", cpu)

    def test_continuity_sums_the_net_flux_of_opposite_faces(self):
        # a uniform flux has exactly zero divergence in any precision
        for precision, target, extension in (("Single", "", "cpp"),
                                             ("Double", "CUDA", "cu")):
            if precision == "Single":
                text = (ek_generated.generated("double")
                        / "ContinuityKernelDoublePrecision.cpp").read_text()
            else:
                text = self.read("continuity", precision, target, extension)
            for k in range(13):
                self.assertRegex(
                    text, rf"net_flux_{k} = -1\.0f? \* _data_flux\[[^;]*\+ "
                    rf"_data_flux\[[^;]*;")
            self.assertEqual(len(re.findall(r"_data_flux\[", text)), 26)

    def test_continuity_updates_the_density_from_26_faces(self):
        cpu = self.read("continuity", extension="cpp")
        # 13 faces of the cell and 13 of its neighbors
        self.assertEqual(len(re.findall(r"_data_flux\[", cpu)), 26)
        self.assertEqual(len(re.findall(r"_data_density\[[^\]]*\] = ", cpu)),
                         1)

    def test_friction_coupling_updates_the_force(self):
        cpu = self.read("friction_coupling", extension="cpp")
        self.assertEqual(len(re.findall(r"_data_force\[[^\]]*\] = ", cpu)), 3)

    def test_precisions_do_not_mix(self):
        for kernel in ek_generated.SINGLE_PRECISION:
            for target, extension in (("", "cpp"), ("CUDA", "cu")):
                single = self.read(kernel, "Single", target, extension)
                self.assertNotIn("double", single, (kernel, target))
                # the global functions of the C math library are the double
                # versions
                self.assertNotRegex(
                    single, r"(?<![:\w.>])(sqrt|fabs|fmin|fmax|pow)\(")
        for kernel in ek_generated.SWEEPS:
            double = self.read(kernel, "Double", extension="cpp")
            self.assertIsNone(re.search(r"\bfloat\b|\d\.\d*f\b", double),
                              kernel)

    def test_cuda_kernels(self):
        for kernel in ek_generated.SWEEPS:
            cuda = self.read(kernel, target="CUDA", extension="cu")
            self.assertIn("__global__", cuda)
            self.assertNotIn("#pragma omp", cuda)

    def test_header_of_the_generator(self):
        for kernel in ek_generated.SWEEPS:
            text = self.read(kernel)
            self.assertIn("Generated by "
                          "maintainer/walberla_kernels/generate_ek_kernels.py",
                          text)
            self.assertIn("namespace walberla::pystencils {", text)


class SweepRegistryTest(unittest.TestCase):

    def test_registered_sweeps(self):
        kernels = generate_ek_kernels.KERNELS
        expected = {
            "advective_flux": "AdvectiveFluxKernel",
            "continuity": "ContinuityKernel",
            "diffusive_flux": "DiffusiveFluxKernel",
            "diffusive_flux_thermalized": "DiffusiveFluxKernelThermalized",
            "diffusive_flux_electrostatic":
                "DiffusiveFluxKernelWithElectrostatic",
            "diffusive_flux_electrostatic_thermalized":
                "DiffusiveFluxKernelWithElectrostaticThermalized",
            "friction_coupling": "FrictionCouplingKernel",
        }
        for key, stem in expected.items():
            kernel = kernels[key]
            self.assertIsInstance(kernel, sfg_driver.Kernel)
            self.assertEqual(kernel.stem, stem)
            self.assertEqual(kernel.targets, ("cpu", "cuda"))
            self.assertFalse(kernel.header_only)
        self.assertEqual(sorted(expected), sorted(ek_generated.SWEEPS))
        ctx = sfg_driver.CodegenContext("single", "cuda")
        self.assertEqual(
            ctx.class_name(kernels["advective_flux"].stem),
            "AdvectiveFluxKernelSinglePrecisionCUDA")


if __name__ == "__main__":
    unittest.main()
