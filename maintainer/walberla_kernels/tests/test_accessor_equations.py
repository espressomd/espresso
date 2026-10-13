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

import dataclasses
import inspect
import re
import unittest

import numpy
import pystencils as ps
import sympy as sp

import accessor_equations
import generate_lb_kernels
import lbmpy
import pystencils_espresso
import relaxation_rates
import sfg_driver

PDFS = {f"f_{i}" for i in range(19)}
FORCES = {f"force_{i}" for i in range(3)}


def make_equations(lb_method=None):
    if lb_method is None:
        fields = pystencils_espresso.generate_fields(
            generate_lb_kernels.STENCIL, ps.create_type("float64"))
        lb_method = generate_lb_kernels.create_lb_method(fields)
    return accessor_equations.AccessorEquations(lb_method)


class EquationInputsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.eqs = make_equations()

    def test_dimensions(self):
        self.assertEqual((self.eqs.D, self.eqs.Q), (3, 19))
        self.assertEqual(set(self.eqs.pdf_names), PDFS)
        self.assertEqual(self.eqs.force_names,
                         ["force_0", "force_1", "force_2"])

    def test_density(self):
        eqs = self.eqs.density
        self.assertEqual(eqs.inputs, PDFS)
        self.assertIn("rho", eqs.outputs)
        self.assertIn("delta_rho", eqs.outputs)

    def test_density_scaled(self):
        eqs = self.eqs.density_scaled
        self.assertEqual(eqs.inputs, PDFS | {"density"})
        rho = dict((lhs.name, rhs) for lhs, rhs in eqs.assignments)["rho"]
        self.assertEqual(rho, (1 + sp.Symbol("delta_rho"))
                         * sp.Symbol("density"))

    def test_momentum_density(self):
        eqs = self.eqs.momentum_density
        self.assertEqual(eqs.inputs, PDFS | FORCES)
        self.assertTrue({"rho", "md_0", "md_1", "md_2"} <= set(eqs.outputs))
        self.assertEqual(self.eqs.used_pdfs(eqs), list(range(19)))
        self.assertEqual(self.eqs.used_forces(eqs), [0, 1, 2])

    def test_momentum_density_force_setter(self):
        eqs = self.eqs.momentum_density_force_setter
        self.assertEqual(eqs.inputs, PDFS | FORCES | {"density"})

    def test_unshifted_momentum_density(self):
        eqs = self.eqs.unshifted_momentum_density
        self.assertEqual(eqs.inputs, PDFS)
        self.assertTrue({"rho", "md_0", "md_1", "md_2"} <= set(eqs.outputs))
        self.assertEqual(self.eqs.used_forces(eqs), [])

    def test_reduce_pruned_momentum_density(self):
        eqs = self.eqs.momentum_density.pruned(["md_0", "md_1", "md_2"])
        self.assertEqual(eqs.inputs, (PDFS - {"f_0"}) | FORCES)
        self.assertEqual(self.eqs.used_pdfs(eqs), list(range(1, 19)))
        for name in ("rho", "delta_rho"):
            self.assertNotIn(name, eqs.outputs)
            self.assertNotIn(name, eqs.inputs)
        self.assertEqual(eqs.outputs[-3:], ["md_0", "md_1", "md_2"])

    def test_second_moment(self):
        eqs = self.eqs.second_moment
        self.assertEqual(eqs.inputs, PDFS - {"f_0"})
        self.assertEqual(eqs.outputs, [f"p_{i}" for i in range(9)])

    def test_velocity_setter(self):
        eqs = self.eqs.velocity_setter
        self.assertEqual(eqs.inputs,
                         FORCES | {"rho", "u_in_0", "u_in_1", "u_in_2"})
        self.assertEqual(eqs.outputs, ["u_0", "u_1", "u_2"])

    def test_equilibrium(self):
        allowed = {"rho", "delta_rho", "u_0", "u_1", "u_2"}
        exprs = self.eqs.equilibrium
        self.assertEqual(len(exprs), 19)
        for expr in exprs:
            self.assertLessEqual({s.name for s in expr.free_symbols}, allowed)

    def test_equilibrium_from_direction(self):
        pairs = self.eqs.equilibrium_from_direction
        self.assertEqual(len(pairs), 19)
        self.assertEqual(pairs[0][0], "C")
        self.assertEqual(len({name for name, _ in pairs}), 19)
        for _, expr in pairs:
            self.assertNotIn("delta_rho", {s.name for s in expr.free_symbols})

    def test_inverse_neighbor_offsets(self):
        offsets = self.eqs.inverse_neighbor_offsets
        names = [name for name, _ in self.eqs.equilibrium_from_direction]
        self.assertEqual(len(offsets), 19)
        known = {"C": (0, 0, 0), "N": (0, -1, 0), "S": (0, 1, 0),
                 "E": (-1, 0, 0), "W": (1, 0, 0), "T": (0, 0, -1),
                 "B": (0, 0, 1), "NE": (-1, -1, 0), "BS": (0, 1, 1)}
        for name, offset in known.items():
            self.assertEqual(offsets[names.index(name)], offset, name)

    def test_pruned_rejects_unknown_targets(self):
        with self.assertRaises(ValueError):
            self.eqs.momentum_density.pruned(["md_7"])

    def test_pruned_keeps_order_and_is_closed(self):
        full = self.eqs.momentum_density
        pruned = full.pruned(["md_1"])
        names = [lhs.name for lhs, _ in pruned.assignments]
        self.assertIn("md_1", names)
        self.assertNotIn("md_0", names)
        original = [lhs.name for lhs, _ in full.assignments]
        self.assertEqual(names, [n for n in original if n in names])
        defined = set()
        for lhs, rhs in pruned.assignments:
            for symbol in rhs.free_symbols:
                self.assertTrue(symbol.name in pruned.inputs
                                or symbol.name in defined, symbol.name)
            defined.add(lhs.name)


def evaluate(eqs, values):
    """Evaluate an equation set sequentially, return all symbol values."""
    values = dict(values)
    for lhs, rhs in eqs.assignments:
        values[lhs.name] = float(rhs.subs(
            {sp.Symbol(k): v for k, v in values.items()}))
    return values


class NumericTest(unittest.TestCase):

    def test_equations_match_lbmpy_numerically(self):
        eqs = make_equations()
        rng = numpy.random.default_rng(42)
        stencil = generate_lb_kernels.STENCIL
        pdfs = rng.uniform(-0.01, 0.01, 19)
        force = rng.uniform(-0.01, 0.01, 3)
        density = float(rng.uniform(0.5, 2.0))
        u_in = rng.uniform(-0.05, 0.05, 3)
        base = {f"f_{i}": pdfs[i] for i in range(19)}
        base.update({f"force_{i}": force[i] for i in range(3)})
        base["density"] = density
        base.update({f"u_in_{i}": u_in[i] for i in range(3)})

        rho = 1.0 + pdfs.sum()
        momentum = numpy.array(stencil).T @ pdfs

        res = evaluate(eqs.density, base)
        self.assertAlmostEqual(res["rho"], rho, places=12)
        res = evaluate(eqs.density_scaled, base)
        self.assertAlmostEqual(res["rho"], rho * density, places=12)
        res = evaluate(eqs.momentum_density, base)
        self.assertAlmostEqual(res["rho"], rho, places=12)
        for i in range(3):
            self.assertAlmostEqual(
                res[f"md_{i}"], momentum[i] + force[i] / 2, places=12)
        res = evaluate(eqs.momentum_density_force_setter, base)
        for i in range(3):
            self.assertAlmostEqual(
                res[f"md_{i}"], momentum[i] + force[i] / density / 2,
                places=12)
        res = evaluate(eqs.unshifted_momentum_density, base)
        for i in range(3):
            self.assertAlmostEqual(res[f"md_{i}"], momentum[i], places=12)
        res = evaluate(eqs.second_moment, base)
        for i in range(3):
            for j in range(3):
                # lbmpy returns the moment of the full populations
                # f_q + w_q, the weights add delta_ij / 3 (lattice sound
                # speed squared)
                reference = sum(
                    stencil[q][i] * stencil[q][j] * pdfs[q]
                    for q in range(19)) + (1 / 3 if i == j else 0)
                self.assertAlmostEqual(res[f"p_{3 * i + j}"], reference,
                                       places=12)

        # velocity setter: requested velocity shifted by half the force
        res = evaluate(eqs.velocity_setter, {**base, "rho": rho})
        for i in range(3):
            self.assertAlmostEqual(
                res[f"u_{i}"], u_in[i] - force[i] / (2 * rho), places=12)

        # equilibrium
        u = rng.uniform(-0.05, 0.05, 3)
        subs = {"rho": rho, "delta_rho": rho - 1,
                **{f"u_{i}": u[i] for i in range(3)}}
        sym = {sp.Symbol(k): v for k, v in subs.items()}
        reference = [float(a.rhs.subs(sym))
                     for a in eqs.lb_method.get_equilibrium().main_assignments]
        values = [float(e.subs(sym)) for e in eqs.equilibrium]
        numpy.testing.assert_allclose(values, reference, atol=1e-14)
        self.assertAlmostEqual(sum(values), rho - 1, places=12)
        numpy.testing.assert_allclose(
            numpy.array(stencil).T @ numpy.array(values), rho * u, atol=1e-12)
        full = [float(e.subs(sym)) for _, e in eqs.equilibrium_from_direction]
        self.assertAlmostEqual(sum(full), rho, places=12)


class PublicApiTest(unittest.TestCase):
    """The equations are sympy expressions: only the expressions that the
    CPU backend switches over are printed."""

    def test_only_expressions_are_printed(self):
        self.assertEqual(accessor_equations.__all__, [
            "EquationSet", "AccessorEquations", "print_expression"])
        self.assertFalse(hasattr(accessor_equations, "print_assignments"))
        self.assertFalse(hasattr(accessor_equations.EquationSet, "lines"))
        self.assertFalse(hasattr(accessor_equations.AccessorEquations,
                                 "delta_rho"))

    def test_equations_do_not_depend_on_the_precision(self):
        self.assertEqual(
            [f.name for f in dataclasses.fields(
                accessor_equations.EquationSet)], ["assignments"])
        self.assertEqual(
            list(inspect.signature(
                accessor_equations.AccessorEquations).parameters),
            ["lb_method"])


class PrintExpressionTest(unittest.TestCase):

    def test_single_precision_has_no_double_literals(self):
        eqs = make_equations()
        single = ps.create_type("float32")
        text = "\n".join(
            accessor_equations.print_expression(e, single)
            for e in [*eqs.equilibrium,
                      *(e for _, e in eqs.equilibrium_from_direction)])
        self.assertNotIn("double", text)
        self.assertNotIn("pow", text)
        # every floating-point literal carries the float suffix
        literals = re.findall(r"(?<![\w.])\d+\.\d*(?:e[+-]?\d+)?f?(?![\w.])",
                              text)
        self.assertTrue(literals)
        self.assertTrue(all(lit.endswith("f") for lit in literals), literals)

    def test_double_precision_has_no_float_literals(self):
        eqs = make_equations()
        text = "\n".join(
            accessor_equations.print_expression(e, ps.create_type("float64"))
            for e in eqs.equilibrium)
        self.assertNotRegex(text, r"\d\.\d*f\b")
        self.assertRegex(text, r"\d\.\d")

    def test_integer_powers_are_expanded(self):
        expr = sp.Symbol("u_0") ** 2
        printed = accessor_equations.print_expression(
            expr, ps.create_type("float32"))
        self.assertEqual(printed, "u_0 * u_0")

    def test_rational_is_folded_in_dtype(self):
        expr = sp.Rational(1, 2) * sp.Symbol("rho")
        self.assertEqual(accessor_equations.print_expression(
            expr, ps.create_type("float32")), "0.5f * rho")
        self.assertEqual(accessor_equations.print_expression(
            expr, ps.create_type("float64")), "0.5 * rho")


class ForceSetterTest(unittest.TestCase):

    def test_force_setter_divides_force_by_density(self):
        eqs = make_equations()
        self.assertNotIn("density", eqs.momentum_density.inputs)
        setter = eqs.momentum_density_force_setter
        self.assertIn("density", setter.inputs)
        plain = dict(eqs.momentum_density.assignments)
        shifted = dict(setter.assignments)
        for i in range(3):
            md, force = sp.Symbol(f"md_{i}"), sp.Symbol(f"force_{i}")
            self.assertEqual(
                sp.simplify(shifted[md] - plain[md]
                            - force / 2 * (1 / sp.Symbol("density") - 1)), 0)


class UnsupportedMethodTest(unittest.TestCase):

    def create_method(self, force=None, **kwargs):
        ctx = sfg_driver.CodegenContext("double", "cpu")
        fields = pystencils_espresso.generate_fields(
            generate_lb_kernels.STENCIL, ctx.data_type)
        config = dict(generate_lb_kernels.LBM_CONFIG_KWARGS)
        config.update(kwargs)
        return lbmpy.creationfunctions.create_mrt_orthogonal(
            stencil=generate_lb_kernels.STENCIL,
            relaxation_rates=relaxation_rates.rr_getter,
            force_model=lbmpy.forcemodels.Schiller(
                fields["force"].center_vector if force is None else force),
            **config)

    def test_incompressible_method_is_rejected(self):
        method = self.create_method(compressible=False)
        with self.assertRaises(NotImplementedError):
            make_equations(method)

    def test_pdfs_that_are_not_zero_centered_are_rejected(self):
        method = self.create_method(zero_centered=False,
                                    delta_equilibrium=False)
        with self.assertRaises(NotImplementedError):
            make_equations(method)

    def test_missing_force_model_is_rejected(self):
        config = dict(generate_lb_kernels.LBM_CONFIG_KWARGS)
        method = lbmpy.creationfunctions.create_mrt_orthogonal(
            stencil=generate_lb_kernels.STENCIL,
            relaxation_rates=relaxation_rates.rr_getter, **config)
        self.assertIsNone(method.force_model)
        with self.assertRaisesRegex(NotImplementedError, "force model"):
            accessor_equations.AccessorEquations(method)

    def test_force_that_is_not_a_field_access_is_rejected(self):
        method = self.create_method(force=(0.1, 0.0, 0.0))
        with self.assertRaisesRegex(NotImplementedError, "field accesses"):
            make_equations(method)


if __name__ == "__main__":
    unittest.main()
