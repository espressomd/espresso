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

import re
import unittest

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


def make_equations(precision="double", lb_method=None):
    ctx = sfg_driver.CodegenContext(precision, "cpu")
    if lb_method is None:
        fields = pystencils_espresso.generate_fields(
            generate_lb_kernels.STENCIL, ctx.data_type)
        lb_method = generate_lb_kernels.create_lb_method(fields)
    return accessor_equations.AccessorEquations(lb_method, ctx.data_type)


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

    def test_delta_rho(self):
        eqs = self.eqs.delta_rho
        self.assertEqual(eqs.inputs, {"rho"})
        self.assertEqual(eqs.outputs, ["delta_rho"])
        (lhs, rhs), = eqs.assignments
        self.assertEqual(rhs, sp.Symbol("rho") - 1)

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
        stencil = generate_lb_kernels.STENCIL
        self.assertEqual(len(offsets), 19)
        for offset, direction in zip(offsets, stencil):
            self.assertEqual(offset, tuple(-c for c in direction))

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


class EquationPrintingTest(unittest.TestCase):

    def test_assignments_are_const_declarations(self):
        eqs = make_equations("double")
        lines = eqs.density.lines()
        self.assertEqual(len(lines), len(eqs.density.assignments))
        for line, name in zip(lines, eqs.density.outputs):
            self.assertTrue(line.startswith(f"const double {name} = "), line)
            self.assertTrue(line.endswith(";"), line)

    def test_single_precision_has_no_double_literals(self):
        eqs = make_equations("single")
        lines = []
        for name in ("density", "density_scaled", "momentum_density",
                     "momentum_density_force_setter",
                     "unshifted_momentum_density", "second_moment",
                     "velocity_setter", "delta_rho"):
            lines += getattr(eqs, name).lines()
        lines += [accessor_equations.print_expression(e, eqs.dtype)
                  for e in eqs.equilibrium]
        text = "\n".join(lines)
        self.assertNotIn("double", text)
        self.assertNotIn("pow", text)
        self.assertTrue(text.startswith("const float "))
        # every floating-point literal carries the float suffix
        literals = re.findall(r"(?<![\w.])\d+\.\d*(?:e[+-]?\d+)?f?(?![\w.])",
                              text)
        self.assertTrue(literals)
        self.assertTrue(all(lit.endswith("f") for lit in literals), literals)

    def test_double_precision_has_no_float_literals(self):
        eqs = make_equations("double")
        text = "\n".join(eqs.momentum_density.lines())
        self.assertNotRegex(text, r"\d\.\d*f\b")
        self.assertIn("0.5 *", text)

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

    def test_force_setter_divides_force_by_density(self):
        eqs = make_equations("single")
        plain = "\n".join(eqs.momentum_density.lines())
        setter = "\n".join(eqs.momentum_density_force_setter.lines())
        for i in range(3):
            self.assertRegex(
                setter, rf"const float md_{i} = momdensity_{i} \+ "
                rf"0\.5f \* force_{i} \* \(1\.0f / density\);")
        self.assertNotRegex(plain, r"\bdensity\b")


class UnsupportedMethodTest(unittest.TestCase):

    def create_method(self, **kwargs):
        ctx = sfg_driver.CodegenContext("double", "cpu")
        fields = pystencils_espresso.generate_fields(
            generate_lb_kernels.STENCIL, ctx.data_type)
        config = dict(generate_lb_kernels.LBM_CONFIG_KWARGS)
        config.update(kwargs)
        return lbmpy.creationfunctions.create_mrt_orthogonal(
            stencil=generate_lb_kernels.STENCIL,
            relaxation_rates=relaxation_rates.rr_getter,
            force_model=lbmpy.forcemodels.Schiller(
                fields["force"].center_vector),
            **config)

    def test_incompressible_method_is_rejected(self):
        method = self.create_method(compressible=False)
        with self.assertRaises(NotImplementedError):
            make_equations(lb_method=method)

    def test_pdfs_that_are_not_zero_centered_are_rejected(self):
        method = self.create_method(zero_centered=False,
                                    delta_equilibrium=False)
        with self.assertRaises(NotImplementedError):
            make_equations(lb_method=method)


if __name__ == "__main__":
    unittest.main()
