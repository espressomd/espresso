#
# Copyright (C) 2021-2023 The ESPResSo project
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

import os
import sympy as sp
import pystencils as ps
import lbmpy_walberla
from pystencils.typing.typed_sympy import TypedSymbol
from pystencils.typing import BasicType, CastFunc, TypedSymbol

# File derived from lbmpy_walberla.walberla_lbm_generation in the
# walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d


def __type_equilibrium_assignments(assignments, config, subs_dict):
    # Function derived from lbmpy_walberla.walberla_lbm_generation.__type_equilibrium_assignments()
    # in the walberla project, commit 9dcd0dd90f50f7b64b0a38bb06327854463fdafd
    from pystencils.node_collection import NodeCollection
    from pystencils.typing.transformations import add_types
    result = assignments.new_with_substitutions(subs_dict)
    result = NodeCollection(result.main_assignments)
    result.evaluate_terms()
    result = add_types(result.all_assignments, config)
    return result


def type_expr(eq, dtype):
    # manually cast floats to dtype since this is not done automatically
    repl = ((rational := sp.Rational(1, i), CastFunc(rational, dtype))
            for i in (2, 3, 4, 6, 8, 9, 12, 24, 18, 36, 72))
    eq = eq.subs(repl)
    return eq.subs({s: TypedSymbol(s.name, dtype)
                   for s in eq.atoms(sp.Symbol)})


def pow_to_mul(eq):
    keep_processing = True
    while keep_processing:
        for expr in sp.preorder_traversal(eq):
            if expr.is_Pow:
                if expr.args[0].is_Symbol and expr.args[1].is_Integer:
                    power = expr.args[1].p
                    if power >= 1:
                        chained_product = expr.args[1].p * [expr.args[0]]
                        expr_mul = sp.Mul(*chained_product, evaluate=False)
                        print(f"folding '{expr}' to '{expr_mul}'")
                        eq = eq.subs(expr, sp.UnevaluatedExpr(expr_mul))
                        break
        else:
            keep_processing = False
    return eq


def make_velocity_getters(cqc, rho_sym, vel_arr_symbols):
    velocity_getter = cqc.equilibrium_input_equations_from_init_values(
        rho_sym, vel_arr_symbols)
    eq = velocity_getter.main_assignments.pop(0)
    assert eq.lhs == rho_sym and eq.rhs == rho_sym
    eq = velocity_getter.main_assignments.pop(0)
    assert eq.lhs.name == f"delta_{rho_sym.name}"
    return velocity_getter


def equations_to_code(equations, variable_prefix="",
                      variables_without_prefix=None, dtype=None, backend=None):
    if dtype is None:
        dtype = BasicType("float64")

    if variables_without_prefix is None:
        variables_without_prefix = []
    if isinstance(equations, ps.AssignmentCollection):
        equations = equations.all_assignments

    variables_without_prefix = list(variables_without_prefix)

    result = []
    left_hand_side_names = [eq.lhs.name for eq in equations]
    for eq in equations:
        lhs, rhs = eq.lhs, eq.rhs
        rhs = lbmpy_walberla.walberla_lbm_generation.field_and_symbol_substitute(
            rhs, variable_prefix, variables_without_prefix + left_hand_side_names)
        lhs = type_expr(lhs, dtype=dtype)
        rhs = type_expr(rhs, dtype=dtype)
        rhs = pow_to_mul(rhs)
        assignment = ps.astnodes.SympyAssignment(lhs, rhs)
        result.append(backend(assignment))
    return "\n".join(result)


def substitute_force_getter_cpp(code):
    field_getter = "force->"
    assert field_getter in code is not None, f"pattern '{field_getter} not found in '''\n{code}\n'''"
    return code.replace(field_getter, "force_field->")


def add_espresso_filters_to_jinja_env(jinja_env):
    jinja_env.filters["substitute_force_getter_cpp"] = substitute_force_getter_cpp


def generate_macroscopic_values_accessors(ctx, config, lb_method, templates):

    # Function derived from lbmpy_walberla.walberla_lbm_generation.__lattice_model()
    # in the walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d
    # with backports from commit de6b00071233a9a1f45d7a6773988363e058f1a0

    from jinja2 import Environment, FileSystemLoader, StrictUndefined
    from sympy.tensor import IndexedBase
    from pystencils.backends.cbackend import CustomSympyPrinter
    from pystencils.backends.cbackend import CBackend
    from pystencils.backends.cuda_backend import CudaBackend
    from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
    from lbmpy_walberla.walberla_lbm_generation import stencil_switch_statement

    cpp_printer = CustomSympyPrinter()
    stencil_name = lb_method.stencil.name
    if not stencil_name:
        raise ValueError(
            "lb_method uses a stencil that is not supported in waLBerla")

    default_dtype = config.data_type.default_factory()
    if config.target == ps.Target.GPU:
        backend = CudaBackend()
    else:
        backend = CBackend()
    kwargs = {
        "backend": backend,
        "variable_prefix": "",
        "dtype": default_dtype}

    cqc = lb_method.conserved_quantity_computation
    vel_symbols = cqc.velocity_symbols
    rho_sym = sp.Symbol("rho")
    pdfs_sym = sp.symbols(f"f_:{lb_method.stencil.Q}")
    vel_arr_symbols = [
        IndexedBase(TypedSymbol("u", dtype=default_dtype), shape=(1,))[i]
        for i in range(len(vel_symbols))]
    momentum_density_symbols = sp.symbols(f"md_:{len(vel_symbols)}")
    second_momentum_symbols = sp.symbols(f"p_:{len(vel_symbols)**2}")

    equilibrium_subs_dict = dict(zip(vel_symbols, vel_arr_symbols))
    equilibrium = lb_method.get_equilibrium()
    lhs_list = [a.lhs for a in equilibrium.main_assignments]
    equilibrium_matrix = sp.Matrix(
        [e.rhs for e in equilibrium.main_assignments])
    equilibrium = ps.AssignmentCollection([ps.Assignment(lhs, rhs)
                                           for lhs, rhs in zip(lhs_list, equilibrium_matrix)])
    equilibrium = __type_equilibrium_assignments(
        equilibrium, config, equilibrium_subs_dict)

    velocity_getters = make_velocity_getters(cqc, rho_sym, vel_arr_symbols)
    density_velocity_setter_macroscopic_values = equations_to_code(
        velocity_getters, variables_without_prefix=["rho", "u"], **kwargs)
    momentum_density_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {"density": rho_sym, "momentum_density": momentum_density_symbols})
    unshifted_momentum_density_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {"density": rho_sym, "momentum_density": momentum_density_symbols})
    for i, eq in reversed(
            list(enumerate(unshifted_momentum_density_getter.main_assignments))):
        if eq.lhs.name.startswith("md_"):
            del unshifted_momentum_density_getter.main_assignments[i]
    second_momentum_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {"moment2": second_momentum_symbols})

    jinja_context = {
        "stencil_name": stencil_name,
        "D": lb_method.stencil.D,
        "Q": lb_method.stencil.Q,
        "compressible": cqc.compressible,
        "zero_centered": cqc.zero_centered_pdfs,
        "dtype": default_dtype,

        "equilibrium_from_direction": stencil_switch_statement(lb_method.stencil, equilibrium),
        "equilibrium": [cpp_printer.doprint(e.rhs) for e in equilibrium],

        "density_getters": equations_to_code(
            cqc.output_equations_from_pdfs(pdfs_sym, {"density": rho_sym}),
            variables_without_prefix=[e.name for e in pdfs_sym], **kwargs),
        "momentum_density_getter": equations_to_code(
            momentum_density_getter, variables_without_prefix=pdfs_sym, **kwargs),
        "second_momentum_getter": equations_to_code(
            second_momentum_getter, variables_without_prefix=pdfs_sym, **kwargs),
        "density_velocity_setter_macroscopic_values": density_velocity_setter_macroscopic_values,
        "unshifted_momentum_density_getter": equations_to_code(unshifted_momentum_density_getter, variables_without_prefix=pdfs_sym, **kwargs),

        "namespace": "lbm",
    }

    env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)),
                      undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)
    add_espresso_filters_to_jinja_env(env)

    for filename, template in templates.items():
        source = env.get_template(template).render(**jinja_context)
        ctx.write_file(filename, source)
