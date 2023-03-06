#
# Copyright (C) 2020-2022 The ESPResSo project
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
import setuptools

import numpy as np
import sympy as sp

import pystencils as ps
import pystencils.astnodes

from pystencils.typing import BasicType, CastFunc, TypedSymbol
import pystencils_walberla
import pystencils_espresso

import lbmpy
import lbmpy.creationfunctions
import lbmpy.forcemodels
import lbmpy.stencils
import lbmpy.enums

import lbmpy_walberla
import lbmpy_espresso

import relaxation_rates
import code_generation_context


target = ps.Target.GPU

# Make sure we have the correct versions of the required dependencies
SpecifierSet = setuptools.version.pkg_resources.packaging.specifiers.SpecifierSet
for module, requirement in [(ps, "==1.1.1"), (lbmpy, "==1.1.1")]:
    assert SpecifierSet(requirement).contains(module.__version__), \
        f"{module.__name__} version {module.__version__} doesn't match requirement {requirement}"





# taken form walberla_lbm_generation.py with modifications
def field_and_symbol_substitute(expr, variable_prefix="lm.", variables_without_prefix=None):
    if variables_without_prefix is None:
        variables_without_prefix = []
    variables_without_prefix = [v.name if isinstance(v, sp.Symbol) else v for v in variables_without_prefix]
    substitutions = {}
    # check for member access
    if variable_prefix.endswith("."):
        postfix = "_"
    else:
        postfix = ""
    for sym in expr.atoms(sp.Symbol):
        if isinstance(sym, ps.Field.Access):
            fa = sym
            prefix = "" if fa.field.name in variables_without_prefix else variable_prefix
            if prefix.endswith("."):
                postfix2 = "_"
            else:
                postfix2 = ""
            if fa.field.index_dimensions == 0:
                substitutions[fa] = sp.Symbol(f"{prefix}{fa.field.name + postfix2}->get(x,y,z)")
            else:
                assert fa.field.index_dimensions == 1, "walberla supports only 0 or 1 index dimensions"
                substitutions[fa] = sp.Symbol(f"{prefix}{fa.field.name + postfix2}->get(x,y,z,{fa.index[0]})")
        else:
            if sym.name not in variables_without_prefix:
                substitutions[sym] = sp.Symbol(variable_prefix + sym.name + postfix)
    return expr.subs(substitutions)


def type_expr(eq, dtype):
    # manually cast floats to dtype since this is somehow not done automatically
    repl = ((rational := sp.Rational(1, i), CastFunc(rational, dtype))
            for i in (2, 3, 4, 6, 8, 9, 12, 24, 18, 36, 72))
    eq = eq.subs(repl)
    return eq.subs({s: TypedSymbol(s.name, dtype) for s in eq.atoms(sp.Symbol)})


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


def equations_to_code(equations, variable_prefix="", variables_without_prefix=None, dtype=None, backend=None):
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
        rhs = field_and_symbol_substitute(
            rhs, variable_prefix, variables_without_prefix + left_hand_side_names)
        lhs = type_expr(lhs, dtype=dtype)
        rhs = type_expr(rhs, dtype=dtype)
        rhs = pow_to_mul(rhs)
        assignment = ps.astnodes.SympyAssignment(lhs, rhs)
        result.append(backend(assignment))
    return "\n".join(result)


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


def generate_macroscopic_values_accessors(
        gen_context, config, lb_method, filename):

    from jinja2 import Environment, FileSystemLoader, StrictUndefined
    from sympy.tensor import IndexedBase
    from pystencils.backends.cbackend import CustomSympyPrinter
    from pystencils.backends.cuda_backend import CudaBackend
    from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
    from lbmpy_walberla.walberla_lbm_generation import stencil_switch_statement

    cpp_printer = CustomSympyPrinter()
    backend = CudaBackend()
    stencil_name = lb_method.stencil.name
    if not stencil_name:
        raise ValueError(
            "lb_method uses a stencil that is not supported in waLBerla")

    default_dtype = config.data_type.default_factory()

    cqc = lb_method.conserved_quantity_computation
    vel_symbols = cqc.velocity_symbols
    rho_sym = sp.Symbol("rho")
    pdfs_sym = sp.symbols(f"f[:{lb_method.stencil.Q}]")
    vel_arr_symbols = [
        IndexedBase(TypedSymbol("u", dtype=default_dtype), shape=(1,))[i]
        for i in range(len(vel_symbols))]
    momentum_density_symbols = sp.symbols(f"md_:{len(vel_symbols)}")
    second_momentum_symbols = sp.symbols(f"p_:{len(vel_symbols)**2}")
    vel_in_symbols = sp.symbols(f"vel_in_[:{lb_method.stencil.D}]")

    equilibrium_subs_dict = {a: b for a, b in zip(vel_symbols, vel_arr_symbols)}
    equilibrium = lb_method.get_equilibrium()
    lhs_list = [a.lhs for a in equilibrium.main_assignments]
    equilibrium_matrix = sp.Matrix(
        [e.rhs for e in equilibrium.main_assignments])
    equilibrium = ps.AssignmentCollection([ps.Assignment(lhs, rhs)
                                           for lhs, rhs in zip(lhs_list, equilibrium_matrix)])
    equilibrium = __type_equilibrium_assignments(
        equilibrium, config, equilibrium_subs_dict)

    eq_input_from_input_eqs = cqc.equilibrium_input_equations_from_init_values(
        sp.Symbol("rho_in"), vel_arr_symbols)
    density_velocity_setter_macroscopic_values = equations_to_code(
        eq_input_from_input_eqs, dtype=default_dtype, variables_without_prefix=["rho_in", "u"], backend=backend)
    momentum_density_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {"density": rho_sym, "momentum_density": momentum_density_symbols})
    second_momentum_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {"moment2": second_momentum_symbols})
    constant_suffix = "f"

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
            variables_without_prefix=[e.name for e in pdfs_sym],
            dtype=default_dtype, backend=backend),
        "momentum_density_getter": equations_to_code(
            cqc.output_equations_from_pdfs(pdfs_sym,
                {"density": rho_sym, "momentum_density": momentum_density_symbols}),
            variables_without_prefix=[e.name for e in pdfs_sym],
            dtype=default_dtype, backend=backend),
        "pdfs_setter": equations_to_code(
            lbmpy.macroscopic_value_kernels.macroscopic_values_setter(
                lb_method, sp.Symbol("rho_in"), vel_in_symbols, pdfs_sym),
            variables_without_prefix=[e.name for e in pdfs_sym],
            dtype=default_dtype, backend=backend),

        "second_momentum_getter": equations_to_code(
            second_momentum_getter, variables_without_prefix=pdfs_sym,
            dtype=default_dtype, backend=backend),
        "density_velocity_setter_macroscopic_values": density_velocity_setter_macroscopic_values,

        "namespace": "lbm",
       }
    for k,v in jinja_context.items():
        print(f"{k}: \n{v}\n")

    env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)),
                      undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template(
        "macroscopic_values_accessors.tmpl.h").render(**jinja_context)

    gen_context.write_file(filename, header)


with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = False
    ctx.cuda = True
    params = {"target": target}

    # codegen configuration
    config = pystencils_espresso.generate_config(ctx, params)

    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    precision_suffix = pystencils_espresso.precision_suffix[ctx.double_accuracy]
    precision_rng = pystencils_espresso.precision_rng[ctx.double_accuracy]
    kT = sp.symbols("kT")
    stencil = lbmpy.stencils.LBStencil(lbmpy.enums.Stencil.D3Q19)
    fields = pystencils_espresso.generate_fields(config, stencil)
    force_field = fields["force"]

    # LB Method definition
    method = lbmpy.creationfunctions.create_mrt_orthogonal(
        stencil=stencil,
        compressible=True,
        weighted=True,
        relaxation_rates=relaxation_rates.rr_getter,
        force_model=lbmpy.forcemodels.Schiller(force_field.center_vector)
    )

    # generate accessors
    if target == ps.Target.GPU:
        #walberla_lbm_generation.generate_macroscopic_values_accessors(
        generate_macroscopic_values_accessors(
            ctx,
            config,
            method,
            f"macroscopic_values_accessors_{precision_suffix}CUDA.h")
