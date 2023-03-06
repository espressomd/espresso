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
from pystencils.node_collection import NodeCollection
from pystencils.typing.transformations import add_types
from pystencils.typing.typed_sympy import TypedSymbol

# File derived from lbmpy_walberla.walberla_lbm_generation in the
# walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d


def __type_equilibrium_assignments(assignments, config, subs_dict):
    # Function derived from lbmpy_walberla.walberla_lbm_generation.__type_equilibrium_assignments()
    # in the walberla project, commit 9dcd0dd90f50f7b64b0a38bb06327854463fdafd
    result = assignments.new_with_substitutions(subs_dict)
    result = NodeCollection(result.main_assignments)
    result.evaluate_terms()
    result = add_types(result.all_assignments, config)
    return result


def generate_macroscopic_values_accessors(
        gen_context, config, lb_method, filename):

    # Function derived from lbmpy_walberla.walberla_lbm_generation.__lattice_model()
    # in the walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d
    # with backports from commit de6b00071233a9a1f45d7a6773988363e058f1a0

    from jinja2 import Environment, FileSystemLoader, StrictUndefined
    from sympy.tensor import IndexedBase
    from pystencils.backends.cbackend import CustomSympyPrinter
    from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
    from lbmpy_walberla.walberla_lbm_generation import equations_to_code, stencil_switch_statement

    cpp_printer = CustomSympyPrinter()
    stencil_name = lb_method.stencil.name
    if not stencil_name:
        raise ValueError(
            "lb_method uses a stencil that is not supported in waLBerla")

    default_dtype = config.data_type.default_factory()

    cqc = lb_method.conserved_quantity_computation
    vel_symbols = cqc.velocity_symbols
    rho_sym = sp.Symbol("rho")
    pdfs_sym = sp.symbols(f"f_:{lb_method.stencil.Q}")
    vel_arr_symbols = [
        IndexedBase(TypedSymbol("u", dtype=default_dtype), shape=(1,))[i]
        for i in range(len(vel_symbols))]
    momentum_density_symbols = sp.symbols(f"md_:{len(vel_symbols)}")
    second_momentum_symbols = sp.symbols(f"p_:{len(vel_symbols)**2}")

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
        eq_input_from_input_eqs, dtype=default_dtype, variables_without_prefix=["rho_in", "u"])
    momentum_density_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {"density": rho_sym, "momentum_density": momentum_density_symbols})
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
            variables_without_prefix=[e.name for e in pdfs_sym],
            dtype=default_dtype),
        "momentum_density_getter": equations_to_code(
            momentum_density_getter, variables_without_prefix=pdfs_sym,
            dtype=default_dtype),
        "second_momentum_getter": equations_to_code(
            second_momentum_getter, variables_without_prefix=pdfs_sym,
            dtype=default_dtype),
        "density_velocity_setter_macroscopic_values": density_velocity_setter_macroscopic_values,

        "namespace": "lbm",
    }

    env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)),
                      undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template(
        "macroscopic_values_accessors.tmpl.h").render(**jinja_context)

    gen_context.write_file(filename, header)
    with open(filename, "r+") as f:
        content = f.read()
        f.seek(0)
        f.truncate(0)
        # remove lattice model
        content = content.replace("lm.force_->", "force_field.")
        f.write(content)
