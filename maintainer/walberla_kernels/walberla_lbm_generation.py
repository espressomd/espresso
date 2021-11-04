#
# Copyright (C) 2021 The ESPResSo project
# Copyright (C) 2020-2021 The waLBerla project
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
import numpy as np
import sympy as sp
import pystencils as ps

# File derived from lbmpy_walberla.walberla_lbm_generation in the
# walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d


def __macroscopic_values_accessors(generation_context, lb_method):

    # Function derived from lbmpy_walberla.walberla_lbm_generation.__lattice_model()
    # in the walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d

    from jinja2 import Environment, FileSystemLoader, StrictUndefined
    from sympy.tensor import IndexedBase
    from pystencils.backends.cbackend import CustomSympyPrinter
    from pystencils.transformations import add_types
    from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
    from lbmpy_walberla.walberla_lbm_generation import get_stencil_name, equations_to_code, stencil_switch_statement

    cpp_printer = CustomSympyPrinter()
    stencil_name = get_stencil_name(lb_method.stencil)
    if not stencil_name:
        raise ValueError(
            "lb_method uses a stencil that is not supported in waLBerla")

    is_float = not generation_context.double_accuracy
    dtype_string = "float32" if is_float else "float64"

    vel_symbols = lb_method.conserved_quantity_computation.first_order_moment_symbols
    rho_sym = sp.Symbol('rho')
    pdfs_sym = sp.symbols(f'f_:{len(lb_method.stencil)}')
    vel_arr_symbols = [
        IndexedBase(sp.Symbol('u'), shape=(1,))[i]
        for i in range(len(vel_symbols))]
    momentum_density_symbols = sp.symbols(f'md_:{len(vel_symbols)}')
    second_momentum_symbols = sp.symbols(f'p_:{len(vel_symbols)**2}')

    equilibrium = lb_method.get_equilibrium()
    equilibrium = equilibrium.new_with_substitutions(
        {a: b for a, b in zip(vel_symbols, vel_arr_symbols)})
    equilibrium = add_types(
        equilibrium.main_assignments, dtype_string, False)[2]
    equilibrium = sp.Matrix([e.rhs for e in equilibrium])

    cqc = lb_method.conserved_quantity_computation

    eq_input_from_input_eqs = cqc.equilibrium_input_equations_from_init_values(
        sp.Symbol('rho_in'), vel_arr_symbols)
    density_velocity_setter_macroscopic_values = equations_to_code(
        eq_input_from_input_eqs, dtype=dtype_string, variables_without_prefix=['rho_in', 'u'])
    momentum_density_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {'density': rho_sym, 'momentum_density': momentum_density_symbols})
    second_momentum_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {'moment2': second_momentum_symbols})

    jinja_context = {
        'stencil_name': stencil_name,
        'D': lb_method.dim,
        'Q': len(lb_method.stencil),
        'compressible': lb_method.conserved_quantity_computation.compressible,

        'equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, equilibrium),
        'equilibrium': [cpp_printer.doprint(e) for e in equilibrium],

        'density_getters': equations_to_code(cqc.output_equations_from_pdfs(pdfs_sym, {"density": rho_sym}),
                                             variables_without_prefix=[e.name for e in pdfs_sym], dtype=dtype_string),
        'momentum_density_getter': equations_to_code(momentum_density_getter, variables_without_prefix=pdfs_sym,
                                                     dtype=dtype_string),
        'second_momentum_getter': equations_to_code(second_momentum_getter, variables_without_prefix=pdfs_sym,
                                                    dtype=dtype_string),
        'density_velocity_setter_macroscopic_values': density_velocity_setter_macroscopic_values,

        'namespace': 'lbm',
    }

    env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)),
                      undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template(
        'macroscopic_values_accessors.tmpl.h').render(**jinja_context)

    generation_context.write_file('macroscopic_values_accessors.h', header)
    with open('macroscopic_values_accessors.h', 'r+') as f:
        content = f.read()
        f.seek(0)
        f.truncate()
        f.write(content.replace('lm.force_->', 'force_field.'))


def generate_macroscopic_values_accessors(
        generation_context, collision_rule, field_layout='zyxf',
        **create_kernel_params):

    # Function derived from lbmpy_walberla.walberla_lbm_generation.generate_lattice_model()
    # in the walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d

    from lbmpy.fieldaccess import StreamPullTwoFieldsAccessor
    from lbmpy.updatekernels import create_lbm_kernel
    from pystencils_walberla.codegen import default_create_kernel_parameters

    # usually a numpy layout is chosen by default i.e. xyzf - which is bad for waLBerla where at least the spatial
    # coordinates should be ordered in reverse direction i.e. zyx
    is_float = not generation_context.double_accuracy
    dtype = np.float32 if is_float else np.float64
    lb_method = collision_rule.method

    q = len(lb_method.stencil)
    dim = lb_method.dim

    create_kernel_params = default_create_kernel_parameters(
        generation_context, create_kernel_params)

    if field_layout == 'fzyx':
        create_kernel_params['cpu_vectorize_info']['assume_inner_stride_one'] = True
    elif field_layout == 'zyxf':
        create_kernel_params['cpu_vectorize_info']['assume_inner_stride_one'] = False

    src_field = ps.Field.create_generic(
        'pdfs',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,))
    dst_field = ps.Field.create_generic(
        'pdfs_tmp', dim, dtype, index_dimensions=1, layout=field_layout, index_shape=(q,))

    stream_collide_update_rule = create_lbm_kernel(
        collision_rule, src_field, dst_field, StreamPullTwoFieldsAccessor())
    stream_collide_ast = ps.create_kernel(
        stream_collide_update_rule,
        **create_kernel_params)
    stream_collide_ast.function_name = 'kernel_streamCollide'
    stream_collide_ast.assumed_inner_stride_one = create_kernel_params[
        'cpu_vectorize_info']['assume_inner_stride_one']

    __macroscopic_values_accessors(generation_context, lb_method)
