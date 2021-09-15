#
# Copyright (C) 2020-2021 The ESPResSo project
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
import jinja2
import sympy as sp
import pystencils as ps
from pystencils.field import Field
from lbmpy.creationfunctions import create_lb_collision_rule, create_mrt_orthogonal, force_model_from_string
from lbmpy.stencils import get_stencil
from pystencils_walberla import CodeGeneration, codegen

from lbmpy_walberla import generate_boundary, generate_lb_pack_info

from lbmpy.boundaries import UBB
from lbmpy_walberla.additional_data_handler import UBBAdditionalDataHandler
import relaxation_rates
from lbmpy.fieldaccess import CollideOnlyInplaceAccessor
from lbmpy.stencils import get_stencil
from lbmpy.updatekernels import create_lbm_kernel, create_stream_pull_with_output_kernel
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter

import lees_edwards
# for collide-push from lbmpy.fieldaccess import StreamPushTwoFieldsAccessor


def adapt_pystencils():
    '''
    Adapt pystencils to the SFINAE method (add the block offset lambda
    callback and the time_step increment).
    '''
    old_add_pystencils_filters_to_jinja_env = codegen.add_pystencils_filters_to_jinja_env

    def new_add_pystencils_filters_to_jinja_env(jinja_env):
        # save original pystencils to adapt
        old_add_pystencils_filters_to_jinja_env(jinja_env)
        old_generate_members = jinja_env.filters['generate_members']
        old_generate_refs_for_kernel_parameters = jinja_env.filters[
            'generate_refs_for_kernel_parameters']

        @jinja2.contextfilter
        def new_generate_members(*args, **kwargs):
            output = old_generate_members(*args, **kwargs)
            token = ' block_offset_0_;'
            if token in output:
                i = output.index(token)
                vartype = output[:i].split('\n')[-1].strip()
                output += f'\nstd::function<void(IBlock *, {vartype}&, {vartype}&, {vartype}&)> block_offset_generator = [](IBlock * const, {vartype}&, {vartype}&, {vartype}&) {{ }};'
            return output

        def new_generate_refs_for_kernel_parameters(*args, **kwargs):
            output = old_generate_refs_for_kernel_parameters(*args, **kwargs)
            if 'block_offset_0' in output:
                old_token = 'auto & block_offset_'
                new_token = 'auto block_offset_'
                assert output.count(old_token) == 3, \
                    f'could not find "{old_token}" in """\n{output}\n"""'
                output = output.replace(old_token, new_token)
                output += '\nblock_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);'
            return output

        # replace pystencils
        jinja_env.filters['generate_members'] = new_generate_members
        jinja_env.filters['generate_refs_for_kernel_parameters'] = new_generate_refs_for_kernel_parameters

    codegen.add_pystencils_filters_to_jinja_env = new_add_pystencils_filters_to_jinja_env


def earmark_generated_kernels():
    '''
    Add an earmark at the beginning of generated kernels to document the
    pystencils/lbmpy toolchain that was used to create them.
    '''
    import lbmpy
    import lbmpy_walberla
    walberla_root = lbmpy_walberla.__file__.split('/python/lbmpy_walberla/')[0]
    with open(os.path.join(walberla_root, '.git/HEAD')) as f:
        walberla_commit = f.read()
    token = '// kernel generated with'
    earmark = (
        f'{token} pystencils v{ps.__version__}, lbmpy v{lbmpy.__version__}, '
        f'lbmpy_walberla/pystencils_walberla from commit {walberla_commit}\n'
    )
    for filename in os.listdir('.'):
        if filename.endswith(('.h', '.cpp')):
            with open(filename, 'r+') as f:
                content = f.read()
                if not content.startswith(token):
                    f.seek(0)
                    f.write(earmark + content)


def generate_fields(lb_method):
    dtype = "float64" 
    field_layout = "fzyx"
    q = len(lb_method.stencil)
    dim = len(lb_method.stencil[0])

    # Symbols for PDF (twice, due to double buffering)
    src_field = ps.Field.create_generic(
        'pdfs',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    dst_field = ps.Field.create_generic(
        'pdfs_tmp',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    vel_field = ps.Field.create_generic(
        'velocity',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )

    return src_field, dst_field, vel_field


def generate_collision_sweep(
        ctx, lb_method, collision_rule, class_name, params):

    # Symbols for PDF (twice, due to double buffering)
    src_field, dst_field, _ = generate_fields(lb_method)

    # Generate collision kernel
    collide_update_rule = create_lbm_kernel(
        collision_rule,
        src_field,
        dst_field,
        CollideOnlyInplaceAccessor())
    collide_ast = ps.create_kernel(collide_update_rule, **params)
    collide_ast.function_name = 'kernel_collide'
    collide_ast.assumed_inner_stride_one = True
    codegen.generate_sweep(ctx, class_name, collide_ast)


def generate_stream_sweep(ctx, lb_method, class_name, params):
    #dtype = "float64"
    #field_layout = "fzyx"

    # Symbols for PDF (twice, due to double buffering)
    src_field, dst_field, velocity_field = generate_fields(lb_method)

    # Generate stream kernel
    stream_update_rule = create_stream_pull_with_output_kernel(lb_method, src_field, dst_field,
                                                               output={'velocity': velocity_field})
    # stream_update_rule = create_stream_pull_only_kernel(
    #    lb_method.stencil, None, 'pdfs', 'pdfs_tmp', field_layout, dtype)
    stream_ast = ps.create_kernel(stream_update_rule, **params)
    stream_ast.function_name = 'kernel_stream'
    stream_ast.assumed_inner_stride_one = True
    codegen.generate_sweep(
        ctx, class_name, stream_ast, field_swaps=[(src_field, dst_field)])


def generate_setters(lb_method):
    pdf_field, _, vel_field = generate_fields(lb_method)

    initial_rho = sp.Symbol('rho_0')
    pdfs_setter = macroscopic_values_setter(lb_method,
                                            initial_rho,
                                            vel_field.center_vector,
                                            pdf_field.center_vector)
    return pdfs_setter


def patch_accessors(classname, name):
    with open(f'{classname}.h') as f:
        with open(f'{name}.h', 'w') as g:
            g.write(f.read().replace('lm.force_->', 'force_field.'))
    os.remove(f'{classname}.h')


def __lattice_model(generation_context, class_name, lb_method,
                    stream_collide_ast, refinement_scaling):

    # Function derived from lbmpy_walberla.walberla_lbm_generation.__lattice_model()
    # in the walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d

    import sympy as sp
    from jinja2 import Environment, FileSystemLoader, StrictUndefined
    from sympy.tensor import IndexedBase
    from pystencils.backends.cbackend import CustomSympyPrinter, get_headers
    from pystencils.data_types import cast_func
    from pystencils.sympyextensions import get_symmetric_part
    from pystencils.transformations import add_types
    from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
    from lbmpy_walberla.walberla_lbm_generation import get_stencil_name, equations_to_code, expression_to_code, stencil_switch_statement

    cpp_printer = CustomSympyPrinter()
    stencil_name = get_stencil_name(lb_method.stencil)
    if not stencil_name:
        raise ValueError(
            "lb_method uses a stencil that is not supported in waLBerla")

    communication_stencil_name = stencil_name if stencil_name != "D3Q15" else "D3Q27"
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

    symmetric_equilibrium = get_symmetric_part(equilibrium, vel_arr_symbols)
    symmetric_equilibrium = symmetric_equilibrium.subs(
        sp.Rational(1, 2), cast_func(sp.Rational(1, 2), dtype_string))
    asymmetric_equilibrium = sp.expand(equilibrium - symmetric_equilibrium)

    force_model = lb_method.force_model
    macroscopic_velocity_shift = None
    if force_model and hasattr(force_model, 'macroscopic_velocity_shift'):
        es = (sp.Rational(1, 2), cast_func(sp.Rational(1, 2), dtype_string))
        macroscopic_velocity_shift = [
            expression_to_code(e.subs(*es), "lm.", ['rho'], dtype=dtype_string)
            for e in force_model.macroscopic_velocity_shift(rho_sym)]

    cqc = lb_method.conserved_quantity_computation

    eq_input_from_input_eqs = cqc.equilibrium_input_equations_from_init_values(
        sp.Symbol('rho_in'), vel_arr_symbols)
    density_velocity_setter_macroscopic_values = equations_to_code(
        eq_input_from_input_eqs, dtype=dtype_string, variables_without_prefix=['rho_in', 'u'])
    momentum_density_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {'density': rho_sym, 'momentum_density': momentum_density_symbols})
    second_momentum_getter = cqc.output_equations_from_pdfs(
        pdfs_sym, {'moment2': second_momentum_symbols})
    constant_suffix = "f" if is_float else ""

    required_headers = get_headers(stream_collide_ast)

    if refinement_scaling:
        refinement_scaling_info = [(e0, e1, expression_to_code(e2, '', dtype=dtype_string)) for e0, e1, e2 in
                                   refinement_scaling.scaling_info]
        # append '_' to entries since they are used as members
        for i in range(len(refinement_scaling_info)):
            updated_entry = (refinement_scaling_info[i][0],
                             refinement_scaling_info[i][1].replace(refinement_scaling_info[i][1],
                                                                   refinement_scaling_info[i][1] + '_'),
                             refinement_scaling_info[i][2].replace(refinement_scaling_info[i][1],
                                                                   refinement_scaling_info[i][1] + '_'),
                             )
            refinement_scaling_info[i] = updated_entry
    else:
        refinement_scaling_info = None

    jinja_context = {
        'class_name': class_name,
        'stencil_name': stencil_name,
        'communication_stencil_name': communication_stencil_name,
        'D': lb_method.dim,
        'Q': len(lb_method.stencil),
        'compressible': lb_method.conserved_quantity_computation.compressible,
        'weights': ",".join(str(w.evalf()) + constant_suffix for w in lb_method.weights),
        'inverse_weights': ",".join(str((1 / w).evalf()) + constant_suffix for w in lb_method.weights),

        'equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, equilibrium),
        'symmetric_equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, symmetric_equilibrium),
        'asymmetric_equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, asymmetric_equilibrium),
        'equilibrium': [cpp_printer.doprint(e) for e in equilibrium],

        'macroscopic_velocity_shift': macroscopic_velocity_shift,
        'density_getters': equations_to_code(cqc.output_equations_from_pdfs(pdfs_sym, {"density": rho_sym}),
                                             variables_without_prefix=[e.name for e in pdfs_sym], dtype=dtype_string),
        'momentum_density_getter': equations_to_code(momentum_density_getter, variables_without_prefix=pdfs_sym,
                                                     dtype=dtype_string),
        'second_momentum_getter': equations_to_code(second_momentum_getter, variables_without_prefix=pdfs_sym,
                                                    dtype=dtype_string),
        'density_velocity_setter_macroscopic_values': density_velocity_setter_macroscopic_values,

        'refinement_scaling_info': refinement_scaling_info,

        'target': 'cpu',
        'namespace': 'lbm',
        'headers': required_headers,
    }

    env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)),
                      undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template(
        'MacroscopicValuesAccessors.tmpl.h').render(**jinja_context)

    generation_context.write_file(f"{class_name}.h", header)


def generate_macroscopic_values_accessors(
        generation_context, class_name, collision_rule, field_layout='zyxf',
        refinement_scaling=None, **create_kernel_params):

    # Function derived from lbmpy_walberla.walberla_lbm_generation.generate_lattice_model()
    # in the walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d

    import numpy as np
    import pystencils as ps
    from lbmpy.fieldaccess import StreamPullTwoFieldsAccessor
    from lbmpy.updatekernels import create_lbm_kernel
    from pystencils import create_kernel
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
    if create_kernel_params['target'] == 'gpu':
        raise ValueError(
            "Lattice Models can only be generated for CPUs. To generate LBM on GPUs use sweeps directly")

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
    stream_collide_ast = create_kernel(
        stream_collide_update_rule,
        **create_kernel_params)
    stream_collide_ast.function_name = 'kernel_streamCollide'
    stream_collide_ast.assumed_inner_stride_one = create_kernel_params[
        'cpu_vectorize_info']['assume_inner_stride_one']

    __lattice_model(
        generation_context,
        class_name,
        lb_method,
        stream_collide_ast,
        refinement_scaling)


adapt_pystencils()

with CodeGeneration() as ctx:
    kT = sp.symbols("kT")
    force_field = ps.fields("force(3): [3D]", layout='fzyx')

    # Vectorization parameters
    cpu_vectorize_info = {
        "instruction_set": "avx",
        "assume_inner_stride_one": True,
        "assume_aligned": True,
        "assume_sufficient_line_padding": True}
    params = {"target": "cpu"}
    params_vec = {"target": "cpu", "cpu_vectorize_info": cpu_vectorize_info}

    # LB Method definition
    method = create_mrt_orthogonal(
        stencil=get_stencil('D3Q19'),
        compressible=True,
        weighted=True,
        relaxation_rate_getter=relaxation_rates.rr_getter,
        force_model=force_model_from_string(
            'schiller', force_field.center_vector)
    )
    generate_stream_sweep(ctx, method, "StreamSweep", params)
    generate_stream_sweep(ctx, method, "StreamSweepAVX", params_vec)

    # generate initial densities
    pdfs_setter = generate_setters(method)
    codegen.generate_sweep(ctx, "InitialPDFsSetter", pdfs_setter)

    # generate unthermalized collision rule
    u_p = lees_edwards.velocity_shift(1, 64)
    velocity_field = ps.fields("velocitye(3): [3D]", layout='fzyx')
    density_field = ps.fields("density(1): [3D]", layout='fzyx')
    collision_rule_unthermalized = create_lb_collision_rule(
        method,
        velocity_input=velocity_field.center_vector + sp.Matrix([u_p, 0, 0]),
        density_input=density_field,
        optimization={'cse_global': True}
    )
    lees_edwards.modify_method(collision_rule_unthermalized, 0, 1)
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_unthermalized,
        "CollideSweep",
        {}
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_unthermalized,
        "CollideSweepAVX",
        {"cpu_vectorize_info": cpu_vectorize_info}
    )

    # generate thermalized LB
    collision_rule_thermalized = create_lb_collision_rule(
        method,
        fluctuating={
            'temperature': kT,
            'block_offsets': 'walberla',
            'rng_node': ps.rng.PhiloxTwoDoubles
        },
        optimization={'cse_global': True}
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_thermalized,
        "CollideSweepThermalized",
        params
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_thermalized,
        "CollideSweepThermalizedAVX",
        params_vec
    )

    # generate accessors
    generate_macroscopic_values_accessors(
        ctx,
        'LBWalberlaImpl',
        collision_rule_unthermalized,
        field_layout="fzyx")
    patch_accessors('LBWalberlaImpl', 'macroscopic_values_accessors')

    # Boundary conditions
    ubb_dynamic = UBB(lambda *args: None, dim=3)
    ubb_data_handler = UBBAdditionalDataHandler(method.stencil, ubb_dynamic)

    generate_boundary(ctx, 'Dynamic_UBB', ubb_dynamic, method,
                      additional_data_handler=ubb_data_handler,
                      streaming_pattern="push")

    # communication
    pdfs = Field.create_generic(
        'pdfs', 3, index_shape=(len(method.stencil),), layout='fzyx')
    generate_lb_pack_info(
        ctx,
        'PushPackInfo',
        method.stencil,
        pdfs,
        streaming_pattern='push'
    )

    # Info header containing correct template definitions for stencil and field
    #ctx.write_file("InfoHeader.h", info_header)

earmark_generated_kernels()
