#
# Copyright (C) 2020 The ESPResSo project
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
from lbmpy.updatekernels import create_lbm_kernel, create_stream_pull_only_kernel
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
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
        old_generate_call = jinja_env.filters['generate_call']
        old_generate_members = jinja_env.filters['generate_members']
        old_generate_refs_for_kernel_parameters = jinja_env.filters[
            'generate_refs_for_kernel_parameters']

        @jinja2.contextfilter
        def new_generate_call(*args, **kwargs):
            output = old_generate_call(*args, **kwargs)
            function_call = output.split('\n')[-1]
            if 'block_offset_0' in function_call:
                output += '\nthis->time_step_++;'
            return output

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
        jinja_env.filters['generate_call'] = new_generate_call
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

    return src_field, dst_field


def generate_collision_sweep(
        ctx, lb_method, collision_rule, class_name, params):

    # Symbols for PDF (twice, due to double buffering)
    src_field, dst_field = generate_fields(lb_method)

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
    dtype = "float64"
    field_layout = "fzyx"

    # Symbols for PDF (twice, due to double buffering)
    src_field, dst_field = generate_fields(lb_method)

    # Generate stream kernel
    stream_update_rule = create_stream_pull_only_kernel(
        lb_method.stencil, None, 'pdfs', 'pdfs_tmp', field_layout, dtype)
    stream_ast = ps.create_kernel(stream_update_rule, **params)
    stream_ast.function_name = 'kernel_stream'
    stream_ast.assumed_inner_stride_one = True
    codegen.generate_sweep(
        ctx, class_name, stream_ast, field_swaps=[(src_field, dst_field)])


def generate_setters(lb_method):
    dtype = "float64"
    field_layout = "fzyx"
    dim = len(lb_method.stencil[0])

    vel_field = ps.Field.create_generic(
        'velocity',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )

    pdfs = generate_fields(lb_method)[0]

    initial_rho = sp.Symbol('rho_0')
    pdfs_setter = macroscopic_values_setter(lb_method,
                                            initial_rho,
                                            vel_field.center_vector,
                                            pdfs.center_vector)
    return pdfs_setter


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
    codegen.generate_sweep(ctx, "InitialPDFsSetter", pdfs_setter, params)
    codegen.generate_sweep(
        ctx,
        "InitialPDFsSetterAVX",
        pdfs_setter,
        params_vec)

    # generate unthermalized collision rule
    collision_rule_unthermalized = create_lb_collision_rule(
        method,
        optimization={'cse_global': True}
    )
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
