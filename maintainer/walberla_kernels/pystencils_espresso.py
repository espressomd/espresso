#
# Copyright (C) 2021-2023 The ESPResSo project
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

import sympy as sp
import lbmpy.fieldaccess
import lbmpy.macroscopic_value_kernels
import lbmpy.updatekernels
import pystencils as ps
import pystencils_walberla


def skip_philox_unthermalized(code, result_symbols, rng_name):
    for r in result_symbols:
        statement = f" {r.name};"
        assert statement in code, f"no declaration for variable '{r.name}' in '{code}'"
        code = code.replace(statement, f" {r.name}{{}};", 1)
    statement = f"{rng_name}("
    assert code.count(statement) == 1, f"need 1 '{rng_name}' call in '{code}'"
    lines = code.rstrip().split("\n")
    assert lines[-1].startswith(rng_name), f"'{rng_name}' not in '{lines[-1]}'"
    lines[-1] = f"if (kT > 0.) {{  \n{lines[-1]}\n}}"
    return "\n".join(lines)


class PhiloxTwoDoubles(ps.rng.PhiloxTwoDoubles):
    def get_code(self, *args, **kwargs):
        code = super().get_code(*args, **kwargs)
        return skip_philox_unthermalized(code, self.result_symbols, self._name)


class PhiloxFourFloats(ps.rng.PhiloxFourFloats):
    def get_code(self, *args, **kwargs):
        code = super().get_code(*args, **kwargs)
        return skip_philox_unthermalized(code, self.result_symbols, self._name)


precision_prefix = {
    True: 'DoublePrecision',
    False: 'SinglePrecision'}
precision_suffix = {
    True: 'double_precision',
    False: 'single_precision'}
precision_rng = {
    True: PhiloxTwoDoubles,
    False: PhiloxFourFloats}
data_type_np = {'double': 'float64', 'float': 'float32'}


def generate_fields(config, stencil):
    dtype = data_type_np[config.data_type.default_factory().c_name]
    field_layout = 'fzyx'
    q = len(stencil)
    dim = len(stencil[0])

    fields = {}
    # Symbols for PDF (twice, due to double buffering)
    fields['pdfs'] = ps.Field.create_generic(
        'pdfs',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['pdfs_tmp'] = ps.Field.create_generic(
        'pdfs_tmp',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['velocity'] = ps.Field.create_generic(
        'velocity',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )
    fields['force'] = ps.Field.create_generic(
        'force',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )

    return fields


def generate_config(ctx, params):
    return pystencils_walberla.codegen.config_from_context(ctx, **params)


def generate_collision_sweep(
        ctx, lb_method, collision_rule, class_name, params):
    config = generate_config(ctx, params)

    # Symbols for PDF (twice, due to double buffering)
    fields = generate_fields(config, lb_method.stencil)

    # Generate collision kernel
    collide_update_rule = lbmpy.updatekernels.create_lbm_kernel(
        collision_rule,
        fields['pdfs'],
        fields['pdfs_tmp'],
        lbmpy.fieldaccess.CollideOnlyInplaceAccessor())
    collide_ast = ps.create_kernel(
        collide_update_rule, config=config, **params)
    collide_ast.function_name = 'kernel_collide'
    collide_ast.assumed_inner_stride_one = True
    pystencils_walberla.codegen.generate_sweep(
        ctx, class_name, collide_ast, **params)


def generate_stream_sweep(ctx, lb_method, class_name, params):
    config = generate_config(ctx, params)

    # Symbols for PDF (twice, due to double buffering)
    fields = generate_fields(config, lb_method.stencil)

    # Generate stream kernel
    stream_update_rule = lbmpy.updatekernels.create_stream_pull_with_output_kernel(
        lb_method, fields['pdfs'], fields['pdfs_tmp'],
        output={'velocity': fields['velocity']})
    stream_ast = ps.create_kernel(stream_update_rule, config=config, **params)
    stream_ast.function_name = 'kernel_stream'
    stream_ast.assumed_inner_stride_one = True
    pystencils_walberla.codegen.generate_sweep(
        ctx, class_name, stream_ast,
        field_swaps=[(fields['pdfs'], fields['pdfs_tmp'])], **params)


def generate_setters(ctx, lb_method, params):
    config = generate_config(ctx, params)
    fields = generate_fields(config, lb_method.stencil)

    initial_rho = sp.Symbol('rho_0')
    pdfs_setter = lbmpy.macroscopic_value_kernels.macroscopic_values_setter(
        lb_method,
        initial_rho,
        fields['velocity'].center_vector,
        fields['pdfs'].center_vector)
    return pdfs_setter
