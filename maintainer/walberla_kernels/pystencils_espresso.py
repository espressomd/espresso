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

import numpy as np
import sympy as sp
import lbmpy.fieldaccess
import lbmpy.macroscopic_value_kernels
import lbmpy.updatekernels
import pystencils as ps
import pystencils_walberla
import pystencils_walberla.utility

from pystencils.astnodes import LoopOverCoordinate
from pystencils.backends.cbackend import CustomCodeNode


def skip_philox_unthermalized(code, result_symbols, rng_name):
    for r in result_symbols:
        statement = f" {r.name};"
        assert statement in code, f"no declaration for variable '{r.name}' in '{code}'"  # nopep8
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


class PhiloxTwoDoublesModulo(ps.rng.PhiloxTwoDoubles):
    def __init__(self, dim, time_step=ps.typing.TypedSymbol(
            "time_step", np.uint32), *args, **kwargs):
        super().__init__(dim, time_step=time_step, *args, **kwargs)

        if kwargs.get("keys") is None:
            keys = (0,) * self._num_keys
        else:
            keys = kwargs.get("keys")

        if kwargs.get("offsets") is None:
            offsets = (0,) * dim
        else:
            offsets = kwargs.get("offsets")

        coordinates = [
            LoopOverCoordinate.get_loop_counter_symbol(i) +
            offsets[i] for i in range(dim)]
        if dim < 3:
            coordinates.append(0)

        # add folding to coordinates with symbols
        field_sizes = [
            ps.typing.TypedSymbol(
                f"field_size_{i}",
                np.uint32) for i in range(dim)]

        new_coordinates = [
            (coord) %
            field_size for coord,
            field_size in zip(
                coordinates,
                field_sizes)]
        self._args = sp.sympify([time_step, *new_coordinates, *keys])
        symbols_read = set.union(*[s.atoms(sp.Symbol) for s in self.args])

        headers = self.headers
        # set headers again, since the constructor of CustomCodeNode resets the
        # header-list
        CustomCodeNode.__init__(
            self,
            "",
            symbols_read=symbols_read,
            symbols_defined=self.result_symbols)
        self.headers = headers


class PhiloxFourFloats(ps.rng.PhiloxFourFloats):
    def get_code(self, *args, **kwargs):
        code = super().get_code(*args, **kwargs)
        return skip_philox_unthermalized(code, self.result_symbols, self._name)


class PhiloxFourFloatsModulo(ps.rng.PhiloxFourFloats):
    def __init__(self, dim, time_step=ps.typing.TypedSymbol(
            "time_step", np.uint32), *args, **kwargs):
        super().__init__(dim, time_step=time_step, *args, **kwargs)

        if kwargs.get("keys") is None:
            keys = (0,) * self._num_keys
        else:
            keys = kwargs.get("keys")

        if kwargs.get("offsets") is None:
            offsets = (0,) * dim
        else:
            offsets = kwargs.get("offsets")

        coordinates = [
            LoopOverCoordinate.get_loop_counter_symbol(i) +
            offsets[i] for i in range(dim)]
        if dim < 3:
            coordinates.append(0)

        # add folding to coordinates with symbols
        field_sizes = [
            ps.typing.TypedSymbol(
                f"field_size_{i}",
                np.uint32) for i in range(dim)]

        new_coordinates = [
            (coord) %
            field_size for coord,
            field_size in zip(
                coordinates,
                field_sizes)]
        self._args = sp.sympify([time_step, *new_coordinates, *keys])
        symbols_read = set.union(*[s.atoms(sp.Symbol) for s in self.args])

        headers = self.headers
        # set headers again, since the constructor of CustomCodeNode resets the
        # header-list
        CustomCodeNode.__init__(
            self,
            "",
            symbols_read=symbols_read,
            symbols_defined=self.result_symbols)
        self.headers = headers


precision_prefix = {
    True: 'DoublePrecision',
    False: 'SinglePrecision'}
precision_suffix = {
    True: 'double_precision',
    False: 'single_precision'}
precision_rng = {
    True: PhiloxTwoDoubles,
    False: PhiloxFourFloats}
precision_rng_modulo = {
    True: PhiloxTwoDoublesModulo,
    False: PhiloxFourFloatsModulo}
data_type_np = {'double': 'float64', 'float': 'float32'}


def generate_fields(config, stencil, field_layout='fzyx'):
    dtype = data_type_np[config.data_type.default_factory().c_name]
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


def generate_pack_info_pdfs_field_assignments(fields, streaming_pattern):
    """
    Visualize the stencil directions with::

       import lbmpy
       import matplotlib.pyplot as plt
       stencil = lbmpy.LBStencil(lbmpy.Stencil.D3Q19)
       stencil.plot(data=[i for i in range(19)])
       plt.show()

    """
    stencil = lbmpy.enums.Stencil.D3Q19
    lbm_config = lbmpy.LBMConfig(stencil=stencil,
                                 method=lbmpy.Method.CUMULANT,
                                 compressible=True,
                                 zero_centered=False,
                                 weighted=True,
                                 streaming_pattern=streaming_pattern,
                                 relaxation_rate=sp.Symbol("omega_shear"),
                                 )
    lbm_opt = lbmpy.LBMOptimisation(
        symbolic_field=fields["pdfs" if streaming_pattern ==
                              "pull" else "pdfs_tmp"],
        symbolic_temporary_field=fields["pdfs" if streaming_pattern ==
                                        "push" else "pdfs_tmp"],
        field_layout=fields['pdfs'].layout)
    lbm_update_rule = lbmpy.create_lb_update_rule(
        lbm_config=lbm_config,
        lbm_optimisation=lbm_opt)
    return lbm_update_rule.all_assignments


def generate_pack_info_vector_field_specifications(config, stencil, layout):
    import collections
    import itertools
    field = ps.Field.create_generic(
        "field",
        3,
        data_type_np[config.data_type.default_factory().c_name],
        index_dimensions=1,
        layout=layout,
        index_shape=(3,)
    )
    q = len(stencil)
    coord = itertools.product(*[(-1, 0, 1)] * 3)
    if q == 19:
        dirs = tuple((i, j, k) for i, j, k in coord if i**2 + j**2 + k**2 != 3)
    else:
        dirs = tuple((i, j, k) for i, j, k in coord)
    spec = collections.defaultdict(set)
    spec[dirs] = {field[0, 0, 0](i) for i in range(3)}
    return spec


def generate_config(ctx, params):
    return pystencils_walberla.utility.config_from_context(ctx, **params)


def generate_collision_sweep(
        ctx, lb_method, collision_rule, class_name, params, **kwargs):
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
    pystencils_walberla.generate_sweep(
        ctx, class_name, collide_ast, **params, **kwargs)


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
    pystencils_walberla.generate_sweep(
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


# this can be removed when https://i10git.cs.fau.de/walberla/walberla/-/issues/247 is fixed
def generate_staggered_flux_boundary(generation_context, class_name, boundary_object,
                                     dim, neighbor_stencil, index_shape, target=ps.Target.CPU, **kwargs):
    assert dim == len(neighbor_stencil[0])
    pystencils_walberla.boundary.generate_boundary(
        generation_context=generation_context,
        class_name=class_name,
        boundary_object=boundary_object,
        field_name='flux',
        neighbor_stencil=neighbor_stencil,
        index_shape=index_shape,
        field_type=ps.FieldType.STAGGERED_FLUX,
        kernel_creation_function=None,
        target=target,
        **kwargs)
