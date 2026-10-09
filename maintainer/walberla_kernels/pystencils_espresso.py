#
# Copyright (C) 2021-2026 The ESPResSo project
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

numpy_types_to_cpp_types = {"float64": "double", "float32": "float"}
precision_prefix = {
    True: 'DoublePrecision',
    False: 'SinglePrecision'}
precision_suffix = {
    True: 'double_precision',
    False: 'single_precision'}


def generate_fields(stencil, data_type, field_layout='fzyx'):
    q = len(stencil)
    dim = len(stencil[0])

    fields = {}
    # Symbols for PDF (twice, due to double buffering)
    fields['pdfs'] = ps.Field.create_generic(
        'pdfs',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['pdfs_tmp'] = ps.Field.create_generic(
        'pdfs_tmp',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['velocity'] = ps.Field.create_generic(
        'velocity',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )
    fields['force'] = ps.Field.create_generic(
        'force',
        dim,
        data_type,
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
                                 method=lbmpy.Method.TRT,
                                 compressible=True,
                                 delta_equilibrium=True,
                                 zero_centered=True,
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


def generate_pack_info_field_specifications(
        stencil, data_type, layout, vec_len=3):
    import collections
    import itertools
    field = ps.Field.create_generic(
        "field",
        3,
        data_type,
        index_dimensions=1,
        layout=layout,
        index_shape=(vec_len,)
    )
    q = len(stencil)
    coord = itertools.product(*[(-1, 0, 1)] * 3)
    if q == 19:
        dirs = tuple((i, j, k) for i, j, k in coord if i**2 + j**2 + k**2 != 3)
    else:
        dirs = tuple((i, j, k) for i, j, k in coord)
    spec = collections.defaultdict(set)
    spec[dirs] = {field[0, 0, 0](i) for i in range(vec_len)}
    return spec


def generate_setters(lb_method, data_type):
    fields = generate_fields(lb_method.stencil, data_type)

    initial_rho = sp.Symbol('rho_0')
    pdfs_setter = lbmpy.macroscopic_value_kernels.macroscopic_values_setter(
        lb_method,
        initial_rho,
        fields['velocity'].center_vector,
        fields['pdfs'].center_vector)
    return pdfs_setter
