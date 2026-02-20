#
# Copyright (C) 2021-2025 The ESPResSo project
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
import re
import lbmpy
import pystencils as ps


# File derived from lbmpy_walberla.walberla_lbm_generation in the
# walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d


def remove_intermediate_variable(code, name):
    re_pat = re.compile(rf"\n *const (float|double|u?int(?:(?:8|16|32|64)?_t)?) {name}(?: *\[\d*\])? = .*?;(?=\n)")  # nopep8
    assert re_pat.search(code) is not None, \
        f"pattern '{re_pat}' not found in '''\n{code}\n'''"
    return re_pat.sub("", code)


def add_espresso_filters_to_jinja_env(jinja_env):
    jinja_env.filters["remove_intermediate_variable"] = remove_intermediate_variable


def generate_accessors(ctx, config, templates):

    # Function derived from lbmpy_walberla.walberla_lbm_generation.__lattice_model()
    # in the walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d
    # with backports from commit de6b00071233a9a1f45d7a6773988363e058f1a0

    from jinja2 import Environment, FileSystemLoader, StrictUndefined
    from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env

    default_dtype = config.data_type.default_factory()

    stencils = list(lbmpy.stencils.LBStencil("D3Q27").stencil_entries)
    stencils_helper = list(map(ps.stencil.offset_to_direction_string,
                           lbmpy.stencils.LBStencil("D3Q27")))
    stencils_helper = dict(zip(stencils_helper, stencils))

    staggered_stencils_helper = ps.field.Field.create_generic("tmp",
                                                              spatial_dimensions=3, index_shape=(13,),
                                                              field_type=ps.field.FieldType.STAGGERED_FLUX).staggered_stencil
    staggered_stencils = dict([(stencils_helper[dir], i)
                              for i, dir in enumerate(staggered_stencils_helper)])

    inverse_stencils_helper = list(
        lbmpy.stencils.LBStencil("D3Q27").inverse_stencil_entries)
    inverse_staggered_stencils = {}
    for i, dir in enumerate(stencils):
        if dir in staggered_stencils:
            inverse_staggered_stencils[inverse_stencils_helper[i]
                                       ] = staggered_stencils[dir]

    jinja_context = {
        "dtype": default_dtype,
        "namespace": "ek",
        "D": 3,
        "FluxCount": 13,
        "Stencils": stencils,
        "InverseStencils": inverse_staggered_stencils,
        "StaggeredStencils": staggered_stencils,
    }

    env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)),
                      undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)
    add_espresso_filters_to_jinja_env(env)

    for filename, template in templates.items():
        source = env.get_template(template).render(**jinja_context)
        ctx.write_file(filename, source)
