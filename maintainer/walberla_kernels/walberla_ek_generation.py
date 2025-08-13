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


# File derived from lbmpy_walberla.walberla_lbm_generation in the
# walberla project, commit 3455bf3eebc64efa9beaecd74ebde3459b98991d


def remove_intermediate_variable(code, name):
    re_pat = re.compile(f"const (float|double) {name} = .*?;\n")
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

    # TODO Find a better way to create these stencils lists/dirs
    # Directions according to walberla directions.h
    stencils = ["C", "N", "S", "W", "E", "T", "B",
                "NW", "NE", "SW", "SE", "TN", "TS",
                "TW", "TE", "BN", "BS", "BW", "BE",
                "TNE", "TNW", "TSE", "TSW", "BNE", "BNW", "BSE", "BSW"]
    # Staggered directions according to pystencils field.py
    staggeredStencils = {"W": 0, "S": 1, "B": 2,
                         "SW": 3, "NW": 4, "BW": 5,
                         "TW": 6, "BS": 7, "TS": 8,
                         "BSW": 9, "TSW": 10, "BNW": 11, "TNW": 12}
    # Inverse directions to the staggered access
    invStencils = {"C": -1, "N": 1, "S": 1, "W": 0, "E": 0, "T": 2, "B": 2,
                   "NW": 4, "NE": 3, "SW": 3, "SE": 4, "TN": 7, "TS": 8,
                   "TW": 6, "TE": 5, "BN": 8, "BS": 7, "BW": 5, "BE": 6,
                   "TNE": 9, "TNW": 12, "TSE": 11, "TSW": 10, "BNE": 10,
                   "BNW": 11, "BSE": 12, "BSW": 9}

    jinja_context = {
        "dtype": default_dtype,
        "namespace": "ek",
        "D": 3,
        "FluxCount": 13,
        "Stencils": stencils,
        "InverseStencils": invStencils,
        "StaggeredStencils": staggeredStencils,
    }

    env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)),
                      undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)
    add_espresso_filters_to_jinja_env(env)

    for filename, template in templates.items():
        source = env.get_template(template).render(**jinja_context)
        ctx.write_file(filename, source)
