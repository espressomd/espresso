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
import re
import sympy as sp
import pystencils as ps
import lbmpy_walberla
from pystencils import TypedSymbol
try:
    from pystencils.typing import CastFunc
except ImportError:
    from pystencils.sympyextensions import CastFunc
try:
    from pystencils.typing import BasicType as PsScalarType
except ImportError:
    from pystencils.types import PsScalarType


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
    from sympy.tensor import IndexedBase
    from pystencils.backends.cbackend import CustomSympyPrinter
    from pystencils.backends.cbackend import CBackend
    from pystencils.backends.cuda_backend import CudaBackend
    from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env

    cpp_printer = CustomSympyPrinter()

    default_dtype = config.data_type.default_factory()
    if config.target == ps.Target.GPU:
        backend = CudaBackend()
    else:
        backend = CBackend()

    jinja_context = {
        "dtype": default_dtype,
        "namespace": "ek",
        "D": 3,
        "FluxCount": 13,
    }

    env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)),
                      undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)
    add_espresso_filters_to_jinja_env(env)

    for filename, template in templates.items():
        source = env.get_template(template).render(**jinja_context)
        ctx.write_file(filename, source)
