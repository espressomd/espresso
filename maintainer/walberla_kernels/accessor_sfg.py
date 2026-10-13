#
# Copyright (C) 2021-2026 The ESPResSo project
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

"""
pystencils-sfg helpers shared by the generators of the field accessors.

The generators add free functions to a namespace of the generated file with
the stock pystencils-sfg composer, which lets the functions of a namespace
overload each other. Relies on pystencils-sfg 0.3 internals; re-check them
whenever ``TOOLCHAIN_REQUIREMENTS`` in `sfg_driver` changes:

* `AllowOverloads`: ``sfg._cursor``, ``SfgFunctionSequencer`` and the
  ``add_entity`` method of the cursor, which registers a function by name
  and rejects a second function of this name;
* `LinearBuffer`: the ``SupportsFieldExtraction`` protocol (``_extract_ptr``,
  ``_extract_size``, ``_extract_stride``) that ``sfg.map_field`` uses to bind
  the buffer to a pystencils field;
* `statement`: ``SfgStatements``.
"""

from pystencils.types import PsCustomType
from pystencilssfg.composer.basic_composer import SfgFunctionSequencer
from pystencilssfg.ir import SfgStatements
from pystencilssfg.lang import (
    AugExpr, SupportsFieldExtraction, asvar, cpptype, depends, includes)

ARRAY = cpptype("std::array< {T}, {N}u >", "<array>")
MATRIX3 = cpptype("walberla::Matrix3< {T} >", "core/math/Matrix3.h")
AUTO = PsCustomType("auto")
UINT = PsCustomType("uint")


class AllowOverloads:
    """Cursor proxy that doesn't register functions by name, so that a
    namespace can hold overloads (pystencils-sfg 0.3 rejects a second
    function of the same name)."""

    def __init__(self, cursor):
        self._cursor = cursor

    def add_entity(self, entity):
        pass

    def __getattr__(self, name):
        return getattr(self._cursor, name)


def function(sfg, name, *params, returns=None, inline=False):
    """
    Start a free function of the current namespace of ``sfg`` that may
    overload another one. The returned sequencer is called with the body.
    """
    sequencer = SfgFunctionSequencer(AllowOverloads(sfg._cursor), name)
    if inline:
        sequencer = sequencer.inline()
    if returns is not None:
        sequencer = sequencer.returns(returns)
    return sequencer.params(*params)


def statement(code, defines=(), uses=()):
    """C++ statement with the variables it defines and uses."""
    deps = set()
    incls = set()
    for expr in uses:
        deps |= depends(expr)
        incls |= includes(expr)
    return SfgStatements(code, [asvar(d) for d in defines], deps, incls)


class LinearBuffer(AugExpr, SupportsFieldExtraction):
    """
    Pointer to a dense array with ``fsize`` values per element, viewed as a
    pystencils field with spatial extents ``sizes`` (1 to 3 expressions):
    the first spatial coordinate is the slowest and the value index the
    fastest, the layout of the buffers of the legacy accessors.
    """

    def __init__(self, ptr_type, sizes, fsize):
        super().__init__(ptr_type)
        self._sizes = [AugExpr.format("{}", s) for s in sizes]
        self._fsize = fsize

    def _extract_ptr(self):
        return AugExpr.format("{}", self)

    def _extract_size(self, coordinate):
        if coordinate < len(self._sizes):
            return self._sizes[coordinate]
        return AugExpr.format("{}", self._fsize)

    def _extract_stride(self, coordinate):
        rank = len(self._sizes)
        if coordinate == rank:
            return AugExpr.format("1")
        stride = AugExpr.format("{}", self._fsize)
        for size in self._sizes[coordinate + 1:]:
            stride = AugExpr.format("{} * {}", stride, size)
        return stride


def kernel_param(kernel, name):
    """The scalar parameter ``name`` of a kernel, or None."""
    for param in kernel.scalar_parameters:
        if param.name == name:
            return param
    return None
