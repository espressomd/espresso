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
pystencils-sfg helpers shared by the backends of the field accessors
(`accessor_backend_cpu`, `accessor_backend_cuda`) and the base class of the
backends, `AccessorBackend`.

The backends add free functions to a namespace of the generated file with
the stock pystencils-sfg composer, which lets the functions of a namespace
overload each other and takes the C++ types from the catalogue
(`accessor_catalogue`). Relies on pystencils-sfg 0.3 internals; re-check
them whenever ``TESTED_TOOLCHAIN`` in `sfg_driver` changes:

* `AllowOverloads`: ``sfg._cursor``, ``SfgFunctionSequencer`` and the
  ``add_entity`` method of the cursor, which registers a function by name
  and rejects a second function of this name;
* `LinearBuffer`: the ``SupportsFieldExtraction`` protocol (``_extract_ptr``,
  ``_extract_size``, ``_extract_stride``) that ``sfg.map_field`` uses to bind
  the buffer to a pystencils field;
* `statement`: ``SfgStatements``.
"""

import dataclasses

import pystencils as ps
from pystencils.types import PsCustomType, constify
from pystencilssfg.composer.basic_composer import SfgFunctionSequencer
from pystencilssfg.ir import SfgStatements
from pystencilssfg.lang import (
    AugExpr, SupportsFieldExtraction, asvar, cpptype, depends, includes)
from pystencilssfg.lang.cpp import std
from sweepgen.api import Cell, CellInterval, Vector3

import accessor_catalogue as catalogue

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


# ---- C++ spelling of the containers of the catalogue

def n_cells(n):
    """Number of values of ``n`` per cell of the interval ``ci``."""
    return "ci.numCells()" if n == 1 else f"ci.numCells() * {n}u"


def value_var(container, dtype, name):
    """Function parameter of the type ``container`` by const reference."""
    kind = container.kind
    if kind == "array":
        return AugExpr(ARRAY(T=dtype, N=container.size, const=True,
                             ref=True)).var(name)
    if kind == "vector3":
        return Vector3(dtype, const=True, ref=True).var(name)
    if kind == "std_vector":
        return std.vector(dtype, const=True, ref=True).var(name)
    raise ValueError(f"no parameter of the container {kind}")


def return_type(container, dtype):
    """Type of a function that returns the container."""
    kind = container.kind
    if kind == "scalar":
        return dtype
    if kind == "array":
        return ARRAY(T=dtype, N=container.size)
    if kind == "vector3":
        return Vector3(dtype).dtype
    if kind == "matrix3":
        return MATRIX3(T=dtype)
    return std.vector(dtype).dtype


def cxx_type(container, T):
    """C++ spelling of the type of a container."""
    kind = container.kind
    if kind == "scalar":
        return T
    if kind == "array":
        return f"std::array<{T}, {container.size}u>"
    if kind == "vector3":
        return f"walberla::Vector3<{T}>"
    if kind == "matrix3":
        return f"walberla::Matrix3<{T}>"
    return f"std::vector<{T}>"


def local_declaration(container, T, name, n=1):
    """C++ declaration of the local variable ``name`` of the container (the
    interval ``ci`` has ``n`` values per cell for a ``std::vector``)."""
    kind = container.kind
    if kind == "scalar":
        return f"{T} {name}{{}};"
    if kind == "std_vector":
        return f"std::vector<{T}> {name}({n_cells(n)});"
    return f"{cxx_type(container, T)} {name};"


def data_pointer(container, name):
    """C++ expression of the pointer to the values of a variable."""
    return f"&{name}" if container.kind == "scalar" else f"{name}.data()"


def construction(container, T, args):
    """C++ expression that builds a vector or matrix from its values."""
    assert container.kind in ("vector3", "matrix3"), container
    return f"{cxx_type(container, T)}({', '.join(args)})"


@dataclasses.dataclass
class Params:
    """Function parameters of an accessor, in the order of the signature."""
    fields: list
    inp: object
    scalars: dict
    domain: object

    @property
    def all(self):
        return [*self.fields, *([self.inp] if self.inp is not None else []),
                *self.scalars.values(),
                *([self.domain] if self.domain is not None else [])]


class AccessorBackend:
    """
    Base class of the backends: generates the wrappers of the accessors of
    a `catalogue.Catalogue` for one target. A subclass sets ``target`` and
    ``inline``, defines the headers, the field parameter types and the kernel
    target, and implements ``emit_<kind>`` for each kind of accessor it
    supports (the name of the `catalogue.Kind`).
    """

    target = None
    # True if the backend generates device code; checked against the
    # target of the code generation context
    gpu = False
    inline = False
    # name of the namespace nested in each accessor namespace that holds the
    # pystencils kernels. The names must differ between the backends: the
    # kernels of the CPU and CUDA backends have equal signatures, and the
    # linker would merge the inline CPU kernels with the host stubs of the
    # CUDA kernels in a build with both.
    kernel_namespace = "kernels"
    # headers of the generated file, and those only the source file needs
    headers = ()
    private_headers = ()
    ps_target = None

    def __init__(self, sfg, ctx, cat):
        if ctx.is_gpu != self.gpu:
            raise ValueError(
                f"the {self.target} accessors don't support the target "
                f"{ctx.target}")
        self.sfg = sfg
        self.cat = cat
        self.dtype = cat.dtype
        self.T = self.dtype.c_string().strip()

    # ---- to define in the subclasses

    def field_param(self, role, const, name):
        """Parameter that is the field ``role`` of the lattice."""
        raise NotImplementedError

    # ---- helpers

    def lattice_fields(self, acc, params):
        """Pairs of the pystencils field and the parameter of the fields of
        an accessor."""
        return [(self.cat.fields[f.role].field, p)
                for f, p in zip(acc.fields, params.fields)]

    def kernel(self, name, assignments):
        config = ps.CreateKernelConfig(
            target=self.ps_target, default_dtype=self.dtype, ghost_layers=0)
        kernels = self.sfg.kernel_namespace(self.kernel_namespace)
        return kernels.create(assignments, name, config)

    def fn(self, name, *params, returns=None):
        """Start a function; the sequencer is called with the body."""
        return function(self.sfg, name, *params, returns=returns,
                        inline=self.inline)

    def scalar_param(self, kernel, name):
        """The scalar ``name`` of an accessor: the kernel parameter, or a
        plain parameter if the kernel doesn't take it."""
        param = kernel_param(kernel, name) if kernel is not None else None
        return param or self.sfg.var(name, constify(self.dtype))

    def params(self, acc, kernel=None):
        """The function parameters of an accessor of the kind ``CELL``,
        ``INTERVAL``, ``ALL`` or ``REDUCE`` (the catalogue only allows to
        rename the fields and the input value: the scalars and the domain
        keep their default names)."""
        def rename(name):
            return acc.param_name(name, self.target)

        fields = [self.field_param(f.role, f.const, rename(f.name))
                  for f in acc.fields]
        inp = None
        if acc.inp is not None:
            name = rename(acc.inp.name)
            if acc.inp.container.kind == "scalar":
                inp = self.sfg.var(name, constify(self.dtype))
            else:
                inp = value_var(acc.inp.container, self.dtype, name)
        scalars = {name: self.scalar_param(kernel, name)
                   for name in acc.scalars}
        domain = None
        if acc.kind == catalogue.Kind.CELL:
            domain = Cell(const=True, ref=True).var("cell")
        elif acc.kind == catalogue.Kind.INTERVAL:
            domain = CellInterval(const=True, ref=True).var("ci")
        return Params(fields, inp, scalars, domain)

    def returns(self, acc):
        """Return type of an accessor, None for ``void``."""
        variant = acc.targets[self.target]
        if variant.returns == "auto":
            return AUTO
        container = acc.out.container if acc.out is not None else (
            acc.reduction.container if acc.reduction is not None else None)
        return None if container is None else return_type(
            container, self.dtype)

    def cell_interval(self, cell):
        """Interval with the single cell ``cell``."""
        ci = CellInterval(const=True).var("ci")
        return ci, statement("walberla::CellInterval ci(cell, cell);",
                             (ci,), (cell,))

    def whole_interval(self, size, field_ptr):
        """Interval over a whole field: ``size`` is the member function of
        the field (e.g. ``xyzSize``)."""
        ci = CellInterval().var("ci")
        return ci, statement(
            f"walberla::CellInterval ci = {field_ptr}->{size}();",
            (ci,), (field_ptr,))

    def broadcast(self, acc, spec, params):
        """
        Nodes that set the scalar parameters of the kernel of a broadcast
        (kind ``ALL``) from the components of the input value (the input
        itself if it is a scalar), in the order of the components.
        """
        symbols = self.cat.component_symbols(*spec.components)
        if acc.inp.container.kind == "scalar":
            values = [AugExpr.format("{}", params.inp)]
        else:
            values = [AugExpr.format("{}[{}u]", params.inp, i)
                      for i in range(len(symbols))]
        assert len(symbols) == len(values), acc.name
        return [self.sfg.set_param(symbol, value)
                for symbol, value in zip(symbols, values)]

    # ---- generate

    def generate(self):
        sfg = self.sfg
        for header in (*self.headers, *self.cat.includes.get(self.target, ())):
            sfg.include(header)
        for header in self.private_headers:
            sfg.include(header, private=True)
        sfg.prelude(self.cat.prelude)
        sfg.namespace(self.cat.namespace)
        for ns in self.cat.for_target(self.target):
            self.generate_namespace(ns)

    def create_kernel(self, spec):
        """Create a kernel of the catalogue."""
        if isinstance(spec, catalogue.ReductionKernel):
            return self.reduction_kernel(spec)
        return self.kernel(spec.name, spec.build())

    def reduction_kernel(self, spec):
        raise NotImplementedError(
            f"the {self.target} backend must implement reduction_kernel")

    def generate_namespace(self, ns):
        kernels, accessors = ns.for_target(self.target)
        with self.sfg.namespace(ns.name):
            shared = {k.name: (k, self.create_kernel(k)) for k in kernels}
            for acc in accessors:
                if isinstance(acc.kernel, str):
                    spec, kernel = shared[acc.kernel]
                elif acc.kernel is None:
                    spec = kernel = None
                else:
                    spec, kernel = acc.kernel, self.create_kernel(acc.kernel)
                emit = getattr(self, f"emit_{acc.kind.value}", None)
                if emit is None:
                    raise ValueError(
                        f"the {self.target} backend doesn't support the "
                        f"accessor {ns.name}::{acc.name} of the kind "
                        f"{acc.kind.value}")
                emit(acc, spec, kernel)
