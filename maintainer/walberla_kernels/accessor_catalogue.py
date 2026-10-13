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
Format of the function catalogue of the lattice field accessors.

The accessors of the waLBerla bridge (``walberla::lbm::accessor``, later
also those of the electrokinetics) read and write the fields of a lattice
for single cells, cell intervals or whole fields. Each accessor is
described once, as an `Accessor` in a `Namespace` of a `Catalogue`; the
backends `accessor_backend_cpu` and `accessor_backend_cuda` turn the
catalogue into the C++ wrappers of one target. This module defines the
format only: it knows neither pystencils-sfg nor a lattice model. The
content of the catalogue (fields, equations, kernels) is supplied by
a specification module such as `lb_accessors`.

An accessor is a C++ function of a namespace that runs one kernel on a cell,
on an interval or on a whole field. The kernel is either shared by all
accessors of the namespace that name it (the kernels of a namespace are
created before its first function) or belongs to a single accessor (it is
created right before the function, which fixes its position in the
generated code). The kind of the accessor (`Kind`) decides how the wrapper
is built:

* ``CELL``, ``INTERVAL``: the kernel runs on the single cell ``cell`` or on
  the interval ``ci``. The kernel exchanges ``n`` values per cell with the
  caller through a buffer: ``out`` is returned to the caller, ``inp`` is
  an argument of the function;
* ``ALL``: the kernel runs on the whole field, ghost layers included, and
  sets the components of ``inp`` (``inp`` itself if it is a scalar) as scalar
  kernel parameters (broadcast);
* ``REDUCE``: the kernel runs on the interior of the first field and the
  values of the cells are summed in the order of the cells;
* ``POINT``: a function of its arguments without field access;
* ``LIST``, ``GATHER``, ``SCATTER``: kernels over lists of cells or of
  particle positions (set the values of a list of cells, interpolate field
  values at positions, spread values onto the cells around positions).

Per-target differences are expressed by `Accessor.targets`: an accessor
exists on the targets that are keys of the mapping, and a `Variant` renames
parameters (the API of each target is frozen: ``tests/data/*.json``) and
overrides the return type.
"""

import dataclasses
import enum
import typing

import pystencils as ps

__all__ = [
    "Kind", "Container", "SCALAR", "VECTOR3", "MATRIX3", "STD_VECTOR",
    "array", "Value", "FieldSpec", "FieldArg", "Variant", "BOTH", "CPU_ONLY",
    "CUDA_ONLY", "Point", "Kernel", "ReductionKernel", "Reduction", "Accessor",
    "Namespace", "Catalogue", "TARGETS"]

TARGETS = ("cpu", "cuda")


class Kind(enum.Enum):
    CELL = "cell"
    INTERVAL = "interval"
    ALL = "all"
    REDUCE = "reduce"
    POINT = "point"
    LIST = "list"
    GATHER = "gather"
    SCATTER = "scatter"


@dataclasses.dataclass(frozen=True)
class Container:
    """
    C++ value type exchanged with the caller of an accessor: ``kind`` is
    one of ``scalar``, ``array`` (``std::array<T, size>``), ``vector3``
    (``walberla::Vector3<T>``), ``matrix3`` (``walberla::Matrix3<T>``) and
    ``std_vector`` (``std::vector<T>``, any number of values).
    """
    kind: str
    size: typing.Optional[int] = None


SCALAR = Container("scalar", 1)
VECTOR3 = Container("vector3", 3)
MATRIX3 = Container("matrix3", 9)
STD_VECTOR = Container("std_vector")


def array(size):
    """The container ``std::array<T, size>``."""
    return Container("array", size)


@dataclasses.dataclass(frozen=True)
class Value:
    """
    Values exchanged with the caller: the C++ ``container`` holds ``n``
    values per cell, ``name`` is the variable name (of the parameter or
    of the local variable that is returned).
    """
    container: Container
    n: int
    name: str


@dataclasses.dataclass(frozen=True)
class FieldSpec:
    """
    A field of the lattice: the pystencils field (with ``fsize`` values per
    cell) and, for the kernels over lists, the same field with
    absolute indexing.
    """
    field: ps.Field
    fsize: int
    absolute: typing.Optional[ps.Field] = None


@dataclasses.dataclass(frozen=True)
class FieldArg:
    """
    Parameter of an accessor that is a field of the lattice: ``role`` is a
    key of `Catalogue.fields`, ``name`` the default name of the parameter.
    """
    role: str
    const: bool
    name: str


@dataclasses.dataclass(frozen=True)
class Variant:
    """
    Spelling of an accessor on one target: ``names`` maps the default name
    of a parameter to the name it has on the target, ``returns`` is the
    spelling of the return type if it differs from the type of the
    returned container (``"auto"``).
    """
    names: typing.Mapping[str, str] = dataclasses.field(default_factory=dict)
    returns: typing.Optional[str] = None


BOTH = {"cpu": Variant(), "cuda": Variant()}
CPU_ONLY = {"cpu": Variant()}
CUDA_ONLY = {"cuda": Variant()}


@dataclasses.dataclass(frozen=True)
class Kernel:
    """
    A kernel of a namespace. ``build()`` returns the pystencils assignments.
    For the kernels of the accessors of kind ``ALL``, ``components`` is
    ``(prefix, n)``: the kernel has the scalar parameters ``prefix_0`` ...
    ``prefix_{n-1}`` (see `Catalogue.component_symbols`), which the wrapper
    sets from the components of the input value.
    """
    name: str
    build: typing.Callable
    components: typing.Optional[tuple] = None
    targets: tuple = TARGETS


@dataclasses.dataclass(frozen=True)
class ReductionKernel:
    """
    Kernel that sums a value per cell over a field. ``build()`` returns the
    pair of the assignments that precede the summation and the list of the
    values that are summed: a backend chooses how to sum them.
    """
    name: str
    build: typing.Callable
    targets: tuple = TARGETS


@dataclasses.dataclass(frozen=True)
class Reduction:
    """
    Result of a reduction: the sums are the components of ``container``
    (named ``name`` in the wrapper), multiplied by the scalar parameter
    ``scale`` of the function if it is not None.
    """
    container: Container
    name: str
    scale: typing.Optional[str] = None


@dataclasses.dataclass(frozen=True)
class Point:
    """
    Data of an accessor of the kind ``POINT``: a function of the stencil
    direction ``direction``, the velocity ``u`` (components ``u_names``) and
    the density ``rho`` that evaluates the expression of the direction
    (``cases``: pairs of the name of the direction and an expression).
    """
    u_names: tuple
    cases: tuple


@dataclasses.dataclass(frozen=True)
class Accessor:
    """
    One accessor function (an overload): see the module documentation.

    The function parameters are, in this order: the lattice fields
    ``fields``, the value ``inp``, the scalars ``scalars`` and the cell
    ``cell`` or the interval ``ci`` (kinds ``CELL`` and ``INTERVAL``).
    A scalar is a parameter of the kernel if the kernel uses a symbol of this
    name, else it is a plain by-value parameter. ``point`` is the `Point` of
    an accessor of kind ``POINT``, ``dims`` the number of values per item
    of the list kinds.
    """
    name: str
    kind: Kind
    kernel: typing.Union[str, Kernel, ReductionKernel, None] = None
    fields: tuple = ()
    inp: typing.Optional[Value] = None
    out: typing.Optional[Value] = None
    scalars: tuple = ()
    reduction: typing.Optional[Reduction] = None
    point: typing.Any = None
    dims: typing.Optional[int] = None
    targets: typing.Mapping[str, Variant] = dataclasses.field(
        default_factory=lambda: dict(BOTH))

    def default_names(self):
        """
        Default names of the parameters that the catalogue defines, in the
        order of the signature. The accessors of the kinds ``LIST``,
        ``GATHER`` and ``SCATTER`` have more parameters (lists, the number
        of ghost layers) that their backend adds.
        """
        domain = {Kind.CELL: ("cell",), Kind.INTERVAL: ("ci",)}
        return (*(f.name for f in self.fields),
                *((self.inp.name,) if self.inp is not None else ()),
                *self.scalars, *domain.get(self.kind, ()))

    def param_name(self, default, target):
        """Name of the parameter ``default`` on ``target``."""
        return self.targets[target].names.get(default, default)

    def parameter_names(self, target):
        """The names of the parameters `default_names` on ``target``."""
        return tuple(self.param_name(n, target) for n in self.default_names())


@dataclasses.dataclass(frozen=True)
class Namespace:
    """Accessors and their kernels, in the order of the generated code."""
    name: str
    kernels: tuple
    accessors: tuple

    def for_target(self, target):
        """Kernels and accessors that exist on ``target``."""
        kernels = tuple(k for k in self.kernels if target in k.targets)
        accessors = tuple(a for a in self.accessors if target in a.targets)
        return kernels, accessors


@dataclasses.dataclass(frozen=True)
class Catalogue:
    """
    The accessors of one lattice model.

    Attributes:
        dtype: floating-point type of the fields
        namespace: C++ namespace that contains the namespaces of the
            accessors
        prelude: comment at the top of the generated files
        includes: waLBerla headers the accessors need on each target, in
            addition to those of the backend
        fields: the `FieldSpec` of the fields of the lattice by role; the
            keys ``pos``, ``indices``, ``values`` and ``out<n>`` are the
            fields of the buffers of the kernels over lists
        dim: number of space dimensions
        namespaces: tuple of `Namespace`
        component_symbols: function ``(prefix, n) -> symbols`` that returns
            the scalar kernel parameters of the kernels of kind ``ALL``
        buffer: function ``n -> field`` that returns the pystencils field
            of the buffer with ``n`` values per cell through which the
            kernels of the kinds ``CELL`` and ``INTERVAL`` exchange values
            with the caller
    """
    dtype: typing.Any
    namespace: str
    prelude: str
    includes: typing.Mapping[str, tuple]
    fields: typing.Mapping[str, FieldSpec]
    dim: int
    namespaces: tuple
    component_symbols: typing.Callable
    buffer: typing.Callable

    def __post_init__(self):
        """Check that the catalogue is consistent."""
        for ns in self.namespaces:
            names = [k.name for k in ns.kernels]
            assert len(names) == len(set(names)), (ns.name, names)
            signatures = set()
            for acc in ns.accessors:
                self._check(ns, acc)
                key = (acc.name, acc.kind, acc.fields, acc.inp and
                       acc.inp.container, acc.scalars, tuple(acc.targets))
                assert key not in signatures, (
                    f"{ns.name}::{acc.name} is defined twice")
                signatures.add(key)

    def _check(self, ns, acc):
        where = f"{ns.name}::{acc.name}"
        assert set(acc.targets) <= set(TARGETS), where
        for field in acc.fields:
            assert field.role in self.fields, (where, field.role)
        for target, variant in acc.targets.items():
            assert set(variant.names) <= set(acc.default_names()), (
                where, target, variant.names)
        if isinstance(acc.kernel, str):
            kernel = [k for k in ns.kernels if k.name == acc.kernel]
            assert len(kernel) == 1, (where, acc.kernel)
            assert set(acc.targets) <= set(kernel[0].targets), (
                where, acc.kernel)
        elif acc.kernel is not None:
            assert set(acc.targets) <= set(acc.kernel.targets), where

    def for_target(self, target):
        """
        The namespaces that have accessors on ``target``, each reduced to
        its kernels and accessors of ``target``.
        """
        assert target in TARGETS, target
        result = []
        for ns in self.namespaces:
            kernels, accessors = ns.for_target(target)
            if accessors:
                result.append(Namespace(ns.name, kernels, accessors))
        return tuple(result)

    def accessors(self, target):
        """All ``(namespace name, accessor)`` of ``target``."""
        return [(ns.name, a) for ns in self.for_target(target)
                for a in ns.accessors]
