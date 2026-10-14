#
# Copyright (C) 2026 The ESPResSo project
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
Boundary conditions on the links between fluid cells and boundary cells, for
the CPU and CUDA targets: the class ``DynamicUBB`` of ``walberla::lbm``,
which has the interface of the class that ``pystencils_walberla`` used to
generate (the value callback of the constructor, ``run``, ``inner``,
``outer``, ``getSweep``, ``fillFromFlagField``, ``getIndexVector``,
``getForceVector``, ``getForce``, ``neighborOffset``, ...).

Design
------
The class is split in a part that doesn't depend on the boundary condition
and a part that does:

* The hand-written CRTP base class ``walberla::IndexVectorBoundary`` of
  ``src/walberla_bridge/src/utils/IndexVectorBoundary.hpp`` collects the links
  from the flag field into index vectors, keeps them on the device if needed,
  and calls the kernel on a vector of links (see its documentation for the
  semantics of the ``ALL``, ``INNER`` and ``OUTER`` vectors and of the force
  vectors). The generated header includes it by its relative path, because
  the include path of the bridge targets doesn't contain ``src/``.
* The generated class derives from it and provides the kernel
  ``linkKernel`` that processes an array of links, the function
  ``assignValue`` that stores the value returned by the value callback in a
  link, the field type ``PdfField_T`` and the stencil offsets.
* `generate_link_boundary` generates this class for assignments over an
  index field and a force field. It is independent of the boundary condition,
  but not fully generic: it requires a force field (and emits a
  ``<name>ForceStruct``), the base class hard-codes the members ``F_0`` to
  ``F_2`` of the force struct and a mandatory force vector, and the class of
  the target field is built from the data type of the code generation
  context and the number of directions of the stencil (LB-shaped: one
  population for each direction in each cell). The boundaries of the
  electrokinetics (a scalar value for each link, no force) need a change of
  the header and a generator with an optional force field; both are deferred
  to sub-project 3. `generate_dynamic_ubb` is the thin layer that defines the
  dynamic velocity bounce-back boundary (a lattice Boltzmann boundary, with a
  velocity for each link) by means of lbmpy.

Index and force structs
-----------------------
A link is a struct ``<name>IndexInfo`` of ``walberla::lbm`` with the
members ``x``, ``y``, ``z`` and ``dir`` and the additional data of the
boundary (``vel_0``, ``vel_1``, ``vel_2`` of the UBB); the force on a link is
a struct ``<name>ForceStruct`` with the members ``F_0``, ``F_1`` and ``F_2``
(for the force vector of lbmpy's ``calculate_force_on_boundary``). Both are
declared with a pystencils ``PsStructType`` that has the qualified name of
the C++ struct, and are emitted into the header together with the
constructor and ``operator==`` that the old code had. The class defines the
aliases ``IndexInfo`` and ``ForceStruct`` of the structs, so that the old
names ``DynamicUBBDoublePrecision::IndexInfo`` still work.
The coordinates ``x``, ``y``, ``z`` and ``dir`` are ``int32_t``, as in the
old code, which keeps the struct small (40 bytes for the UBB in double
precision). pystencils requires the coordinates to have the index type of
the kernel, so the strides and all the index arithmetic of the kernel are
``int32_t`` too, while the old kernel computed the offsets in the PDF field
with 64 bits: the PDF field of a block (with ghost layers and all
populations) must have fewer than 2^31 elements. The generated code doesn't
check this; the bridge has to check the size of the blocks (the check
will be added to ``LBWalberlaImpl`` in sub-project 2c-3).

Kernels
-------
The link kernel is a sparse pystencils kernel: pystencils 2.1 turns the
assignments over a field with ``FieldType.INDEXED`` and the
``CreateKernelConfig.index_field`` option into a loop over the elements of
the index field. The CPU kernel is parallelized with OpenMP and the CUDA
kernel uses one thread for each link. `lbm_link_assignments` builds the
assignments exactly like ``create_lattice_boltzmann_boundary_kernel`` of
lbmpy 2.1.1, which compiles the kernel itself and doesn't return the
assignments, while pystencils-sfg needs the assignments to create a kernel
handle. The code relies on the following internals; re-check them whenever
the versions in ``sfg_driver.TOOLCHAIN_REQUIREMENTS`` or ``WALBERLA_COMMIT``
change:

* lbmpy: ``BetweenTimestepsIndexing`` (``proxy_fields``, ``dir_symbol``,
  ``inverse_dir_symbol``, ``substitute_proxies``, ``create_code_node``),
  the call signature of the boundary classes of ``lbmpy.boundaries`` and
  their ``additional_data``, ``get_additional_code_nodes`` and the
  ``all_assignments`` of their result; the branch
  ``add_subexpressions_for_field_reads`` of
  ``create_lattice_boltzmann_boundary_kernel`` (taken when the dtype of the
  PDF field differs from the default dtype) is omitted, which is why
  `generate_link_boundary` checks the dtype of the PDF field;
* pystencils: the ``index_field``, ``index_dtype`` and
  ``skip_independence_check`` options of the kernel configuration, and the
  ``PsStructType`` of the elements of an indexed field;
* pystencils-sfg: ``sfg.kernel_namespace``, ``sfg.struct``, ``sfg.klass``,
  and ``sfg.map_field`` with an object that implements
  ``SupportsFieldExtraction`` (`ArrayExtraction`);
* sweepgen: ``sweepgen.build_config.get_build_config`` and the field
  pointer types ``GhostLayerFieldPtr`` and ``GpuFieldPtr``.

Kernel namespaces
-----------------
The kernels of a class live in the namespace ``<ClassName>_kernels`` nested
in ``walberla::lbm``. The CPU and CUDA classes have different names, so
their kernels have different qualified names, which the linker needs in
builds with CUDA (``tests/test_boundary_symbols.py``).
"""

import dataclasses

import numpy as np
from lbmpy.advanced_streaming import Timestep
from lbmpy.advanced_streaming.indexing import BetweenTimestepsIndexing
from lbmpy.boundaries import UBB
from pystencils import Assignment, Field, FieldType
from pystencils.types import PsCustomType, PsPointerType, PsStructType
from pystencils.types import create_type
from pystencilssfg.lang import AugExpr, Ref
from pystencilssfg.lang.gpu import CudaAPI
from sweepgen.api import GhostLayerFieldPtr, GpuFieldPtr
from sweepgen.build_config import get_build_config


BASE_CLASS = "walberla::IndexVectorBoundary"
# relative to the directory of the generated files, as seen by the compiler:
# the include path of the bridge targets doesn't contain ``src/``
BASE_HEADER = "../../utils/IndexVectorBoundary.hpp"
NAMESPACE = "walberla::lbm"
INDEX_TYPE = create_type("int32")


@dataclasses.dataclass(frozen=True)
class LinkValue:
    """
    The value that the value callback returns for a boundary cell, and how it
    is stored in a link: ``members`` are pairs of a member of the index
    struct and a C++ expression of the value ``v``.
    """
    cpp_type: str
    members: tuple


def index_struct_type(name, additional_data):
    """Type of a link: the coordinates, the direction and ``additional_data``,
    a list of pairs of a member name and its type."""
    members = [(member, INDEX_TYPE) for member in ("x", "y", "z", "dir")]
    return PsStructType(members + list(additional_data), name)


def force_struct_type(name, dtype, dim):
    """Type of the force on a link: one component for each dimension."""
    return PsStructType([(f"F_{i}", dtype) for i in range(dim)], name)


def link_fields(name, dtype, dim, additional_data):
    """
    The index field and the force field of a boundary class ``name``: indexed
    fields of the structs ``<NAMESPACE>::<name>IndexInfo`` and
    ``<NAMESPACE>::<name>ForceStruct``.
    """
    index_type = index_struct_type(f"{NAMESPACE}::{name}IndexInfo",
                                   additional_data)
    force_type = force_struct_type(f"{NAMESPACE}::{name}ForceStruct", dtype,
                                   dim)
    return (Field.create_generic("indexVector", 1, index_type,
                                 field_type=FieldType.INDEXED),
            Field.create_generic("forceVector", 1, force_type,
                                 field_type=FieldType.INDEXED))


def lbm_link_assignments(boundary, lb_method, pdf_field, index_field,
                         force_field, streaming_pattern="pull"):
    """
    Assignments of the kernel of the lbmpy boundary condition ``boundary``,
    as built by ``create_lattice_boltzmann_boundary_kernel`` of lbmpy 2.1.1.

    The ``add_subexpressions_for_field_reads`` branch of lbmpy, which is taken
    when the dtype of the PDF field differs from the default dtype, is
    intentionally omitted; `generate_link_boundary` rejects PDF fields whose
    dtype differs from the data type of the code generation context.
    """
    indexing = BetweenTimestepsIndexing(
        pdf_field, lb_method.stencil, Timestep.BOTH, streaming_pattern,
        np.int32, np.int32)
    f_out, f_in = indexing.proxy_fields
    dir_symbol = indexing.dir_symbol
    boundary_assignments = boundary(
        f_out, f_in, dir_symbol, indexing.inverse_dir_symbol, lb_method,
        index_field, force_field)
    boundary_assignments = indexing.substitute_proxies(boundary_assignments)
    elements = list(indexing.create_code_node().get_array_declarations())
    for node in boundary.get_additional_code_nodes(lb_method)[::-1]:
        elements += node.get_array_declarations()
    elements.append(Assignment(dir_symbol, index_field[0]("dir")))
    return elements + boundary_assignments.all_assignments


class ArrayExtraction:
    """
    Field extraction (``SupportsFieldExtraction`` of pystencils-sfg) of a
    one-dimensional array given by a pointer variable and a size variable,
    for ``sfg.map_field``: the stride is one element.
    """

    def __init__(self, pointer, size):
        self._pointer = pointer
        self._size = size

    def _extract_ptr(self):
        return AugExpr.format("{}", self._pointer)

    def _extract_size(self, coordinate):
        return AugExpr.format("{}", self._size) if coordinate == 0 else None

    def _extract_stride(self, coordinate):
        return AugExpr.format("1") if coordinate == 0 else None


def _emit_struct(sfg, struct_type, constructor_members):
    """
    Emit the struct of ``struct_type`` with the constructor that takes the
    first ``constructor_members`` members as arguments and initializes the
    others to zero, and an ``operator==`` that compares the floating-point
    members with ``walberla::floatIsEqual``.
    """
    members = [(m.name, m.dtype) for m in struct_type.members]
    name = struct_type.name.rsplit("::", 1)[-1]
    arguments = [sfg.var(f"{member}_", dtype)
                 for member, dtype in members[:constructor_members]]
    constructor = sfg.constructor(*arguments)
    for (member, dtype), argument in zip(members, arguments):
        constructor = constructor.init(sfg.var(member, dtype))(argument)
    for member, dtype in members[constructor_members:]:
        constructor = constructor.init(sfg.var(member, dtype))(
            f"{dtype.c_string()}(0)")
    comparisons = " && ".join(
        f"walberla::floatIsEqual({member}, o.{member})" if dtype.is_float()
        else f"{member} == o.{member}" for member, dtype in members)
    sfg.struct(name)(
        *[sfg.var(member, dtype) for member, dtype in members],
        constructor,
        f"bool operator==(const {name} &o) const {{ return {comparisons}; }}")


def _offsets_table(stencil):
    """C++ initializer of the array ``neighborOffset[dim][q]``."""
    dim = len(stencil[0])
    rows = (", ".join(str(int(velocity[i])) for velocity in stencil)
            for i in range(dim))
    return ",\n".join("{" + row + "}" for row in rows)


def generate_link_boundary(sfg, ctx, name, *, assignments, index_field,
                           force_field, pdf_field, stencil, value):
    """
    Generate the structs of a link and of the force on it, the link kernel
    and the boundary class ``name`` that derives from
    ``walberla::IndexVectorBoundary`` (target CPU or CUDA, see
    ``ctx``), all in ``NAMESPACE``.

    :param assignments: the kernel over ``index_field`` (the links) and
        ``force_field``, that also reads and writes the field ``pdf_field``
    :param index_field: indexed field of the structs of the links, as returned
        by `link_fields`
    :param force_field: indexed field of the structs of the forces
    :param pdf_field: the field the links refer to, a pystencils field of
        ``len(stencil)`` components for each cell, of the data type of
        ``ctx``; the class of the target field is derived from it
    :param stencil: the velocities of the lattice, a sequence of tuples
    :param value: the `LinkValue` of the value callback
    :raises ValueError: if ``pdf_field`` doesn't have ``len(stencil)``
        components or its dtype differs from ``ctx.data_type``
    """
    dim = len(stencil[0])
    n_directions = len(stencil)
    if tuple(pdf_field.index_shape) != (n_directions,):
        raise ValueError(
            f"the PDF field '{pdf_field.name}' has the index shape "
            f"{tuple(pdf_field.index_shape)}, but the stencil has "
            f"{n_directions} directions")
    if pdf_field.dtype.c_string() != ctx.data_type.c_string():
        raise ValueError(
            f"the dtype of the PDF field '{pdf_field.name}' "
            f"({pdf_field.dtype.c_string()}) differs from the data type of "
            f"the code generation context ({ctx.data_type.c_string()})")
    dtype = pdf_field.dtype.c_string()
    index_type, force_type = index_field.dtype, force_field.dtype
    index_name = index_type.name.rsplit("::", 1)[-1]
    force_name = force_type.name.rsplit("::", 1)[-1]

    config = get_build_config(sfg).get_pystencils_config()
    config.index_field = index_field
    config.index_dtype = INDEX_TYPE
    config.skip_independence_check = True

    # core/DataTypes.h defines walberla::floatIsEqual, used by the operator==
    # of the structs of the links and of the forces
    for header in ("<array>", "<cstddef>", "<cstdint>", "<functional>",
                   "<memory>", "core/DataTypes.h", "core/cell/Cell.h",
                   "core/math/Vector3.h",
                   "blockforest/StructuredBlockForest.h",
                   f"stencil/D{dim}Q{n_directions}.h", BASE_HEADER):
        sfg.include(header)

    with sfg.namespace(NAMESPACE):
        _emit_struct(sfg, index_type, 4)
        _emit_struct(sfg, force_type, 0)

        kernel = sfg.kernel_namespace(f"{name}_kernels").create(
            assignments, "linkKernel", config)

        if ctx.is_gpu:
            pdf_ptr = GpuFieldPtr.create(pdf_field)
            memory = "walberla::IndexVectorMemory::Device"
            pdf_class = f"walberla::gpu::GPUField< {dtype} >"
        else:
            pdf_ptr = GhostLayerFieldPtr.create(pdf_field)
            memory = "walberla::IndexVectorMemory::Host"
            pdf_class = (f"walberla::field::GhostLayerField< {dtype}, "
                         f"{n_directions} >")
        index_ptr = sfg.var("indexVector_", PsPointerType(
            PsCustomType(index_name)))
        force_ptr = sfg.var("forceVector_", PsPointerType(
            PsCustomType(force_name)))
        size = sfg.var("size", PsCustomType("std::size_t"))
        parameters = [pdf_ptr, index_ptr, force_ptr, size]
        if ctx.is_gpu:
            stream = CudaAPI.stream_t().var("stream")
            parameters.append(stream)
            launch = sfg.gpu_invoke(kernel, stream=stream)
        else:
            launch = sfg.call(kernel)

        base = (f"{BASE_CLASS}< {name}, {index_name}, {force_name}, "
                f"walberla::stencil::D{dim}Q{n_directions}, {value.cpp_type}, "
                f"{memory} >")
        callback = (
            f"std::function< {value.cpp_type}(const walberla::Cell &, "
            "const walberla::shared_ptr< walberla::StructuredBlockForest > &, "
            "walberla::IBlock &) >")
        assign_value = "\n".join(f"e.{member} = {expr};"
                                 for member, expr in value.members)
        sfg.klass(name, bases=[f"public {base}"])(
            sfg.public(
                f"using Base = {base};",
                f"using IndexInfo = {index_name};",
                f"using ForceStruct = {force_name};",
                f"using PdfField_T = {pdf_class};",
                sfg.constructor(
                    sfg.var("blocks", Ref(PsCustomType(
                        "std::shared_ptr< walberla::StructuredBlockForest >",
                        const=True))),
                    sfg.var("pdfsID_", PsCustomType("walberla::BlockDataID")),
                    sfg.var("valueCallback", Ref(PsCustomType(
                        callback, const=True))),
                ).init("Base")("blocks", "pdfsID_", "valueCallback",
                               f'"{name}"'),
                f"static void assignValue(IndexInfo &e, {value.cpp_type} "
                f"const &v) {{\n{assign_value}\n}}",
                sfg.method("linkKernel").params(*parameters)(
                    sfg.map_field(pdf_field, pdf_ptr.with_cell_interval(None)),
                    sfg.map_field(index_field,
                                  ArrayExtraction(index_ptr, size)),
                    sfg.map_field(force_field,
                                  ArrayExtraction(force_ptr, size)),
                    launch),
                f"static constexpr std::array< std::array< int, "
                f"{n_directions}u >, {dim}u > neighborOffset = {{{{\n"
                f"{_offsets_table(stencil)}\n}}}};",
            ))


def generate_dynamic_ubb(sfg, ctx, name, lb_method, pdf_field):
    """
    Generate the dynamic velocity bounce-back boundary ``name`` of the
    lattice Boltzmann method ``lb_method`` acting on ``pdf_field``: the
    velocity of each boundary link is returned by the value callback, and the
    force on each link is stored in the force vector.
    """
    dim = lb_method.dim
    # the velocity is only used to tell lbmpy that it is taken from the
    # index field: the callback is not called by the generated code
    ubb = UBB(lambda *args: None, dim=dim, data_type=ctx.data_type,
              calculate_force_on_boundary=True, name=name)
    index_field, force_field = link_fields(
        name, ctx.data_type, dim, ubb.additional_data)
    stencil = [tuple(int(c) for c in velocity) for velocity in
               lb_method.stencil]
    generate_link_boundary(
        sfg, ctx, name,
        assignments=lbm_link_assignments(
            ubb, lb_method, pdf_field, index_field, force_field),
        index_field=index_field, force_field=force_field,
        pdf_field=pdf_field, stencil=stencil,
        value=LinkValue(
            f"walberla::Vector3< {ctx.data_type.c_string()} >",
            tuple((f"vel_{i}", f"v[{i}]") for i in range(dim))))
