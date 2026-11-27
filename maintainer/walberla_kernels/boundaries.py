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
Boundary conditions on the links between cells and boundary cells, for the
CPU and CUDA targets: the class ``DynamicUBB`` of ``walberla::lbm``, which has
the interface of the class that ``pystencils_walberla`` used to generate (the
value callback of the constructor, ``run``, ``inner``, ``outer``,
``getSweep``, ``fillFromFlagField``, ``getIndexVector``, ``getForceVector``,
``getForce``, ``neighborOffset``, ...), and the classes ``Dirichlet`` and
``FixedFlux`` of ``walberla::ek``, which have the interface of the old classes
of the electrokinetics (the same, without the force).

Design
------
The class is split in a part that doesn't depend on the boundary condition
and a part that does:

* The hand-written CRTP base class ``walberla::IndexVectorBoundary`` of
  ``src/walberla_bridge/src/utils/IndexVectorBoundary.hpp`` collects the links
  from the flag field into index vectors, keeps them on the device if needed,
  and calls the kernel on a vector of links (see its documentation for the
  semantics of the ``ALL``, ``INNER`` and ``OUTER`` vectors, of the force
  vectors and of the link policies). The generated header includes it by its
  relative path, because the include path of the bridge targets doesn't
  contain ``src/``.
* The generated class derives from it and provides the kernel
  ``linkKernel`` that processes an array of links, the function
  ``assignValue`` that stores the value returned by the value callback in a
  link, the field type ``PdfField_T`` and the stencil offsets.
* `generate_link_boundary` generates this class for assignments over an
  index field and, optionally, a force field. It is independent of the
  boundary condition. Without a force field (``ForceStruct = void`` in the
  base class), there are no force vectors, no ``<name>ForceStruct`` and no
  force argument of the kernel. The link policy of the base class is a
  parameter. The class of the target field is built from the data type of
  the field and its number of components.
  `generate_dynamic_ubb` is the thin layer that defines the dynamic velocity
  bounce-back boundary (a lattice Boltzmann boundary, with a velocity for each
  link and the force on the boundary, one population for each direction in
  each cell) by means of lbmpy, `generate_dirichlet` the boundary of the
  electrokinetics that sets the density of the boundary cells (a scalar value
  for each link, no force) and `generate_fixed_flux` the one that sets the
  flux through the faces between a cell and a boundary cell (a vector for
  each link, no force).

Index and force structs
-----------------------
A link is a struct ``<name>IndexInfo`` of the namespace of the class
(``walberla::lbm`` or ``walberla::ek``) with the members ``x``, ``y``, ``z``
and ``dir`` and the additional data of the boundary (``vel_0``, ``vel_1``,
``vel_2`` of the UBB, ``value`` of the Dirichlet boundary, ``flux_0``,
``flux_1``, ``flux_2`` of the fixed flux boundary); the force on a
link is a struct ``<name>ForceStruct`` with the members ``F_0``, ``F_1`` and
``F_2`` (for the force vector of lbmpy's ``calculate_force_on_boundary``).
Both are declared with a pystencils ``PsStructType`` that has the qualified
name of the C++ struct, and are emitted into the header together with the
constructor and ``operator==`` that the old code had. The class defines the
aliases ``IndexInfo`` and ``ForceStruct`` of the structs, so that the old
names ``DynamicUBBDoublePrecision::IndexInfo`` still work.
The coordinates ``x``, ``y``, ``z`` and ``dir`` are ``int32_t``, as in the
old code, which keeps the struct small (40 bytes for the UBB in double
precision). pystencils requires the coordinates to have the index type of
the kernel, so the strides and all the index arithmetic of the kernel are
``int32_t`` too, while the old kernel computed the offsets in the field
with 64 bits: the field of a block (with ghost layers and all components)
must have fewer than 2^31 elements. The generated code doesn't check this;
the bridge has to check the size of the blocks. The LB bridge checks the PDF
field (``LBWalberlaImpl``, see the Readme), and the flux field of the
electrokinetics, which has 13 components, with the density field
(``EKinWalberlaImpl``).

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
the versions in ``sfg_driver.TESTED_TOOLCHAIN`` or the waLBerla commit
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

Fixed flux
----------
The class ``FixedFlux`` of ``walberla::ek`` sets the flux through the faces
between a domain cell ``c`` and its neighbors in the boundary: the flux
through the face between ``c`` and ``c + d`` along ``d`` is ``v . d / 9``
for the value ``v`` that the callback returns for the boundary cell
``c + d`` (``9`` is the square of the number of dimensions, as in the old
class). The link policy is ``DomainLinksNoCenter`` (the links start in all
domain cells, also in boundary cells) and the stencil is the ``D3Q27``
stencil of waLBerla. The flux is stored in the flux field of the
electrokinetics (see `flux_layout`): the cell that owns the face, the
component and the sign of the stored flux depend on the direction ``d`` of
the link, which the kernel takes from tables that are indexed with ``dir``.
The face between two boundary cells is written by two links, one from each
cell, with the values of the two boundary cells: the kernels may write the
same element from two links (see ``IndexVectorBoundary.hpp``), as the old
kernels did. With the same value in both cells, which is the case for a
boundary that sets the same flux on all its cells (e.g. no flux), the result
doesn't depend on the order of the links; with different values, the value of
one of the two cells is stored, and it depends on the order of the links
(on the CPU, on the schedule of the OpenMP threads; on the GPU, on the
threads) which one.

Kernel namespaces
-----------------
The kernels of a class live in the namespace ``<ClassName>_kernels`` nested
in the namespace of the class (``walberla::lbm`` or ``walberla::ek``). The
CPU and CUDA classes have different names, so their kernels have different
qualified names, which the linker needs in builds with CUDA
(``tests/test_boundary_symbols.py``).
"""

import dataclasses

import lbmpy
import numpy as np
import sympy as sp
from lbmpy.advanced_streaming import Timestep
from lbmpy.advanced_streaming.indexing import BetweenTimestepsIndexing
from lbmpy.boundaries import UBB
from pystencils import Assignment, Field, FieldType, TypedSymbol
from pystencils.defaults import DEFAULTS
from pystencils.types import PsCustomType, PsPointerType, PsStructType
from pystencils.types import create_type
from pystencils.types.quick import Arr
from pystencilssfg.lang import AugExpr, Ref
from pystencilssfg.lang.gpu import CudaAPI
from sweepgen.api import GhostLayerFieldPtr, GpuFieldPtr
from sweepgen.build_config import get_build_config

import flux_layout


BASE_CLASS = "walberla::IndexVectorBoundary"
# relative to the directory of the generated files, as seen by the compiler:
# the include path of the bridge targets doesn't contain ``src/``
BASE_HEADER = "../../utils/IndexVectorBoundary.hpp"
NAMESPACE = "walberla::lbm"  # of the lattice Boltzmann boundaries
EK_NAMESPACE = "walberla::ek"
INDEX_TYPE = create_type("int32")
# the link policies of IndexVectorBoundary.hpp; the default one is not named
# in the generated base class, so that the generated files of the LB
# boundaries don't depend on the policies
FLUID_LINKS = "FluidLinks"
DOMAIN_LINKS_NO_CENTER = "DomainLinksNoCenter"
DOMAIN_LINKS_WITH_GHOST_LAYERS = "DomainLinksWithGhostLayers"
LINK_POLICIES = (FLUID_LINKS, DOMAIN_LINKS_NO_CENTER,
                 DOMAIN_LINKS_WITH_GHOST_LAYERS)


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


def link_index_field(name, additional_data, namespace=NAMESPACE):
    """
    The index field of a boundary class ``name``: an indexed field of the
    struct ``<namespace>::<name>IndexInfo``.
    """
    index_type = index_struct_type(f"{namespace}::{name}IndexInfo",
                                   additional_data)
    return Field.create_generic("indexVector", 1, index_type,
                                field_type=FieldType.INDEXED)


def lbm_link_fields(name, dtype, dim, additional_data):
    """
    The index field and the force field of a lattice Boltzmann boundary
    class ``name`` of ``NAMESPACE``: indexed fields of the structs
    ``<NAMESPACE>::<name>IndexInfo`` and ``<NAMESPACE>::<name>ForceStruct``.
    """
    force_type = force_struct_type(f"{NAMESPACE}::{name}ForceStruct", dtype,
                                   dim)
    return (link_index_field(name, additional_data),
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
                           force_field=None, field, stencil, value,
                           links=FLUID_LINKS):
    """
    Generate the structs of a link and (with a force field) of the force on
    it, the link kernel and the boundary class ``name`` that derives from
    ``walberla::IndexVectorBoundary`` (target CPU or CUDA, see ``ctx``), all
    in the namespace of the struct of the links (see `link_index_field`).

    :param assignments: the kernel over ``index_field`` (the links) and
        ``force_field``, that also reads and writes the field ``field``
    :param index_field: indexed field of the structs of the links
    :param force_field: indexed field of the structs of the forces, or
        ``None`` for a boundary without force: no struct, no force vectors,
        and the kernel has no force argument
    :param field: the field the links refer to, a pystencils field of the
        data type of ``ctx``; the class of the target field is derived from
        it (``GhostLayerField`` or ``GPUField``, with the number of its
        components if it has index dimensions)
    :param stencil: the velocities of the links, a sequence of tuples: the
        stencil of waLBerla with this number of directions and dimensions
        (e.g. ``D3Q19``), or ``walberla::CellStencil`` for the single
        velocity zero, the link of a cell to itself. Only the number of
        directions and dimensions is validated: the order of the velocities
        must be the direction order of waLBerla's ``D{d}Q{q}`` stencil, which
        the generated ``neighborOffset`` and the direction index of the links
        refer to
    :param value: the `LinkValue` of the value callback
    :param links: the link policy of ``IndexVectorBoundary.hpp``, one of
        `LINK_POLICIES`; ``DomainLinksWithGhostLayers`` requires the
        single velocity zero, and the single velocity zero requires
        ``DomainLinksWithGhostLayers`` (the other policies never produce a
        link with it)
    :raises ValueError: if the dtype of ``field`` differs from
        ``ctx.data_type``, the policy is unknown, the stencil of one
        velocity isn't zero, or the policy doesn't match the stencil
    """
    dim = len(stencil[0])
    n_directions = len(stencil)
    if field.dtype.c_string() != ctx.data_type.c_string():
        raise ValueError(
            f"the dtype of the field '{field.name}' "
            f"({field.dtype.c_string()}) differs from the data type of "
            f"the code generation context ({ctx.data_type.c_string()})")
    if links not in LINK_POLICIES:
        raise ValueError(f"unknown link policy {links}, expected one of "
                         f"{', '.join(LINK_POLICIES)}")
    cell_stencil = n_directions == 1
    if cell_stencil and any(stencil[0]):
        raise ValueError(
            f"a stencil of one direction must be the velocity zero, not "
            f"{stencil[0]}")
    if links == DOMAIN_LINKS_WITH_GHOST_LAYERS and not cell_stencil:
        raise ValueError(
            f"the link policy {links} requires a stencil of one direction")
    if cell_stencil and links != DOMAIN_LINKS_WITH_GHOST_LAYERS:
        # a link to the cell itself needs a boundary cell that is also a
        # domain cell: fluid cells are never boundary cells, and the policy
        # without the center skips the only direction
        raise ValueError(
            f"the link policy {links} never produces a link with a stencil "
            f"of one direction, use {DOMAIN_LINKS_WITH_GHOST_LAYERS}")
    dtype = field.dtype.c_string()
    n_components = int(np.prod(field.index_shape)) if field.index_shape else 1
    index_type = index_field.dtype
    index_name = index_type.name.rsplit("::", 1)[-1]
    namespace = index_type.name.rsplit("::", 1)[0]
    has_force = force_field is not None
    if has_force:
        force_type = force_field.dtype
        force_name = force_type.name.rsplit("::", 1)[-1]
    stencil_class = ("walberla::CellStencil" if cell_stencil
                     else f"walberla::stencil::D{dim}Q{n_directions}")

    config = get_build_config(sfg).get_pystencils_config()
    config.index_field = index_field
    config.index_dtype = INDEX_TYPE
    config.skip_independence_check = True

    # core/DataTypes.h defines walberla::floatIsEqual, used by the operator==
    # of the structs of the links and of the forces
    headers = ["<array>", "<cstddef>", "<cstdint>", "<functional>",
               "<memory>", "core/DataTypes.h", "core/cell/Cell.h",
               "core/math/Vector3.h", "blockforest/StructuredBlockForest.h"]
    if not cell_stencil:
        headers.append(f"stencil/D{dim}Q{n_directions}.h")
    for header in headers + [BASE_HEADER]:
        sfg.include(header)

    with sfg.namespace(namespace):
        _emit_struct(sfg, index_type, 4)
        if has_force:
            _emit_struct(sfg, force_type, 0)

        kernel = sfg.kernel_namespace(f"{name}_kernels").create(
            assignments, "linkKernel", config)

        if ctx.is_gpu:
            field_ptr = GpuFieldPtr.create(field)
            memory = "walberla::IndexVectorMemory::Device"
            field_class = f"walberla::gpu::GPUField< {dtype} >"
        else:
            field_ptr = GhostLayerFieldPtr.create(field)
            memory = "walberla::IndexVectorMemory::Host"
            field_class = (f"walberla::field::GhostLayerField< {dtype}, "
                           f"{n_components} >")
        index_ptr = sfg.var("indexVector_", PsPointerType(
            PsCustomType(index_name)))
        size = sfg.var("size", PsCustomType("std::size_t"))
        parameters = [field_ptr, index_ptr]
        mappings = [
            sfg.map_field(field, field_ptr.with_cell_interval(None)),
            sfg.map_field(index_field, ArrayExtraction(index_ptr, size))]
        if has_force:
            force_ptr = sfg.var("forceVector_", PsPointerType(
                PsCustomType(force_name)))
            parameters.append(force_ptr)
            mappings.append(sfg.map_field(
                force_field, ArrayExtraction(force_ptr, size)))
        parameters.append(size)
        if ctx.is_gpu:
            stream = CudaAPI.stream_t().var("stream")
            parameters.append(stream)
            launch = sfg.gpu_invoke(kernel, stream=stream)
        else:
            launch = sfg.call(kernel)

        policy = "" if links == FLUID_LINKS else f", walberla::{links}"
        base = (f"{BASE_CLASS}< {name}, {index_name}, "
                f"{force_name if has_force else 'void'}, {stencil_class}, "
                f"{value.cpp_type}, {memory}{policy} >")
        callback = (
            f"std::function< {value.cpp_type}(const walberla::Cell &, "
            "const walberla::shared_ptr< walberla::StructuredBlockForest > &, "
            "walberla::IBlock &) >")
        assign_value = "\n".join(f"e.{member} = {expr};"
                                 for member, expr in value.members)
        field_id = f"{field.name}ID_"
        sfg.klass(name, bases=[f"public {base}"])(
            sfg.public(
                f"using Base = {base};",
                f"using IndexInfo = {index_name};",
                *((f"using ForceStruct = {force_name};",) if has_force
                  else ()),
                f"using PdfField_T = {field_class};",
                sfg.constructor(
                    sfg.var("blocks", Ref(PsCustomType(
                        "std::shared_ptr< walberla::StructuredBlockForest >",
                        const=True))),
                    sfg.var(field_id, PsCustomType("walberla::BlockDataID")),
                    sfg.var("valueCallback", Ref(PsCustomType(
                        callback, const=True))),
                ).init("Base")("blocks", field_id, "valueCallback",
                               f'"{name}"'),
                f"static void assignValue(IndexInfo &e, {value.cpp_type} "
                f"const &v) {{\n{assign_value}\n}}",
                sfg.method("linkKernel").params(*parameters)(
                    *mappings, launch),
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

    :raises ValueError: if ``pdf_field`` doesn't have one component for each
        direction of the stencil of ``lb_method``
    """
    dim = lb_method.dim
    stencil = [tuple(int(c) for c in velocity) for velocity in
               lb_method.stencil]
    if tuple(pdf_field.index_shape) != (len(stencil),):
        raise ValueError(
            f"the PDF field '{pdf_field.name}' has the index shape "
            f"{tuple(pdf_field.index_shape)}, but the stencil has "
            f"{len(stencil)} directions")
    # the velocity is only used to tell lbmpy that it is taken from the
    # index field: the callback is not called by the generated code
    ubb = UBB(lambda *args: None, dim=dim, data_type=ctx.data_type,
              calculate_force_on_boundary=True, name=name)
    index_field, force_field = lbm_link_fields(
        name, ctx.data_type, dim, ubb.additional_data)
    generate_link_boundary(
        sfg, ctx, name,
        assignments=lbm_link_assignments(
            ubb, lb_method, pdf_field, index_field, force_field),
        index_field=index_field, force_field=force_field,
        field=pdf_field, stencil=stencil,
        value=LinkValue(
            f"walberla::Vector3< {ctx.data_type.c_string()} >",
            tuple((f"vel_{i}", f"v[{i}]") for i in range(dim))))


def generate_dirichlet(sfg, ctx, name):
    """
    Generate the Dirichlet boundary ``name`` of the electrokinetics in
    ``EK_NAMESPACE``: the value that the callback returns for a boundary cell
    (the density) is stored in the field of that cell, in the interior and in
    the ghost layers. The links of the boundary cells point to themselves
    (``walberla::CellStencil``, ``DomainLinksWithGhostLayers``), and there is
    no force.
    """
    field = Field.create_generic(
        "field", 3, ctx.data_type, index_dimensions=0, layout="fzyx")
    index_field = link_index_field(
        name, [("value", ctx.data_type)], EK_NAMESPACE)
    generate_link_boundary(
        sfg, ctx, name,
        assignments=[Assignment(field.center, index_field[0]("value"))],
        index_field=index_field, field=field, stencil=[(0, 0, 0)],
        links=DOMAIN_LINKS_WITH_GHOST_LAYERS,
        value=LinkValue(ctx.data_type.c_string(), (("value", "v"),)))


def fixed_flux_assignments(dtype, flux_field, index_field, stencil):
    """
    The assignments of the kernel of the fixed flux boundary over the links
    of ``index_field`` (the struct with the members ``flux_0``, ``flux_1``
    and ``flux_2``) for the velocities ``stencil`` (the stencil of the links
    in the order of waLBerla), which store ``v . d / 9`` as the flux through
    the face between the cell of a link and its neighbor in the direction
    ``d``. The cell, the component and the sign of the flux field are given
    by `flux_layout.face_location`; the kernel looks them up in tables that
    are indexed with the direction of the link. There is no link in the
    direction zero (the tables have zeros).
    """
    dim = flux_layout.DIM
    index_type = create_type("int32")
    n_directions = len(stencil)
    locations = [flux_layout.face_location(d) if any(d) else ((0,) * dim, 0, 0)
                 for d in stencil]
    arrays = []
    owners, coefficients = [], []
    for i, axis in enumerate("xyz"):
        owner = TypedSymbol(f"owner_{axis}", Arr(index_type, n_directions))
        arrays.append(Assignment(owner, tuple(o[i] for o, _, _ in locations)))
        owners.append(sp.IndexedBase(owner, shape=(1,)))
        coefficient = TypedSymbol(
            f"coefficient_{axis}", Arr(dtype, n_directions))
        arrays.append(Assignment(coefficient, tuple(
            sp.Rational(sign * d[i], dim ** 2)
            for d, (_, _, sign) in zip(stencil, locations))))
        coefficients.append(sp.IndexedBase(coefficient, shape=(1,)))
    component = TypedSymbol("component", Arr(index_type, n_directions))
    arrays.append(Assignment(component, tuple(c for _, c, _ in locations)))
    component = sp.IndexedBase(component, shape=(1,))
    link = index_field[0]
    direction = link("dir")
    # the counters of a sparse kernel are the coordinates of the link
    cell = DEFAULTS.spatial_counters
    value = [link(f"flux_{i}") for i in range(dim)]
    target = flux_field.absolute_access(
        tuple(cell[i] + owners[i][direction] for i in range(dim)),
        (component[direction],))
    flux = sum(coefficients[i][direction] * value[i] for i in range(dim))
    return arrays + [Assignment(target, flux)]


def generate_fixed_flux(sfg, ctx, name):
    """
    Generate the fixed flux boundary ``name`` of the electrokinetics in
    ``EK_NAMESPACE`` (see the section "Fixed flux" of the module
    documentation): the value that the callback returns for a boundary cell
    (a vector) sets the flux through the faces between the boundary cell and
    its neighbors in the domain, in the 13-component flux field. The links
    have all the directions of the ``D3Q27`` stencil but the center
    (``DomainLinksNoCenter``), and there is no force.
    """
    dim = flux_layout.DIM
    stencil = [tuple(int(c) for c in velocity)
               for velocity in lbmpy.LBStencil(lbmpy.Stencil.D3Q27)]
    flux_field = Field.create_generic(
        "flux", dim, ctx.data_type, index_dimensions=1,
        index_shape=(flux_layout.FLUX_COUNT,), layout="fzyx",
        field_type=FieldType.CUSTOM)
    index_field = link_index_field(
        name, [(f"flux_{i}", ctx.data_type) for i in range(dim)],
        EK_NAMESPACE)
    generate_link_boundary(
        sfg, ctx, name,
        assignments=fixed_flux_assignments(
            ctx.data_type, flux_field, index_field, stencil),
        index_field=index_field, field=flux_field, stencil=stencil,
        links=DOMAIN_LINKS_NO_CENTER,
        value=LinkValue(
            f"walberla::Vector3< {ctx.data_type.c_string()} >",
            tuple((f"flux_{i}", f"v[{i}]") for i in range(dim))))
