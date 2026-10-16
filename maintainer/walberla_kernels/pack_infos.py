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
Pack infos of the waLBerla bridge for the CPU and CUDA targets: the classes
``PackInfoPdf`` and ``PackInfoVec`` of namespace ``walberla::pystencils``,
which have the interface of the classes that ``pystencils_walberla`` used to
generate (constructor from a ``BlockDataID``, ``pack``, ``unpack``, ``size``
and ``communicateLocal``).

A pack info is described by the field components that it communicates in each
direction, a dictionary ``{direction: (component, ...)}``:

* `pdf_streaming_sets`: the populations of a PDF field that cross the face,
  edge or corner of a direction during pull streaming (streaming-aware pack
  info), so that a face sends 5 of the 19 populations of D3Q19;
* `vector_sets`: all components of a vector field, for each direction of
  the stencil.

`FieldPackInfo` generates, for each distinct set of components, a pack, an
unpack and a local copy kernel, and a class that dispatches on the direction.
The class derives from a CRTP base that implements the waLBerla pack info
interface and calls ``doPack``, ``doUnpack``, ``doLocalCopy`` and
``elementsPerCell`` of the generated class:

* CPU: ``walberla::UniformCpuFieldPackInfoBase`` of
  ``src/walberla_bridge/src/utils/UniformCpuFieldPackInfoBase.hpp``.
* CUDA: ``walberla::sweepgen::UniformGpuFieldPackInfoBase`` of the waLBerla
  fork that sweepgen comes from (commit ``WALBERLA_COMMIT`` of `sfg_driver`).

The code follows ``sweepgen.communication.pack_infos.GpuPdfFieldPackInfo``,
which only supports PDFs on CUDA. It relies on the following sweepgen and
pystencils-sfg 0.3 internals; re-check them whenever ``WALBERLA_COMMIT`` or the
pystencils-sfg version change:

* `sweepgen.communication.pack_infos.PackingKernelsContext` and
  `CaptureBufferPointer`, which bind the buffer field of the kernels to the
  ``std::span`` that the pack info methods receive;
* `sweepgen.api.GhostLayerFieldExtraction`, whose ``_extract_size`` and
  ``_extract_stride`` methods `ReversedAxesExtraction` overrides;
* ``sfg.kernel_namespace``, ``sfg.klass`` and ``sfg.switch`` of the stock
  pystencils-sfg composer.

Buffer layout
-------------
pystencils linearizes the index of a buffer field in the order of the loops:
the first axis of the symbolic field varies slowest in the buffer. The
waLBerla fields of ESPResSo are ``fzyx`` fields, whose first axis would be
``x``, so that the buffer would be x-slowest and the loops would run over
``x`` first. The old pack infos stored the buffer x-fastest, with the
components of a cell next to each other, which matches the memory order of
the field. Therefore the symbolic fields are declared with the axes
``(z, y, x)``, and `ReversedAxesExtraction` maps these to the sizes and
strides of the waLBerla field in reverse order (the component axis is not
reversed). The loops of the generated kernels then run over ``z``, ``y`` and
``x`` from the outside to the inside, and the buffer is x-fastest. The layout
only matters within one run, because the same pack info class fills and
empties the buffer.

Kernel namespaces
-----------------
The kernels live in the namespace ``<ClassName>_kernels`` nested in
``walberla::pystencils``. The CPU and CUDA classes have different names, so
their kernels have different qualified names, which the linker needs in
builds with CUDA (``tests/test_pack_info_symbols.py``).

The CPU kernels are not parallelized with OpenMP, like the old ones: they
pack thin slices of a block.
"""

import functools
import itertools
import warnings

from pystencils import Assignment, Field, FieldType
from pystencils.stencil import offset_to_direction_string
from pystencilssfg.composer.custom import CustomGenerator
from pystencilssfg.lang import strip_ptr_ref
from pystencilssfg.lang.cpp import std
from pystencilssfg.lang.gpu import CudaAPI
from sweepgen.api import (
    CellInterval, Direction, GhostLayerFieldExtraction, GhostLayerFieldPtr,
    GpuFieldPtr, uint_t)
from sweepgen.build_config import get_build_config
from sweepgen.communication.pack_infos import (
    CaptureBufferPointer, PackingKernelsContext)

import pystencils_espresso
import sfg_driver

CPU_BASE_CLASS = "walberla::UniformCpuFieldPackInfoBase"
# relative to the directory of the generated files, as seen by the compiler:
# the include path of the bridge targets doesn't contain ``src/``
CPU_BASE_HEADER = "../../utils/UniformCpuFieldPackInfoBase.hpp"
GPU_BASE_CLASS = "walberla::sweepgen::UniformGpuFieldPackInfoBase"
GPU_BASE_HEADER = (
    "walberla/sweepgen/communication/UniformGpuFieldPackInfoBase.hpp")

# all directions of a 3D block, in the order of the cases of the generated code
COMMUNICATION_DIRECTIONS = tuple(
    d for d in itertools.product((-1, 0, 1), repeat=3) if d != (0, 0, 0))


def pdf_streaming_sets(stencil):
    """
    The populations that cross the face, edge or corner of each communication
    direction during pull streaming.

    A population with velocity ``e`` is advected out of the block through
    the boundary in direction ``d`` if ``e`` and ``d`` agree in all non-zero
    components of ``d``. Directions without such populations are not
    communicated (the corners of D3Q19).

    :returns: dict from the direction to the sorted indices in ``stencil``
    """
    sets = {}
    for direction in COMMUNICATION_DIRECTIONS:
        indices = tuple(
            index for index, velocity in enumerate(stencil)
            if all(v == d for v, d in zip(velocity, direction) if d != 0))
        if indices:
            sets[direction] = indices
    return sets


def vector_sets(stencil, n_components):
    """
    All components of a field, in all communication directions that are
    velocities of the stencil.

    :returns: dict from the direction to the component indices
    """
    velocities = set(map(tuple, stencil))
    return {d: tuple(range(n_components)) for d in COMMUNICATION_DIRECTIONS
            if d in velocities}


class ReversedAxesExtraction(GhostLayerFieldExtraction):
    """
    Binds the symbolic field axes ``(a0, a1, a2)`` to ``(z, y, x)`` of a
    waLBerla field: sizes and strides of the spatial axes are taken in
    reverse order. The axis of the components is not changed.
    """

    @staticmethod
    def _reverse(coordinate):
        return 2 - coordinate if coordinate < 3 else coordinate

    def _extract_size(self, coordinate):
        return super()._extract_size(self._reverse(coordinate))

    def _extract_stride(self, coordinate):
        return super()._extract_stride(self._reverse(coordinate))


class FieldPackInfo(CustomGenerator):
    """
    Pack info class ``name`` for the 3D field ``field`` (a pystencils field
    with one index dimension, or a scalar field without index dimension,
    whose only component is 0), which communicates the components
    ``comm_sets[direction]`` in each direction (see `pdf_streaming_sets`).
    The target (CPU or CUDA) is the one of the build configuration of the
    ``sfg`` the generator is run with.
    """

    def __init__(self, name, field, comm_sets):
        self._name = name
        self._dtype = field.dtype
        self._comm_sets = comm_sets
        self._scalar = field.index_dimensions == 0
        if field.index_dimensions > 1:
            raise ValueError(
                f"Field {field.name} has {field.index_dimensions} index "
                "dimensions, at most one is supported")
        n_components = 1 if self._scalar else field.index_shape[0]
        # axes (z, y, x), see the module docstring
        self._src = self._symbolic_field(f"{field.name}_src", n_components)
        self._dst = self._symbolic_field(f"{field.name}_dst", n_components)

    def _symbolic_field(self, name, n_components):
        if self._scalar:
            return Field.create_generic(
                name, 3, self._dtype, index_dimensions=0, layout="c")
        return Field.create_generic(
            name, 3, self._dtype, index_dimensions=1, layout="c",
            index_shape=(n_components,))

    def _access(self, symbolic_field, component):
        """Center access of ``component`` (none for a scalar field)."""
        if self._scalar:
            assert component == 0
            return symbolic_field.center
        return symbolic_field.center(component)

    def _buffer(self, n_components):
        return Field.create_generic(
            "buffer", 1, field_type=FieldType.BUFFER, dtype=self._dtype,
            index_shape=(n_components,))

    def _create_kernels(self, kernels_context, directions):
        """
        Pack, unpack and local copy kernels of all ``directions`` that
        communicate the same components.
        """
        components = self._comm_sets[directions[0]]
        buffer = self._buffer(len(components))
        name = "_".join(offset_to_direction_string(d) for d in directions)
        create = kernels_context.kns.create
        config = kernels_context.cfg
        src = functools.partial(self._access, self._src)
        dst = functools.partial(self._access, self._dst)
        pack = [Assignment(buffer(i), src(c))
                for i, c in enumerate(components)]
        unpack = [Assignment(dst(c), buffer(i))
                  for i, c in enumerate(components)]
        local_copy = [Assignment(dst(c), src(c)) for c in components]
        return (create(pack, f"pack{name}", config),
                create(unpack, f"unpack{name}", config),
                create(local_copy, f"localCopy{name}", config))

    def generate(self, sfg):
        config = get_build_config(sfg).get_pystencils_config()
        is_gpu = config.get_target().is_gpu()
        if is_gpu:
            base = f"{GPU_BASE_CLASS}< {self._name} >"
            sfg.include(GPU_BASE_HEADER)
            field_ptr = GpuFieldPtr.create
        else:
            config.cpu.openmp.enable = False
            base = f"{CPU_BASE_CLASS}< {self._name} >"
            sfg.include(CPU_BASE_HEADER)
            field_ptr = GhostLayerFieldPtr.create
        kernels_context = PackingKernelsContext(
            sfg, kns=sfg.kernel_namespace(f"{self._name}_kernels"),
            cfg=config)

        # directions with the same components share their kernels
        by_components = {}
        for direction, components in self._comm_sets.items():
            by_components.setdefault(components, []).append(direction)
        kernels = {}
        for directions in by_components.values():
            handles = self._create_kernels(kernels_context, directions)
            for direction in directions:
                kernels[direction] = handles

        src_ptr = field_ptr(self._src)
        dst_ptr = field_ptr(self._dst)
        buffer = std.span(self._dtype).var("buffer")
        dir_var = Direction().var("dir")
        src_interval = CellInterval(const=True, ref=True).var("srcInterval")
        dst_interval = CellInterval(const=True, ref=True).var("dstInterval")
        stream = CudaAPI.stream_t().var("stream") if is_gpu else None
        stream_params = [] if stream is None else [stream]

        def bind(symbolic_field, ptr, interval):
            return sfg.map_field(
                symbolic_field, ReversedAxesExtraction(ptr, interval))

        def dispatch(index):
            """Switch over the directions that calls kernel ``index``."""
            if is_gpu:
                def call(handle):
                    return sfg.gpu_invoke(handle, stream=stream)
            else:
                call = sfg.call
            return sfg.switch(dir_var).cases({
                Direction.from_offset(d): call(kernels[d][index])
                for d in self._comm_sets}).default("/* unreachable */")

        elements_per_cell = sfg.switch(dir_var, autobreak=False).cases({
            Direction.from_offset(d): f"return {len(components)};"
            for d, components in self._comm_sets.items()
        }).default("return 0;")

        sfg.klass(self._name, bases=[f"public {base}"])(
            sfg.public(
                f"using Base = {base};",
                "using Base::Base;",
                f"using Field_T = "
                f"{strip_ptr_ref(src_ptr.get_dtype()).c_string()};",
                sfg.method("doPack").const().params(
                    src_ptr, buffer, dir_var, src_interval, *stream_params)(
                    bind(self._src, src_ptr, src_interval),
                    CaptureBufferPointer(sfg, "buffer", buffer),
                    dispatch(0)),
                sfg.method("doUnpack").const().params(
                    dst_ptr, buffer, dir_var, dst_interval, *stream_params)(
                    bind(self._dst, dst_ptr, dst_interval),
                    CaptureBufferPointer(sfg, "buffer", buffer),
                    dispatch(1)),
                sfg.method("doLocalCopy").const().params(
                    src_ptr, src_interval, dst_ptr, dst_interval, dir_var,
                    *stream_params)(
                    bind(self._src, src_ptr, src_interval),
                    bind(self._dst, dst_ptr, dst_interval),
                    dispatch(2)),
                sfg.method("elementsPerCell").inline().const().params(
                    dir_var).returns(uint_t)(elements_per_cell),
            ))


def _generate(sfg, ctx, name, field, comm_sets):
    assert ctx.target in ("cpu", "cuda"), ctx.target
    with warnings.catch_warnings():
        # pystencils-sfg 0.3 warns about the kernel parameters that carry
        # the pointer property of the buffer, which is intended here
        warnings.filterwarnings(
            "ignore", message="Encountered two non-identical variables",
            category=UserWarning)
        with sfg.namespace(sfg_driver.KERNEL_NAMESPACE):
            sfg.generate(FieldPackInfo(name, field, comm_sets))


def generate_pdf_pack_info(sfg, ctx, name, stencil):
    """
    Generate the streaming-aware pack info ``name`` of the PDF field of the
    lattice of ``stencil``, in ``walberla::pystencils``.
    """
    fields = pystencils_espresso.generate_fields(stencil, ctx.data_type)
    _generate(sfg, ctx, name, fields["pdfs"], pdf_streaming_sets(stencil))


def generate_vector_pack_info(sfg, ctx, name, stencil):
    """
    Generate the pack info ``name`` of the vector fields (force, velocity)
    of the lattice of ``stencil``, in ``walberla::pystencils``.
    """
    fields = pystencils_espresso.generate_fields(stencil, ctx.data_type)
    _generate(sfg, ctx, name, fields["force"],
              vector_sets(stencil, fields["force"].index_shape[0]))
