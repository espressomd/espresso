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
CUDA backend of the field accessors: generates the wrappers of a
`accessor_catalogue.Catalogue` as functions declared in a header ``.cuh``
and defined in a source file ``.cu``, with pystencils-sfg.

Device work:

* the accessors of the kinds ``CELL``, ``INTERVAL`` and ``ALL`` run
  pystencils CUDA kernels, launched from the generated host functions
  (``sfg.gpu_invoke``). The fields are `walberla::gpu::GPUField` views over
  the requested `CellInterval`; host data is exchanged through a linear
  device buffer that is mapped to a pystencils field (`LinearBuffer`,
  x-major, f-fastest);
* the reductions (``REDUCE``) are per-cell kernels that fill a device
  buffer, which is summed on the host in the order of the cells
  (deterministic);
* the accessors over lists (``LIST``, ``GATHER``) are kernels over a 1D list
  of particle positions (or indices) that address the lattice fields with
  ``Field.absolute_access``;
* ``SCATTER`` spreads values onto the cells with ``atomicAdd``, which
  pystencils can't express: its ``__global__`` kernel is a C++ string.
  It is written for three space dimensions and uses the same
  trilinear weights as the kernels of ``GATHER``
  (``LbAccessors.interpolation_weights``).

Relies on pystencils-sfg 0.3 internals, see `accessor_sfg`.
"""

import pystencils as ps
from pystencils.types import PsCustomType, PsPointerType
from pystencilssfg.lang import AugExpr, cpptype
from pystencilssfg.lang.cpp import std
from sweepgen.api import CellInterval, GpuFieldPtr

import accessor_catalogue as catalogue
from accessor_sfg import (
    AccessorBackend, LinearBuffer, UINT, cxx_type, data_pointer,
    local_declaration, n_cells, statement)

_DeviceVector = cpptype("thrust::device_vector< {T} >",
                        "<thrust/device_vector.h>")


class CudaBackend(AccessorBackend):
    target = "cuda"
    ps_target = ps.Target.CUDA
    headers = ("core/DataTypes.h", "core/cell/Cell.h",
               "core/cell/CellInterval.h", "core/math/Matrix3.h",
               "core/math/Vector3.h", "gpu/GPUField.h",
               "<thrust/device_vector.h>", "<array>", "<vector>")
    private_headers = ("gpu/FieldAccessor.h", "gpu/FieldIndexing.h",
                       "<cassert>", "<cstdint>", "<thrust/copy.h>")

    def field_param(self, role, const, name):
        spec = self.cat.fields[role]
        return GpuFieldPtr(self.dtype, self.cat.dim, spec.fsize,
                           const=const).var(name)

    def reduction_kernel(self, spec):
        """The values of the cells are written to the buffer."""
        assignments, values = spec.build()
        buf = self.cat.buffer(len(values))
        return self.kernel(spec.name, [
            *assignments,
            *(ps.Assignment(buf(i), v) for i, v in enumerate(values))])

    # ---- helpers

    def device_vector(self, name, size_args, uses=(), shape=None, fsize=None):
        """
        Device vector ``name`` (of the working precision) and the pointer
        ``name_ptr`` to its data. With ``shape`` (spatial extents) the
        pointer is a `LinearBuffer` with ``fsize`` values per element.
        """
        T = self.T
        vector = AugExpr(_DeviceVector(T=self.dtype)).var(name)
        ptr_type = PsPointerType(self.dtype)
        if shape is None:
            ptr = AugExpr(ptr_type).var(f"{name}_ptr")
        else:
            ptr = LinearBuffer(ptr_type, shape, fsize).var(f"{name}_ptr")
        decl = [
            statement(f"thrust::device_vector<{T}> {name}({size_args});",
                      (vector,), uses),
            statement(f"{T} *{name}_ptr = "
                      f"thrust::raw_pointer_cast({name}.data());",
                      (ptr,), (vector,))]
        return ptr, decl

    def slice_body(self, kernel, fields, *, ci, buffer=None, params=(),
                   post=()):
        """
        Body of a function that runs ``kernel`` on the interval ``ci``.

        * ``fields``: pairs of a pystencils field and a GPU field,
        * ``buffer``: ``(field, size_args, uses)`` of the linear device
          buffer that is mapped to the pystencils field ``field``;
          ``size_args`` is the argument list of the constructor of the
          device vector,
        * ``params``: nodes that set the scalar parameters of the kernel,
        * ``post``: nodes that follow the launch.
        """
        nodes = []
        ptr = None
        if buffer is not None:
            field, size_args, uses = buffer
            shape = [AugExpr.format(f"{{}}.{c}Size()", ci) for c in "xyz"]
            ptr, decl = self.device_vector(
                "dev_data", size_args, (*uses, ci), shape,
                field.index_shape[0])
            nodes += decl
        nodes += [self.sfg.map_field(f, g.with_cell_interval(ci))
                  for f, g in fields]
        if buffer is not None:
            nodes.append(self.sfg.map_field(buffer[0], ptr))
        return [*nodes, *params, self.sfg.gpu_invoke(kernel), *post]

    def lattice_fields(self, acc, params):
        """Pairs of the pystencils field and the parameter of the fields of
        an accessor."""
        return [(self.cat.fields[f.role].field, p)
                for f, p in zip(acc.fields, params.fields)]

    @staticmethod
    def copy_to_host(target, vector="dev_data"):
        return f"thrust::copy({vector}.begin(), {vector}.end(), {target});"

    @staticmethod
    def constructor_arguments(value, name):
        """Arguments of the constructor of a device vector that copies the
        host value ``name``."""
        kind = value.container.kind
        if kind == "array" or kind == "std_vector":
            return f"{name}.begin(), {name}.end()"
        if kind == "vector3":
            return f"{name}.data(), {name}.data() + {value.container.size}u"
        return f"1u, {name}"

    # ---- kinds CELL and INTERVAL

    def emit_cell(self, acc, spec, kernel):
        self.emit_slice(acc, kernel)

    def emit_interval(self, acc, spec, kernel):
        self.emit_slice(acc, kernel)

    def emit_slice(self, acc, kernel):
        """Run the kernel on a cell or an interval; the values are
        exchanged through a device buffer."""
        T = self.T
        params = self.params(acc, kernel)
        cell = acc.kind == catalogue.Kind.CELL
        nodes = []
        if cell:
            ci, decl = self.cell_interval(params.domain)
            nodes.append(decl)
        else:
            ci = params.domain
        value = acc.inp or acc.out
        buffer_field = self.cat.buffer(value.n)
        post = ()
        if acc.inp is not None:
            name = acc.param_name(acc.inp.name, self.target)
            if not cell:
                nodes.append(f"assert(uint_c({name}.size()) == "
                             f"{n_cells(value.n)});")
            size_args = self.constructor_arguments(acc.inp, name)
            uses = (params.inp,)
        else:
            name = acc.out.name
            size_args = f"{value.n}u" if cell else n_cells(value.n)
            uses = ()
            post = (local_declaration(acc.out.container, T, name, value.n),
                    self.copy_to_host(data_pointer(acc.out.container, name)),
                    f"return {name};")
        self.fn(acc.name, *params.all, returns=self.returns(acc))(
            *nodes,
            *self.slice_body(
                kernel, self.lattice_fields(acc, params), ci=ci,
                buffer=(buffer_field, size_args, uses), post=post))

    # ---- kind ALL

    def emit_all(self, acc, spec, kernel):
        """Run the kernel on a whole field, ghost layers included; the
        components of the input are kernel parameters."""
        params = self.params(acc, kernel)
        ci, decl = self.whole_interval("xyzSizeWithGhostLayer",
                                       params.fields[0])
        self.fn(acc.name, *params.all)(
            decl, *self.slice_body(
                kernel, self.lattice_fields(acc, params), ci=ci,
                params=self.broadcast(acc, spec, params)))

    # ---- kind REDUCE

    def emit_reduce(self, acc, spec, kernel):
        """The values of all cells of the interior are summed on the host,
        in a fixed order: the reduction is deterministic."""
        T = self.T
        params = self.params(acc, kernel)
        result = acc.reduction
        n = result.container.size
        ci, decl = self.whole_interval("xyzSize", params.fields[0])
        name = result.name
        result_type = cxx_type(result.container, T)
        value = f"{name} * {result.scale}" if result.scale else name
        self.fn(acc.name, *params.all, returns=self.returns(acc))(
            decl, *self.slice_body(
                kernel, self.lattice_fields(acc, params), ci=ci,
                buffer=(self.cat.buffer(n), n_cells(n), ()),
                post=(local_declaration(catalogue.STD_VECTOR, T, "out", n),
                      self.copy_to_host("out.data()"),
                      f"{result_type} {name}({T}{{0}});",
                      f"for (auto i = uint_t{{0u}}; i < out.size(); "
                      f"i += {n}u) {{",
                      *(f"{name}[{i}u] += out[i + {i}u];" for i in range(n)),
                      "}", f"return {value};")))

    # ---- kind LIST

    def emit_list(self, acc, spec, kernel):
        """Set field values at the cells of a list (COO format) with
        indices relative to the origin of the halo of ``gl`` ghost layers."""
        T, D = self.T, acc.dims
        params = self.params(acc, kernel)
        field = params.fields[0]
        fname = acc.param_name(acc.fields[0].name, self.target)
        field_mut = self.field_param(acc.fields[0].role, False, "field_mut")
        gl = self.sfg.var("gl", UINT)
        idx_t = cpptype("thrust::device_vector< int >",
                        "<thrust/device_vector.h>")
        indices = AugExpr(idx_t(const=True, ref=True)).var("indices")
        values = AugExpr(_DeviceVector(
            T=self.dtype, const=True, ref=True)).var("values")
        n = AugExpr(PsCustomType("uint")).var("n_entries")
        n_decl = statement(
            f"uint n_entries = static_cast<uint>(indices.size() / {D}ul);",
            (n,), (indices,))
        ci = CellInterval().var("ci")
        ci_decl = statement(
            f"walberla::CellInterval ci = {fname}->xyzSize();\n"
            "ci.expand(walberla::cell_idx_c(gl));", (ci,), (field, gl))
        # the kernel writes the field, which is passed by const pointer
        field_decl = statement(
            f"auto *field_mut = const_cast<walberla::gpu::GPUField<{T}> *>"
            f"({fname});", (field_mut,), (field,))
        f_idx = self.cat.fields["indices"].field
        f_val = self.cat.fields["values"].field
        ptr_idx = LinearBuffer(PsPointerType(ps.create_type("int32")),
                               [n], D).var("indices_ptr")
        ptr_val = LinearBuffer(PsPointerType(self.dtype), [n],
                               D).var("values_ptr")
        decl = [
            statement("int const *indices_ptr = "
                      "thrust::raw_pointer_cast(indices.data());",
                      (ptr_idx,), (indices,)),
            statement(f"{T} const *values_ptr = "
                      f"thrust::raw_pointer_cast(values.data());",
                      (ptr_val,), (values,))]
        absolute = self.cat.fields[acc.fields[0].role].absolute
        self.fn(acc.name, field, indices, values, gl)(
            "assert(values.size() == indices.size());",
            "if (indices.empty()) {", "  return;", "}",
            n_decl, ci_decl, field_decl, *decl,
            self.sfg.map_field(absolute, field_mut.with_cell_interval(ci)),
            self.sfg.map_field(f_idx, ptr_idx),
            self.sfg.map_field(f_val, ptr_val),
            self.sfg.gpu_invoke(kernel))

    # ---- kind GATHER

    def emit_gather(self, acc, spec, kernel):
        """
        Function that runs ``kernel`` over the flat list ``pos`` of particle
        positions (``dim`` values per particle) and returns the ``dims``
        values per particle that it writes to the output buffer.
        """
        T, D = self.T, self.cat.dim
        params = self.params(acc, kernel)
        field = params.fields[0]
        fname = acc.param_name(acc.fields[0].name, self.target)
        gl = self.sfg.var("gl", UINT)
        pos = std.vector(self.dtype, const=True, ref=True).var("pos")
        n = AugExpr(PsCustomType("uint")).var("n_pos")
        ci = CellInterval().var("ci")
        n_out = "n_pos" if acc.dims == 1 else f"{acc.dims}u * n_pos"
        nodes = [
            "if (pos.empty()) {", "  return {};", "}",
            statement(
                f"uint n_pos = static_cast<uint>(pos.size() / {D}ul);",
                (n,), (pos,)),
            statement(f"walberla::CellInterval ci = {fname}->xyzSize();\n"
                      "ci.expand(walberla::cell_idx_c(gl));", (ci,),
                      (field, gl))]
        absolute = self.cat.fields[acc.fields[0].role].absolute
        maps = [self.sfg.map_field(absolute, field.with_cell_interval(ci))]
        uploads = (
            ("dev_pos", "pos.begin(), pos.end()", (pos,),
             self.cat.fields["pos"].field),
            ("dev_out", n_out, (), self.cat.fields[f"out{acc.dims}"].field))
        for vec_name, args, uses, f_buffer in uploads:
            ptr, decl = self.device_vector(
                vec_name, args, (*uses, n), [n], f_buffer.index_shape[0])
            nodes += decl
            maps.append(self.sfg.map_field(f_buffer, ptr))
        ghost_layers = ps.TypedSymbol("ghost_layers", "int64")
        set_gl = self.sfg.set_param(
            ghost_layers, AugExpr.format("int64_t({})", gl))
        post = (f"std::vector<{T}> out({n_out});",
                self.copy_to_host("out.data()", "dev_out"), "return out;")
        self.fn(acc.name, field, pos, *params.scalars.values(), gl,
                returns=self.returns(acc))(
            *nodes, *maps, set_gl, self.sfg.gpu_invoke(kernel), *post)

    # ---- kind SCATTER

    def emit_scatter(self, acc, spec, kernel):
        """
        Spread forces onto the lattice nodes. The nodes are shared by
        several particles, which needs ``atomicAdd``; pystencils kernels
        have no atomic updates of field entries: this kernel is written as
        a C++ string.
        """
        # the weights of ``calculate_weights`` are the same as those of the
        # pystencils kernels (``LbAccessors.interpolation_weights``): keep
        # both in sync
        T, D = self.T, acc.dims
        assert D == 3, "the scatter kernel is written for 3 dimensions"
        params = self.params(acc, kernel)
        pos = std.vector(self.dtype, const=True, ref=True).var("pos")
        forces = std.vector(self.dtype, const=True, ref=True).var("forces")
        gl = self.sfg.var("gl", UINT)
        self.sfg.code(_ADD_FORCE_KERNEL.format(T=T), impl=True)
        self.fn(acc.name, params.fields[0], pos, forces, gl)(
            "assert(forces.size() == pos.size());",
            "if (pos.empty()) {", "  return;", "}",
            f"thrust::device_vector<{T}> dev_pos(pos.begin(), pos.end());",
            f"thrust::device_vector<{T}> dev_for(forces.begin(), "
            "forces.end());",
            f"auto const *dev_pos_ptr = thrust::raw_pointer_cast("
            "dev_pos.data());",
            f"auto const *dev_for_ptr = thrust::raw_pointer_cast("
            "dev_for.data());",
            "auto const threads_per_block = uint(64u);",
            f"auto const n_pos = static_cast<uint>(pos.size() / {D}ul);",
            "auto const dim_grid = (n_pos + threads_per_block - 1u) / "
            "threads_per_block;",
            f"auto const accessor = walberla::gpu::FieldIndexing<{T}>::"
            "withGhostLayerXYZ(*field, gl).gpuAccess();",
            "/* clang-format off */",
            "kernel_add_force<<<dim_grid, threads_per_block, 0u, nullptr>>>("
            "accessor, dev_pos_ptr, dev_for_ptr, n_pos, gl);",
            "/* clang-format on */")


# the double braces are literal braces: the text is formatted with ``T``
_ADD_FORCE_KERNEL = """
/** @brief Calculate interpolation weights. */
static __forceinline__ __device__ void calculate_weights(
    {T} const *RESTRICT const pos, int *RESTRICT const corner,
    {T} *RESTRICT const weights, uint gl) {{
#pragma unroll
  for (int dim = 0; dim < 3; ++dim) {{
    auto const fractional_index = pos[dim] - {T}{{0.5}};
    auto const nmp = floor(fractional_index);
    auto const distance = fractional_index - nmp - {T}{{0.5}};
    corner[dim] = static_cast<int>(nmp) + static_cast<int>(gl);
    weights[dim * 2 + 0] = {T}{{0.5}} - distance;
    weights[dim * 2 + 1] = {T}{{0.5}} + distance;
  }}
}}

__global__ void kernel_add_force(
    walberla::gpu::FieldAccessor<{T}> force,
    {T} const *RESTRICT const pos, {T} const *RESTRICT const forces,
    uint n_pos, uint gl) {{
  uint const pos_index = blockDim.x * blockIdx.x + threadIdx.x;
  force.set({{0u, 0u, 0u}}, {{0u, 0u, 0u}});
  if (force.isValidPosition() and pos_index < n_pos) {{
    auto const array_offset = pos_index * uint(3u);
    int corner[3];
    {T} weights[3][2];
    calculate_weights(pos + array_offset, corner, &weights[0][0], gl);
#pragma unroll
    for (int i = 0; i < 2; i++) {{
      auto const cx = corner[0] + i;
      auto const wx = weights[0][i];
#pragma unroll
      for (int j = 0; j < 2; j++) {{
        auto const cy = corner[1] + j;
        auto const wxy = wx * weights[1][j];
#pragma unroll
        for (int k = 0; k < 2; k++) {{
          auto const cz = corner[2] + k;
          auto const weight = wxy * weights[2][k];
#pragma unroll
          for (uint cf = 0u; cf < 3u; ++cf) {{
            atomicAdd(&force.getNeighbor(cx, cy, cz, cf),
                      weight * forces[array_offset + cf]);
          }}
        }}
      }}
    }}
  }}
}}
"""
