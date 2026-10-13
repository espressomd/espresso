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
CUDA lattice field accessors of the waLBerla bridge, generated with
pystencils-sfg. Replaces the Jinja template ``FieldAccessors.tmpl.{cuh,cu}``.

The header ``FieldAccessors<Precision>PrecisionCUDA.cuh`` declares the
functions of the namespaces ``walberla::lbm::accessor::{Population, Vector,
Interpolation, Density, Velocity, Force, MomentumDensity, PressureTensor}``,
the source file ``FieldAccessors<Precision>PrecisionCUDA.cu`` defines them.

Device work:

* all per-cell and per-interval operations (the slice operations) are
  pystencils CUDA kernels, created from the equations of
  `accessor_equations` and launched from the generated host functions
  (``sfg.gpu_invoke``). The fields are `walberla::gpu::GPUField` views over
  the requested `CellInterval`; host data is exchanged through a linear
  device buffer that is mapped to a pystencils field
  (`LinearBuffer`, x-major, f-fastest);
* the reductions are per-cell pystencils kernels that fill a device buffer,
  which is summed on the host in the order of the cells (deterministic);
* the trilinear interpolation of particle positions and
  ``Vector::set_from_list`` are pystencils kernels over a 1D particle
  (or index) list that address the lattice fields with
  ``Field.absolute_access``;
* ``Interpolation::add_force`` scatters into the lattice field with
  ``atomicAdd``, which pystencils can't express: its ``__global__`` kernel
  is a C++ string.

Relies on pystencils-sfg 0.3 internals, see `accessor_sfg`.
"""

import itertools

import pystencils as ps
import sympy as sp
from pystencils.sympyextensions.typed_sympy import TypeCast
from pystencils.types import PsCustomType, PsPointerType, constify
from pystencilssfg.lang import AugExpr, cpptype
from pystencilssfg.lang.cpp import std
from sweepgen.api import Cell, CellInterval, GpuFieldPtr, Vector3

import accessor_equations
import accessor_sfg
from accessor_sfg import ARRAY, MATRIX3, UINT, LinearBuffer, statement

_DIM = 3
_DeviceVector = cpptype("thrust::device_vector< {T} >",
                        "<thrust/device_vector.h>")


def _instantiate(eqs, inputs, suffix=""):
    """
    Pystencils assignments of the equation set ``eqs``: the input symbols
    (by name) are replaced by ``inputs``, the assigned symbols get ``suffix``
    appended to their names.
    """
    subs = {sp.Symbol(name): value for name, value in inputs.items()}
    assignments = []
    for lhs, rhs in eqs.assignments:
        new = sp.Symbol(lhs.name + suffix)
        assignments.append(ps.Assignment(new, rhs.xreplace(subs)))
        subs[lhs] = new
    return assignments, subs


class _CudaFieldAccessors:
    def __init__(self, sfg, ctx, lb_method):
        if not ctx.is_gpu:
            raise ValueError("the CUDA accessors require a GPU target")
        self.sfg = sfg
        self.dtype = ctx.data_type
        self.T = self.dtype.c_string().strip()
        self.equations = accessor_equations.AccessorEquations(
            lb_method, self.dtype)
        self.Q = self.equations.Q
        self.D = self.equations.D
        if self.D != _DIM:
            raise ValueError("the field accessors require a "
                             "three-dimensional lattice")
        T, D, Q = self.dtype, self.D, self.Q
        # parameter types, as in the legacy header
        self.pdf_const = GpuFieldPtr(T, D, Q, const=True).var("pdf_field")
        self.pdf_mut = GpuFieldPtr(T, D, Q).var("pdf_field")
        self.vec_const = lambda name: GpuFieldPtr(
            T, D, D, const=True).var(name)
        self.vec_mut = lambda name: GpuFieldPtr(T, D, D).var(name)
        self.cell = Cell(const=True, ref=True).var("cell")
        self.ci = CellInterval(const=True, ref=True).var("ci")
        self.vector = lambda name: Vector3(T, const=True, ref=True).var(name)
        self.density_param = sfg.var("density", constify(T))
        self.values = std.vector(T, const=True, ref=True).var("values")
        self.pop = AugExpr(ARRAY(T=T, N=Q, const=True, ref=True)).var("pop")
        self.gl = sfg.var("gl", UINT)
        self.vector_return = std.vector(T).dtype
        self.vector3_return = Vector3(T).dtype
        self.matrix3_return = MATRIX3(T=T)
        self.array_return = ARRAY(T=T, N=Q)

        # pystencils fields; sizes and strides are set at runtime
        def field(name, rank, fsize, dtype=None, **kwargs):
            return ps.Field.create_generic(
                name, rank, dtype or self.dtype, index_dimensions=1,
                layout="fzyx" if rank == 3 else "zyxf",
                index_shape=(fsize,), **kwargs)
        self.f_pdf = field("pdfs", 3, Q)
        self.f_vel = field("velocity", 3, D)
        self.f_force = field("force", 3, D)
        self.f_vec = field("vec", 3, D)
        self.f_buf = {n: field(f"buf{n}", 3, n) for n in (1, D, D * D, Q)}
        # fields of the kernels over particle lists
        custom = dict(field_type=ps.FieldType.CUSTOM)
        self.f_pdf_abs = field("pdfs", 3, Q, **custom)
        self.f_vec_abs = field("vec", 3, D, **custom)
        self.f_pos = field("pos", 1, D)
        self.f_out = {n: field(f"out{n}", 1, n) for n in (1, D)}
        self.f_indices = field("indices", 1, D, dtype=ps.create_type("int32"))
        self.f_values = field("values", 1, D)

    # ---- helpers

    def sym(self, name):
        return ps.TypedSymbol(name, self.dtype)

    def kernel(self, name, assignments):
        config = ps.CreateKernelConfig(
            target=ps.Target.CUDA, default_dtype=self.dtype, ghost_layers=0)
        return self.sfg.kernels.create(assignments, name, config)

    def density_of(self, kernel):
        """The parameter ``density`` of a kernel, or the plain variable if
        the kernel doesn't take it."""
        return (accessor_sfg.kernel_param(kernel, "density")
                or self.density_param)

    def fn(self, name, *params, returns=None):
        """Start a (non-inline) function; the sequencer is called with the
        body."""
        return accessor_sfg.function(self.sfg, name, *params, returns=returns)

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

    def slice_body(self, kernel, fields, *, ci, buffer=None, scalars=None,
                   extra=(), post=()):
        """
        Body of a function that runs ``kernel`` on the interval ``ci``.

        * ``fields``: pairs of a pystencils field and a GPU field,
        * ``buffer``: ``(field, size_args, uses)`` of the linear device
          buffer that is mapped to the pystencils field ``field``;
          ``size_args`` is the argument list of the constructor of the
          device vector,
        * ``scalars``: ``(symbols, vector)``: the kernel parameters are the
          components of the array or vector,
        * ``extra``: nodes that precede the launch,
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
        if scalars is not None:
            # in the order of the components: the generated code must not
            # depend on the iteration order of a set
            symbols, vector = scalars
            nodes += [self.sfg.set_param(
                sym, AugExpr.format("{}[{}u]", vector, i))
                for i, sym in enumerate(symbols)]
        return [*nodes, *extra, self.sfg.gpu_invoke(kernel), *post]

    def cell_interval(self):
        """Interval with the single cell ``cell``."""
        ci = CellInterval(const=True).var("ci")
        return ci, statement("walberla::CellInterval ci(cell, cell);",
                             (ci,), (self.cell,))

    def whole_interval(self, size, field_ptr):
        """Interval over a whole field: ``size`` is the member function of
        the field (e.g. ``xyzSize``)."""
        ci = CellInterval().var("ci")
        return ci, statement(
            f"walberla::CellInterval ci = {field_ptr}->{size}();",
            (ci,), (field_ptr,))

    def copy_to_host(self, target, vector="dev_data"):
        return f"thrust::copy({vector}.begin(), {vector}.end(), {target});"

    # ---- generate

    def generate(self):
        sfg = self.sfg
        for header in ("core/DataTypes.h", "core/cell/Cell.h",
                       "core/cell/CellInterval.h", "core/math/Matrix3.h",
                       "core/math/Vector3.h", "gpu/GPUField.h",
                       "<thrust/device_vector.h>", "<array>", "<vector>"):
            sfg.include(header)
        for header in ("gpu/FieldAccessor.h", "gpu/FieldIndexing.h",
                       "<cassert>", "<cstdint>", "<thrust/copy.h>"):
            sfg.include(header, private=True)
        sfg.prelude("\nLattice field accessors.\n"
                    "Adapted from the waLBerla source file "
                    "python/lbmpy_walberla/templates/LatticeModel.tmpl.h "
                    "(commit a16141524c58ab88386e2a0f8fdd7c63c5edd704).")
        sfg.namespace("walberla::lbm::accessor")
        for name in ("population", "vector", "interpolation", "density",
                     "velocity", "force", "momentum_density",
                     "pressure_tensor"):
            getattr(self, f"gen_{name}")()

    # ---- equations as pystencils assignments

    def reads(self, eqs, pdf=None, force=None, **others):
        """Replacement of the input symbols of an equation set: the
        populations ``f_i`` by ``pdf(i)``, the forces by ``force(i)``."""
        E = self.equations
        inputs = dict(others)
        if pdf is not None:
            inputs.update({E.pdf_names[i]: pdf(i) for i in E.used_pdfs(eqs)})
        if force is not None:
            inputs.update({E.force_names[i]: force(i)
                           for i in E.used_forces(eqs)})
        return inputs

    def equilibrium(self, pdf, rho, u):
        """Assignments that set all populations of ``pdf`` to the
        equilibrium of the density ``rho`` and the velocity ``u``."""
        E = self.equations
        delta_rho = sp.Symbol("delta_rho_eq")
        subs = {sp.Symbol("rho"): rho, sp.Symbol("delta_rho"): delta_rho,
                **{sp.Symbol(n): v for n, v in zip(E.u_names, u)}}
        return [ps.Assignment(delta_rho, rho - 1)] + [
            ps.Assignment(pdf(i), expr.xreplace(subs))
            for i, expr in enumerate(E.equilibrium)]

    def momenta_over_density(self, rho_inv):
        """Velocity components ``md_j * rho_inv``."""
        return [sp.Symbol(f"md_{j}") * rho_inv for j in range(self.D)]

    # ---- namespace Population

    def gen_population(self):
        T, Q = self.T, self.Q
        E = self.equations
        pdf, pdf_mut, cell, ci = (self.pdf_const, self.pdf_mut, self.cell,
                                  self.ci)
        pop, values = self.pop, self.values
        vel = self.vec_mut("velocity_field")
        frc = self.vec_const("force_field")
        f_pdf, f_vel, f_force, buf = (self.f_pdf, self.f_vel, self.f_force,
                                      self.f_buf[Q])
        scalars = [self.sym(f"pop_{i}") for i in range(Q)]
        md = E.momentum_density
        rho, rho_inv = sp.Symbol("rho"), sp.Symbol("rho_inv")

        with self.sfg.namespace("Population"):
            k_get = self.kernel("get", [
                ps.Assignment(buf(i), f_pdf(i)) for i in range(Q)])
            k_set = self.kernel("set", [
                ps.Assignment(f_pdf(i), buf(i)) for i in range(Q)])
            k_broadcast = self.kernel("broadcast", [
                ps.Assignment(f_pdf(i), s) for i, s in enumerate(scalars)])
            # the populations are read from the buffer
            temporaries, _ = _instantiate(
                md, self.reads(md, pdf=buf, force=f_force.center))
            k_set_vel = self.kernel("set_vel", [
                *temporaries, ps.Assignment(rho_inv, 1 / rho),
                *(ps.Assignment(f_pdf(i), buf(i)) for i in range(Q)),
                *(ps.Assignment(f_vel(j), m) for j, m in enumerate(
                    self.momenta_over_density(rho_inv)))])

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("get", pdf, cell, returns=self.array_return)(
                ci_cell_decl,
                *self.slice_body(
                    k_get, [(f_pdf, pdf)], ci=ci_cell,
                    buffer=(buf, f"{Q}u", ()),
                    post=(f"std::array<{T}, {Q}u> pop;",
                          self.copy_to_host("pop.data()"), "return pop;")))

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("set", pdf_mut, pop, cell)(
                ci_cell_decl,
                *self.slice_body(
                    k_set, [(f_pdf, pdf_mut)], ci=ci_cell,
                    buffer=(buf, "pop.begin(), pop.end()", (pop,))))

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("set", pdf_mut, vel, frc, pop, cell)(
                ci_cell_decl,
                *self.slice_body(
                    k_set_vel, [(f_pdf, pdf_mut), (f_vel, vel),
                                (f_force, frc)], ci=ci_cell,
                    buffer=(buf, "pop.begin(), pop.end()", (pop,))))

            ci_all, ci_all_decl = self.whole_interval(
                "xyzSizeWithGhostLayer", pdf_mut)
            self.fn("initialize", pdf_mut, pop)(
                ci_all_decl,
                *self.slice_body(
                    k_broadcast, [(f_pdf, pdf_mut)], ci=ci_all,
                    scalars=(scalars, pop)))

            self.fn("get", pdf, ci, returns=self.vector_return)(
                *self.slice_body(
                    k_get, [(f_pdf, pdf)], ci=ci,
                    buffer=(buf, f"ci.numCells() * {Q}u", ()),
                    post=(f"std::vector<{T}> out(ci.numCells() * {Q}u);",
                          self.copy_to_host("out.data()"), "return out;")))

            self.fn("set", pdf_mut, values, ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * {Q}u);",
                *self.slice_body(
                    k_set, [(f_pdf, pdf_mut)], ci=ci,
                    buffer=(buf, "values.begin(), values.end()", (values,))))

            self.fn("set", pdf_mut, vel, frc, values, ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * {Q}u);",
                *self.slice_body(
                    k_set_vel, [(f_pdf, pdf_mut), (f_vel, vel),
                                (f_force, frc)], ci=ci,
                    buffer=(buf, "values.begin(), values.end()", (values,))))

    # ---- namespace Vector

    def gen_vector(self):
        T, D = self.T, self.D
        vec_c, vec_m = self.vec_const("field"), self.vec_mut("field")
        vec_c_ci = self.vec_const("vec_field")
        vec = self.vector("vec")
        cell, ci, values = self.cell, self.ci, self.values
        f_vec, buf = self.f_vec, self.f_buf[D]
        scalars = [self.sym(f"vec_{j}") for j in range(D)]
        comps = range(D)

        with self.sfg.namespace("Vector"):
            k_get = self.kernel("get", [
                ps.Assignment(buf(j), f_vec(j)) for j in comps])
            k_set = self.kernel("set", [
                ps.Assignment(f_vec(j), buf(j)) for j in comps])
            k_add = self.kernel("add", [
                ps.Assignment(f_vec(j), f_vec(j) + buf(j)) for j in comps])
            k_broadcast = self.kernel("broadcast", [
                ps.Assignment(f_vec(j), s) for j, s in enumerate(scalars)])
            k_broadcast_add = self.kernel("broadcast_add", [
                ps.Assignment(f_vec(j), f_vec(j) + s)
                for j, s in enumerate(scalars)])

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("get", vec_c, cell, returns=self.vector3_return)(
                ci_cell_decl,
                *self.slice_body(
                    k_get, [(f_vec, vec_c)], ci=ci_cell,
                    buffer=(buf, f"{D}u", ()),
                    post=(f"walberla::Vector{D}<{T}> vec;",
                          self.copy_to_host("vec.data()"), "return vec;")))

            for name, kernel in (("set", k_set), ("add", k_add)):
                ci_cell, ci_cell_decl = self.cell_interval()
                self.fn(name, vec_m, vec, cell)(
                    ci_cell_decl,
                    *self.slice_body(
                        kernel, [(f_vec, vec_m)], ci=ci_cell,
                        buffer=(buf, f"vec.data(), vec.data() + {D}u",
                                (vec,))))

            for name, kernel in (("initialize", k_broadcast),
                                 ("add_to_all", k_broadcast_add)):
                ci_all, ci_all_decl = self.whole_interval(
                    "xyzSizeWithGhostLayer", vec_m)
                self.fn(name, vec_m, vec)(
                    ci_all_decl,
                    *self.slice_body(
                        kernel, [(f_vec, vec_m)], ci=ci_all,
                        scalars=(scalars, vec)))

            self.fn("get", vec_c_ci, ci, returns=self.vector_return)(
                *self.slice_body(
                    k_get, [(f_vec, vec_c_ci)], ci=ci,
                    buffer=(buf, f"ci.numCells() * {D}u", ()),
                    post=(f"std::vector<{T}> out(ci.numCells() * {D}u);",
                          self.copy_to_host("out.data()"), "return out;")))

            vec_m_ci_set = self.vec_mut("vec_field")
            self.fn("set", vec_m_ci_set, values, ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * {D}u);",
                *self.slice_body(
                    k_set, [(f_vec, vec_m_ci_set)], ci=ci,
                    buffer=(buf, "values.begin(), values.end()", (values,))))

            self.gen_set_from_list()

    def gen_set_from_list(self):
        """Set field values at the cells of a list (COO format) with
        indices relative to the origin of the halo of ``gl`` ghost layers."""
        T, D = self.T, self.D
        field = self.vec_const("field")
        field_mut = self.vec_mut("field_mut")
        idx_t = cpptype("thrust::device_vector< int >",
                        "<thrust/device_vector.h>")
        indices = AugExpr(idx_t(const=True, ref=True)).var("indices")
        values = AugExpr(_DeviceVector(
            T=self.dtype, const=True, ref=True)).var("values")
        f_idx, f_val, f_vec = self.f_indices, self.f_values, self.f_vec_abs
        gl = self.gl
        idx = tuple(f_idx(j) for j in range(D))
        kernel = self.kernel("set_from_list", [
            ps.Assignment(f_vec.absolute_access(idx, (j,)), f_val(j))
            for j in range(D)])
        n = AugExpr(PsCustomType("uint")).var("n_entries")
        n_decl = statement(
            f"uint n_entries = static_cast<uint>(indices.size() / {D}ul);",
            (n,), (indices,))
        ci = CellInterval().var("ci")
        ci_decl = statement(
            "walberla::CellInterval ci = field->xyzSize();\n"
            "ci.expand(walberla::cell_idx_c(gl));", (ci,), (field, gl))
        # the kernel writes the field, which is passed by const pointer
        field_decl = statement(
            f"auto *field_mut = const_cast<walberla::gpu::GPUField<{T}> *>"
            f"(field);", (field_mut,), (field,))
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
        self.fn("set_from_list", field, indices, values, gl)(
            f"assert(values.size() == indices.size());",
            "if (indices.empty()) {", "  return;", "}",
            n_decl, ci_decl, field_decl, *decl,
            self.sfg.map_field(f_vec, field_mut.with_cell_interval(ci)),
            self.sfg.map_field(f_idx, ptr_idx),
            self.sfg.map_field(f_val, ptr_val),
            self.sfg.gpu_invoke(kernel))

    # ---- namespace Interpolation

    def interpolation_weights(self):
        """
        Assignments of the trilinear interpolation weights of the particle
        position ``pos``, which is relative to the origin of the local
        halo-free domain, with cells centered at half-integer positions.

        The same weights are computed by the C++ function
        ``calculate_weights`` in `gen_add_force`: keep both in sync (they
        are covered by the interpolation tests and the GPU equivalence
        harness).
        Returns the assignments, the corner cell (in ``pdf`` index space,
        which starts ``gl`` cells below the origin), and a function that
        returns the weight of the corner ``(i, j, k)``.
        """
        half = sp.Rational(1, 2)
        int64 = ps.create_type("int64")
        gl = ps.TypedSymbol("ghost_layers", int64)
        assignments, corner, weights = [], [], []
        for d in range(self.D):
            fractional = sp.Symbol(f"fractional_index_{d}")
            nmp = sp.Symbol(f"nmp_{d}")
            distance = sp.Symbol(f"distance_{d}")
            c = ps.TypedSymbol(f"corner_{d}", int64)
            assignments += [
                ps.Assignment(fractional, self.f_pos(d) - half),
                ps.Assignment(nmp, sp.floor(fractional)),
                ps.Assignment(distance, fractional - nmp - half),
                ps.Assignment(c, TypeCast(nmp, int64) + gl)]
            corner.append(c)
            weights.append((half - distance, half + distance))

        def weight(i, j, k):
            return weights[0][i] * weights[1][j] * weights[2][k]
        return assignments, corner, weight

    def interpolate(self, name, fsize, out, cell_values):
        """
        Pystencils kernel that interpolates ``fsize`` values per particle,
        returned by ``cell_values(corner, n)`` for the corner cell with
        absolute indices ``corner`` and the unique name suffix ``n``.
        The corners are accumulated in the order of the nested loops
        ``i, j, k``.
        """
        assignments, corner, weight = self.interpolation_weights()
        accumulated = [None] * fsize
        for n, (i, j, k) in enumerate(itertools.product((0, 1), repeat=3)):
            cell = (corner[0] + i, corner[1] + j, corner[2] + k)
            w = sp.Symbol(f"weight_{n}")
            assignments.append(ps.Assignment(w, weight(i, j, k)))
            prologue, values = cell_values(cell, n)
            assignments += prologue
            for cf in range(fsize):
                total = w * values[cf]
                if accumulated[cf] is not None:
                    total = accumulated[cf] + total
                accumulated[cf] = sp.Symbol(f"accumulated_{cf}_{n}")
                assignments.append(ps.Assignment(accumulated[cf], total))
        assignments += [ps.Assignment(out(cf), accumulated[cf])
                        for cf in range(fsize)]
        return self.kernel(name, assignments)

    def particle_function(self, name, field, f_field, kernel, uploads,
                          post=(), returns=None, params=(), pos=None):
        """
        Function that runs ``kernel`` over the flat list ``pos`` of particle
        positions (``D`` values per particle). ``uploads`` are tuples
        ``(device vector name, constructor arguments, variables used,
        pystencils field)`` of the device buffers, which have ``D`` or 1
        values per particle.
        """
        D = self.D
        gl = self.gl
        pos = pos or std.vector(self.dtype, const=True, ref=True).var("pos")
        n = AugExpr(PsCustomType("uint")).var("n_pos")
        ci = CellInterval().var("ci")
        nodes = [
            "if (pos.empty()) {",
            "  return {};" if returns is not None else "  return;", "}",
            statement(
                f"uint n_pos = static_cast<uint>(pos.size() / {D}ul);",
                (n,), (pos,)),
            statement("walberla::CellInterval ci = field->xyzSize();\n"
                      "ci.expand(walberla::cell_idx_c(gl));", (ci,),
                      (field, gl))]
        maps = [self.sfg.map_field(f_field, field.with_cell_interval(ci))]
        for vec_name, args, uses, f_buffer in uploads:
            ptr, decl = self.device_vector(
                vec_name, args, (*uses, n), [n], f_buffer.index_shape[0])
            nodes += decl
            maps.append(self.sfg.map_field(f_buffer, ptr))
        ghost_layers = ps.TypedSymbol("ghost_layers", "int64")
        set_gl = self.sfg.set_param(
            ghost_layers, AugExpr.format("int64_t({})", gl))
        return self.fn(name, field, *params, gl, returns=returns)(
            *nodes, *maps, set_gl, self.sfg.gpu_invoke(kernel), *post)

    def gen_interpolation(self):
        T, D = self.T, self.D
        E = self.equations
        pos = std.vector(self.dtype, const=True, ref=True).var("pos")
        forces = std.vector(self.dtype, const=True, ref=True).var("forces")
        density_sym = self.sym("density")
        scaled = E.density_scaled

        def cell_density(cell, n):
            reads = self.reads(
                scaled, pdf=lambda q: self.f_pdf_abs.absolute_access(
                    cell, (q,)), density=density_sym)
            temporaries, subs = _instantiate(scaled, reads, f"_{n}")
            return temporaries, [subs[sp.Symbol("rho")]]

        def cell_velocity(cell, n):
            return [], [self.f_vec_abs.absolute_access(cell, (cf,))
                        for cf in range(D)]

        copy_out = "thrust::copy(dev_out.begin(), dev_out.end(), out.data());"
        with self.sfg.namespace("Interpolation"):
            k_rho = self.interpolate(
                "get_rho", 1, self.f_out[1], cell_density)
            k_vel = self.interpolate(
                "get_vel", D, self.f_out[D], cell_velocity)

            field = GpuFieldPtr(self.dtype, D, self.Q, const=True).var(
                "field")
            self.particle_function(
                "get_rho", field, self.f_pdf_abs, k_rho,
                [("dev_pos", "pos.begin(), pos.end()", (pos,), self.f_pos),
                 ("dev_out", "n_pos", (), self.f_out[1])],
                (f"std::vector<{T}> out(n_pos);", copy_out, "return out;"),
                returns=self.vector_return, params=(
                    pos, self.density_of(k_rho)),
                pos=pos)

            field = GpuFieldPtr(self.dtype, D, D, const=True).var("field")
            self.particle_function(
                "get_vel", field, self.f_vec_abs, k_vel,
                [("dev_pos", "pos.begin(), pos.end()", (pos,), self.f_pos),
                 ("dev_out", f"{D}u * n_pos", (), self.f_out[D])],
                (f"std::vector<{T}> out({D}u * n_pos);", copy_out,
                 "return out;"),
                returns=self.vector_return, params=(pos,), pos=pos)
            self.gen_add_force(pos, forces)

    def gen_add_force(self, pos, forces):
        """
        Spread forces onto the lattice nodes. The nodes are shared by
        several particles, which needs ``atomicAdd``; pystencils kernels
        have no atomic updates of field entries: this kernel is written as
        a C++ string.
        """
        # the weights of ``calculate_weights`` are the same as those of the
        # pystencils kernels (`interpolation_weights`): keep both in sync
        T, D = self.T, self.D
        field = GpuFieldPtr(self.dtype, D, D, const=True).var("field")
        self.sfg.code(f"""
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
""", impl=True)
        self.fn("add_force", field, pos, forces, self.gl)(
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

    # ---- namespace Density

    def gen_density(self):
        T, D = self.T, self.D
        E = self.equations
        pdf, pdf_mut, cell, ci, values = (self.pdf_const, self.pdf_mut,
                                          self.cell, self.ci, self.values)
        rho_param = self.sfg.var("rho", constify(self.dtype))
        f_pdf, buf = self.f_pdf, self.f_buf[1]
        density_sym = self.sym("density")
        rho, rho_inv = sp.Symbol("rho"), sp.Symbol("rho_inv")
        scaled = E.density_scaled
        unshifted = E.unshifted_momentum_density

        with self.sfg.namespace("Density"):
            temporaries, _ = _instantiate(
                scaled, self.reads(scaled, pdf=f_pdf, density=density_sym))
            k_get = self.kernel("get", [
                *temporaries, ps.Assignment(buf(0), rho)])
            # the density is changed at constant velocity (without the half
            # force shift)
            temporaries, _ = _instantiate(
                unshifted, self.reads(unshifted, pdf=f_pdf))
            rho_new = sp.Symbol("rho_new")
            u_old = [sp.Symbol(f"u_old_{j}") for j in range(D)]
            k_set = self.kernel("set", [
                *temporaries, ps.Assignment(rho_inv, 1 / rho),
                *(ps.Assignment(u, m) for u, m in zip(
                    u_old, self.momenta_over_density(rho_inv))),
                ps.Assignment(rho_new, buf(0) / density_sym),
                *self.equilibrium(f_pdf, rho_new, u_old)])

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("get", pdf, self.density_of(k_get), cell,
                    returns=self.dtype)(
                ci_cell_decl,
                *self.slice_body(
                    k_get, [(f_pdf, pdf)], ci=ci_cell,
                    buffer=(buf, "1u", ()),
                    post=(f"{T} rho{{}};",
                          self.copy_to_host("&rho"), "return rho;")))

            self.fn("get", pdf, self.density_of(k_get), ci,
                    returns=self.vector_return)(
                *self.slice_body(
                    k_get, [(f_pdf, pdf)], ci=ci,
                    buffer=(buf, "ci.numCells()", ()),
                    post=(f"std::vector<{T}> out(ci.numCells());",
                          self.copy_to_host("out.data()"), "return out;")))

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("set", pdf_mut, rho_param, self.density_of(k_set), cell)(
                ci_cell_decl,
                *self.slice_body(
                    k_set, [(f_pdf, pdf_mut)], ci=ci_cell,
                    buffer=(buf, "1u, rho", (rho_param,))))

            self.fn("set", pdf_mut, values, self.density_of(k_set), ci)(
                "assert(uint_c(values.size()) == ci.numCells());",
                *self.slice_body(
                    k_set, [(f_pdf, pdf_mut)], ci=ci,
                    buffer=(buf, "values.begin(), values.end()", (values,))))

    # ---- namespace Velocity

    def gen_velocity(self):
        T, D = self.T, self.D
        E = self.equations
        pdf, pdf_mut, cell, ci, values = (self.pdf_const, self.pdf_mut,
                                          self.cell, self.ci, self.values)
        vel = self.vec_mut("velocity_field")
        frc = self.vec_const("force_field")
        u = self.vector("u")
        f_pdf, f_vel, f_force, buf = (self.f_pdf, self.f_vel, self.f_force,
                                      self.f_buf[D])
        rho, rho_inv = sp.Symbol("rho"), sp.Symbol("rho_inv")
        md = E.momentum_density
        dens = E.density
        setter = E.velocity_setter

        with self.sfg.namespace("Velocity"):
            temporaries, _ = _instantiate(
                md, self.reads(md, pdf=f_pdf, force=f_force.center))
            k_get = self.kernel("get", [
                *temporaries, ps.Assignment(rho_inv, 1 / rho),
                *(ps.Assignment(buf(j), m) for j, m in enumerate(
                    self.momenta_over_density(rho_inv)))])
            # the velocity is shifted by half the force, the stored velocity
            # is not
            temporaries, _ = _instantiate(
                dens, self.reads(dens, pdf=f_pdf))
            shift, _ = _instantiate(
                setter, self.reads(setter, force=f_force.center,
                                   **{n: buf(j) for j, n in
                                      enumerate(E.u_in_names)}))
            k_set = self.kernel("set", [
                *temporaries, *shift,
                *(ps.Assignment(f_vel(j), buf(j)) for j in range(D)),
                *self.equilibrium(
                    f_pdf, rho, [sp.Symbol(n) for n in E.u_names])])

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("get", pdf, frc, cell, returns=self.vector3_return)(
                ci_cell_decl,
                *self.slice_body(
                    k_get, [(f_pdf, pdf), (f_force, frc)], ci=ci_cell,
                    buffer=(buf, f"{D}u", ()),
                    post=(f"walberla::Vector{D}<{T}> vec;",
                          self.copy_to_host("vec.data()"), "return vec;")))

            # the buffer holds the vectors of all cells of the interval
            self.fn("get", pdf, frc, ci, returns=self.vector_return)(
                *self.slice_body(
                    k_get, [(f_pdf, pdf), (f_force, frc)], ci=ci,
                    buffer=(buf, f"ci.numCells() * {D}u", ()),
                    post=(f"std::vector<{T}> out(ci.numCells() * {D}u);",
                          self.copy_to_host("out.data()"), "return out;")))

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("set", pdf_mut, vel, frc, u, cell)(
                ci_cell_decl,
                *self.slice_body(
                    k_set, [(f_pdf, pdf_mut), (f_vel, vel), (f_force, frc)],
                    ci=ci_cell,
                    buffer=(buf, f"u.data(), u.data() + {D}u", (u,))))

            self.fn("set", pdf_mut, vel, frc, values, ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * {D}u);",
                *self.slice_body(
                    k_set, [(f_pdf, pdf_mut), (f_vel, vel), (f_force, frc)],
                    ci=ci,
                    buffer=(buf, "values.begin(), values.end()", (values,))))

    # ---- namespace Force

    def gen_force(self):
        D = self.D
        E = self.equations
        pdf = self.pdf_const
        vel = self.vec_mut("velocity_field")
        frc = self.vec_mut("force_field")
        u = self.vector("u")
        cell, ci, values = self.cell, self.ci, self.values
        f_pdf, f_vel, f_force, buf = (self.f_pdf, self.f_vel, self.f_force,
                                      self.f_buf[D])
        density_sym = self.sym("density")
        rho, rho_inv = sp.Symbol("rho"), sp.Symbol("rho_inv")
        eqs = E.momentum_density_force_setter

        with self.sfg.namespace("Force"):
            temporaries, _ = _instantiate(
                eqs, self.reads(
                    eqs, pdf=f_pdf, density=density_sym,
                    **{n: buf(j) for j, n in enumerate(E.force_names)}))
            density_inv = sp.Symbol("density_inv")
            k_set = self.kernel("set", [
                *temporaries, ps.Assignment(rho_inv, 1 / rho),
                ps.Assignment(density_inv, 1 / density_sym),
                *(ps.Assignment(f_force(j), buf(j) * density_inv)
                  for j in range(D)),
                *(ps.Assignment(f_vel(j), m) for j, m in enumerate(
                    self.momenta_over_density(rho_inv)))])

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("set", pdf, vel, frc, u, self.density_of(k_set), cell)(
                ci_cell_decl,
                *self.slice_body(
                    k_set, [(f_pdf, pdf), (f_vel, vel), (f_force, frc)],
                    ci=ci_cell,
                    buffer=(buf, f"u.data(), u.data() + {D}u", (u,))))

            self.fn("set", pdf, vel, frc, values, self.density_of(k_set), ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * {D}u);",
                *self.slice_body(
                    k_set, [(f_pdf, pdf), (f_vel, vel), (f_force, frc)],
                    ci=ci,
                    buffer=(buf, "values.begin(), values.end()", (values,))))

    # ---- reductions

    def host_sum(self, name, n, value_type):
        """
        Statements that sum, on the host and in the order of the cells, the
        ``n`` values per cell of the device buffer copy ``out`` into
        ``name``.
        """
        components = (f"{name}[{i}u] += out[i + {i}u];" for i in range(n))
        return (f"{value_type} {name}({self.T}{{0}});",
                f"for (auto i = uint_t{{0u}}; i < out.size(); i += {n}u) {{",
                *components, "}")

    # ---- namespace MomentumDensity

    def gen_momentum_density(self):
        """The momentum density of all cells is summed on the host, in a
        fixed order: the reduction is deterministic."""
        T, D = self.T, self.D
        E = self.equations
        pdf = self.pdf_const
        frc = self.vec_const("force_field")
        density_sym = self.sym("density")
        f_pdf, f_force, buf = self.f_pdf, self.f_force, self.f_buf[D]
        # only the shifted momentum density is needed: this skips the
        # density and the populations it depends on
        eqs = E.momentum_density.pruned([f"md_{j}" for j in range(D)])
        with self.sfg.namespace("MomentumDensity"):
            temporaries, _ = _instantiate(
                eqs, self.reads(eqs, pdf=f_pdf, force=f_force.center))
            kernel = self.kernel("get", [
                *temporaries,
                *(ps.Assignment(buf(j), sp.Symbol(f"md_{j}") * density_sym)
                  for j in range(D))])
            ci, ci_decl = self.whole_interval("xyzSize", pdf)
            self.fn("reduce", pdf, frc, self.density_of(kernel),
                    returns=self.vector3_return)(
                ci_decl,
                *self.slice_body(
                    kernel, [(f_pdf, pdf), (f_force, frc)], ci=ci,
                    buffer=(buf, f"ci.numCells() * {D}u", ()),
                    post=(f"std::vector<{T}> out(ci.numCells() * {D}u);",
                          self.copy_to_host("out.data()"),
                          *self.host_sum("momentum_density", D,
                                         f"walberla::Vector{D}<{T}>"),
                          "return momentum_density;")))

    # ---- namespace PressureTensor

    def gen_pressure_tensor(self):
        T, D = self.T, self.D
        E = self.equations
        pdf, cell, ci = self.pdf_const, self.cell, self.ci
        density = self.density_param
        density_sym = self.sym("density")
        f_pdf, buf = self.f_pdf, self.f_buf[D * D]
        eqs = E.second_moment
        offsets = E.inverse_neighbor_offsets
        comps = range(D * D)
        # the populations are pulled from the neighbor cells
        used = {E.pdf_names[i]: f_pdf[offsets[i]](i)
                for i in E.used_pdfs(eqs)}
        with self.sfg.namespace("PressureTensor"):
            temporaries, _ = _instantiate(eqs, used)
            p = [sp.Symbol(f"p_{i}") for i in comps]
            k_get = self.kernel("get", [
                *temporaries,
                *(ps.Assignment(buf(i), p[i] * density_sym) for i in comps)])
            k_reduce = self.kernel("reduce", [
                *temporaries,
                *(ps.Assignment(buf(i), p[i]) for i in comps)])

            ci_cell, ci_cell_decl = self.cell_interval()
            self.fn("get", pdf, self.density_of(k_get), cell,
                    returns=self.matrix3_return)(
                ci_cell_decl,
                *self.slice_body(
                    k_get, [(f_pdf, pdf)], ci=ci_cell,
                    buffer=(buf, f"{D * D}u", ()),
                    post=(f"walberla::Matrix{D}<{T}> out;",
                          self.copy_to_host("out.data()"), "return out;")))

            self.fn("get", pdf, self.density_of(k_get), ci,
                    returns=self.vector_return)(
                *self.slice_body(
                    k_get, [(f_pdf, pdf)], ci=ci,
                    buffer=(buf, f"ci.numCells() * {D * D}u", ()),
                    post=(f"std::vector<{T}> out(ci.numCells() * "
                          f"{D * D}u);", self.copy_to_host("out.data()"),
                          "return out;")))

            ci_all, ci_all_decl = self.whole_interval("xyzSize", pdf)
            self.fn("reduce", pdf, density, returns=self.matrix3_return)(
                ci_all_decl,
                *self.slice_body(
                    k_reduce, [(f_pdf, pdf)], ci=ci_all,
                    buffer=(buf, f"ci.numCells() * {D * D}u", ()),
                    post=(f"std::vector<{T}> out(ci.numCells() * {D * D}u);",
                          self.copy_to_host("out.data()"),
                          *self.host_sum("pressure_tensor", D * D,
                                         f"walberla::Matrix{D}<{T}>"),
                          "return pressure_tensor * density;")))


def generate_field_accessors_cuda(sfg, ctx, lb_method):
    """
    Add the complete set of CUDA field accessors of ``lb_method`` to ``sfg``
    (header ``.cuh`` and implementation ``.cu``, see ``KernelFileGenerator``).
    """
    _CudaFieldAccessors(sfg, ctx, lb_method).generate()
