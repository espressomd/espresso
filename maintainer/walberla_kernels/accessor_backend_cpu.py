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
CPU backend of the field accessors: generates the wrappers of a
`accessor_catalogue.Catalogue` as ``inline`` functions of a header-only
file, with pystencils-sfg. The header is self-contained: the kernels
(pystencils ``GenericCPU`` kernels over the strides of the fields,
``inline`` functions in the namespace ``cpu_kernels`` of each namespace of
accessors) are defined in it, so the waLBerla bridge needs no source file
and no CMake change for the accessors.

* the accessors of the kinds ``CELL``, ``INTERVAL`` and ``ALL`` run a kernel
  on the `walberla::field::GhostLayerField` views over the requested
  `CellInterval`; host data is exchanged through a linear host buffer
  (the container that is returned or passed, `LinearBuffer`) that is
  mapped to a pystencils field: x-major, f-fastest, as the legacy
  accessors;
* a reduction (``REDUCE``) is a kernel with ``AddReductionAssignment``
  that sums sequentially in the order of the cells (z, y, x), which is
  deterministic and equal to the order of the legacy loops;
* ``POINT`` evaluates an expression per stencil direction in a ``switch``.

The kernels of the CPU target need no `sfg_driver.patch_kernel_source`:
it patches ``.cpp`` and ``.cu`` files, and the accessor kernels don't use
the functions it replaces (the RNG, ``sqrt``, ``fmin``, ``fmax``);
``tests/test_field_accessors.py`` guards this.

Relies on pystencils-sfg 0.3 internals, see `accessor_sfg`.
"""

import pystencils as ps
from pystencils.types import PsPointerType, constify
from pystencilssfg.lang import AugExpr
from sweepgen.api import Direction, GhostLayerFieldPtr, Vector3

import accessor_catalogue as catalogue
import accessor_equations
from accessor_sfg import (
    AccessorBackend, LinearBuffer, construction, data_pointer,
    local_declaration, n_cells, statement)


class CpuBackend(AccessorBackend):
    target = "cpu"
    gpu = False
    inline = True
    kernel_namespace = "cpu_kernels"
    ps_target = ps.Target.GenericCPU
    headers = ("core/DataTypes.h", "core/Abort.h", "core/cell/Cell.h",
               "core/cell/CellInterval.h", "core/math/Matrix3.h",
               "core/math/Vector3.h", "field/GhostLayerField.h", "<array>",
               "<cassert>", "<vector>")

    def field_param(self, role, const, name):
        spec = self.cat.fields[role]
        return GhostLayerFieldPtr(self.dtype, self.cat.dim, spec.fsize,
                                  const=const).var(name)

    def reduction_kernel(self, spec):
        """The values of the cells are added to the scalars ``sum_i``,
        which the kernel takes by pointer."""
        assignments, values = spec.build()
        return self.kernel(spec.name, [
            *assignments,
            *(ps.AddReductionAssignment(self.sum_symbol(i), v)
              for i, v in enumerate(values))])

    # ---- helpers

    def sum_symbol(self, i):
        return ps.TypedSymbol(f"sum_{i}", self.dtype)

    def map_fields(self, acc, params, ci):
        return [self.sfg.map_field(f, g.with_cell_interval(ci))
                for f, g in self.lattice_fields(acc, params)]

    # ---- kinds CELL and INTERVAL

    def emit_cell(self, acc, _spec, kernel):
        self.emit_slice(acc, kernel)

    def emit_interval(self, acc, _spec, kernel):
        self.emit_slice(acc, kernel)

    def emit_slice(self, acc, kernel):
        """Run the kernel on a cell or an interval; the values are
        exchanged through a host buffer."""
        # The ``statement(code, defines, uses)`` tuples tell pystencils-sfg
        # which variables a line of code defines and uses, which orders the
        # lines (e.g. ``ci`` before the buffer that is sized with it). ``ci``
        # is a fixed local name that is shared with ``n_cells()`` and
        # ``cell_interval()`` in accessor_sfg.py. The buffer pointer must be
        # a named variable (``buffer_ptr``) because a ``LinearBuffer`` can
        # only wrap a variable.
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
        shape = [AugExpr.format(f"{{}}.{c}Size()", ci) for c in "xyz"]
        post = ()
        if acc.inp is not None:
            name = acc.param_name(acc.inp.name, self.target)
            if not cell:
                nodes.append(f"assert(uint_c({name}.size()) == "
                             f"{n_cells(value.n)});")
            ptr = LinearBuffer(PsPointerType(constify(self.dtype)), shape,
                               value.n).var("buffer_ptr")
            buffer = [statement(
                f"{T} const *buffer_ptr = "
                f"{data_pointer(acc.inp.container, name)};",
                (ptr,), (ci, params.inp))]
        else:
            name = acc.out.name
            ptr = LinearBuffer(PsPointerType(self.dtype), shape,
                               value.n).var("buffer_ptr")
            buffer = [
                statement(local_declaration(acc.out.container, T, name,
                                            value.n), (), (ci,)),
                statement(f"{T} *buffer_ptr = "
                          f"{data_pointer(acc.out.container, name)};",
                          (ptr,), ())]
            post = (f"return {name};",)
        self.fn(acc.name, *params.all, returns=self.returns(acc))(
            *nodes, *buffer, *self.map_fields(acc, params, ci),
            self.sfg.map_field(self.cat.buffer(value.n), ptr),
            self.sfg.call(kernel), *post)

    # ---- kind ALL

    def emit_all(self, acc, spec, kernel):
        """Run the kernel on a whole field, ghost layers included; the
        components of the input are kernel parameters."""
        params = self.params(acc, kernel)
        ci, decl = self.whole_interval("xyzSizeWithGhostLayer",
                                       params.fields[0])
        self.fn(acc.name, *params.all)(
            decl, *self.map_fields(acc, params, ci),
            *self.broadcast(acc, spec, params), self.sfg.call(kernel))

    # ---- kind REDUCE

    def emit_reduce(self, acc, _spec, kernel):
        """The kernel adds the values of the cells of the interior to local
        accumulators, in the order of the cells."""
        T = self.T
        params = self.params(acc, kernel)
        result = acc.reduction
        n = result.container.size
        ci, decl = self.whole_interval("xyzSize", params.fields[0])
        # The local accumulators ``acc_<i>`` must be named differently from
        # the kernel parameters ``sum_<i>`` (they are passed by pointer):
        # otherwise pystencils emits ``sum = &sum``. The tests pin this
        # naming.
        names = [f"acc_{i}" for i in range(n)]
        accumulators = [AugExpr(self.dtype).var(name) for name in names]
        value = construction(result.container, T, names)
        if result.scale:
            value = f"{value} * {result.scale}"
        self.fn(acc.name, *params.all, returns=self.returns(acc))(
            decl,
            *(statement(f"{T} {name}{{0}};", (a,), ())
              for name, a in zip(names, accumulators)),
            *self.map_fields(acc, params, ci),
            *(self.sfg.set_param(self.sum_symbol(i),
                                 AugExpr.format("&{}", a))
              for i, a in enumerate(accumulators)),
            self.sfg.call(kernel), f"return {value};")

    # ---- kind POINT

    def emit_point(self, acc, _spec, _kernel):
        """A function of a stencil direction, the velocity and the density
        that returns an expression per direction.

        ``POINT`` is LB-specific (``EquilibriumDistribution``): the
        function takes a ``stencil::Direction``, the velocity ``u`` as a
        ``walberla::Vector3`` (default 0) and the density ``rho``
        (default 1), the cases are named after the ``stencil::`` direction
        constants, and ``point.u_names`` are the components of ``u``.
        """
        T = self.T
        point = acc.point
        direction = Direction(const=True).var("direction")
        u = Vector3(self.dtype, const=True, ref=True).var("u")
        rho = AugExpr(self.dtype).var("rho")
        switch = self.sfg.switch("direction", autobreak=False)
        for name, expr in point.cases:
            value = accessor_equations.print_expression(expr, self.dtype)
            switch.case(f"stencil::{name}")(f"return {value};")
        switch.default('WALBERLA_ABORT("Invalid Direction")')
        self.fn(acc.name, direction, (u, f"walberla::Vector3<{T}>({T}{{0}})"),
                (rho, f"{T}{{1}}"), returns=self.dtype)(
            *(f"const {T} {name} = u[{i}u];"
              for i, name in enumerate(point.u_names)),
            switch)
