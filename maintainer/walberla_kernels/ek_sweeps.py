#
# Copyright (C) 2022-2026 The ESPResSo project
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
The sweeps of the electrokinetics (EK) that are not reactions, as plain
sweepgen classes of ``walberla::pystencils``: ``DiffusiveFluxKernel``
(variants with an electrostatic potential and with thermal fluctuations),
``AdvectiveFluxKernel``, ``ContinuityKernel`` and ``FrictionCouplingKernel``.

The density ``rho`` of a species in a cell changes by the divergence of the
fluxes through the faces of the cell (the finite-volume method of pystencils
2.1, ``ps.fd.FVM1stOrder``). The flux field is staggered: the flux through
the face between two neighbors, along the direction of the 26 neighbors of a
cell, is stored in the cell that owns the face (see `flux_layout` for the
layout, which has 13 components for each cell). The kernels of one time step
run in this order:

1. ``DiffusiveFluxKernel``: ``flux = -D grad(rho)``, with a term
   ``-D / kT * rho * z * (grad(phi) - f_ext)`` with an electrostatic potential
   ``phi`` and the external field ``f_ext``, and with the noise of the
   thermal fluctuations in the thermalized variants, in each face;
2. ``AdvectiveFluxKernel`` (optional): ``flux += rho u`` through each face by
   the volume-of-fluid method of pystencils (``ps.fd.VOF``), with the
   velocity ``u`` of the fluid;
3. ``FrictionCouplingKernel`` (optional): the force that the fluxes exert on
   the fluid, ``f += kT / (2 D rho_lb) * sum over the 26 faces of the
   flux through the face times the direction``;
4. ``ContinuityKernel``: ``rho -= divergence of the flux``.

Domain of the sweeps
--------------------
The continuity and friction coupling kernels read the fluxes through all the
faces of a cell, so the flux kernels also compute the faces between an
interior cell and a ghost cell: the flux field is never communicated, each
rank computes its own ghost layer. The flux kernels run on the interior of a
block expanded by one cell in each direction (``runOnCellInterval`` with
``field->xyzSize()`` expanded by one, which needs at least one ghost layer),
and `flux_blocks` guards the store of each flux component with
`face_condition`, so that every face that touches an interior cell is
computed exactly once, by the cell that owns it. The density (and the
potential and velocity) must be up to date in the first ghost layer, so the
ghost layers have to be communicated before the flux kernels run. The
continuity and friction coupling kernels run on the interior of the block
(``operator()``).

Gradient on the body diagonals
------------------------------
``ps.fd.FVM1stOrder.discrete_flux`` of pystencils 2.1.0 chooses the
stencil of the gradient along the four body diagonals, ``(+-1, +-1, +-1)``,
from a family of stencils with free weights, and the choice depends on
``PYTHONHASHSEED``: the generated code differed between runs. `discrete_flux`
uses the same stencils as pystencils for the faces along the axes and the
face diagonals and the two-point difference
``d_i f = e_i (f(c + e) - f(c)) / |e|^2`` for the body diagonals ``e``, which
is exact for linear fields like the other stencils.

Thermal fluctuations
--------------------
The noise through a face is ``sqrt(2 D rho_face / |e| / A0) * w`` with a
random number ``w`` of zero mean and unit variance from the Philox RNG of the
pystencils 2.1 runtime (``rho_face`` is the mean density of the two cells,
``A0`` the normalization of the fluxes, the sum of the lengths of the 13
stored directions divided by 3). The random numbers of a cell depend on its
global coordinates and not on the block that it belongs to, so that the noise
doesn't depend on the domain decomposition and the two blocks that share a
face compute the same flux through it. The bridge therefore sets the property
``block_offset`` to the global coordinates of the first cell of the interval
that it runs the sweep on (the lower corner of the block, minus one in each
direction), so that the sweep sees the cell with the counters ``ctr`` at the
global coordinates ``block_offset + ctr``. ``grid_size`` is the size of the
periodic global grid: the random numbers of a cell outside of the global grid
are those of its periodic image. The offsets of the RNG are
``block_offset + grid_size``, which makes all coordinates positive for the
remainder of the periodicity (the first cell of the interval of the block at
the lower end of the global grid has the coordinate -1). ``seed`` and
``time_step`` (the counter that the bridge advances after each time step) are
the key and the counter of the RNG.
"""

import functools

import pystencils as ps
import sympy as sp
from pystencils import DynamicType, Field, TypedSymbol
from pystencils.defaults import DEFAULTS
from pystencils.fd.derivation import FiniteDifferenceStaggeredStencilDerivation
from pystencils.fd.finitevolumes import get_access_and_direction
from pystencils.flow import Store
from pystencils.sympyextensions.random import Philox

import flux_layout
import sfg_driver

DIM = flux_layout.DIM


def scalar_field(name, dtype):
    """Scalar field of one value for each cell (``fzyx``)."""
    return Field.create_generic(name, DIM, dtype, index_dimensions=0,
                                layout="fzyx")


def vector_field(name, dtype):
    """Field of a three-dimensional vector for each cell (``fzyx``)."""
    return Field.create_generic(name, DIM, dtype, index_dimensions=1,
                                index_shape=(DIM,), layout="fzyx")


def face_condition(edge, shape):
    """
    The condition on the counters of the cell loop under which the flux
    through the face between the cell ``c`` and ``c + edge`` is computed,
    for a sweep over the interior of a block of the shape ``shape`` expanded
    by one cell: the counters of the interior cells are 1 to ``shape - 2``,
    and the cells of the ghost layers have the counters 0 and ``shape - 1``.
    The face is computed if ``c`` and ``c + edge`` are in the extended
    region and, on the axes that the edge doesn't cross, ``c`` is an
    interior cell. This is one range of the counter on each axis, which
    keeps the conditions of the stores cheap. All the faces that touch an
    interior cell are computed, and a few faces between ghost cells too.
    """
    counters = DEFAULTS.spatial_counters
    conditions = []
    for i in range(DIM):
        lower = 1 if edge[i] <= 0 else 0
        upper = shape[i] - 1 if edge[i] < 0 else shape[i] - 2
        conditions += [counters[i] >= lower, counters[i] <= upper]
    return sp.And(*conditions)


def flux_blocks(stores, flux, shape, name, preconditions=()):
    """
    One flow block for each stored flux component ``k``: ``stores`` are the
    ``Store`` of the staggered flux field, in the order of
    ``flux_layout.STORED_DIRECTIONS``; the lowered store is guarded by
    `face_condition` of the direction ``k``. ``preconditions`` are the blocks
    that the stores depend on.
    """
    nodes = []
    for edge, store in zip(flux_layout.STORED_DIRECTIONS, stores):
        lhs = flux_layout.lower(store.lhs, flux)
        rhs = flux_layout.lower(store.rhs, flux)

        @ps.flow.guarded_block(face_condition(edge, shape),
                               preds=preconditions,
                               name=f"{name}_{len(nodes)}")
        def block(let, lhs=lhs, rhs=rhs):
            let.store[lhs] = rhs

        nodes.append(block)
    return nodes


def discrete_flux(flux, staggered_flux):
    """
    The stores of the flux ``flux`` of the scalar ``density`` (a vector of
    expressions of the derivatives and values of fields) to the components of
    the staggered flux field: the flux through the face between ``c`` and
    ``c + e`` along ``e`` is the projection of the flux on ``e``, with the
    values at the center of the face and the gradients at the face, divided
    by the length of ``e`` and normalized by ``A0``. This is
    ``ps.fd.FVM1stOrder.discrete_flux`` with a deterministic stencil on the
    body diagonals (see the module documentation).
    """
    def discretize(term, neighbor):
        if isinstance(term, ps.field.Field.Access):
            return (term.get_shifted(*neighbor) + term) * sp.Rational(1, 2)
        if isinstance(term, ps.fd.Diff):
            access, direction = get_access_and_direction(term)
            if sp.Matrix(neighbor).dot(neighbor) > 2:
                # the two-point difference along the body diagonal
                (i,) = direction
                return sp.Rational(neighbor[i], 3) * (
                    access.get_shifted(*neighbor) - access)
            return FiniteDifferenceStaggeredStencilDerivation(
                neighbor, access.field.spatial_dimensions,
                direction).apply(access)
        if term.args:
            return term.func(*[discretize(a, neighbor) for a in term.args])
        return term

    fluxes = flux.applyfunc(ps.fd.derivative.expand_diff_full)
    edges = flux_layout.STORED_DIRECTIONS
    a0 = sum(sp.Matrix(e).norm() for e in edges) / DIM
    stores = []
    for edge in edges:
        directional = sum((fluxes[i] * int(edge[i]) for i in range(DIM)),
                          sp.Integer(0))
        discretized = sp.simplify(discretize(directional, edge))
        value = sp.simplify(discretized / sp.Matrix(edge).norm()) / a0
        stores.append(Store(
            staggered_flux[staggered_flux.stencil.staggered_index(edge)],
            value))
    return stores


def continuity(dtype):
    """
    ``rho -= divergence of the fluxes``, a flow graph for ``generate_sweep``
    over the interior cells. The net flux through the 13 pairs of opposite
    faces of the cell is computed first (``flux(c)[k] - flux(c - e_k)[k]``,
    the flux out of the face of the cell in the direction ``e_k`` and the
    one into the opposite face) and then summed: a uniform flux gives
    exactly zero in any precision, and the density is exactly stationary.
    (``ps.fd.FVM1stOrder.discrete_continuity`` sums the 26 fluxes in an
    unspecified order, whose rounding errors don't cancel: in single
    precision, the density of a uniform species with a uniform flux drifted
    by about the machine precision in each time step.)
    """
    density = scalar_field("density", dtype)
    flux = flux_layout.flux_field(dtype)
    nets = [sp.Symbol(f"net_flux_{k}")
            for k in range(flux_layout.FLUX_COUNT)]

    @ps.flow.block(name="continuity")
    def block(let):
        for k, (net, edge) in enumerate(
                zip(nets, flux_layout.STORED_DIRECTIONS)):
            let[net] = flux.center(k) - flux[tuple(-e for e in edge)](k)
        let.store[density.center] = density.center - sum(
            nets[1:], nets[0])

    return ps.flow.tie(block, name="ContinuityKernel")


def friction_coupling(dtype):
    """
    ``f += kT / (2 D rho_lb) * sum over the 26 faces of the flux through the
    face times its direction``: the assignments of the sweep over the
    interior cells.
    """
    force = vector_field("force", dtype)
    flux = flux_layout.flux_field(dtype)
    staggered_flux = flux_layout.staggered_flux_field(dtype)
    diffusion = TypedSymbol("diffusion", dtype)
    kT = TypedSymbol("kT", dtype)
    lb_density = TypedSymbol("lb_density", dtype)
    total = sp.Matrix([0, 0, 0])
    for direction in flux_layout.NEIGHBORS:
        total += flux_layout.lower(
            staggered_flux.face[0, 0, 0, direction], flux) * sp.Matrix(
                direction)
    center = force.center_vector
    rhs = center + kT / (2 * diffusion * lb_density) * total
    return ps.AssignmentCollection(
        [ps.Assignment(center[i], rhs[i]) for i in range(DIM)])


@functools.lru_cache(maxsize=None)
def advective_flux(dtype):
    """
    ``flux += rho u`` through each face (volume of fluid), a flow graph for
    ``generate_sweep`` over the interior cells expanded by one cell. The
    result is cached: building the flow graph takes seconds, and it is the
    same for the CPU and the CUDA target.
    """
    density = scalar_field("density", dtype)
    velocity = vector_field("velocity", dtype)
    flux = flux_layout.flux_field(dtype)
    vof = ps.fd.VOF(flux_layout.staggered_flux_field(dtype), velocity,
                    density)
    stores = [Store(st.lhs, st.lhs + st.rhs) for st in vof.assignments]
    nodes = flux_blocks(stores, flux, density.shape, "advection")
    return ps.flow.tie(*nodes, name="AdvectiveFluxKernel")


def fluctuation_noise(dtype, density, diffusion):
    """
    The random numbers and the noise of the faces of the thermalized flux:
    returns the flow block that draws the random numbers, and the list of
    the noise added to the flux in each face, in the order of
    ``flux_layout.STORED_DIRECTIONS``.
    """
    edges = flux_layout.STORED_DIRECTIONS
    index_type = DynamicType.INDEX_TYPE
    seed = TypedSymbol("seed", index_type)
    time_step = TypedSymbol("time_step", index_type)
    offsets = tuple(TypedSymbol(f"block_offset_{i}", index_type)
                    for i in range(DIM))
    sizes = tuple(TypedSymbol(f"grid_size_{i}", index_type)
                  for i in range(DIM))
    # the cells of the ghost layers have the coordinates -1 and N: the
    # offset N makes all counters positive, so that the C remainder wraps
    rng = Philox("fluctuation_rng", dtype, seed=seed,
                 offsets=tuple(o + n for o, n in zip(offsets, sizes)),
                 periodicity=(True,) * DIM, num_cells=sizes)
    n_draws = -(-len(edges) // rng.vector_size)
    draws = [rng.get_random_vector(time_step) for _ in range(n_draws)]
    a0 = sum(sp.Matrix(e).norm() for e in edges) / DIM
    noise = []
    for i, edge in enumerate(edges):
        vector, _ = draws[i // rng.vector_size]
        number = vector[i % rng.vector_size]
        density_face = (density.center + density[tuple(edge)]) / 2
        noise.append(
            sp.sqrt(2 * diffusion * density_face / sp.Matrix(edge).norm()
                    / a0) * (number - sp.Rational(1, 2)) * sp.sqrt(12))

    @ps.flow.block(name="rng")
    def rng_block(let):
        for _, assignment in draws:
            let.export[assignment.lhs] = assignment.rhs

    return rng_block, noise


@functools.lru_cache(maxsize=None)
def diffusive_flux_stores(dtype, electrostatic=False):
    """
    The stores of the diffusive flux to the 13 components of the staggered
    flux field (see `discrete_flux`): ``-D grad(rho)``, with ``- D / kT * rho
    * z * (grad(phi) - f_ext)`` if ``electrostatic``. The result is cached:
    the derivation takes seconds.
    """
    density = scalar_field("density", dtype)
    potential = scalar_field("potential", dtype)
    diffusion = TypedSymbol("diffusion", dtype)
    kT = TypedSymbol("kT", dtype)
    valency = TypedSymbol("valency", dtype)
    ext_efield = [TypedSymbol(f"ext_efield_{i}", dtype) for i in range(DIM)]
    expression = -diffusion * sp.Matrix(
        [ps.fd.diff(density, i) for i in range(DIM)])
    if electrostatic:
        field = sp.Matrix([ps.fd.diff(potential, i) - ext_efield[i]
                           for i in range(DIM)])
        expression += -diffusion / kT * density.center * valency * field
    return tuple(discrete_flux(
        expression, flux_layout.staggered_flux_field(dtype)))


@functools.lru_cache(maxsize=None)
def diffusive_flux(dtype, electrostatic=False, thermalized=False):
    """
    The diffusive flux (see `diffusive_flux_stores`) with the noise of the
    thermal fluctuations if ``thermalized``, a flow graph for
    ``generate_sweep`` over the interior cells expanded by one cell. The
    properties of the sweep are ``diffusion`` and the fields ``density`` and
    ``flux``; with ``electrostatic`` the field ``potential`` and the scalars
    ``kT``, ``valency`` and the external field ``ext_efield`` (three
    components); with ``thermalized`` the scalars of the RNG (``seed``,
    ``time_step``, ``block_offset`` and ``grid_size``, see the module
    documentation). The result is cached (see `advective_flux`).
    """
    density = scalar_field("density", dtype)
    flux = flux_layout.flux_field(dtype)
    stores = diffusive_flux_stores(dtype, electrostatic)
    preconditions = []
    if thermalized:
        rng_block, noise = fluctuation_noise(
            dtype, density, TypedSymbol("diffusion", dtype))
        stores = [Store(st.lhs, st.rhs + n) for st, n in zip(stores, noise)]
        preconditions = [rng_block]
    nodes = flux_blocks(stores, flux, density.shape, "flux", preconditions)
    return ps.flow.tie(*nodes, name="DiffusiveFluxKernel")


def generate_diffusive_flux(sfg, ctx, name, electrostatic=False,
                            thermalized=False):
    """Generate the sweep ``name`` of `diffusive_flux`."""
    sfg_driver.generate_sweep(
        sfg, name, diffusive_flux(ctx.data_type, electrostatic, thermalized))


def generate_advective_flux(sfg, ctx, name):
    """Generate the sweep ``name`` of `advective_flux`."""
    sfg_driver.generate_sweep(sfg, name, advective_flux(ctx.data_type))


def generate_continuity(sfg, ctx, name):
    """Generate the sweep ``name`` of `continuity`."""
    sfg_driver.generate_sweep(sfg, name, continuity(ctx.data_type))


def generate_friction_coupling(sfg, ctx, name):
    """Generate the sweep ``name`` of `friction_coupling`."""
    sfg_driver.generate_sweep(sfg, name, friction_coupling(ctx.data_type))
