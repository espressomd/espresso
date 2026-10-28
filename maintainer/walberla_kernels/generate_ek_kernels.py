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

"""
Generate the electrokinetics kernels of the waLBerla bridge with sweepgen and
pystencils-sfg. Each kernel is written to its own file pair, named after its
class, e.g. ``DirichletDoublePrecisionCUDA.h``. The boundary classes live in
``walberla::ek``, the sweeps in ``walberla::pystencils``, the accessors in
``walberla::ek::accessor``.

* ``diffusive_flux[_electrostatic][_thermalized]``, ``advective_flux``,
  ``continuity`` and ``friction_coupling``: the sweeps
  ``DiffusiveFluxKernel[WithElectrostatic][Thermalized]``,
  ``AdvectiveFluxKernel``, ``ContinuityKernel`` and
  ``FrictionCouplingKernel`` (see ``ek_sweeps.py``), which use the layout of
  the flux field of ``flux_layout.py``.
* ``reaction_bulk_<N>`` and ``reaction_indexed_<N>`` for ``N`` = 1 to 5
  reactants: the classes ``ReactionKernelBulk<N>`` and
  ``ReactionKernelIndexed<N>`` (see ``reactions.py``).
* ``dirichlet`` and ``fixed_flux``: the boundaries ``Dirichlet`` (the density
  of the boundary cells) and ``FixedFlux`` (the flux through the faces of the
  boundary cells), see ``boundaries.py``.
* ``field_accessors`` (CPU) and ``field_accessors_cuda``: the accessors of
  the fields of the electrokinetics, ``EKFieldAccessors``, see
  ``ek_accessors.py``.

There is no AVX variant of any of them.
"""

import functools

import boundaries
import ek_accessors
import ek_sweeps
import reactions
import sfg_driver


def reaction_kernels():
    """The reaction kernels of all numbers of reactants."""
    kernels = {}
    for n in range(1, reactions.MAX_NUM_REACTANTS + 1):
        for kind, generate in (("bulk", reactions.generate_bulk),
                               ("indexed", reactions.generate_indexed)):
            kernels[f"reaction_{kind}_{n}"] = sfg_driver.Kernel(
                f"ReactionKernel{kind.capitalize()}{n}", ("cpu", "cuda"),
                functools.partial(generate, num_reactants=n))
    return kernels


def flux_kernels():
    """The diffusive flux kernels: with or without an electrostatic
    potential, with or without thermal fluctuations."""
    kernels = {}
    for electrostatic in (False, True):
        for thermalized in (False, True):
            key = "diffusive_flux"
            stem = "DiffusiveFluxKernel"
            if electrostatic:
                key += "_electrostatic"
                stem += "WithElectrostatic"
            if thermalized:
                key += "_thermalized"
                stem += "Thermalized"
            kernels[key] = sfg_driver.Kernel(
                stem, ("cpu", "cuda"), functools.partial(
                    ek_sweeps.generate_diffusive_flux,
                    electrostatic=electrostatic, thermalized=thermalized))
    return kernels


def generate_accessors(sfg, ctx, name):
    ek_accessors.generate_field_accessors(sfg, ctx)


KERNELS = {
    "dirichlet": sfg_driver.Kernel("Dirichlet", ("cpu", "cuda"),
                                   boundaries.generate_dirichlet),
    "fixed_flux": sfg_driver.Kernel("FixedFlux", ("cpu", "cuda"),
                                    boundaries.generate_fixed_flux),
    **flux_kernels(),
    "advective_flux": sfg_driver.Kernel(
        "AdvectiveFluxKernel", ("cpu", "cuda"),
        ek_sweeps.generate_advective_flux),
    "field_accessors": sfg_driver.Kernel(
        "EKFieldAccessors", ("cpu",), generate_accessors, header_only=True),
    "field_accessors_cuda": sfg_driver.Kernel(
        "EKFieldAccessors", ("cuda",), generate_accessors,
        header_extension="cuh"),
    "continuity": sfg_driver.Kernel(
        "ContinuityKernel", ("cpu", "cuda"), ek_sweeps.generate_continuity),
    "friction_coupling": sfg_driver.Kernel(
        "FrictionCouplingKernel", ("cpu", "cuda"),
        ek_sweeps.generate_friction_coupling),
    **reaction_kernels(),
}


if __name__ == "__main__":
    sfg_driver.generate_main(
        "Generate the waLBerla EK kernels.", "generate_ek_kernels.py",
        KERNELS)
