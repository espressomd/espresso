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
class, e.g. ``DirichletDoublePrecisionCUDA.h``. The boundary classes live
in ``walberla::ek``, the sweeps in ``walberla::pystencils``.

Only the Dirichlet boundary and the reactions have been migrated; the other
electrokinetics kernels are still the committed output of pystencils 1.4.
The kernels ``reaction_bulk_<N>`` and ``reaction_indexed_<N>`` for
``N`` = 1 to 5 reactants are the classes ``ReactionKernelBulk<N>`` and
``ReactionKernelIndexed<N>``.
"""

import functools

import boundaries
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


KERNELS = {
    "dirichlet": sfg_driver.Kernel("Dirichlet", ("cpu", "cuda"),
                                   boundaries.generate_dirichlet),
    **reaction_kernels(),
}


if __name__ == "__main__":
    sfg_driver.generate_main(
        "Generate the waLBerla EK kernels.", "generate_ek_kernels.py",
        KERNELS)
