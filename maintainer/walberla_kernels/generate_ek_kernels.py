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
class, e.g. ``DirichletDoublePrecisionCUDA.h``. The classes live in
``walberla::ek`` (the boundaries) or ``walberla::pystencils`` (the sweeps).

Only the Dirichlet boundary has been migrated; the other electrokinetics
kernels are still the committed output of pystencils 1.4.
"""

import boundaries
import sfg_driver


def generate_dirichlet(sfg, ctx, name):
    boundaries.generate_dirichlet(sfg, ctx, name)


KERNELS = {
    "dirichlet": sfg_driver.Kernel("Dirichlet", ("cpu", "cuda"),
                                   generate_dirichlet),
}


if __name__ == "__main__":
    sfg_driver.generate_main(
        "Generate the waLBerla EK kernels.", "generate_ek_kernels.py",
        KERNELS)
