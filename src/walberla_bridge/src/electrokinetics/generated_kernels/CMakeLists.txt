#
# Copyright (C) 2022-2023 The ESPResSo project
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

foreach(precision double_precision single_precision)
  target_sources(
    espresso_walberla
    PRIVATE DiffusiveFluxKernel_${precision}.cpp
            DiffusiveFluxKernelThermalized_${precision}.cpp
            DiffusiveFluxKernelWithElectrostatic_${precision}.cpp
            DiffusiveFluxKernelWithElectrostaticThermalized_${precision}.cpp
            ContinuityKernel_${precision}.cpp
            AdvectiveFluxKernel_${precision}.cpp
            FrictionCouplingKernel_${precision}.cpp
            FixedFlux_${precision}.cpp
            Dirichlet_${precision}.cpp)
  # disable UBSAN on large files due to excessive compile times
  set_source_files_properties(
    AdvectiveFluxKernel_${precision}.cpp
    DiffusiveFluxKernelWithElectrostatic_${precision}.cpp
    DiffusiveFluxKernelWithElectrostaticThermalized_${precision}.cpp
    TARGET_DIRECTORY espresso_walberla PROPERTIES COMPILE_FLAGS
                                                  -fno-sanitize=undefined)
endforeach()
