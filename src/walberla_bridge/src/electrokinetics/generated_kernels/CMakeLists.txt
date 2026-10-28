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

# kernels generated with sweepgen (see maintainer/walberla_kernels), named like
# the ones of the lattice Boltzmann method: <Name><Precision>Precision
function(configure_sweepgen_kernel kernel)
  set(sweepgen_backends ".cpp")
  if(WALBERLA_BUILD_WITH_CUDA)
    list(APPEND sweepgen_backends "CUDA.cu")
  endif()
  foreach(precision Single Double)
    foreach(backend ${sweepgen_backends})
      set(source_file ${kernel}${precision}Precision${backend})
      target_sources(espresso_walberla_codegen PRIVATE ${source_file})
      espresso_set_walberla_codegen_kernel_source_properties(${source_file})
    endforeach()
  endforeach()
endfunction()

configure_sweepgen_kernel(AdvectiveFluxKernel)
configure_sweepgen_kernel(ContinuityKernel)
configure_sweepgen_kernel(DiffusiveFluxKernel)
configure_sweepgen_kernel(DiffusiveFluxKernelThermalized)
configure_sweepgen_kernel(DiffusiveFluxKernelWithElectrostatic)
configure_sweepgen_kernel(DiffusiveFluxKernelWithElectrostaticThermalized)
configure_sweepgen_kernel(Dirichlet)
configure_sweepgen_kernel(FixedFlux)
configure_sweepgen_kernel(FrictionCouplingKernel)
configure_sweepgen_kernel(PackInfoScalar)
# keep the range in sync with ReactionKernelSelector::Sizes (in
# walberla_bridge/src/electrokinetics/reactions/ReactionKernelSelector.hpp) and
# with the numbers of reactants in maintainer/walberla_kernels/reactions.py
foreach(num RANGE 1 5)
  configure_sweepgen_kernel(ReactionKernelBulk${num})
  configure_sweepgen_kernel(ReactionKernelIndexed${num})
endforeach()

# the accessors are header-only on the CPU
if(WALBERLA_BUILD_WITH_CUDA)
  foreach(precision Single Double)
    set(source_file EKFieldAccessors${precision}PrecisionCUDA.cu)
    target_sources(espresso_walberla_codegen PRIVATE ${source_file})
    espresso_set_walberla_codegen_kernel_source_properties(${source_file})
  endforeach()
endif()
