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

function(configure_codegen_kernel kernel)
  set(backends ".cpp")
  if(WALBERLA_BUILD_WITH_CUDA)
    list(APPEND backends "_CUDA.cu")
  endif()
  foreach(precision single double)
    foreach(backend ${backends})
      set(source_file ${kernel}_${precision}_precision${backend})
      target_sources(espresso_walberla_codegen PRIVATE ${source_file})
      espresso_set_walberla_codegen_kernel_source_properties(${source_file})
    endforeach()
  endforeach()
endfunction()

# kernels generated with sweepgen (see maintainer/walberla_kernels), named
# like the ones of the lattice Boltzmann method: <Name><Precision>Precision
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

configure_codegen_kernel(DiffusiveFluxKernel)
configure_codegen_kernel(DiffusiveFluxKernelThermalized)
configure_codegen_kernel(DiffusiveFluxKernelWithElectrostatic)
configure_codegen_kernel(DiffusiveFluxKernelWithElectrostaticThermalized)
configure_codegen_kernel(ContinuityKernel)
configure_codegen_kernel(AdvectiveFluxKernel)
configure_codegen_kernel(FrictionCouplingKernel)
configure_codegen_kernel(FixedFlux)
configure_sweepgen_kernel(Dirichlet)
# keep the range in sync with ReactionKernelSelector::Sizes (in
# walberla_bridge/src/electrokinetics/reactions/ReactionKernelSelector.hpp) and
# with the numbers of reactants in maintainer/walberla_kernels/reactions.py
foreach(num RANGE 1 5)
  configure_sweepgen_kernel(ReactionKernelBulk${num})
  configure_sweepgen_kernel(ReactionKernelIndexed${num})
endforeach()

if(WALBERLA_BUILD_WITH_CUDA)
  foreach(precision single double)
    target_sources(espresso_walberla_codegen
                   PRIVATE EK_FieldAccessors_${precision}_precision_CUDA.cu)
  endforeach()
endif()
