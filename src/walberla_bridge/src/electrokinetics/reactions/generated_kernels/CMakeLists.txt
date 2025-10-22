#
# Copyright (C) 2022-2025 The ESPResSo project
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
      foreach(num RANGE 1 5)
        set(source_file ${kernel}_${num}_${precision}_precision${backend})
        target_sources(espresso_walberla_codegen PRIVATE ${source_file})
        espresso_set_walberla_codegen_kernel_source_properties(${source_file})
      endforeach()
    endforeach()
  endforeach()
endfunction()

configure_codegen_kernel(ReactionKernelBulk)
configure_codegen_kernel(ReactionKernelIndexed)
