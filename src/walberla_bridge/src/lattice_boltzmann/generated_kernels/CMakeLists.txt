#
# Copyright (C) 2021-2023 The ESPResSo project
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

foreach(precision DoublePrecision SinglePrecision)
  target_sources(
    espresso_walberla_codegen
    PRIVATE PackInfoPdf${precision}.cpp PackInfoVec${precision}.cpp
            InitialPDFsSetter${precision}.cpp DynamicUBB${precision}.cpp
            UpdateVelFromPDF${precision}.cpp)
  if(ESPRESSO_BUILD_WITH_WALBERLA_AVX)
    target_sources(
      espresso_walberla_codegen
      PRIVATE StreamCollideSweep${precision}LeesEdwardsAVX.cpp
              StreamCollideSweep${precision}ThermalizedAVX.cpp)
  else()
    target_sources(
      espresso_walberla_codegen
      PRIVATE StreamCollideSweep${precision}LeesEdwards.cpp
              StreamCollideSweep${precision}Thermalized.cpp)
  endif()
  if(WALBERLA_BUILD_WITH_CUDA)
    target_sources(
      espresso_walberla_codegen_cuda
      PRIVATE StreamCollideSweep${precision}LeesEdwardsCUDA.cu
              StreamCollideSweep${precision}ThermalizedCUDA.cu
              FieldAccessors${precision}CUDA.cu
              InitialPDFsSetter${precision}CUDA.cu
              DynamicUBB${precision}CUDA.cu PackInfoPdf${precision}CUDA.cu
              UpdateVelFromPDF${precision}CUDA.cu
              PackInfoVec${precision}CUDA.cu)
  endif()
endforeach()
