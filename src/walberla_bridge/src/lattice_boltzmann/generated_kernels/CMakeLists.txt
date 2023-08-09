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

target_sources(
  espresso_walberla
  PRIVATE InitialPDFsSetterSinglePrecision.cpp
          InitialPDFsSetterDoublePrecision.cpp Dynamic_UBB_single_precision.cpp
          Dynamic_UBB_double_precision.cpp VelFieldUpdateSinglePrecision.cpp
          VelFieldUpdateDoublePrecision.cpp)
if(ESPRESSO_BUILD_WITH_WALBERLA_AVX)
  target_sources(
    espresso_walberla
    PRIVATE StreamCollideSweepSinglePrecisionLeesEdwardsAVX.cpp
            StreamCollideSweepDoublePrecisionLeesEdwardsAVX.cpp
            StreamCollideSweepSinglePrecisionThermalizedAVX.cpp
            StreamCollideSweepDoublePrecisionThermalizedAVX.cpp)
else()
  target_sources(
    espresso_walberla
    PRIVATE StreamCollideSweepSinglePrecisionLeesEdwards.cpp
            StreamCollideSweepDoublePrecisionLeesEdwards.cpp
            StreamCollideSweepSinglePrecisionThermalized.cpp
            StreamCollideSweepDoublePrecisionThermalized.cpp)
endif()
