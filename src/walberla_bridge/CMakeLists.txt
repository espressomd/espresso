#
# Copyright (C) 2020 The ESPResSo project
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
set(WalberlaBridge_SRC WalberlaBlockForest.hpp
                       generated_kernels/StreamSweep.cpp generated_kernels/CollideSweep.cpp generated_kernels/CollideSweepThermalized.cpp
                       generated_kernels/InitialPDFsSetter.cpp generated_kernels/Dynamic_UBB.cpp
                       lb_walberla_init.cpp
                       generated_kernels/DiffusiveFluxKernel.cpp generated_kernels/ContinuityKernel.cpp
                       ekin_walberla_init.cpp)
if (WALBERLA_USE_AVX)
  list(APPEND WalberlaBridge_SRC generated_kernels/CollideSweepAVX.cpp generated_kernels/CollideSweepThermalizedAVX.cpp)
endif()
add_library(WalberlaBridge SHARED ${WalberlaBridge_SRC})
target_link_libraries(WalberlaBridge PUBLIC MPI::MPI_CXX EspressoUtils
                      PRIVATE Boost::boost cxx_interface ${WALBERLA_LIBS})
target_include_directories(
        WalberlaBridge PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${WALBERLA_INCLUDE_DIRS}
        PRIVATE ${walberla_BINARY_DIR}/src)

install(TARGETS WalberlaBridge LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)
