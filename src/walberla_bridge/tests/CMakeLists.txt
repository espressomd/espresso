# Copyright (C) 2020-2021 The ESPResSo project
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

include(unit_test)

unit_test(NAME LBWalberlaImpl_unit_tests SRC LBWalberlaImpl_unit_tests.cpp DEPENDS WalberlaBridge
          EspressoUtils Boost::mpi ${WALBERLA_LIBS})
target_include_directories(LBWalberlaImpl_unit_tests PRIVATE ${WALBERLA_INCLUDE_DIRS}
                                                   ${walberla_BINARY_DIR}/src)

unit_test(NAME LBWalberlaImpl_statistical_tests SRC LBWalberlaImpl_statistical_tests.cpp DEPENDS WalberlaBridge
          EspressoUtils Boost::mpi ${WALBERLA_LIBS})
target_include_directories(LBWalberlaImpl_statistical_tests PRIVATE ${WALBERLA_INCLUDE_DIRS}
                                                   ${walberla_BINARY_DIR}/src)

unit_test(NAME LBWalberlaImpl_flow_tests SRC LBWalberlaImpl_flow_tests.cpp DEPENDS WalberlaBridge
          EspressoUtils Boost::mpi ${WALBERLA_LIBS})
target_include_directories(LBWalberlaImpl_flow_tests PRIVATE ${WALBERLA_INCLUDE_DIRS}
                                                   ${walberla_BINARY_DIR}/src)