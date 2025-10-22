#
# Copyright (C) 2021-2024 The ESPResSo project
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

include(espresso_unit_test)

espresso_unit_test(SRC SingleReaction_test.cpp DEPENDS espresso::core)
espresso_unit_test(SRC ReactionAlgorithm_test.cpp DEPENDS espresso::core
                   Boost::mpi MPI::MPI_CXX NUM_PROC 2)
espresso_unit_test(SRC particle_tracking_test.cpp DEPENDS espresso::core
                   Boost::mpi MPI::MPI_CXX)
espresso_unit_test(SRC reaction_methods_utils_test.cpp DEPENDS espresso::core)
