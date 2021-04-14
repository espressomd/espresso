# Copyright (C) 2021 The ESPResSo project
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

unit_test(NAME SingleReaction_test SRC SingleReaction_test.cpp DEPENDS
          EspressoCore)
unit_test(NAME ConstantpHEnsemble_test SRC ConstantpHEnsemble_test.cpp DEPENDS
          EspressoCore Boost::mpi MPI::MPI_CXX)
unit_test(NAME ReactionAlgorithm_test SRC ReactionAlgorithm_test.cpp DEPENDS
          EspressoCore Boost::mpi MPI::MPI_CXX)
unit_test(NAME ReactionEnsemble_test SRC ReactionEnsemble_test.cpp DEPENDS
          EspressoCore Boost::mpi MPI::MPI_CXX)
unit_test(NAME WangLandauReactionEnsemble_test SRC
          WangLandauReactionEnsemble_test.cpp DEPENDS EspressoCore Boost::mpi
          MPI::MPI_CXX)
unit_test(NAME particle_tracking_test SRC particle_tracking_test.cpp DEPENDS
          EspressoCore Boost::mpi MPI::MPI_CXX)
unit_test(NAME reaction_methods_utils_test SRC reaction_methods_utils_test.cpp
          DEPENDS EspressoCore)
