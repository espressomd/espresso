#
# Copyright (C) 2020-2026 The ESPResSo project
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
from libcpp.memory cimport shared_ptr
from libcpp.memory cimport unique_ptr

cdef extern from "MpiCallbacks.hpp" namespace "Communication":
    cppclass MpiCallbacks:
        pass

cdef extern from "communication.hpp":
    cppclass CommunicationEnvironment:
        CommunicationEnvironment()
        shared_ptr[MpiCallbacks] mpiCallbacksHandle()
    unique_ptr[CommunicationEnvironment] communication_environment
    void mpi_loop()
    int this_node
