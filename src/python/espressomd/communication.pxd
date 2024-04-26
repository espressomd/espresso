#
# Copyright (C) 2020-2022 The ESPResSo project
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
from boost cimport environment

cdef extern from "MpiCallbacks.hpp" namespace "Communication":
    cppclass MpiCallbacks:
        pass

cdef extern from "communication.hpp":
    shared_ptr[environment] mpi_init()
    void mpi_loop()
    int this_node

cdef extern from "communication.hpp" namespace "Communication":
    MpiCallbacks & mpiCallbacks()
    shared_ptr[MpiCallbacks] mpiCallbacksHandle()
    void init(shared_ptr[environment])
    void deinit()
