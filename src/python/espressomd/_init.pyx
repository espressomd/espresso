#
# Copyright (C) 2013-2019 The ESPResSo project
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
import sys
from . cimport script_interface
from . cimport communication
from libcpp.memory cimport shared_ptr
from boost cimport environment

# Main code
cdef shared_ptr[environment] mpi_env = communication.mpi_init()
communication.init(mpi_env)

# Initialize script interface
# Has to be _after_ mpi_init
script_interface.init(communication.mpiCallbacks())

# Block the worker nodes in the callback loop.
# The head node is just returning to the user script.
if communication.this_node != 0:
    communication.mpi_loop()
    sys.exit(0)
