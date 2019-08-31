#
# Copyright (C) 2013-2018 The ESPResSo project
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
from . import script_interface

cdef extern from "communication.hpp":
    void mpi_init()
    void mpi_loop()
    int this_node

# Main code
mpi_init()

# Initialize script interface
# Has to be _after_ mpi_init
script_interface.init()

# Block the slaves in the callback loop
# The master is just returning to the user script
if this_node != 0:
    mpi_loop()
    sys.exit(0)
