#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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

cdef extern from "communication.hpp":
    void mpi_init(int * argc, char ** *argv)
    int this_node

cdef extern from "initialize.hpp":
    void on_program_start()
    void mpi_loop()

# Main code
mpi_init(NULL, NULL)

# Main slave loop
if this_node != 0:
    on_program_start()
    mpi_loop()
    sys.exit()


def setup():
    on_program_start()
