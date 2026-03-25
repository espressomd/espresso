#
# Copyright (C) 2013-2026 The ESPResSo project
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
import atexit
from . cimport script_interface
from . cimport communication
from .communication cimport CommunicationEnvironment
from .communication cimport communication_environment
from libcpp.memory cimport make_unique

# Main code
communication_environment = make_unique[CommunicationEnvironment]()

# Initialize script interface
# Has to be _after_ mpi_init
script_interface.init(communication_environment.get().mpiCallbacksHandle())


def session_shutdown():
    script_interface.deinit()
    communication_environment.reset()


atexit.register(session_shutdown)

# Block the worker nodes in the callback loop.
# The head node is just returning to the user script.
if communication.this_node != 0:
    communication.mpi_loop()
    sys.exit(0)
