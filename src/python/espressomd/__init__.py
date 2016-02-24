# Copyright (C) 2016 The ESPResSo project
# Copyright (C) 2014 Olaf Lenz
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
#
# Define the espressomd package

# Initialize MPI, start the main loop on the slaves
import espressomd._init

# Initialize the Tcl Interpreter if available
try:
    import espressomd._tcl
    tcl = espressomd._tcl.TclInterpreter()
except ImportError:
    pass

espressomd._init.setup()

from espressomd.system import System
from espressomd.code_info import features
