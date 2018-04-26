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

from espressomd.system import System
from espressomd.code_info import features

def has_features(*args):
    """Tests whether a list of features is a subset of the compiled-in features"""

    if len(args) == 1 and type(args[0]) is not str and hasattr(args[0], "__iter__"):
        return set(args[0]) < set(features())

    return set(args) < set(features())


def missing_features(*args):
    """Returns a list of the missing features in the argument"""

    if len(args) == 1 and type(args[0]) is not str and hasattr(args[0], "__iter__"):
            return set(args[0]) - set(features())

    return set(args) - set(features())


def assert_features(ExceptionType = Exception, *args):
    """Raises an excpetion when a list of features is not a subset of the compiled-in features"""

    if not has_features(*args):
        raise ExceptionType("Missing features " + ", ".join(missing_features(*args)))
