# Copyright (C) 2016-2019 The ESPResSo project
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
from . import _init

from .system import System
from .code_info import features, all_features
from .cuda_init import gpu_available


class FeaturesError(Exception):

    def __init__(self, missing_features_list):
        message = "Missing features " + ", ".join(missing_features_list)
        super().__init__(message)


def has_features(*args):
    """Tests whether a list of features is a subset of the compiled-in features"""

    if len(args) == 1 and not isinstance(
            args[0], str) and hasattr(args[0], "__iter__"):
        check_set = set(args[0])
    else:
        check_set = set(args)

    if not check_set < all_features():
        raise RuntimeError(
            "'{}' is not a feature".format(','.join(check_set - all_features())))

    return check_set <= set(features())


def missing_features(*args):
    """Returns a list of the missing features in the argument"""

    if len(args) == 1 and not isinstance(
            args[0], str) and hasattr(args[0], "__iter__"):
        return set(args[0]) - set(features())

    return set(args) - set(features())


def assert_features(*args):
    """Raises an exception when a list of features is not a subset of the compiled-in features"""

    if not has_features(*args):
        raise FeaturesError(missing_features(*args))
