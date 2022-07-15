#
# Copyright (C) 2016-2022 The ESPResSo project
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

from . import code_info
from . import utils


class FeaturesError(Exception):

    def __init__(self, missing_features):
        super().__init__(f"Missing features {', '.join(missing_features)}")


def has_features(*args):
    """
    Check whether a list of features is a subset of the compiled-in features.
    """

    lvl = utils.nesting_level(args)
    assert lvl in [1, 2], "has_features() takes strings as argument"
    if lvl == 2:
        check_set = set(args[0])
    else:
        check_set = set(args)

    if not check_set <= code_info.all_features():
        unknown_features = check_set - code_info.all_features()
        raise RuntimeError(f"unknown features {','.join(unknown_features)}")

    return check_set <= set(code_info.features())


def missing_features(*args):
    """
    Return a list of the missing features in the arguments.
    """

    lvl = utils.nesting_level(args)
    assert lvl in [1, 2], "missing_features() takes strings as argument"
    if lvl == 2:
        features = set(args[0])
    else:
        features = set(args)

    return sorted(features - set(code_info.features()))


def assert_features(*args):
    """
    Raise an exception when a list of features is not a subset of the
    compiled-in features.
    """

    if not has_features(*args):
        raise FeaturesError(missing_features(*args))
