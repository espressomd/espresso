# Copyright (C) 2019 The ESPResSo project
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
import sys
import unittest
import espressomd  # pylint: disable=import-error
from espressomd.utils import to_str


def _id(x):
    return x


def skipIfMissingFeatures(*args):
    """Unittest skipIf decorator for missing Espresso features."""

    if not espressomd.has_features(*args):
        missing_features = espressomd.missing_features(*args)
        return unittest.skip("Skipping test: missing feature{} {}".format(
            's' if len(missing_features) else '', ', '.join(missing_features)))
    return _id


def skipIfMissingModules(*args):
    """Unittest skipIf decorator for missing Python modules."""

    if len(args) == 1 and not isinstance(
            args[0], str) and hasattr(args[0], "__iter__"):
        args = set(args[0])
    else:
        args = set(args)
    missing_modules = set(args) - set(sys.modules.keys())
    if missing_modules:
        return unittest.skip("Skipping test: missing python module{} {}".format(
            's' if len(missing_modules) else '', ', '.join(missing_modules)))
    return _id


def skipIfMissingGPU():
    """Unittest skipIf decorator for missing GPU."""

    if not espressomd.gpu_available():
        return unittest.skip("Skipping test: no GPU available")
    devices = espressomd.cuda_init.CudaInitHandle().device_list
    current_device_id = espressomd.cuda_init.CudaInitHandle().device
    return _id
