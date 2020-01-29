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
import importlib
import setuptools
import unittest

import espressomd


def _id(x):
    return x


def skipIfMissingFeatures(*args):
    """Unittest skipIf decorator for missing Espresso features."""
    if not espressomd.has_features(*args):
        missing_features = espressomd.missing_features(*args)
        return unittest.skip("Skipping test: missing feature{} {}".format(
            's' if missing_features else '', ', '.join(missing_features)))
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
            's' if missing_modules else '', ', '.join(missing_modules)))
    return _id


def skipIfMissingGPU():
    """Unittest skipIf decorator for missing GPU."""
    if not espressomd.gpu_available():
        return unittest.skip("Skipping test: no GPU available")
    return _id


def skipIfUnmetModuleVersionRequirement(module, version_requirement):
    """Unittest skipIf decorator for unmet module version requirement."""
    try:
        _module = importlib.import_module(module)
    except ImportError:
        return skipIfMissingModules(module)
    if not setuptools.version.pkg_resources.packaging.specifiers.SpecifierSet(
            version_requirement).contains(_module.__version__):
        return unittest.skip(
            "Skipping test: version requirement not met for module {}".format(module))
    return _id
