#
# Copyright (C) 2019-2022 The ESPResSo project
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
import importlib
import pkg_resources
import unittest

import espressomd
import espressomd.code_info
import espressomd.code_features


def no_skip(x):
    return x


def skipIfMissingFeatures(*args):
    """Unittest skipIf decorator for missing Espresso features."""
    if not espressomd.has_features(*args):
        missing_features = espressomd.code_features.missing_features(*args)
        return unittest.skip("Skipping test: missing feature{} {}".format(
            's' if missing_features else '', ', '.join(missing_features)))
    return no_skip


def skipIfMissingScafacosMethod(method_name):
    """Unittest skipIf decorator for missing ScaFaCoS features."""
    if method_name not in espressomd.code_info.scafacos_methods():
        return unittest.skip(f"ScaFaCoS method '{method_name}' not available")
    return no_skip


def skipIfMissingModules(*args):
    """Unittest skipIf decorator for missing Python modules."""
    if len(args) == 1 and not isinstance(
            args[0], str) and hasattr(args[0], "__iter__"):
        args = set(args[0])
    else:
        args = set(args)
    missing_modules = sorted(set(args) - set(sys.modules.keys()))
    if missing_modules:
        return unittest.skip("Skipping test: missing python module{} {}".format(
            's' if len(missing_modules) > 1 else '', ', '.join(missing_modules)))
    return no_skip


def skipIfMissingGPU():
    """Unittest skipIf decorator for missing GPU."""
    if not espressomd.gpu_available():
        return unittest.skip("Skipping test: no GPU available")
    return no_skip


def skipIfUnmetModuleVersionRequirement(module, version_requirement):
    """Unittest skipIf decorator for unmet module version requirement."""
    try:
        _module = importlib.import_module(module)
    except ImportError:
        return skipIfMissingModules(module)
    if not pkg_resources.packaging.specifiers.SpecifierSet(
            version_requirement).contains(_module.__version__):
        return unittest.skip(
            "Skipping test: version requirement not met for module {}".format(module))
    return no_skip
