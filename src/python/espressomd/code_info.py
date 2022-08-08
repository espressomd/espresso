#
# Copyright (C) 2022 The ESPResSo project
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

from .script_interface import script_interface_register, ScriptInterfaceHelper


@script_interface_register
class _CodeInfo(ScriptInterfaceHelper):
    _so_name = "CodeInfo::CodeInfo"
    _so_creation_policy = "LOCAL"
    _so_bind_methods = (
        "build_type", "features", "all_features", "scafacos_methods"
    )


def build_type():
    """
    Get the CMake build type of this build of ESPResSo.
    Can be e.g. Debug, Release, RelWithAssert, RelWithDebInfo, Coverage, etc.
    """
    return _CodeInfo().build_type()


def features():
    """Get the list of features available in this build of ESPResSo."""
    return _CodeInfo().features()


def all_features():
    """Get the list of all features that can be activated in ESPResSo."""
    return _CodeInfo().all_features()


def scafacos_methods():
    """Lists long-range methods available in the ScaFaCoS library."""
    return _CodeInfo().scafacos_methods()
