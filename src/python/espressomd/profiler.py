#
# Copyright (C) 2023 The ESPResSo project
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
from .code_features import assert_features
from .script_interface import script_interface_register, ScriptInterfaceHelper


@script_interface_register
class Caliper(ScriptInterfaceHelper):
    """
    Add Caliper section markers at runtime.

    Methods
    -------
    begin_section()
        Start named section.

        Parameters
        ----------
        label : :obj:`str`
            Name of the section.

    end_section()
        End named section.

        Parameters
        ----------
        label : :obj:`str`
            Name of the section.

    """
    _so_name = "ScriptInterface::Profiler::Caliper"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("begin_section", "end_section")

    def __init__(self, *args, **kwargs):
        assert_features(["CALIPER"])
        super().__init__(*args, **kwargs)
