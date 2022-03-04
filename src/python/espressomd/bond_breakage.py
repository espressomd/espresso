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
from .script_interface import script_interface_register, ScriptObjectMap, ScriptInterfaceHelper
from . import interactions


@script_interface_register
class BreakageSpec(ScriptInterfaceHelper):
    """
    Specifications for bond breakage.
    See :ref:`Deleting bonds when particles are pulled apart` for more details.

    Parameters
    ----------
    breakage_length: :obj:`float`
        Maximal bond extension until the bond breaks.
    action_type: :obj:`str`, \{'revert_center_bond', 'revert_vs_bond', 'none'\}
        Action triggered when the bond reaches its maximal extension.

    """

    _so_name = "BondBreakage::BreakageSpec"


@script_interface_register
class BreakageSpecs(ScriptObjectMap):
    _so_name = "BondBreakage::BreakageSpecs"

    def _get_key(self, key):
        """Convert a bond object to a bond id."""
        if isinstance(key, interactions.BondedInteraction):
            key = key._bond_id
            if key == -1:
                raise ValueError("Bond needs to be added to the system first")
        return key

    def __getitem__(self, key):
        return super().__getitem__(self._get_key(key))

    def __setitem__(self, key, value):
        return super().__setitem__(self._get_key(key), value)
