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
    action_type: :obj:`str`, {'delete_bond', 'revert_bind_at_point_of_collision', 'none'}
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


@script_interface_register
class BondBreakage(ScriptInterfaceHelper):
    """Bond breakage interface.

    This class provides methods to manage and manipulate bond breakage.

    Methods
    -------
    execute(parameters)
        Executes the bond breakage with the given parameters.
    """

    _so_name = "BondBreakage::BondBreakage"
    _so_bind_methods = ("execute",)

    def execute(self, parameters):
        """
        Execute the bond breakage process.

        Parameters
        ----------
        parameters : dict
            Dictionary containing parameters for the bond breakage process.

        Returns
        -------
        result : any
            Result of the bond breakage execution.
        """
        return self.call_method("execute", parameters=parameters)


@script_interface_register
class BondBreakageInterface(ScriptObjectMap):
    """Interface for managing multiple BondBreakage instances."""

    _so_name = "BondBreakage::BondBreakageInterface"
    _so_bind_methods = BondBreakage._so_bind_methods + ("size", "clear")

    def add(self, bond_breakage):
        """
        Add a BondBreakage instance to the interface.

        Parameters
        ----------
        bond_breakage : BondBreakage
            An instance of BondBreakage to add.
        """

        def _add(self, bond_breakage):
            if isinstance(bond_breakage, BondBreakage):
                self.call_method("add", object=bond_breakage)
            else:
                raise ValueError("Only BondBreakage instances can be added.")

        if isinstance(bond_breakage, collections.abc.Iterable):
            for bb in bond_breakage:
                _add(self, bb)
        else:
            _add(self, bond_breakage)

    def remove(self, bond_breakage):
        """
        Remove a BondBreakage instance from the interface.

        Parameters
        ----------
        bond_breakage : BondBreakage
            An instance of BondBreakage to remove.
        """

        def _remove(self, bond_breakage):
            if isinstance(bond_breakage, BondBreakage):
                self.call_method("remove", object=bond_breakage)
            else:
                raise ValueError("Only BondBreakage instances can be removed.")

        if isinstance(bond_breakage, collections.abc.Iterable):
            for bb in bond_breakage:
                _remove(self, bb)
        else:
            _remove(bond_breakage)
