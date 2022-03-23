# Copyright (C) 2010-2019 The ESPResSo project
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

from ...script_interface import script_interface_register, ScriptInterfaceHelper  # pylint: disable=import
from ...__init__ import assert_features


class UnitSystem:
    """
    Data class for writing H5MD trajectories with
    `physical units <https://nongnu.org/h5md/modules/units.html>`__.
    There are four settable units: 'mass', 'length', 'time', 'charge'.
    Units should be written as strings following the specifications defined
    `here <https://nongnu.org/h5md/modules/units.html#unit-string>`__,
    e.g. ``UnitSystem(time='ps', mass='u', length='nm', charge='e')``.
    """

    def __init__(self, **kwargs):
        self.mass = ''
        self.time = ''
        self.length = ''
        self.charge = ''
        for key, value in kwargs.items():
            assert hasattr(self, key), 'unknown dimension ' + key
            setattr(self, key, value or '')

        if self.length and self.mass and self.time:
            self.force = f'{self.length} {self.mass} {self.time}-2'
        else:
            self.force = ''
        if self.length and self.time:
            self.velocity = f'{self.length} {self.time}-1'
        else:
            self.velocity = ''


@script_interface_register
class H5md(ScriptInterfaceHelper):

    """
    H5md file object.

    .. note::
       Bonds will be written to the file automatically if they exist.
       The pypresso script will be written in the metadata.

    Parameters
    ----------
    file_path : :obj:`str`
        Path to the trajectory file.
    unit_system : :obj:`UnitSystem`, optional
        Physical units for the data.

    Methods
    -------
    get_params()
        Get the parameters from the script interface.

    write()
        Call the H5md write method.

    flush()
        Call the H5md flush method.

    close()
        Close the H5md file.

    Attributes
    ----------
    file_path: :obj:`str`
        Path to the trajectory file.
    script_path: :obj:`str`
        Path to the pypresso script, or empty string for interactive sessions.
    mass_unit: :obj:`str`
    length_unit: :obj:`str`
    time_unit: :obj:`str`
    force_unit: :obj:`str`
    velocity_unit: :obj:`str`
    charge_unit: :obj:`str`

    """
    _so_name = "ScriptInterface::Writer::H5md"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("write", "flush", "close")

    def __init__(self, file_path, unit_system=UnitSystem()):
        assert_features("H5MD")
        super().__init__(
            file_path=file_path,
            script_path=sys.argv[0],
            mass_unit=unit_system.mass,
            length_unit=unit_system.length,
            time_unit=unit_system.time,
            force_unit=unit_system.force,
            velocity_unit=unit_system.velocity,
            charge_unit=unit_system.charge
        )

    def __reduce__(self):
        raise RuntimeError("H5md doesn't support checkpointing")
