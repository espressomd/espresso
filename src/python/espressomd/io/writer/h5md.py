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
"""Interface module for the H5md core implementation."""

import sys

from ...script_interface import PScriptInterface  # pylint: disable=import
from ...code_info import features


class UnitSystem:
    """
    Data class for writing H5MD trajectories with
    `physical units <https://nongnu.org/h5md/modules/units.html>`_.
    There are four settable units: 'mass', 'length', 'time', 'charge'.
    Units should be written as strings following the specifications defined
    `here <https://nongnu.org/h5md/modules/units.html#unit-string>`_,
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


if 'H5MD' not in features():
    class H5md:
        def __init__(self, *args, **kwargs):
            raise RuntimeError("H5md not available.")
else:
    class H5md:

        """H5md file object.

        Used for accessing the H5MD core implementation.

        .. note::
           Bonds will be written to the file automatically if they exist.

        Parameters
        ----------
        file_path : :obj:`str`
            Path to the trajectory file.
        unit_system : :obj:`UnitSystem`, optional
            Physical units for the data.

        """

        def __init__(self, file_path, unit_system=UnitSystem()):
            self.h5md_instance = PScriptInterface(
                "ScriptInterface::Writer::H5md", file_path=file_path, script_path=sys.argv[0],
                mass_unit=unit_system.mass, length_unit=unit_system.length, 
                time_unit=unit_system.time,
                force_unit=unit_system.force,
                velocity_unit=unit_system.velocity,
                charge_unit=unit_system.charge
            )

        def get_params(self):
            """Get the parameters from the script interface."""
            return self.h5md_instance.get_params()

        def write(self):
            """Call the H5md write method."""
            self.h5md_instance.call_method("write")

        def flush(self):
            """Call the H5md flush method."""
            self.h5md_instance.call_method("flush")

        def close(self):
            """Close the H5md file."""
            self.h5md_instance.call_method("close")
