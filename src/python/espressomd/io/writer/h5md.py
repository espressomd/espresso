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
from ... import utils


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
            assert hasattr(self, key), f'unknown dimension {key}'
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
       Bonds will be written to the file if they exist.
       The pypresso script will be written in the metadata.

    Parameters
    ----------
    file_path : :obj:`str`
        Path to the trajectory file, or an existing file to append data to
        (it must have the same specifications).
    unit_system : :obj:`UnitSystem`, optional
        Physical units for the data.
    fields : :obj:`set` or :obj:`str`, optional
        List of fields to write to the trajectory file. Defaults to ``'all'``.
        See :meth:`~espressomd.io.writer.h5md.H5md.valid_fields()` for the
        list of valid fields. This list defines the H5MD specifications.
        If the file in ``file_path`` already exists but has different
        specifications, an exception is raised.

    Methods
    -------
    get_params()
        Get the parameters from the script interface.

    valid_fields()
        Get the list of valid fields.

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
    fields: :obj:`list`
        List of fields to write to the trajectory file.
    mass_unit: :obj:`str`
    length_unit: :obj:`str`
    time_unit: :obj:`str`
    force_unit: :obj:`str`
    velocity_unit: :obj:`str`
    charge_unit: :obj:`str`

    """
    _so_name = "ScriptInterface::Writer::H5md"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("valid_fields", "write", "flush", "close")

    def __init__(self, **kwargs):
        assert_features("H5MD")

        if "sip" in kwargs:
            super().__init__(**kwargs)
            return

        params = self.default_params()
        params.update(kwargs)
        unit_system = params["unit_system"]
        fields = params["fields"]
        fields = [fields] if isinstance(fields, str) else list(fields)
        params["fields"] = fields
        self.validate_params(params)
        super().__init__(
            file_path=params["file_path"],
            script_path=sys.argv[0],
            fields=fields,
            mass_unit=unit_system.mass,
            length_unit=unit_system.length,
            time_unit=unit_system.time,
            force_unit=unit_system.force,
            velocity_unit=unit_system.velocity,
            charge_unit=unit_system.charge
        )

    def default_params(self):
        return {"unit_system": UnitSystem(), "fields": "all"}

    def required_keys(self):
        return {"file_path"}

    def valid_keys(self):
        return {"file_path", "unit_system", "fields"}

    def validate_params(self, params):
        """Check validity of given parameters.
        """
        utils.check_type_or_throw_except(
            params["file_path"], 1, str, "'file_path' should be a string")
        for item in params["fields"]:
            utils.check_type_or_throw_except(
                item, 1, str, "'fields' should be a string or a list of strings")
