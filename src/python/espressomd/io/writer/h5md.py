#
# Copyright (C) 2010-2026 The ESPResSo project
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
import pathlib

from ...script_interface import script_interface_register, ScriptInterfaceHelper  # pylint: disable=import


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
    file_path : :obj:`str` or :obj:`pathlib.Path`
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
    chunk_size : :obj:`int`
        The chunk size for hdf5 write operations. Must be greater than 0.

    Methods
    -------
    get_params()
        Get the parameters from the script interface.

    valid_fields()
        Get the list of valid fields.

    write()
        Call the H5md write method.

    flush()
        Call the H5md flush method. Only the dataset buffers will be flushed.
        The file superblock, which contains the datasets and metadata address
        space information (End of Allocation, EOA), won't be flushed. Only the
        H5md close method can flush the superblock. A file that was manually
        flushed via this method but not properly closed might have a stale EOA
        and thus be unreadable. This can happen when a simulation receives a
        termination signal and doesn't have time to release all file handles.

    close()
        Close the H5md file and flush both the data and superblock.

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
    chunk_size: :obj:`int`

    """
    _so_name = "ScriptInterface::Writer::H5md"
    _so_features = ("H5MD",)
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("valid_fields", "write", "flush", "close")

    def __init__(self, **kwargs):
        if "sip" in kwargs:
            super().__init__(**kwargs)
            return

        params = self.default_params()
        params.update(kwargs)
        unit_system = params.pop("unit_system")
        if isinstance(params["fields"], str):
            params["fields"] = [params["fields"]]
        script_path = ""
        if sys.argv and sys.argv[0]:
            script_path = str(pathlib.Path(sys.argv[0]).resolve())
        if not isinstance(params["file_path"], (str, pathlib.Path)):
            raise TypeError(
                "Parameter 'file_path' should be a string or a path")
        params["file_path"] = str(pathlib.Path(params["file_path"]).resolve())
        super().__init__(
            script_path=script_path,
            mass_unit=unit_system.mass,
            length_unit=unit_system.length,
            time_unit=unit_system.time,
            force_unit=unit_system.force,
            velocity_unit=unit_system.velocity,
            charge_unit=unit_system.charge,
            **params
        )

    def default_params(self):
        return {"unit_system": UnitSystem(), "fields": "all",
                "chunk_size": 1000}

    def required_keys(self):
        return {"file_path"}

    def valid_keys(self):
        return {"file_path", "unit_system", "fields", "chunk_size"}
