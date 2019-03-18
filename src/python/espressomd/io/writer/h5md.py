# Copyright (C) 2010-2018 The ESPResSo project
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
"""Interface module for the H5md core implementation."""


import sys

from six import iteritems

from ...script_interface import PScriptInterface  # pylint: disable=import
from ...code_info import features

if not 'H5MD' in features():
    class H5md(object):

        def __init__(self, *args, **kwargs):
            raise RuntimeError("H5md not available.")
else:
    class H5md(object):

        """H5md file object.

        Used for accessing the H5MD core implementation via the
        PScriptInterface.

        .. note::
           Bonds will be written to the file automatically if they exist.

        Parameters
        ----------
        filename : :obj:`str`
                   Name of the trajectory file.
        write_pos : :obj:`bool`, optional
                    If positions should be written.
        write_vel : :obj:`bool`, optional
                    If velocities should be written.
        write_force : :obj:`bool`, optional
                      If forces should be written.
        write_species : :obj:`bool`, optional
                     If types (called 'species' in the H5MD specification) should be written.
        write_mass : :obj:`bool`, optional
                     If masses should be written.
        write_charge : :obj:`bool`, optional
                       If charges should be written.
        write_ordered : :obj:`bool`, optional
                        If particle properties should be ordered according to
                        ids.

        """

        def __init__(self, write_ordered=True, **kwargs):
            self.valid_params = ['filename', "write_ordered"]
            if 'filename' not in kwargs:
                raise ValueError("'filename' parameter missing.")
            self.what = {'write_pos': 1 << 0,
                         'write_vel': 1 << 1,
                         'write_force': 1 << 2,
                         'write_species': 1 << 3,
                         'write_mass': 1 << 4,
                         'write_charge': 1 << 5}
            self.valid_params.append(self.what.keys())
            self.what_bin = 0
            for i, j in iteritems(kwargs):
                if i in self.what.keys():
                    if isinstance(j, bool):
                        if j:
                            self.what_bin += self.what[i]
                    else:
                        raise ValueError(
                            "{} has to be a bool value.".format(i))
                elif i not in self.valid_params:
                    raise ValueError(
                        "Unknown parameter {} for H5MD writer.".format(i))

            self.h5md_instance = PScriptInterface(
                "ScriptInterface::Writer::H5mdScript")
            self.h5md_instance.set_params(filename=kwargs['filename'],
                                          what=self.what_bin,
                                          scriptname=sys.argv[0],
                                          write_ordered=write_ordered)
            self.h5md_instance.call_method("init_file")

        def get_params(self):
            """Get the parameters from the scriptinterface."""
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
