# Copyright (C) 2010-2020 The ESPResSo project
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

from . import script_interface


class PyMetric:
    def __init__(self, desc):
        """
        Constructor.

        Parameters
        ----------

        'desc' : string, description of the metric
        """
        self.__instance = script_interface.PScriptInterface(
            "ScriptInterface::GenericDD::Metric")
        self.__instance.set_params(metric=desc)

    def get_metric(self):
        return self.__instance.get_parameter("metric")

    def average(self):
        return self.__instance.call_method("average")

    def maximum(self):
        return self.__instance.call_method("maximum")

    def imbalance(self):
        return self.__instance.call_method("imbalance")

    def _get_instance(self):
        return self.__instance


class GenericDD:
    def __init__(self, grid_type):
        if grid_type not in supported_grid_types():
            raise RuntimeError(
                "Grid type {} not supported by librepa.".format(grid_type))
        self._grid_type = grid_type
        mpi_bcast_generic_dd_grid(grid_type.encode())
        mpi_bcast_cell_structure(CELL_STRUCTURE_GENERIC_DD)
        self.__instance = script_interface.PScriptInterface(
            "ScriptInterface::GenericDD::GenericDD")

    def metric(self, desc):
        """
        Constructs a metric.

        Parameters
        ----------

        'desc' : string, description of the metric
        """
        return PyMetric(desc)

    def repart(self, m, args=None):
        """
        Repartitions the underlying grid.

        Parameters
        ----------

        'm' : :obj:`PyMetric`, metric defining the cell weights
        'args' : :obj:`string`, optional
                 command args processed before the call to repart
        """
        if args is not None:
            self.command(args)
        self.__instance.call_method("repart", metric=m.get_metric())

    def command(self, cmd):
        """
        Deliver an implementation-defined command to the repartitioner.

        Parameters
        ----------

        'cmd' : string, command
        """
        self.__instance.call_method("command", cmd=cmd)

    def grid_type(self):
        """
        Returns the type of grid currently in use.
        """
        return self._grid_type


def supported_grid_types():
    l = librepa_supported_grid_types()
    return [t.decode() for t in l]
