#
# Copyright (C) 2013-2022 The ESPResSo project
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
class CudaInitHandle(ScriptInterfaceHelper):
    """
    Attributes
    ----------
    device: :obj:`int`
        Device id to use.

    Methods
    -------
    list_devices()
        List devices.

        Returns
        -------
        :obj:`dict` :
            Available CUDA devices sorted by device id.

    """
    _so_name = "System::CudaInitHandle"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("list_devices",)

    def list_devices_properties(self):
        """
        List devices with their properties on each host machine.

        Returns
        -------
        :obj:`dict` :
            Available CUDA devices with their properties sorted by hostname
            and device id.

        """
        out = self.call_method("list_devices_properties")
        for listing in out.values():
            for dev in listing.values():
                dev["compute_capability"] = tuple(dev["compute_capability"])
        return out


def gpu_available():
    """
    Check whether there is at least one compatible GPU available.
    """
    n_compatible_gpus = CudaInitHandle().call_method("get_n_gpus")
    return n_compatible_gpus > 0
