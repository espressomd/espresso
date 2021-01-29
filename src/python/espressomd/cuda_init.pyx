#
# Copyright (C) 2013-2019 The ESPResSo project
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
include "myconfig.pxi"
from . cimport cuda_init
from . import utils

cdef class CudaInitHandle:
    def __init__(self):
        IF CUDA != 1:
            raise Exception("CUDA is not compiled in")

    IF CUDA == 1:
        @property
        def device(self):
            """
            Get device.

            Returns
            -------
            :obj:`int` :
                Id of current set device.

            """
            dev = cuda_get_device()
            return dev

        @device.setter
        def device(self, int dev):
            """
            Specify which device to use.

            Parameters
            ----------
            dev : :obj:`int`
                Set the device id of the graphics card to use.

            """
            cuda_set_device(dev)

        def list_devices(self):
            """
            List devices.

            Returns
            -------
            :obj:`dict` :
                List of available CUDA devices.

            """
            cdef char gpu_name_buffer[4 + 64]
            n_gpus = 0
            try:
                n_gpus = cuda_get_n_gpus()
            except RuntimeError:
                pass
            devices = dict()
            for i in range(n_gpus):
                try:
                    cuda_get_gpu_name(i, gpu_name_buffer)
                except RuntimeError:
                    continue
                devices[i] = utils.to_str(gpu_name_buffer)
            return devices

        def list_devices_properties(self):
            """
            List devices with their properties on each host machine.

            Returns
            -------
            :obj:`dict` :
                List of available CUDA devices with their properties.

            """
            cdef vector[EspressoGpuDevice] devices
            cdef EspressoGpuDevice dev
            try:
                devices = cuda_gather_gpus()
            except RuntimeError:
                pass
            resources = dict()
            for i in range(devices.size()):
                dev = devices[i]
                hostname = utils.to_str(dev.proc_name)
                if hostname not in resources:
                    resources[hostname] = {}
                resources[hostname][dev.id] = {
                    'name': utils.to_str(dev.name),
                    'compute_capability': (
                        dev.compute_capability_major,
                        dev.compute_capability_minor
                    ),
                    'cores': dev.n_cores,
                    'total_memory': dev.total_memory,
                }
            return resources

IF CUDA:
    def gpu_available():
        try:
            return cuda_get_n_gpus() > 0
        except RuntimeError:
            return False
ELSE:
    def gpu_available():
        return False
