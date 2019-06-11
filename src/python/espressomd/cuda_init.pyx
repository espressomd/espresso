#
# Copyright (C) 2013-2018 The ESPResSo project
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
from __future__ import print_function, absolute_import
include "myconfig.pxi"
from . cimport cuda_init

cdef class CudaInitHandle(object):
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
            if dev == -1:
                raise Exception("cuda device get error")
            return dev

        @device.setter
        def device(self, int _dev):
            """
            Specify which device to use.

            Parameters
            ----------
            dev : :obj:`int`
                Set the device id of the graphics card to use.

            """
            cuda_set_device(_dev)

    IF CUDA == 1:
        @property
        def device_list(self):
            """
            List devices.

            Returns
            -------
            :obj:`list` :
                List of available CUDA devices.

            """
            cdef char gpu_name_buffer[4 + 64]
            devices = dict()
            for i in range(cuda_get_n_gpus()):
                cuda_get_gpu_name(i, gpu_name_buffer)
                devices[i] = gpu_name_buffer
            return devices

        @device_list.setter
        def device_list(self, dict _dev_dict):
            raise Exception("cuda device list is read only")


IF CUDA:
    def gpu_available():
        return cuda_get_n_gpus() > 0
ELSE:
    def gpu_available():
        return False
