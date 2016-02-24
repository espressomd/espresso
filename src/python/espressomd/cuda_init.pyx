#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
cimport cuda_init

cdef class CudaInitHandle:
    def __init__(self):
        IF CUDA != 1:
            raise Exception("Cuda is not compiled in")

    property device:
        """cuda device to use"""

        IF CUDA == 1:
            def __set__(self, int _dev):
                if cuda_set_device(_dev):
                    raise Exception("cuda device set error")

            def __get__(self):
                dev = cuda_get_device()
                if dev == -1:
                    raise Exception("cuda device get error")
                return dev

    # property device_list:
    #   IF CUDA == 1:
    #     def __set__(self, int _dev):
    #       raise Exception("cuda device list is read only")
    #     def __get__(self):
    #       cdef int _p_devl
    #       cdef char _devname[4+64]
    #       if getdevicelist(&_p_devl, _devname):
    #         raise Exception("cuda devicelist error")
    #       return _devname
