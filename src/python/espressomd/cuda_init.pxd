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

from libcpp.vector cimport vector

cdef extern from "cuda_init.hpp":
    cdef struct EspressoGpuDevice:
        int id
        char name[64]
        char proc_name[64]
        int node
        int compute_capability_major
        int compute_capability_minor
        size_t total_memory
        int n_cores

    void cuda_set_device(int dev) except +
    int cuda_get_device() except +
    int cuda_get_n_gpus() except +
    void cuda_get_gpu_name(int dev, char name[64]) except +
    vector[EspressoGpuDevice] cuda_gather_gpus()
