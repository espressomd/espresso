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
import numpy as np

from .grid cimport box_geo
from .globals cimport time_step
from .globals cimport mpi_set_time_step
from .globals cimport min_global_cut
from .globals cimport sim_time
from .globals cimport timing_samples
from .globals cimport forcecap_set
from .globals cimport forcecap_get
from .utils import array_locked
from .utils cimport Vector3d, make_array_locked

cdef class Globals:
    property box_l:
        def __set__(self, _box_l):
            cdef Vector3d temp_box_l
            if len(_box_l) != 3:
                raise ValueError("Box length must be of length 3")
            for i in range(3):
                if _box_l[i] <= 0:
                    raise ValueError(
                        "Box length must be > 0  in all directions")
                temp_box_l[i] = _box_l[i]
            box_geo.set_length(temp_box_l)
            mpi_bcast_parameter(FIELD_BOXL)

        def __get__(self):
            return make_array_locked(< Vector3d > box_geo.length())

    property time_step:
        def __set__(self, time_step):
            mpi_set_time_step(time_step)

        def __get__(self):
            global time_step
            return time_step

    property min_global_cut:
        def __set__(self, _min_global_cut):
            global min_global_cut
            min_global_cut = _min_global_cut
            mpi_bcast_parameter(FIELD_MIN_GLOBAL_CUT)

        def __get__(self):
            global min_global_cut
            return min_global_cut

    property periodicity:
        def __set__(self, _periodic):
            for i in range(3):
                box_geo.set_periodic(i, _periodic[i])

            mpi_bcast_parameter(FIELD_PERIODIC)

        def __get__(self):
            periodicity = np.zeros(3)

            for i in range(3):
                periodicity[i] = box_geo.periodic(i)

            return array_locked(periodicity)

    property time:
        def __set__(self, double _time):
            global sim_time
            sim_time = _time
            mpi_bcast_parameter(FIELD_SIMTIME)

        def __get__(self):
            global sim_time
            return sim_time

    property timings:
        def __set__(self, int _timings):
            global timing_samples
            if _timings <= 0:
                timing_samples = 0
            else:
                timing_samples = _timings

        def __get__(self):
            global timing_samples
            return timing_samples

    property force_cap:
        def __set__(self, cap):
            forcecap_set(cap)

        def __get__(self):
            return forcecap_get()

    def __getstate__(self):
        state = {'box_l': self.box_l,
                 'time_step': self.time_step,
                 'min_global_cut': self.min_global_cut,
                 'periodicity': self.periodicity,
                 'time': self.time,
                 'timings': self.timings,
                 'force_cap': self.force_cap}
        return state

    def __setstate__(self, state):
        self.box_l = state['box_l']
        self.time_step = state['time_step']
        self.min_global_cut = state['min_global_cut']
        self.periodicity = state['periodicity']
        self.time = state['time']
        self.timings = state['timings']
        self.force_cap = state['force_cap']
