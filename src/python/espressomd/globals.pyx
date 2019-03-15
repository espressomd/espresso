import numpy as np

cimport grid
from globals cimport time_step
from globals cimport mpi_set_time_step
from globals cimport min_global_cut
from grid cimport periodic
from globals cimport sim_time
from globals cimport timing_samples
from globals cimport forcecap_set
from globals cimport forcecap_get
from espressomd.utils import array_locked, is_valid_type
from utils cimport Vector3d

cdef class Globals(object):
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
            grid.box_l = temp_box_l
            mpi_bcast_parameter(FIELD_BOXL)

        def __get__(self):
            return array_locked(np.array([grid.box_l[0], grid.box_l[1], grid.box_l[2]]))

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
            global periodic
            periodic = 4 * _periodic[2] + 2 * _periodic[1] + _periodic[0]
            # first 3 bits of periodic determine the periodicity
            mpi_bcast_parameter(FIELD_PERIODIC)

        def __get__(self):
            global periodic
            periodicity = np.zeros(3)
            periodicity[0] = periodic % 2
            periodicity[1] = int(periodic / 2) % 2
            periodicity[2] = int(periodic / 4) % 2
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
