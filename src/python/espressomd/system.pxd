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
from __future__ import print_function, absolute_import
include "myconfig.pxi"
import particle_data
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "grid.hpp":
    cdef void rescale_boxl(int dir, double d_new)

cdef extern from "communication.hpp" namespace "Random":
    void mpi_random_seed(int cnt, vector[int] & seed)

cdef extern from "forcecap.hpp":
    double forcecap_get()
    void forcecap_set(double forcecap)

from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector
ctypedef vector[string] string_vec
cdef extern from "random.hpp" namespace "Random":
    string mpi_random_get_stat()
    void mpi_random_set_stat(const vector[string] &stat)
    int get_state_size_of_generator()

cdef extern from "utils.hpp":
    void get_mi_vector(double* res,double* a, double* b)

cdef extern from "rotate_system.hpp":
    void rotate_system(double phi, double theta, double alpha)

IF EXCLUSIONS:
    cdef extern from "particle_data.hpp":
        void auto_exclusions(int distance)

cdef bool skin_set

cdef extern from "particle_data.hpp":
    int init_type_map(int type) except +
    int get_random_p_id(int type)  except +
    int number_of_particles_with_type(int type)  except +
