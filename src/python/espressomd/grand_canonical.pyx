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
from utils cimport *
cimport globals
from globals cimport max_seen_particle


def is_valid_type(type_id):
    return (not isinstance(type_id, int) or type_id < 0 or type_id > globals.n_particle_types)


def check_valid_type(type_id):
    if is_valid_type(type_id):
        raise ValueError("type", type_id, "does not exist!")


def setup(type_list=None):
    if not hasattr(type_list, "__iter__"):
        raise ValueError("type_list has to be iterable.")

    for type_id in type_list:
        if (type_id < 0):
            raise ValueError("type", type_id, "is invalid!")
        if (max_seen_particle < 0):
            raise ValueError(
                "The system contains no particles. Create one particle with arbitrary type first!")
        init_type_array(type_id)


def delete_particles(type_id=None):
    check_valid_type(type_id)
    delete_particle_of_type(type_id)


def find_particle(type_id=None):
    check_valid_type(type_id)
    cdef int pid
    find_particle_type(type_id, & pid)
    return int(pid)


def number_of_particles(type_id=None):
    check_valid_type(type_id)
    cdef int number
    number_of_particles_with_type(type_id, & number)
    return int(number)
