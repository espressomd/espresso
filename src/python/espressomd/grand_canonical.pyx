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
from espressomd.utils cimport *
cimport globals
from globals cimport max_seen_particle

es_error=1
def is_valid_type(current_type):
    return (not (isinstance(current_type, int) or current_type < 0 or current_type > globals.n_particle_types))


def check_valid_type(current_type):
    if is_valid_type(current_type):
        raise ValueError("type", current_type, "does not exist!")


def setup(type_list=None):
    """
    For using Espresso conveniently for simulations in the grand canonical
    ensemble, or other purposes, when particles of certain types are created
    and deleted frequently. Particle ids can be stored in lists for each
    individual type and so random ids of particles of a certain type can be
    drawn. If you want Espresso to keep track of particle ids of a certain type
    you have to initialize the method by calling the setup function. After that
    Espresso will keep track of particle ids of that type.

    """
    if not hasattr(type_list, "__iter__"):
        raise ValueError("type_list has to be iterable.")

    for current_type in type_list:
        if (current_type < 0):
            raise ValueError("type", current_type, "is invalid!")
        if (max_seen_particle < 0):
            raise ValueError(
                "The system contains no particles. Create one particle with arbitrary type first!")
        status=init_type_array(current_type)
        if status==es_error:
            raise Exception("gc init failed")
        handle_errors("init_type_array -> updatePartCfg failed")


def number_of_particles(current_type=None):
    """
    Parameters
    ----------
    current_type : :obj:`int` (:attr:`espressomd.particle_data.ParticleHandle.type`)
                   Particle type to count the number for. 

    Returns
    -------
    :obj:`int`
        The number of particles which share the given type.

    """
    check_valid_type(current_type)
    cdef int number
    if ( number_of_particles_with_type(current_type, &number) == -3 ):
        raise Exception("no list for particle type ", current_type)
    number_of_particles_with_type(current_type, & number)
    return int(number)

def find_particle(current_type=None):
    """
    The command will return a randomly chosen particle id, for a particle of
    the given type.
    
    """
    check_valid_type(current_type)
    cdef int pid
    status=find_particle_type(current_type, & pid)
    if(status== es_error):
        print("error no particle found")
        return -1
    else:
        return int(pid)

def status(current_type=None):
    """
    Parameters
    ----------
    current_type : :obj:`int` (:attr:`espressomd.particle_data.ParticleHandle.type`)
                   Particle type to get the ids for. 

    Returns
    -------
    list of :obj:`int`
        The id list for particles which share the given type.

    """
    check_valid_type(current_type)
    if ( (type_array!=NULL) and type_array[Index.type[current_type]].max_entry!= 0 ):
        indexed=0;
        for i in range(Type.max_entry):
            if (current_type==Type.index[i]):
                indexed=1;
                break;
    if ( indexed==1 ):
        id_list=[]
        for i in range( type_array[Index.type[current_type]].max_entry):
            id_list.append(type_array[Index.type[current_type]].id_list[i])
        return id_list
    else:
        print("no list for particle")
        return []
