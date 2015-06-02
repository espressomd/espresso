#
# Copyright (C) 2013,2014 The ESPResSo project
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
# For C-extern Analysis
include "myconfig.pxi"
cimport c_analyze
cimport utils
cimport particle_data
import utils
import code_info
import particle_data
from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector
from interactions import *
from interactions cimport *
import numpy as np
cimport numpy as np
#
# Minimal distance between particles
#


def mindist(system=None, p1='default', p2='default'):
    """Minimal distance between particles
      mindist(p1="default",p2="default")

      p1, p2: lists of particle types
    """

    cdef IntList * set1
    cdef IntList * set2

    if p1 == 'default' and p2 == 'default':
        result = c_analyze.mindist(NULL, NULL)
    elif (p1 == 'default' and not p2 == 'default') or \
         (not p1 == 'default' and p2 == 'default'):
        raise Exception("Both, p1 and p2 have to be specified\n" + __doc__)
    else:
        for i in range(len(p1)):
            if not isinstance(p1[i], int):
                raise ValueError(
                    "Particle types in p1 and p2 have to be of type int: " + str(p1[i]))

        for i in range(len(p2)):
            if not isinstance(p2[i], int):
                raise ValueError(
                    "Particle types in p1 and p2 have to be of type int" + str(p2[i]))

        set1 = create_IntList_from_python_object(p1)
        set2 = create_IntList_from_python_object(p2)

        result = c_analyze.mindist(set1, set2)

        realloc_intlist(set1, 0)
        realloc_intlist(set2, 0)

# The following lines are probably not necessary.
#  free (set1)
#  free (set2)

    return result

#
# Distance to particle or point
#


def distto(system=None, id=None, pos=None):
    """"Distance to particle or point
        distto(id=None,pos=None)
        id: id of the particle
        pos: position as array
    """

    if id == None and pos == None:
        raise Exception("Either id or pos have to be specified\n" + __doc__)

    if id != None and pos != None:
        raise Exception("Only one of id or pos may be specified\n" + __doc__)

    cdef double cpos[3]
    if system.n_part == 0:
        raise Exception("no particles")

    # Get position
    # If particle id specified
    if id != None:
        if not isinstance(id, int):
            raise ValueError("Id has to be an integer")
        _pos = system.part[id].pos
        for i in range(3):
            cpos[i] = _pos[i]
    else:
        for i in range(3):
            cpos[i] = pos[i]
            _id = -1
    return c_analyze.distto(cpos, _id)

#
# Analyze Linear Momentum
#


def analyze_linear_momentum(system=None, include_particles=True, include_lbfluid=True):
    """"Analyze the system's linear momentum
        analyze_linear_momentum(system, include_particles=True, include_lbfluid=True)
        one can calculate the linear momentum of particles and the lbfluid with this command
    """
    return c_analyze.calc_linear_momentum(include_particles, include_lbfluid)

# get all particles in neighborhood r_catch of pos and return their ids
# in il. plane can be used to specify the distance in the xy, xz or yz
# plane


def nbhood(system=None, pos=None, r_catch=None, plane='3d'):
    """nbhood(pos=None, r_catch=None, plane = '3d'):
       get all particles in neighborhood r_catch of pos and return their ids
       in il. plane can be used to specify the distance in the xy, xz or yz
       plane
       nbhood(pos=None, r_catch=None, plane = '3d'):
    """

    cdef int planedims[3]
    cdef IntList * il = NULL
    cdef double c_pos[3]

    checkTypeOrExcept(
        pos, 3, float, "_pos=(float,float,float) must be passed to nbhood")
    checkTypeOrExcept(
        r_catch, 1, float, "r_catch=float needs to be passed to nbhood")

    # default 3d takes into account dist in x, y and z
    planedims[0] = 1
    planedims[1] = 1
    planedims[2] = 1
    if plane == 'xy':
        planedims[2] = 0
    elif plane == 'xz':
        planedims[1] = 0
    elif plane == 'yz':
        planedims[0] = 0
    elif plane != '3d':
        raise Exception(
            'Invalid argument for specifying plane, must be xy, xz, or yz plane')

    for i in range(3):
        c_pos[i] = pos[i]

    il = <IntList * > malloc(sizeof(IntList))
    c_analyze.nbhood(c_pos, r_catch, il, planedims)

    result = create_nparray_from_IntList(il)
    realloc_intlist(il, 0)
    free(il)
    return result


#
# Pressure analysis
#
def pressure(system=None, pressure_type='all', id1='default', id2='default', v_comp=False):
    """
       pressure(pressure_type = 'all', id1 = 'default', id2 = 'default', v_comp=False)
    """
    cdef vector[string] pressure_labels
    cdef vector[double] pressures

    checkTypeOrExcept(v_comp, 1, bool, "v_comp must be a boolean")

    if pressure_type == 'all':
        c_analyze.analyze_pressure_all(pressure_labels, pressures, v_comp)
        return dict(zip(pressure_labels, pressures))
    elif id1 == 'default' and id2 == 'default':
        pressure = c_analyze.analyze_pressure(pressure_type, v_comp)
        return pressure
    elif id1 != 'default' and id2 == 'default':
        checkTypeOrExcept(id1, 1, int, "id1 must be an int")
        pressure = c_analyze.analyze_pressure_single(
            pressure_type, id1, v_comp)
        return pressure
    else:
        checkTypeOrExcept(id1, 1, int, "id1 must be an int")
        checkTypeOrExcept(id2, 1, int, "id2 must be an int")
        pressure = c_analyze.analyze_pressure_pair(
            pressure_type, id1, id2, v_comp)
        return pressure


def stress_tensor(system=None, stress_type='all', id1='default', id2='default', v_comp=False):
    """stress_tensor(system, stress_type = 'all', id1 = 'default', id2 = 'default', v_comp=False)"""
    cdef vector[string] stress_labels
    cdef vector[double] stresses

    checkTypeOrExcept(v_comp, 1, bool, "v_comp must be a boolean")

    if stress_type == 'all':
        c_analyze.analyze_stress_tensor_all(stress_labels, stresses, v_comp)
        return dict(zip(stress_labels, stresses))
    elif id1 == 'default' and id2 == 'default':
        if (c_analyze.analyze_stress_tensor(stress_type, v_comp, stresses)):
            raise Exception("Error while calculating stress tensor")
        return stresses
    elif id1 != 'default' and id2 == 'default':
        checkTypeOrExcept(id1, 1, int, "id1 must be an int")
        if (c_analyze.analyze_stress_single(stress_type, id1, v_comp, stresses)):
            raise Exception("Error while calculating stress tensor")
        return stresses
    else:
        checkTypeOrExcept(id1, 1, int, "id1 must be an int")
        checkTypeOrExcept(id2, 1, int, "id2 must be an int")
        if (c_analyze.analyze_stress_pair(stress_type, id1, id2, v_comp, stresses)):
            raise Exception("Error while calculating stress tensor")
        return stresses


def local_stress_tensor(system=None, periodicity=(1, 1, 1), range_start=(0.0, 0.0, 0.0), stress_range=(1.0, 1.0, 1.0), bins=(1, 1, 1)):
    """local_stress_tensor(periodicity=(1, 1, 1), range_start=(0.0, 0.0, 0.0), stress_range=(1.0, 1.0, 1.0), bins=(1, 1, 1))
    """

    cdef DoubleList * local_stress_tensor = NULL
    cdef int[3] c_periodicity, c_bins
    cdef double[3] c_range_start, c_stress_range

    for i in range(3):
        c_bins[i] = bins[i]
        c_periodicity[i] = periodicity[i]
        c_range_start[i] = range_start[i]
        c_stress_range[i] = stress_range[i]

    if c_analyze.analyze_local_stress_tensor(c_periodicity, c_range_start, c_stress_range, c_bins, local_stress_tensor):
        raise Exception("Error while calculating local stress tensor")
    stress_tensor = create_nparray_from_DoubleList(local_stress_tensor)
    free(local_stress_tensor)
    return stress_tensor


def stress_tensor(self, v_comp=0):
    """stress_tensor(v_comp=0)
    """
    checkTypeOrExcept(v_comp, 1, int, "v_comp must be a boolean")

    # Dict to store the results
    p = {}

    # Update in espresso core if necessary
    if (c_analyze.total_p_tensor.init_status != 1 + v_comp):
        c_analyze.update_pressure(v_comp)
#
    # Individual components of the pressure

    # Total pressure
    cdef int i
    cdef double tmp
    tmp = 0
    for i in range(c_analyze.total_p_tensor.data.n):
        tmp += c_analyze.total_p_tensor.data.e[i]

    p["total"] = tmp

    # Ideal
    p["ideal"] = create_nparray_from_double_array(
        c_analyze.total_p_tensor.data.e, 9)

    # Nonbonded
    total_bonded = np.zeros((3, 3))
    for i in range(c_analyze.n_bonded_ia):
        if (bonded_ia_params[i].type != 0):
            p["bonded", i] = np.reshape(create_nparray_from_double_array(c_analyze.obsstat_bonded( & c_analyze.total_p_tensor, i), 9), (3, 3))
            total_bonded += p["bonded", i]
    p["bonded"] = total_bonded

    # Non-Bonded interactions, total as well as intra and inter molecular
    cdef int j
    total_non_bonded = np.zeros((3, 3))
    total_non_bonded_intra = np.zeros((3, 3))
    total_non_bonded_inter = np.zeros((3, 3))

    for i in range(c_analyze.n_particle_types):
        for j in range(c_analyze.n_particle_types):
            #      if checkIfParticlesInteract(i, j):

            p["nonBonded", i, j] = np.reshape(create_nparray_from_double_array(c_analyze.obsstat_nonbonded( & c_analyze.total_p_tensor, i, j), 9), (3, 3))
            total_non_bonded += p["nonBonded", i, j]

            p["nonBondedIntra", i, j] = np.reshape(create_nparray_from_double_array(c_analyze.obsstat_nonbonded_intra( & c_analyze.total_p_tensor_non_bonded, i, j), 9), (3, 3))
            total_non_bonded_intra += p["nonBondedIntra", i, j]

            p["nonBondedInter", i, j] = np.reshape(create_nparray_from_double_array(c_analyze.obsstat_nonbonded_inter( & c_analyze.total_p_tensor_non_bonded, i, j), 9), (3, 3))
            total_non_bonded_inter += p["nonBondedInter", i, j]

    p["nonBondedIntra"] = total_non_bonded_intra
    p["nonBondedInter"] = total_non_bonded_inter
    p["nonBonded"] = total_non_bonded

    # Electrostatics
    IF ELECTROSTATICS == 1:
        total_coulomb = np.zeros(9)
        for i in range(c_analyze.total_p_tensor.n_coulomb):
            p["coulomb", i] = np.reshape(
                create_nparray_from_double_array(c_analyze.total_p_tensor.coulomb, 9), (3, 3))
            total_coulomb = p["coulomb", i]
        p["coulomb"] = total_coulomb

    # Dipoles
    IF DIPOLES == 1:
        total_dipolar = np.zeros(9)
        for i in range(c_analyze.total_p_tensor.n_dipolar):
            p["dipolar", i] = np.reshape(
                create_nparray_from_double_array(c_analyze.total_p_tensor.dipolar, 9), (3, 3))
            total_dipolar = p["dipolar", i]
        p["dipolar"] = total_dipolar

    # virtual sites
    IF VIRTUAL_SITES_RELATIVE == 1:
        p["vs_relative"] = np.reshape(create_nparray_from_double_array(
            c_analyze.total_p_tensor.vs_relative, 9), (3, 3))

    return p

#
# Energy analysis
#


def energy(system=None, etype='all', id1='default', id2='default'):
    """energy(system, etype = 'all', id1 = 'default', id2 = 'default')"""
    if system.n_part == 0:
        raise Exception('no particles')

    if c_analyze.total_energy.init_status == 0:
        c_analyze.init_energies( & c_analyze.total_energy)
        c_analyze.master_energy_calc()
    _value = 0.0

    if etype == 'all':
        _result = energy(system, 'total') + ' ' + energy(system, 'kinetic')
        _result += energy(system, 'nonbonded', 0, 0)
        # todo: check for existing particle and bond types
        # and add those to _result
        return _result

    if etype == 'total':
        if id1 != 'default' or id2 != 'default':
            print ('warning: energy(\'total\') does not need '
                   'further arguments, ignored.')
        for i in range(c_analyze.total_energy.data.n):
            _value += c_analyze.total_energy.data.e[i]
        return '{ energy: %f }' % _value

    if etype == 'kinetic':
        if id1 != 'default' or id2 != 'default':
            print ('warning: energy(\'kinetic\') does not need '
                   'further arguments, ignored.')
        _value = c_analyze.total_energy.data.e[0]
        return '{ kinetic: %f }' % _value

    # coulomb interaction
    if etype == 'coulomb':
        if(code_info.electrostatics_defined()):
            for i in range(c_analyze.total_energy.n_coulomb):
                _value += c_analyze.total_energy.coulomb[i]
            return '{ coulomb: %f }' % _value
        else:
            print 'error: ELECTROSTATICS not compiled'
            return 'error: ELECTROSTATICS not compiled'

    if etype == 'magnetic':
        if(code_info.dipoles_defined()):
            for i in range(c_analyze.total_energy.n_dipolar):
                _value += c_analyze.total_energy.dipolar[i]
            return '{ magnetic: %f }' % _value
        else:
            print 'error: DIPOLES not compiled'
            return 'error: DIPOLES not compiled'

    # bonded interactions
    if etype == 'bonded':
        if not isinstance(id1, int):
            print ('error: analyze.energy(\'bonded\',<bondid>): '
                   '<bondid> must be integer')
            raise TypeError('analyze.energy(\'bonded\',<bondid>): '
                            '<bondid> must be integer')
        else:
            # todo: check if bond type id1 exist
            _value = c_analyze.obsstat_bonded( & c_analyze.total_energy, id1)[0]
            return '{ %d bonded: %f }' % (id1, _value)

    # nonbonded interactions
    if etype == 'nonbonded':
        if not isinstance(id1, int):
            print ('error: analyze.energy(\'bonded\',<bondid>): '
                   '<bondid> must be integer')
            raise TypeError('analyze.energy(\'bonded\',<bondid>): '
                            '<bondid> must be integer')
        if not isinstance(id2, int):
            print ('error: analyze.energy(\'bonded\',<bondid>): '
                   '<bondid> must be integer')
            raise TypeError('analyze.energy(\'bonded\',<bondid>): '
                            '<bondid> must be integer')
        else:
            # todo: check if particle types id1 and id2 exist
            _value = c_analyze.obsstat_nonbonded( & c_analyze.total_energy, id1, id2)[0]
            return '{ %d %d nonbonded: %f }' % (id1, id2, _value)

    return 'error: unknown feature of analyze energy: \'%s\'' % etype


def calc_re(system=None, chain_start=None, number_of_chains=None, chain_length=None):
    cdef double * re = NULL
    check_topology(system, chain_start, number_of_chains, chain_length)
    c_analyze.calc_re( & re)
    tuple_re = (re[0], re[1], re[2])
    free(re)
    return tuple_re


def calc_rg(system=None, chain_start=None, number_of_chains=None, chain_length=None):
    cdef double * rg = NULL
    check_topology(system, chain_start, number_of_chains, chain_length)
    c_analyze.calc_rg( & rg)
    tuple_rg = (rg[0], rg[1], rg[2])
    free(rg)
    return tuple_rg


def calc_rh(system=None, chain_start=None, number_of_chains=None, chain_length=None):
    cdef double * rh = NULL
    check_topology(system, chain_start, number_of_chains, chain_length)
    c_analyze.calc_rh( & rh)
    tuple_rh = (rh[0], rh[1], rh[2])
    free(rh)
    return tuple_rh


def check_topology(system=None, chain_start=None, number_of_chains=None, chain_length=None):
    checkTypeOrExcept(
        chain_start, 1, int, "chain_start=int is a required argument")
    checkTypeOrExcept(
        number_of_chains, 1, int, "number_of_chains=int is a required argument")
    checkTypeOrExcept(
        chain_length, 1, int, "chain_length=int is a required argument")
    if not system:
        raise ValueError('Must pass an instance of ESPResSo to thi function')
    if chain_start < 0:
        raise ValueError('chain_start must be greater than zero')
    if chain_length < 0:
        raise ValueError('chain_length must be greater than zero')
    if number_of_chains < 0:
        raise ValueError('number_of_chains must be greater than zero')
    c_analyze.sortPartCfg()
    if chain_start + chain_length * number_of_chains >= system.n_part:
        raise ValueError(
            'start+number_of_chains*chain_length cannot be greater than the total number of particles.')
    c_analyze.chain_start = chain_start
    c_analyze.chain_n_chains = number_of_chains
    c_analyze.chain_length = chain_length

#
# Structure factor
#


def structure_factor(system=None, sf_type='default', sf_order='default'):
    """Structure Factor
       structure_factor(system = None, sf_type = 'default', sf_order = 'default' )
    """
    cdef double * sf

    checkTypeOrExcept(sf_type, 1, int, "sf_type must be an int")
    checkTypeOrExcept(sf_order, 1, int, "sf_order must be an int")

    # Used to take the WITHOUT_BONDS define
    c_analyze.updatePartCfg(0)
    c_analyze.calc_structurefactor(sf_type, sf_order, & sf)

    return c_analyze.modify_stucturefactor(sf_order, sf)
