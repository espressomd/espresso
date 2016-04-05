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
from libcpp.map cimport map  # import std::map as map
from interactions import *
from interactions cimport *
import numpy as np
cimport numpy as np
from globals cimport n_configs, min_box_l


#
# Minimal distance between particles
#


def mindist(system=None, p1='default', p2='default'):
    """Minimal distance between particles
      mindist(p1="default",p2="default")

      p1, p2: lists of particle types
    """

    cdef int_list * set1
    cdef int_list * set2

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

        set1 = create_int_list_from_python_object(p1)
        set2 = create_int_list_from_python_object(p2)

        result = c_analyze.mindist(set1, set2)

        realloc_intlist(set1, 0)
        realloc_intlist(set2, 0)

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


#
# Analyze center of mass
#


def centermass(system=None, part_type=None):
    """Analyze the system's center of mass
    """
    return c_analyze.centerofmass(part_type)


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
    cdef int_list * il = NULL
    cdef double c_pos[3]

    check_type_or_throw_except(
        pos, 3, float, "_pos=(float,float,float) must be passed to nbhood")
    check_type_or_throw_except(
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

    il = <int_list * > malloc(sizeof(int_list))
    c_analyze.nbhood(c_pos, r_catch, il, planedims)

    result = create_nparray_from_int_list(il)
    realloc_intlist(il, 0)
    free(il)
    return result


def cylindrical_average(system=None, center=None, direction=None,
                        length=None, radius=None,
                        bins_axial=None, bins_radial=None,
                        types=[-1]):

    # Check the input types
    check_type_or_throw_except(
        center,      3, float, "center has to be 3 floats")
    check_type_or_throw_except(
        direction,   3, float, "direction has to be 3 floats")
    check_type_or_throw_except(
        length,      1, float, "length has to be a float")
    check_type_or_throw_except(
        radius,      1, float, "radius has to be a floats")
    check_type_or_throw_except(
        bins_axial,  1, int,   "bins_axial has to be an int")
    check_type_or_throw_except(
        bins_radial, 1, int,   "bins_radial has to be an int")

    # Convert Python types to C++ types
    cdef vector[double] c_center = center
    cdef vector[double] c_direction = direction
    cdef double c_length = length
    cdef double c_radius = radius
    cdef int c_bins_axial = bins_axial
    cdef int c_bins_radial = bins_radial
    cdef vector[int] c_types = types

    cdef map[string, vector[vector[vector[double]]]] distribution
    c_analyze.calc_cylindrical_average(c_center, c_direction, c_length,
                                       c_radius, c_bins_axial, c_bins_radial, c_types,
                                       distribution)

    cdef double binwd_axial = c_length / c_bins_axial
    cdef double binwd_radial = c_radius / c_bins_radial
    cdef double binvolume, pos_radial, pos_axial

    cdef vector[string] names = ["density", "v_r", "v_t"]

    buffer = np.empty(
        [bins_radial * bins_axial, 5 + c_types.size() * names.size()])

    cdef int index_radial, index_axial
    cdef unsigned int type_id, name
    for index_radial in range(0, bins_radial):
        for index_axial in range(0, bins_axial):
            pos_radial = (index_radial + .5) * binwd_radial
            pos_axial = (index_axial + .5) * binwd_axial - .5 * c_length

            if (index_radial == 0):
                binvolume = np.pi * binwd_radial * binwd_radial * c_length
            else:
                binvolume = np.pi * \
                    (index_radial * index_radial + 2 * index_radial) * \
                    binwd_radial * binwd_radial * c_length

            buffer[index_axial + bins_axial * index_radial, 0] = index_radial
            buffer[index_axial + bins_axial * index_radial, 1] = index_axial
            buffer[index_axial + bins_axial * index_radial, 2] = pos_radial
            buffer[index_axial + bins_axial * index_radial, 3] = pos_axial
            buffer[index_axial + bins_axial * index_radial, 4] = binvolume
            for type_id in range(0, c_types.size()):
                for name in range(0, names.size()):
                    buffer[index_axial + bins_axial * index_radial, 5 + name + type_id * names.size()] \
                        = distribution[names[name]][type_id][index_radial][index_axial]

    return buffer

#
# Pressure analysis
#


def pressure(self, v_comp=0):
    """Pressure calculation
       pressure(v_comp=False)
       """
    check_type_or_throw_except(v_comp, 1, int, "v_comp must be a boolean")
#
    # Dict to store the results
    p = {}

    # Update in espresso core if necessary
    if (c_analyze.total_pressure.init_status != 1 + v_comp):
        c_analyze.update_pressure(v_comp)
#
    # Individual components of the pressure

    # Total pressure
    cdef int i
    cdef double tmp
    tmp = 0
    for i in range(c_analyze.total_pressure.data.n):
        tmp += c_analyze.total_pressure.data.e[i]

    p["total"] = tmp

    # Ideal
    p["ideal"] = c_analyze.total_pressure.data.e[0]

    # Nonbonded
    cdef double total_bonded
    total_bonded = 0
    for i in range(c_analyze.n_bonded_ia):
        if (bonded_ia_params[i].type != 0):
            p["bonded", i] = c_analyze.obsstat_bonded(& c_analyze.total_pressure, i)[0]
            total_bonded += c_analyze.obsstat_bonded( & c_analyze.total_pressure, i)[0]
    p["bonded"] = total_bonded

    # Non-Bonded interactions, total as well as intra and inter molecular
    cdef int j
    cdef double total_intra
    cdef double total_inter
    cdef double total_non_bonded
    total_inter = 0
    total_intra = 0
    total_non_bonded = 0

    for i in range(c_analyze.n_particle_types):
        for j in range(c_analyze.n_particle_types):
            #      if checkIfParticlesInteract(i, j):
            p["nonBonded", i, j] = c_analyze.obsstat_nonbonded(& c_analyze.total_pressure, i, j)[0]
            total_non_bonded = c_analyze.obsstat_nonbonded(& c_analyze.total_pressure, i, j)[0]
            total_intra += c_analyze.obsstat_nonbonded_intra(& c_analyze.total_pressure_non_bonded, i, j)[0]
            p["nonBondedIntra", i, j] = c_analyze.obsstat_nonbonded_intra(& c_analyze.total_pressure_non_bonded, i, j)[0]
            p["nonBondedInter", i, j] = c_analyze.obsstat_nonbonded_inter(& c_analyze.total_pressure_non_bonded, i, j)[0]
            total_inter += c_analyze.obsstat_nonbonded_inter(& c_analyze.total_pressure_non_bonded, i, j)[0]
    p["nonBondedIntra"] = total_intra
    p["nonBondedInter"] = total_inter
    p["nonBondedInter"] = total_inter
    p["nonBonded"] = total_non_bonded

    # Electrostatics
    IF ELECTROSTATICS == 1:
        cdef double total_coulomb
        total_coulomb = 0
        for i in range(c_analyze.total_pressure.n_coulomb):
            total_coulomb += c_analyze.total_pressure.coulomb[i]
            p["coulomb", i] = c_analyze.total_pressure.coulomb[i]
        p["coulomb"] = total_coulomb

    # Dipoles
    IF DIPOLES == 1:
        cdef double total_dipolar
        total_dipolar = 0
        for i in range(c_analyze.total_pressure.n_dipolar):
            total_dipolar += c_analyze.total_pressure.dipolar[i]
            p["dipolar", i] = c_analyze.total_pressure.coulomb[i]
        p["dipolar"] = total_dipolar

    # virtual sites
    IF VIRTUAL_SITES_RELATIVE == 1:
        p["vs_relative"] = c_analyze.total_pressure.vs_relative[0]

    return p


def stress_tensor(self, v_comp=0):
    """stress_tensor(v_comp=0)
    """
    check_type_or_throw_except(v_comp, 1, int, "v_comp must be a boolean")

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


def local_stress_tensor(system=None, periodicity=(1, 1, 1), range_start=(0.0, 0.0, 0.0), stress_range=(1.0, 1.0, 1.0), bins=(1, 1, 1)):
    """local_stress_tensor(periodicity=(1, 1, 1), range_start=(0.0, 0.0, 0.0), stress_range=(1.0, 1.0, 1.0), bins=(1, 1, 1))
    """

    cdef double_list * local_stress_tensor = NULL
    cdef int[3] c_periodicity, c_bins
    cdef double[3] c_range_start, c_stress_range

    for i in range(3):
        c_bins[i] = bins[i]
        c_periodicity[i] = periodicity[i]
        c_range_start[i] = range_start[i]
        c_stress_range[i] = stress_range[i]

    if c_analyze.analyze_local_stress_tensor(c_periodicity, c_range_start, c_stress_range, c_bins, local_stress_tensor):
        raise Exception("Error while calculating local stress tensor")
    stress_tensor = create_nparray_from_double_list(local_stress_tensor)
    free(local_stress_tensor)
    return stress_tensor

#
# Energy analysis
#


def energy(system, etype='all', id1='default', id2='default'):
    """energy()
    """
#  if system.n_part == 0:
#    raise Exception('no particles')

    e = {}

    if c_analyze.total_energy.init_status == 0:
        c_analyze.init_energies( & c_analyze.total_energy)
        c_analyze.master_energy_calc()

    # Individual components of the pressur

    # Total energy
    cdef int i
    cdef double tmp
    tmp = 0
    for i in range(c_analyze.total_energy.data.n):
        tmp += c_analyze.total_energy.data.e[i]

    e["total"] = tmp

    # Ideal
    e["ideal"] = c_analyze.total_energy.data.e[0]

    # Nonbonded
    cdef double total_bonded
    total_bonded = 0
    for i in range(c_analyze.n_bonded_ia):
        if (bonded_ia_params[i].type != 0):
            e["bonded", i] = c_analyze.obsstat_bonded(& c_analyze.total_energy, i)[0]
            total_bonded += c_analyze.obsstat_bonded( & c_analyze.total_energy, i)[0]
    e["bonded"] = total_bonded

    # Non-Bonded interactions, total as well as intra and inter molecular
    cdef int j
    cdef double total_intra
    cdef double total_inter
    cdef double total_non_bonded
    total_inter = 0
    total_intra = 0
    total_non_bonded = 0

    for i in range(c_analyze.n_particle_types):
        for j in range(c_analyze.n_particle_types):
            #      if checkIfParticlesInteract(i, j):
            e["nonBonded", i, j] = c_analyze.obsstat_nonbonded(& c_analyze.total_energy, i, j)[0]
            total_non_bonded = c_analyze.obsstat_nonbonded(& c_analyze.total_energy, i, j)[0]
#        total_intra +=c_analyze.obsstat_nonbonded_intra(&c_analyze.total_energy_non_bonded, i, j)[0]
#        e["nonBondedIntra",i,j] =c_analyze.obsstat_nonbonded_intra(&c_analyze.total_energy_non_bonded, i, j)[0]
#        e["nonBondedInter",i,j] =c_analyze.obsstat_nonbonded_inter(&c_analyze.total_energy_non_bonded, i, j)[0]
#        total_inter+= c_analyze.obsstat_nonbonded_inter(&c_analyze.total_energy_non_bonded, i, j)[0]
#  e["nonBondedIntra"]=total_intra
#  e["nonBondedInter"]=total_inter
    e["nonBonded"] = total_non_bonded

    # Electrostatics
    IF ELECTROSTATICS == 1:
        cdef double total_coulomb
        total_coulomb = 0
        for i in range(c_analyze.total_energy.n_coulomb):
            total_coulomb += c_analyze.total_energy.coulomb[i]
            e["coulomb", i] = c_analyze.total_energy.coulomb[i]
        e["coulomb"] = total_coulomb

    # Dipoles
    IF DIPOLES == 1:
        cdef double total_dipolar
        total_dipolar = 0
        for i in range(c_analyze.total_energy.n_dipolar):
            total_dipolar += c_analyze.total_energy.dipolar[i]
            e["dipolar", i] = c_analyze.total_energy.coulomb[i]
        e["dipolar"] = total_dipolar

    return e


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
    check_type_or_throw_except(
        chain_start, 1, int, "chain_start=int is a required argument")
    check_type_or_throw_except(
        number_of_chains, 1, int, "number_of_chains=int is a required argument")
    check_type_or_throw_except(
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

    check_type_or_throw_except(sf_type, 1, int, "sf_type must be an int")
    check_type_or_throw_except(sf_order, 1, int, "sf_order must be an int")

    # Used to take the WITHOUT_BONDS define
    c_analyze.updatePartCfg(0)
    c_analyze.calc_structurefactor(sf_type, sf_order, & sf)

    return c_analyze.modify_stucturefactor(sf_order, sf)

#
# RDF
#

def rdf(system=None, rdf_type=None, type_list_a=None, type_list_b=None,
        r_min = 0.0, r_max = None, r_bins = 100, n_conf = None):

    if rdf_type is None:
        raise ValueError("rdf_type must not be empty!")
    if (type_list_a is None) or (not hasattr(type_list_a, '__iter__')):
        raise ValueError("type_list_a has to be a list!")
    if (type_list_b is None) or (not hasattr(type_list_b, '__iter__')):
        raise ValueError("type_list_b has to be a list!")

    if rdf_type != 'rdf':
        if n_configs == 0:
            raise ValueError("No configurations founds!\n",
                             "Use 'analyze append' to save some,",
                             "or 'analyze rdf' to only look at current RDF!""")
        if n_conf is None:
            n_conf = n_configs

    if r_max is None: r_max = min_box_l / 2.0;

    cdef vector[double] rdf
    rdf.resize(r_bins)
    cdef vector[int] p1_types = type_list_a
    cdef vector[int] p2_types = type_list_b

    c_analyze.updatePartCfg(0)
    if rdf_type == 'rdf':
        c_analyze.calc_rdf(p1_types, p2_types, r_min, r_max, r_bins, rdf)
    elif rdf_type == '<rdf>':
        c_analyze.calc_rdf_av(p1_types, p2_types,r_min, r_max, r_bins, rdf, n_conf)
    elif rdf_type == '<rdf-intermol>':
        c_analyze.calc_rdf_intermol_av(p1_types, p2_types, r_min, r_max, r_bins, rdf, n_conf)
    else:
        raise Exception("rdf_type has to be one of 'rdf', '<rdf>', and '<rdf_intermol>'")

    r = np.empty(r_bins)
    bin_width = (r_max - r_min) / r_bins
    rr = r_min + bin_width/2.0
    for i in range(r_bins):
        r[i] = rr
        rr += bin_width

    return np.array([r, rdf])

#
# angularmomentum
#

def angularmomentum(system=None, p_type=None):
    print "p_type = ", p_type
    check_type_or_throw_except(
        p_type, 1, int,   "p_type has to be an int")

    cdef double[3] com
    cdef int p1 = p_type

    c_analyze.angularmomentum(p1, com)

    return np.array([com[0], com[1], com[2]])

#
# gyration_tensor
#

def gyration_tensor(system=None, p_type=None):
    cdef vector[double] gt

    if p_type is not None:
        check_type_or_throw_except(p_type, 1, int, "p_type has to be an int")
        if (p_type < 0 or p_type >= c_analyze.n_particle_types):
            raise ValueError("Particle type",p_type, "does not exist!")
    else:
        p_type = -1

    c_analyze.calc_gyration_tensor(p_type, gt)

    return { "Rg^2" : gt[3],
             "shape" : [ gt[4], gt[5], gt[6] ],
             "eva0" : [ gt[0], [ gt[7], gt[8], gt[9] ] ],
             "eva1" : [ gt[1], [ gt[10], gt[11], gt[12] ] ],
             "eva2" : [ gt[2], [ gt[13], gt[14], gt[15] ] ] }

#
# momentofinertiamatrix
#

def momentofinertiamatrix(system=None, p_type=None):
    cdef double[9] MofImatrix

    if p_type is not None:
        check_type_or_throw_except(p_type, 1, int, "p_type has to be an int")
        if (p_type < 0 or p_type >= c_analyze.n_particle_types):
            raise ValueError("Particle type",p_type, "does not exist!")

        c_analyze.momentofinertiamatrix(p_type, MofImatrix)

        MofImatrix_np = np.empty((9))
        for i in range(9): MofImatrix_np[i] = MofImatrix[i]

        return MofImatrix_np

#
# rdfchain
#

def rdfchain(system=None, r_min=None, r_max=None, r_bins=None,
             chain_start=None, n_chains=None, chain_length=None):
    cdef double *f1
    cdef double *f2
    cdef double *f3

    check_type_or_throw_except(r_min, 1, float, "r_min has to be a float")
    check_type_or_throw_except(r_max, 1, float, "r_max has to be a float")
    check_type_or_throw_except(r_bins, 1, int, "r_bins has to be an int")

    check_topology(system=system, chain_start=chain_start,
                   number_of_chains=number_of_chains, chain_length=chain_length)

    if (c_analyze.chain_n_chains == 0 or chain_length == 0):
        raise Exception("The chain topology has not been set")
    if (r_bins <=0):
        raise Exception("Nothing to be done - choose <r_bins> greater zero!")
    if (r_min <= 0.):
        raise Exception("<r_min> has to be positive")
    if (r_max <= r_min):
        raise Exception("<r_max> has to be larger than <r_min>")

    # Used to take the WITHOUT_BONDS define
    c_analyze.updatePartCfg(0)
    c_analyze.analyze_rdfchain(r_min, r_max, r_bins, &f1, &f2, &f3);

    rdfchain = np.empty((r_bins,4))
    bin_width = (r_max - r_min) / float(r_bins)
    r = r_min + bin_width/2.0;
    for i in range(r_bins):
        rdfchain[i,0] = r
        rdfchain[i,1] = f1[i]
        rdfchain[i,2] = f2[i]
        rdfchain[i,3] = f3[i]
        r += bin_width;

    return rdfchain

#
# Vkappa
#

_Vkappa = {
    "Vk1" : 0.0,
    "Vk2" : 0.0,
    "avk" : 0.0
}
def Vkappa(system=None, mode=None, Vk1=None, Vk2=None, avk=None):
    check_type_or_throw_except(mode, 1, str, "mode has to be a string")

    if ( mode == "reset" ):
        _Vkappa["Vk1"] = 0.0
        _Vkappa["Vk2"] = 0.0
        _Vkappa["avk"] = 0.0
    elif ( mode == "read" ):
        return _Vkappa
    elif ( mode == "set" ):
        check_type_or_throw_except(Vk1, 1, float, "Vk1 has to be a float")
        _Vkappa["Vk1"] = Vk1
        check_type_or_throw_except(Vk2, 1, float, "Vk2 has to be a float")
        _Vkappa["Vk2"] = Vk2
        check_type_or_throw_except(avk, 1, float, "avk has to be a float")
        _Vkappa["avk"] = avk
        if ( _Vkappa["avk"] <= 0.0 ):
            raise Exception("ERROR: # of averages <avk> has to be positive! Resetting values.")
            result = _Vkappa["Vk1"] = _Vkappa["Vk2"] = _Vkappa["avk"] = 0.0
        result = _Vkappa["Vk2"] / _Vkappa["avk"] - (_Vkappa["Vk1"] / _Vkappa["avk"])**2
    else:
        raise Exception("ERROR: Unknown mode.")
    #else: # <- WTF?  Why is there a second `else'?
    #    _Vkappa["Vk1"] += system.box_l[0] * system.box_l[1] * system.box_l[2]
    #    _Vkappa["Vk2"] += (system.box_l[0] * system.box_l[1] * system.box_l[2])**2
    #    _Vkappa["avk"] += 1.0
    #    result = _Vkappa["Vk2"] / _Vkappa["avk"] - (_Vkappa["Vk1"] / _Vkappa["avk"])**2

    return result
