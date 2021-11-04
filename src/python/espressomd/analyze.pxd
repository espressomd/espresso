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

# For C-extern Analysis

cimport numpy as np
import numpy as np
from .utils cimport Vector3i, Vector3d, Vector9d, Span
from .utils cimport create_nparray_from_double_span
from libcpp.vector cimport vector  # import std::vector as vector
from libcpp cimport bool as cbool
from libcpp.memory cimport shared_ptr
from .interactions cimport bonded_ia_params_zero_based_type
from .interactions cimport enum_bonded_interaction
from .interactions cimport bonded_ia_params_next_key
include "myconfig.pxi"

cdef extern from "<array>" namespace "std" nogil:
    cdef cppclass array4 "std::array<double, 4>":
        array4() except+
        double & operator[](size_t)

    cdef cppclass array2 "std::array<double, 2>":
        array2() except+
        double & operator[](size_t)

cdef extern from "PartCfg.hpp":
    cppclass PartCfg:
        pass

cdef extern from "partCfg_global.hpp":
    PartCfg & partCfg()

cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
    int max_seen_particle_type

cdef extern from "Observable_stat.hpp":
    cdef cppclass Observable_stat:
        Span[double] kinetic
        Span[double] coulomb
        Span[double] dipolar
        Span[double] virtual_sites
        Span[double] external_fields
        double accumulate(...)
        double accumulate2 "accumulate"(double acc, size_t column)
        double accumulate1 "accumulate"(double acc)
        Span[double] bonded_contribution(int bond_id)
        Span[double] non_bonded_intra_contribution(int type1, int type2)
        Span[double] non_bonded_inter_contribution(int type1, int type2)
        size_t get_chunk_size()

cdef extern from "statistics.hpp":
    cdef void calc_structurefactor(PartCfg & , const vector[int] & p_types, int order, vector[double] & wavevectors, vector[double] & intensities) except +
    cdef double mindist(PartCfg & , const vector[int] & set1, const vector[int] & set2)
    cdef vector[int] nbhood(PartCfg & , const Vector3d & pos, double r_catch, const Vector3i & planedims)
    cdef vector[double] calc_linear_momentum(int include_particles, int include_lbfluid)
    cdef vector[double] centerofmass(PartCfg & , int part_type)

    Vector3d angularmomentum(PartCfg & , int p_type)

    void momentofinertiamatrix(PartCfg & , int p_type, double * MofImatrix)

    void calc_part_distribution(
        PartCfg &, const vector[int] & p1_types, const vector[int] & p2_types,
        double r_min, double r_max, int r_bins, bint log_flag, double * low,
        double * dist)

cdef extern from "statistics_chain.hpp":
    array4 calc_re(int, int, int)
    array4 calc_rg(int, int, int) except +
    array2 calc_rh(int, int, int)

cdef extern from "pressure.hpp":
    cdef shared_ptr[Observable_stat] calculate_pressure()

cdef extern from "energy.hpp":
    cdef shared_ptr[Observable_stat] calculate_energy()
    double calculate_current_potential_energy_of_system()

cdef extern from "dpd.hpp":
    Vector9d dpd_stress()

cdef inline get_obs_contribs(Span[double] contributions,
                             cbool calc_scalar_pressure):
    """
    Convert an Observable_stat range of contributions into a correctly
    shaped numpy array.

    Parameters
    ----------
    contributions : (N,) array_like of :obj:`float`
        Flattened array of energy/pressure contributions from an observable.
    size : :obj:`int`, \{1, 9\}
        Dimensionality of the data.
    calc_scalar_pressure : :obj:`bool`
        Whether to calculate a scalar pressure (only relevant when
        ``contributions`` is a pressure tensor).

    """
    cdef np.ndarray value
    value = create_nparray_from_double_span(contributions)
    if contributions.size() == 9:
        if calc_scalar_pressure:
            return np.einsum('...ii', value.reshape((-1, 3, 3))) / 3
        else:
            return value.reshape((-1, 3, 3))
    else:
        return value

cdef inline get_obs_contrib(Span[double] contribution,
                            cbool calc_scalar_pressure):
    """
    Convert an Observable_stat contribution into a correctly
    shaped numpy array. If the size is 1, decay to a float.

    Parameters
    ----------
    contributions : (N,) array_like of :obj:`float`
        Flattened array of energy/pressure contributions from an observable.
    calc_scalar_pressure : :obj:`bool`
        Whether to calculate a scalar pressure (only relevant when
        ``contributions`` is a pressure tensor).

    """
    cdef np.ndarray value
    value = get_obs_contribs(contribution, calc_scalar_pressure)
    if value.shape[0] == 1:
        return value[0]
    return value

cdef inline observable_stat_matrix(size_t size, cbool calc_sp):
    if size == 9 and not calc_sp:
        return np.zeros((3, 3), dtype=float)
    else:
        return 0.0

cdef inline Observable_stat_to_dict(Observable_stat * obs, cbool calc_sp):
    """Transform an ``Observable_stat`` object to a python dict.

    Parameters
    ----------
    obs :
        Core observable.
    calc_sp : :obj:`bool`
        Whether to calculate a scalar pressure (only relevant when
        ``obs`` is a pressure tensor observable).

    Returns
    -------
    :obj:`dict`
        A dictionary with the following keys:

        * ``"total"``: total contribution
        * ``"kinetic"``: kinetic contribution
        * ``"bonded"``: total bonded contribution
        * ``"bonded", <bond_type>``: bonded contribution which arises from the given bond_type
        * ``"non_bonded"``: total non-bonded contribution
        * ``"non_bonded", <type_i>, <type_j>``: non-bonded contribution which arises from the interactions between type_i and type_j
        * ``"non_bonded_intra", <type_i>, <type_j>``: non-bonded contribution between short ranged forces between type i and j and with the same mol_id
        * ``"non_bonded_inter", <type_i>, <type_j>``: non-bonded contribution between short ranged forces between type i and j and different mol_ids
        * ``"coulomb"``: Coulomb contribution, how it is calculated depends on the method
        * ``"coulomb", <i>``: Coulomb contribution from particle pairs (``i=0``), electrostatics solvers (``i=1``)
        * ``"dipolar"``: dipolar contribution, how it is calculated depends on the method
        * ``"dipolar", <i>``: dipolar contribution from particle pairs and magnetic field constraints (``i=0``), magnetostatics solvers (``i=1``)
        * ``"virtual_sites"``: virtual sites contribution
        * ``"virtual_sites", <i>``: contribution from virtual site i
        * ``"external_fields"``: external fields contribution

    """

    cdef size_t i
    cdef size_t j
    cdef size_t obs_dim = obs.get_chunk_size()
    cdef size_t n_bonded = bonded_ia_params_next_key()
    cdef size_t n_nonbonded = <size_t > max_seen_particle_type
    cdef double[9] total

    # Dict to store the results
    p = {}

    # Total contribution

    for i in range(obs_dim):
        total[i] = obs.accumulate(0.0, i)
    p["total"] = get_obs_contrib(Span[double](total, obs_dim), calc_sp)

    # Kinetic
    p["kinetic"] = get_obs_contrib(obs.kinetic, calc_sp)

    # External
    p["external_fields"] = get_obs_contrib(obs.external_fields, calc_sp)

    # Bonded
    total_bonded = observable_stat_matrix(obs_dim, calc_sp)
    for i in range(n_bonded):
        if bonded_ia_params_zero_based_type(
                i) != enum_bonded_interaction.BONDED_IA_NONE:
            val = get_obs_contrib(obs.bonded_contribution(i), calc_sp)
            p["bonded", i] = val
            total_bonded += val
    p["bonded"] = total_bonded

    # Non-Bonded interactions, total as well as intra and inter molecular
    total_intra = observable_stat_matrix(obs_dim, calc_sp)
    total_inter = observable_stat_matrix(obs_dim, calc_sp)
    for i in range(n_nonbonded):
        for j in range(i, n_nonbonded):
            intra = get_obs_contrib(
                obs.non_bonded_intra_contribution(i, j), calc_sp)
            total_intra += intra
            p["non_bonded_intra", i, j] = intra
            inter = get_obs_contrib(
                obs.non_bonded_inter_contribution(i, j), calc_sp)
            total_inter += inter
            p["non_bonded_inter", i, j] = inter
            p["non_bonded", i, j] = intra + inter

    p["non_bonded_intra"] = total_intra
    p["non_bonded_inter"] = total_inter
    p["non_bonded"] = total_intra + total_inter

    # Electrostatics
    IF ELECTROSTATICS == 1:
        cdef np.ndarray coulomb
        coulomb = get_obs_contribs(obs.coulomb, calc_sp)
        for i in range(coulomb.shape[0]):
            p["coulomb", i] = coulomb[i]
        p["coulomb"] = np.sum(coulomb, axis=0)

    # Dipoles
    IF DIPOLES == 1:
        cdef np.ndarray dipolar
        dipolar = get_obs_contribs(obs.dipolar, calc_sp)
        for i in range(dipolar.shape[0]):
            p["dipolar", i] = dipolar[i]
        p["dipolar"] = np.sum(dipolar, axis=0)

    # virtual sites
    IF VIRTUAL_SITES == 1:
        cdef np.ndarray virtual_sites
        virtual_sites = get_obs_contribs(obs.virtual_sites, calc_sp)
        for i in range(virtual_sites.shape[0]):
            p["virtual_sites", i] = virtual_sites[i]
        p["virtual_sites"] = np.sum(virtual_sites, axis=0)

    return p

cdef inline get_scalar_pressure():
    return Observable_stat_to_dict(calculate_pressure().get(), True)

cdef inline get_pressure_tensor():
    return Observable_stat_to_dict(calculate_pressure().get(), False)

cdef inline get_energy():
    return Observable_stat_to_dict(calculate_energy().get(), False)
