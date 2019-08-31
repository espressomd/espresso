#
# Copyright (C) 2013-2018 The ESPResSo project
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
from . cimport analyze
from . cimport utils
from . cimport particle_data
from . import utils
from . import code_info
from . import particle_data
from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector
from libcpp.map cimport map  # import std::map as map
from .interactions import *
from espressomd.interactions cimport *
import numpy as np
cimport numpy as np
from globals cimport n_configs

from .utils import array_locked

from collections import OrderedDict
from .system import System
from espressomd.utils import is_valid_type


class Analysis:

    def __init__(self, system):
        if not isinstance(system, System):
            raise TypeError("An instance of System is required as argument")
        self._system = system

    #
    # Append configs
    #

    def append(self):
        """Append configuration for averaged analysis."""
        if analyze.n_part == 0:
            raise Exception("No particles to append!")
        if (analyze.n_configs > 0) and (analyze.n_part_conf != analyze.n_part):
            raise Exception(
                "All configurations stored must have the same length")

        analyze.analyze_append(analyze.partCfg())

    #
    # Minimal distance between particles
    #

    def min_dist2(self, p1, p2):
        """Minimal distance between two three-dimensional coordinates p1 and p2.

        Parameters
        ----------
        p1, p2 : arrays of :obj:`float`

        """
        cdef Vector3d p1c
        cdef Vector3d p2c
        for i in range(3):
            p1c[i] = p1[i]
            p2c[i] = p2[i]
        return analyze.min_distance2(p1c, p2c)

    def min_dist(self, p1='default', p2='default'):
        """Minimal distance between two sets of particles.

        Parameters
        ----------
        p1, p2 : lists of :obj:`int`
            Particle :attr:`~espressomd.particle_data.ParticleHandle.type` in
            both sets. If both are set to ``'default'``, the minimum distance
            of all pairs is returned.

        """

        cdef List[int] set1
        cdef List[int] set2

        if p1 == 'default' and p2 == 'default':
            pass
        elif (p1 == 'default' and not p2 == 'default') or \
             (not p1 == 'default' and p2 == 'default'):
            raise Exception("Both, p1 and p2 have to be specified\n" + __doc__)
        else:
            for i in range(len(p1)):
                if not is_valid_type(p1[i], int):
                    raise ValueError(
                        "Particle types in p1 and p2 have to be of type int: " + str(p1[i]))

            for i in range(len(p2)):
                if not is_valid_type(p2[i], int):
                    raise ValueError(
                        "Particle types in p1 and p2 have to be of type int" + str(p2[i]))

            set1 = create_int_list_from_python_object(p1)
            set2 = create_int_list_from_python_object(p2)
        return analyze.mindist(analyze.partCfg(), set1, set2)

    #
    # Distance to particle or point
    #

    def dist_to(self, id=None, pos=None):
        """
        Calculates the distance to a point or particle.

        Parameters
        ----------
        id : :obj:`int`, optional
            Calculate distance to particle with
            :attr:`~espressomd.particle_data.ParticleHandle.id` `id`.
        pos : array of :obj:`float`, optional
            Calculate distance to position `pos`.

        Returns
        -------
        :obj:`float`
            The calculated distance.

        """

        if id is None and pos is None:
            raise Exception(
                "Either id or pos have to be specified\n" + __doc__)

        if (id is not None) and (pos is not None):
            raise Exception(
                "Only one of id or pos may be specified\n" + __doc__)

        cdef Vector3d cpos
        if len(self._system.part) == 0:
            raise Exception("no particles")

        # Get position
        # If particle id specified
        if id is not None:
            if not is_valid_type(id, int):
                raise ValueError("Id has to be an integer")
            if not id in self._system.part[:].id:
                raise ValueError(
                    "Id has to be an index of an existing particle")
            _pos = self._system.part[id].pos
            for i in range(3):
                cpos[i] = _pos[i]
            _id = id
        else:
            for i in range(3):
                cpos[i] = pos[i]
            _id = -1
        return analyze.distto(analyze.partCfg(), cpos, _id)

    #
    # Analyze Linear Momentum
    #

    def linear_momentum(self, include_particles=True,
                                include_lbfluid=True):
        """
        Calculates the systems linear momentum.

        Parameters
        ----------
        include_particles : :obj:`bool`, optional
            whether to include the particles contribution to the linear
            momentum.
        include_lbfluid : :obj:`bool`, optional
            whether to include the Lattice Boltzmann fluid contribution
            to the linear momentum.

        Returns
        -------
        :obj:`float`
            The linear momentum of the system.

        """
        return analyze.calc_linear_momentum(include_particles, include_lbfluid)

    #
    # Analyze center of mass
    #

    def center_of_mass(self, p_type=None):
        """
        Calculates the systems center of mass.

        Parameters
        ----------
        p_type : :obj:`int`
            Particle :attr:`~espressomd.particle_data.ParticleHandle.type` for
            which to calculate the center of mass.

        Returns
        -------
        array of :obj:`float`
            The center of mass of the system.

        """
        if p_type is None:
            raise ValueError(
                "The p_type keyword argument must be provided (particle type)")
        check_type_or_throw_except(p_type, 1, int, "p_type has to be an int")
        if (p_type < 0 or p_type >= analyze.max_seen_particle_type):
            raise ValueError("Particle type", p_type, "does not exist!")

        return analyze.centerofmass(analyze.partCfg(), p_type)

    # get all particles in neighborhood r_catch of pos and return their ids
    # in il. plane can be used to specify the distance in the xy, xz or yz
    # plane

    def nbhood(self, pos=None, r_catch=None, plane='3d'):
        """
        Get all particles in a defined neighborhood.

        Parameters
        ----------
        pos : array of :obj:`float`
            Reference position for the neighborhood.
        r_catch : :obj:`float`
            Radius of the region.
        plane : :obj:`str`, \{'xy', 'xz', 'yz'\}
            If given, `r_catch` is the distance to the respective plane.

        Returns
        -------
        array of :obj:`int`
            The neighbouring particle ids.

        """

        cdef Vector3i planedims
        cdef List[int] ids
        cdef Vector3d c_pos

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

        ids = analyze.nbhood(analyze.partCfg(), c_pos, r_catch, planedims)

        return create_nparray_from_int_list(ids)

    def cylindrical_average(self, center=None, axis=None,
                            length=None, radius=None,
                            bins_axial=None, bins_radial=None,
                            types=[-1]):
        """
        Calculates the particle distribution using cylindrical binning.

        Parameters
        ----------
        center : (3,) array_like of :obj:`float`
            Coordinates of the centre of the cylinder.
        axis : (3,) array_like of :obj:`float`
            Axis vectory of the cylinder, does not need to be normalized.
        length : :obj:`float`
            Length of the cylinder.
        radius : :obj:`float`
            Radius of the cylinder.
        bins_axial : :obj:`int`
            Number of axial bins.
        bins_radial : :obj:`int`
            Number of radial bins.
        types : lists of :obj:`int`
            Particle :attr:`~espressomd.particle_data.ParticleHandle.type`

        Returns
        -------
        list of lists
            columns indicate `index_radial`, `index_axial`, `pos_radial`, `pos_axial`, `binvolume`, `density`, `v_radial`, `v_axial`, `density`, `v_radial` and `v_axial`.
            Note that the columns `density`, `v_radial` and `v_axial` appear for each type indicated in `types` in the same order.

        """

        # Check the input types
        check_type_or_throw_except(
            center, 3, float, "center has to be 3 floats")
        check_type_or_throw_except(
            axis, 3, float, "direction has to be 3 floats")
        check_type_or_throw_except(
            length, 1, float, "length has to be a float")
        check_type_or_throw_except(
            radius, 1, float, "radius has to be a floats")
        check_type_or_throw_except(
            bins_axial, 1, int, "bins_axial has to be an int")
        check_type_or_throw_except(
            bins_radial, 1, int, "bins_radial has to be an int")

        # Convert Python types to C++ types
        cdef vector[double] c_center = center
        cdef vector[double] c_direction = axis
        cdef double c_length = length
        cdef double c_radius = radius
        cdef int c_bins_axial = bins_axial
        cdef int c_bins_radial = bins_radial
        cdef vector[int] c_types = types

        cdef map[string, vector[vector[vector[double]]]] distribution
        analyze.calc_cylindrical_average(
            analyze.partCfg(), c_center, c_direction, c_length,
                                           c_radius, c_bins_axial, c_bins_radial, c_types,
                                           distribution)

        cdef double binwd_axial = c_length / c_bins_axial
        cdef double binwd_radial = c_radius / c_bins_radial
        cdef double binvolume, pos_radial, pos_axial

        cdef vector[string] names = [b"density", b"v_r", b"v_t"]

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

                buffer[index_axial + bins_axial *
                       index_radial, 0] = index_radial
                buffer[index_axial + bins_axial *
                       index_radial, 1] = index_axial
                buffer[index_axial + bins_axial * index_radial, 2] = pos_radial
                buffer[index_axial + bins_axial * index_radial, 3] = pos_axial
                buffer[index_axial + bins_axial * index_radial, 4] = binvolume
                for type_id in range(0, c_types.size()):
                    for name in range(0, names.size()):
                        buffer[index_axial + bins_axial * index_radial, 5 + name + type_id * names.size()] \
                            = distribution[names[name]][type_id][index_radial][index_axial]

        return buffer

    def pressure(self, v_comp=False):
        """Calculates the instantaneous pressure (in parallel). This is only
        sensible in an isotropic system which is homogeneous (on average)! Do
        not use this in an anisotropic or inhomogeneous system. In order to
        obtain the pressure the ensemble average needs to be calculated.

        Returns
        -------
        dict
            A dictionary with the following keys:

            * ``"total"``: total pressure
            * ``"kinetic"``: kinetic pressure
            * ``"bonded"``: total bonded pressure
            * ``"bonded", <bond_type>``: bonded pressure which arises from the given bond_type
            * ``"nonbonded"``: total nonbonded pressure
            * ``"nonbonded", <type_i>, <type_j>``: nonbonded pressure which arises from the interactions between type_i and type_j
            * ``"nonbonded_intra", <type_i>, <type_j>``: nonbonded pressure between short ranged forces between type i and j and with the same mol_id
            * ``"nonbonded_inter", <type_i>, <type_j>``: nonbonded pressure between short ranged forces between type i and j and different mol_ids
            * ``"coulomb"``: Coulomb pressure, how it is calculated depends on the method. It is equivalent to 1/3 of the trace of the Coulomb stress tensor.
              For how the stress tensor is calculated, see below. The averaged value in an isotropic NVT simulation is equivalent to the average of
              :math:`E^{coulomb}/(3V)`, see :cite:`brown95a`.
            * ``"dipolar"``: TODO
            * ``"virtual_sites"``: Stress contribution due to virtual sites

        """
        v_comp = int(v_comp)

        check_type_or_throw_except(v_comp, 1, int, "v_comp must be a boolean")

        # Dict to store the results
        p = OrderedDict()

        # Update in espresso core if necessary
        if (analyze.total_pressure.init_status != 1 + v_comp):
            analyze.update_pressure(v_comp)

        # Individual components of the pressure

        # Total pressure
        cdef int i
        total = 0
        for i in range(analyze.total_pressure.data.n):
            total += analyze.total_pressure.data.e[i]

        p["total"] = total

        # kinetic
        p["kinetic"] = analyze.total_pressure.data.e[0]

        # Bonded
        cdef double total_bonded
        total_bonded = 0
        for i in range(bonded_ia_params.size()):
            if (bonded_ia_params[i].type != BONDED_IA_NONE):
                p["bonded", i] = analyze.obsstat_bonded(& analyze.total_pressure, i)[0]
                total_bonded += analyze.obsstat_bonded(& analyze.total_pressure, i)[0]
        p["bonded"] = total_bonded

        # Non-Bonded interactions, total as well as intra and inter molecular
        cdef int j
        cdef double total_intra
        cdef double total_inter
        cdef double total_non_bonded
        total_inter = 0
        total_intra = 0
        total_non_bonded = 0

        for i in range(analyze.max_seen_particle_type):
            for j in range(i, analyze.max_seen_particle_type):
                #      if checkIfParticlesInteract(i, j):
                p["non_bonded", i, j] = analyze.obsstat_nonbonded(& analyze.total_pressure, i, j)[0]
                total_non_bonded += analyze.obsstat_nonbonded(& analyze.total_pressure, i, j)[0]
                total_intra += analyze.obsstat_nonbonded_intra(& analyze.total_pressure_non_bonded, i, j)[0]
                p["non_bonded_intra", i, j] = analyze.obsstat_nonbonded_intra(& analyze.total_pressure_non_bonded, i, j)[0]
                p["non_bonded_inter", i, j] = analyze.obsstat_nonbonded_inter(& analyze.total_pressure_non_bonded, i, j)[0]
                total_inter += analyze.obsstat_nonbonded_inter(& analyze.total_pressure_non_bonded, i, j)[0]
        p["non_bonded_intra"] = total_intra
        p["non_bonded_inter"] = total_inter
        p["non_bonded"] = total_non_bonded

        # Electrostatics
        IF ELECTROSTATICS == 1:
            cdef double total_coulomb
            total_coulomb = 0
            for i in range(analyze.total_pressure.n_coulomb):
                total_coulomb += analyze.total_pressure.coulomb[i]
                p["coulomb", i] = analyze.total_pressure.coulomb[i]
            p["coulomb"] = total_coulomb

        # Dipoles
        IF DIPOLES == 1:
            cdef double total_dipolar
            total_dipolar = 0
            for i in range(analyze.total_pressure.n_dipolar):
                total_dipolar += analyze.total_pressure.dipolar[i]
                p["dipolar", i] = analyze.total_pressure.coulomb[i]
            p["dipolar"] = total_dipolar

        # virtual sites
        IF VIRTUAL_SITES == 1:
            p_vs = 0.
            for i in range(analyze.total_pressure.n_virtual_sites):
                p_vs += analyze.total_pressure.virtual_sites[i]
                p["virtual_sites", i] = analyze.total_pressure.virtual_sites[
                    0]
            if analyze.total_pressure.n_virtual_sites:
                p["virtual_sites"] = p_vs

        return p

    def stress_tensor(self, v_comp=False):
        """Calculates the instantaneous stress tensor (in parallel). This is
        sensible in an anisotropic system. Still it assumes that the system is
        homogeneous since the volume averaged stress tensor is used. Do not use
        this stress tensor in an (on average) inhomogeneous system. If the
        system is (on average inhomogeneous) then use a local stress tensor.
        In order to obtain the stress tensor the ensemble average needs to be
        calculated.

        Returns
        -------
        dict
            A dictionary with the following keys:

            * ``"total"``: total stress tensor
            * ``"kinetic"``: kinetic stress tensor
            * ``"bonded"``: total bonded stress tensor
            * ``"bonded", <bond_type>``: bonded stress tensor which arises from the given bond_type
            * ``"nonbonded"``: total nonbonded stress tensor
            * ``"nonbonded", <type_i>, <type_j>``: nonbonded stress tensor which arises from the interactions between type_i and type_j
            * ``"nonbonded_intra" <type_i>, <type_j>``: nonbonded stress tensor between short ranged forces between type i and j and with the same mol_id
            * ``"nonbonded_inter" <type_i>, <type_j>``: nonbonded stress tensor between short ranged forces between type i and j and different mol_ids
            * ``"coulomb"``: Maxwell stress tensor, how it is calculated depends on the method
            * ``"dipolar"``: TODO
            * ``"virtual_sites"``: Stress tensor contribution for virtual sites

        """
        v_comp = int(v_comp)

        check_type_or_throw_except(v_comp, 1, int, "v_comp must be a boolean")

        # Dict to store the results
        p = OrderedDict()

        # Update in espresso core if necessary
        if (analyze.total_p_tensor.init_status != 1 + v_comp):
            analyze.update_pressure(v_comp)

        # Individual components of the pressure

        # Total pressure
        cdef int i
        total = np.zeros(9)
        for i in range(9):
            for k in range(analyze.total_p_tensor.data.n // 9):
                total[i] += analyze.total_p_tensor.data.e[9 * k + i]

        p["total"] = total.reshape((3, 3))

        # kinetic
        p["kinetic"] = create_nparray_from_double_array(
            analyze.total_p_tensor.data.e, 9)
        p["kinetic"] = p["kinetic"].reshape((3, 3))

        # Bonded
        total_bonded = np.zeros((3, 3))
        for i in range(bonded_ia_params.size()):
            if (bonded_ia_params[i].type != BONDED_IA_NONE):
                p["bonded", i] = np.reshape(create_nparray_from_double_array(
                    analyze.obsstat_bonded(& analyze.total_p_tensor, i), 9),
                    (3, 3))
                total_bonded += p["bonded", i]
        p["bonded"] = total_bonded

        # Non-Bonded interactions, total as well as intra and inter molecular
        cdef int j
        total_non_bonded = np.zeros((3, 3))
        total_non_bonded_intra = np.zeros((3, 3))
        total_non_bonded_inter = np.zeros((3, 3))

        for i in range(analyze.max_seen_particle_type):
            for j in range(i, analyze.max_seen_particle_type):
                #      if checkIfParticlesInteract(i, j):
                p["non_bonded", i, j] = np.reshape(
                    create_nparray_from_double_array(analyze.obsstat_nonbonded(
                        & analyze.total_p_tensor, i, j), 9), (3, 3))
                total_non_bonded += p["non_bonded", i, j]

                p["non_bonded_intra", i, j] = np.reshape(
                    create_nparray_from_double_array(
                        analyze.obsstat_nonbonded_intra(
                            & analyze.total_p_tensor_non_bonded, i, j), 9), (3, 3))
                total_non_bonded_intra += p["non_bonded_intra", i, j]

                p["non_bonded_inter", i, j] = np.reshape(
                    create_nparray_from_double_array(
                        analyze.obsstat_nonbonded_inter(
                            & analyze.total_p_tensor_non_bonded, i, j), 9), (3, 3))
                total_non_bonded_inter += p["non_bonded_inter", i, j]

        p["non_bonded_intra"] = total_non_bonded_intra
        p["non_bonded_inter"] = total_non_bonded_inter
        p["non_bonded"] = total_non_bonded

        # Electrostatics
        IF ELECTROSTATICS == 1:
            total_coulomb = np.zeros((3, 3))
            for i in range(analyze.total_p_tensor.n_coulomb):
                p["coulomb", i] = np.reshape(
                    create_nparray_from_double_array(
                        analyze.total_p_tensor.coulomb + 9 * i, 9), (3, 3))
                total_coulomb += p["coulomb", i]
            p["coulomb"] = total_coulomb

        # Dipoles
        IF DIPOLES == 1:
            total_dipolar = np.zeros((3, 3))
            for i in range(analyze.total_p_tensor.n_dipolar):
                p["dipolar", i] = np.reshape(
                    create_nparray_from_double_array(
                        analyze.total_p_tensor.dipolar + 9 * i, 9), (3, 3))
                total_dipolar += p["dipolar", i]
            p["dipolar"] = total_dipolar

        # virtual sites
        IF VIRTUAL_SITES_RELATIVE == 1:
            total_vs = np.zeros((3, 3))
            for i in range(analyze.total_p_tensor.n_virtual_sites):
                p["virtual_sites", i] = np.reshape(
                    create_nparray_from_double_array(
                        analyze.total_p_tensor.virtual_sites + 9 * i, 9), (3, 3))
                total_vs += p["virtual_sites", i]
            if analyze.total_p_tensor.n_virtual_sites:
                p["virtual_sites"] = total_vs

        return p

    IF DPD == 1:
        def dpd_stress(self):
            cdef Vector9d p
            p = dpd_stress()
            return array_locked((
                p[0], p[1], p[2],
                p[3], p[4], p[5],
                p[6], p[7], p[8])).reshape((3, 3))

    #
    # Energy analysis
    #

    def energy(self):
        """Calculate the systems energy.

        Returns
        -------
        :obj:`dict`
            A dictionary with keys `total`, `kinetic`, `bonded`, `nonbonded`,
            `coulomb`, `external_fields`.


        Examples
        --------
        >>> energy = system.analysis.energy()
        >>> print(energy["total"])
        >>> print(energy["kinetic"])
        >>> print(energy["bonded"])
        >>> print(energy["non_bonded"])
        >>> print(energy["external_fields"])

        """
    #  if system.n_part == 0:
    #    raise Exception('no particles')

        e = OrderedDict()

        if analyze.total_energy.init_status == 0:
            analyze.init_energies(& analyze.total_energy)
            analyze.master_energy_calc()
            handle_errors("calc_long_range_energies failed")

        # Individual components of the pressure

        # Total energy
        cdef int i
        total = analyze.total_energy.data.e[0]  # kinetic energy
        total += calculate_current_potential_energy_of_system()

        e["total"] = total
        e["external_fields"] = analyze.total_energy.external_fields[0]

        # Kinetic energy
        e["kinetic"] = analyze.total_energy.data.e[0]

        # Nonbonded
        cdef double total_bonded
        total_bonded = 0
        for i in range(bonded_ia_params.size()):
            if (bonded_ia_params[i].type != BONDED_IA_NONE):
                e["bonded", i] = analyze.obsstat_bonded(& analyze.total_energy, i)[0]
                total_bonded += analyze.obsstat_bonded(& analyze.total_energy, i)[0]
        e["bonded"] = total_bonded

        # Non-Bonded interactions, total as well as intra and inter molecular
        cdef int j
        cdef double total_intra
        cdef double total_inter
        cdef double total_non_bonded
        total_inter = 0
        total_intra = 0
        total_non_bonded = 0.

        for i in range(analyze.max_seen_particle_type):
            for j in range(analyze.max_seen_particle_type):
                #      if checkIfParticlesInteract(i, j):
                e["non_bonded", i, j] = analyze.obsstat_nonbonded(& analyze.total_energy, i, j)[0]
                if i <= j:
                    total_non_bonded += analyze.obsstat_nonbonded(& analyze.total_energy, i, j)[0]
    #        total_intra +=analyze.obsstat_nonbonded_intra(&analyze.total_energy_non_bonded, i, j)[0]
    #        e["non_bonded_intra",i,j] =analyze.obsstat_nonbonded_intra(&analyze.total_energy_non_bonded, i, j)[0]
    #        e["nonBondedInter",i,j] =analyze.obsstat_nonbonded_inter(&analyze.total_energy_non_bonded, i, j)[0]
    #        total_inter+= analyze.obsstat_nonbonded_inter(&analyze.total_energy_non_bonded, i, j)[0]
    #  e["nonBondedIntra"]=total_intra
    #  e["nonBondedInter"]=total_inter
        e["non_bonded"] = total_non_bonded

        # Electrostatics
        IF ELECTROSTATICS == 1:
            cdef double total_coulomb
            total_coulomb = 0
            for i in range(analyze.total_energy.n_coulomb):
                total_coulomb += analyze.total_energy.coulomb[i]
                e["coulomb", i] = analyze.total_energy.coulomb[i]
            e["coulomb"] = total_coulomb

        # Dipoles
        IF DIPOLES == 1:
            cdef double total_dipolar
            total_dipolar = 0
            for i in range(analyze.total_energy.n_dipolar):
                total_dipolar += analyze.total_energy.dipolar[i]
                e["dipolar", i] = analyze.total_energy.dipolar[i]
            e["dipolar"] = total_dipolar

        return e

    def calc_re(self, chain_start=None, number_of_chains=None,
                chain_length=None):
        """
        Calculates the mean end-to-end distance of chains and its
        standard deviation, as well as mean square end-to-end distance of
        chains and its standard deviation.

        This requires that a set of chains of equal length which start with the
        particle with particle number ``chain_start`` and are consecutively
        numbered, the last particle in that topology has id number
        ``chain_start + number_of_chains * chain_length - 1``.

        Parameters
        ----------
        chain_start : :obj:`int`
            The id of the first monomer of the first chain.
        number_of_chains : :obj:`int`
            Number of chains contained in the range.
        chain_length : :obj:`int`
            The length of every chain.

        Returns
        -------
        (4,) array_like of :obj:`float`
            Where [0] is the mean end-to-end distance of chains and [1] its
            standard deviation, [2] the mean square end-to-end distance and
            [3] its standard deviation.

        """
        self.check_topology(chain_start, number_of_chains, chain_length)
        re = analyze.calc_re(analyze.partCfg())
        return np.array([re[0], re[1], re[2], re[3]])

    def calc_rg(self, chain_start=None, number_of_chains=None,
                chain_length=None):
        """
        Calculates the mean radius of gyration of chains and its standard deviation,
        as well as the mean square radius of gyration of chains and its
        standard deviation.

        This requires that a set of chains of equal length which start with the
        particle with particle number ``chain_start`` and are consecutively
        numbered, the last particle in that topology has id number

        Parameters
        ----------
        chain_start : :obj:`int`
            The id of the first monomer of the first chain.
        number_of_chains : :obj:`int`
            Number of chains contained in the range.
        chain_length : :obj:`int`
            The length of every chain.

        Returns
        -------
        (4,) array_like of :obj:`float`
            Where [0] is the mean radius of gyration of the chains and [1] its
            standard deviation, [2] the mean square radius of gyration and [3]
            its standard deviation.

        """
        self.check_topology(chain_start, number_of_chains, chain_length)
        rg = analyze.calc_rg(analyze.partCfg())
        return np.array([rg[0], rg[1], rg[2], rg[3]])

    def calc_rh(self, chain_start=None, number_of_chains=None,
                chain_length=None):
        """
        Calculates the hydrodynamic mean radius of chains and its standard deviation.

        This requires that a set of chains of equal length which start with the
        particle with particle number ``chain_start`` and are consecutively
        numbered (the last particle in that topology has id number :
        ``chain_start + number_of_chains * chain_length - 1``).

        Parameters
        ----------
        chain_start : :obj:`int`
            The id of the first monomer of the first chain
        number_of_chains : :obj:`int`
            Number of chains contained in the range.
        chain_length : :obj:`int`
            The length of every chain.

        Returns
        -------
        (2,) array_like of :obj:`float`:
            Where [0] is the mean hydrodynamic radius of the chains
            and [1] its standard deviation.

        """

        self.check_topology(chain_start, number_of_chains, chain_length)
        rh = analyze.calc_rh(analyze.partCfg())
        return np.array([rh[0], rh[1]])

    def check_topology(self, chain_start=None, number_of_chains=None,
                       chain_length=None):
        check_type_or_throw_except(
            chain_start, 1, int, "chain_start=int is a required argument")
        check_type_or_throw_except(
            number_of_chains, 1, int, "number_of_chains=int is a required argument")
        check_type_or_throw_except(
            chain_length, 1, int, "chain_length=int is a required argument")
        id_min = chain_start
        id_max = chain_start + chain_length * number_of_chains
        for i in range(id_min, id_max):
            if (not self._system.part.exists(i)):
                raise ValueError('particle with id {0:.0f} does not exist\ncannot perform analysis on the range chain_start={1:.0f}, n_chains={2:.0f}, chain_length={3:.0f}\nplease provide a contiguous range of particle ids'.format(
                    i, chain_start, number_of_chains, chain_length));
        analyze.chain_start = chain_start
        analyze.chain_n_chains = number_of_chains
        analyze.chain_length = chain_length

    #
    # Structure factor
    #

    def structure_factor(self, sf_types=None, sf_order=None):
        """
        Calculate the structure factor for given types.  Returns the
        spherically averaged structure factor of particles specified in
        `types`.  The structure factor is calculated for all possible wave
        vectors q up to `order` Do not choose parameter `order` too large
        because the number of calculations grows as `order` to the third power.

        Parameters
        ----------
        sf_types : list of :obj:`int`
            Specifies which particle :attr:`~espressomd.particle_data.ParticleHandle.type`
            should be considered.
        sf_order : :obj:`int`
            Specifies the maximum wavevector.

        Returns
        -------
        :obj:`ndarray`
            Where [0] contains q
            and [1] contains the structure factor s(q)

        """

        if (sf_types is None) or (not hasattr(sf_types, '__iter__')):
            raise ValueError("sf_types has to be a list!")
        check_type_or_throw_except(
            sf_order, 1, int, "sf_order has to be an int!")

        p_types = create_int_list_from_python_object(sf_types)

        sf = analyze.calc_structurefactor(analyze.partCfg(), p_types.e,
                                          p_types.n, sf_order)

        return np.transpose(analyze.modify_stucturefactor(sf_order, sf.data()))

    #
    # RDF
    #

    def rdf(self, rdf_type=None, type_list_a=None, type_list_b=None,
            r_min=0.0, r_max=None, r_bins=100, n_conf=None):
        """
        Calculate a radial distribution function.
        The result is normalized by the spherical bin shell, the total number
        of particle pairs and the system volume.


        Parameters
        ----------
        rdf_type : :obj:`str`
            ``'rdf'`` or ``'<rdf>'``.
        type_list_a : lists of :obj:`int`
            Left :attr:`~espressomd.particle_data.ParticleHandle.type` of the rdf.
        type_list_b : lists of :obj:`int`, optional
            Right :attr:`~espressomd.particle_data.ParticleHandle.type` of the rdf.
        r_min : :obj:`float`
            Minimal distance to consider.
        r_max : :obj:`float`
            Maximal distance to consider
        r_bins : :obj:`int`
            Number of bins.
        n_conf : :obj:`int`, optional
            If rdf_type is ``'<rdf>'`` this determines
            the number of stored configs that are used.

        Returns
        -------
        :obj:`ndarray`
            Where [0] contains the midpoints of the bins,
            and [1] contains the values of the rdf.

        """

        if rdf_type is None:
            raise ValueError("rdf_type must not be empty!")
        if (type_list_a is None) or (not hasattr(type_list_a, '__iter__')):
            raise ValueError("type_list_a has to be a list!")
        if type_list_b and (not hasattr(type_list_b, '__iter__')):
            raise ValueError("type_list_b has to be a list!")
        if type_list_b is None:
            type_list_b = type_list_a

        if rdf_type != 'rdf':
            if n_configs == 0:
                raise ValueError("No configurations founds!\n",
                                 "Use 'analyze append' to save some,",
                                 "or 'analyze rdf' to only look at current RDF!""")
            if n_conf is None:
                n_conf = n_configs

        if r_max is None:
            r_max = min_box_l / 2.0

        cdef vector[double] rdf
        rdf.resize(r_bins)
        cdef vector[int] p1_types = type_list_a
        cdef vector[int] p2_types = type_list_b

        if rdf_type == 'rdf':
            analyze.calc_rdf(analyze.partCfg(), p1_types,
                               p2_types, r_min, r_max, r_bins, rdf)
        elif rdf_type == '<rdf>':
            analyze.calc_rdf_av(
                analyze.partCfg(), p1_types, p2_types, r_min,
                                  r_max, r_bins, rdf, n_conf)
        else:
            raise Exception(
                "rdf_type has to be one of 'rdf', '<rdf>', and '<rdf_intermol>'")

        r = np.empty(r_bins)
        bin_width = (r_max - r_min) / r_bins
        rr = r_min + bin_width / 2.0
        for i in range(r_bins):
            r[i] = rr
            rr += bin_width

        return np.array([r, rdf])

    #
    # distribution
    #

    def distribution(self, type_list_a=None, type_list_b=None,
                     r_min=0.0, r_max=None, r_bins=100, log_flag=0, int_flag=0):
        """
        Calculates the distance distribution of particles (probability of
        finding a particle of type A at a certain distance around a particle of
        type B, disregarding the fact that a spherical shell of a larger radius
        covers a larger volume). The distance is defined as the minimal distance
        between a particle of group ``type_list_a`` to any of the group
        ``type_list_b``. Returns two arrays, the bins and the (normalized)
        distribution.

        Parameters
        ----------
        type_list_a : list of :obj:`int`
            List of particle :attr:`~espressomd.particle_data.ParticleHandle.type`,
            only consider distances from these types.
        type_list_b : list of :obj:`int`
            List of particle :attr:`~espressomd.particle_data.ParticleHandle.type`,
            only consider distances to these types.
        r_min : :obj:`float`
            Minimum distance.
        r_max : :obj:`float`
            Maximum distance.
        r_bins : :obj:`int`
            Number of bins.
        log_flag : :obj:`int`
            When set to 0, the bins are linearly equidistant; when set to 1,
            the bins are logarithmically equidistant.
        int_flag : :obj:`int`
            When set to 1, the result is an integrated distribution.

        Returns
        -------
        :obj:`ndarray`
            Where [0] contains the midpoints of the bins,
            and [1] contains the values of the rdf.

        """

        if (type_list_a is None) or (not hasattr(type_list_a, '__iter__')):
            raise ValueError("type_list_a has to be a list!")
        if (type_list_b is None) or (not hasattr(type_list_b, '__iter__')):
            raise ValueError("type_list_b has to be a list!")

        if r_max is None:
            r_max = min_box_l / 2.0

        if r_min < 0.0 or (log_flag == 1 and r_min == 0.0):
            raise ValueError("r_min was chosen too small!")
        if r_max <= r_min:
            raise ValueError("r_max has to be greater than r_min!")
        if r_bins < 1:
            raise ValueError("r_bins has to be greater than zero!")

        cdef double low
        cdef vector[double] distribution
        distribution.resize(r_bins)

        p1_types = create_int_list_from_python_object(type_list_a)
        p2_types = create_int_list_from_python_object(type_list_b)

        analyze.calc_part_distribution(
            analyze.partCfg(), p1_types.e, p1_types.n, p2_types.e, p2_types.n,
            r_min, r_max, r_bins, log_flag, & low, distribution.data())

        np_distribution = create_nparray_from_double_array(
            distribution.data(), r_bins)

        if int_flag:
            np_distribution[0] += low
            for i in xrange(r_bins - 1):
                np_distribution[i + 1] += np_distribution[i]

        r = np.zeros(r_bins)

        if log_flag:
            log_fac = (float(r_max) / r_min)**(1.0 / r_bins)
            r[0] = r_min * np.sqrt(log_fac)
            for i in xrange(1, r_bins):
                r[i] = r[i - 1] * log_fac
        else:
            bin_width = (r_max - r_min) / float(r_bins)
            r = np.linspace(r_min + bin_width / 2.0,
                            r_max - bin_width / 2.0, r_bins)

        return np.array([r, np_distribution])

    #
    # angularmomentum
    #

    def angular_momentum(self, p_type=None):
        check_type_or_throw_except(
            p_type, 1, int, "p_type has to be an int")

        cdef double[3] com
        cdef int p1 = p_type

        analyze.angularmomentum(analyze.partCfg(), p1, com)

        return np.array([com[0], com[1], com[2]])

    #
    # gyration_tensor
    #

    def gyration_tensor(self, p_type=None):
        """
        Analyze the gyration tensor of particles of a given type or of all
        particles in the system if no type is given.

        Parameters
        ----------
        p_type : list of :obj:`int`, optional
            A particle :attr:`~espressomd.particle_data.ParticleHandle.type`,
            or list of all particle types to be considered.

        Returns
        -------
        dict
            A dictionary with the following keys:

            * ``"Rg^2"``: squared radius of gyration
            * ``"shape"``: three shape descriptors (asphericity, acylindricity, and relative shape anisotropy)
            * ``"eva0"``: eigenvalue 0 of the gyration tensor and its corresponding eigenvector.
            * ``"eva1"``: eigenvalue 1 of the gyration tensor and its corresponding eigenvector.
            * ``"eva2"``: eigenvalue 2 of the gyration tensor and its corresponding eigenvector.

            The eigenvalues are sorted in descending order.

        """
        if p_type is None:
            raise ValueError(
                "The p_type keyword argument must be provided (particle type)")
        if not hasattr(p_type, '__iter__'):
            p_type = [p_type]
        for type in p_type:
            check_type_or_throw_except(
                type, 1, int, "particle type has to be an int")
            if (type < 0 or type >= analyze.max_seen_particle_type):
                raise ValueError("Particle type", type, "does not exist!")
        selection = np.in1d(self._system.part[:].type, p_type)

        cm = np.mean(self._system.part[selection].pos, axis=0)
        mat = np.zeros(shape=(3, 3))
        for i, j in np.ndindex((3, 3)):
            mat[i, j] = np.mean(((self._system.part[selection].pos)[:, i] - cm[i]) * (
                (self._system.part[selection].pos)[:, j] - cm[j]))
        w, v = np.linalg.eig(mat)
        # return eigenvalue/vector tuples in order of increasing eigenvalues
        order = np.argsort(np.abs(w))[::-1]
        rad_gyr_sqr = mat[0, 0] + mat[1, 1] + mat[2, 2]
        aspheric = w[order[0]] - 0.5 * (w[order[1]] + w[order[2]])
        acylindric = w[order[1]] - w[order[2]]
        rel_shape_anis = (aspheric**2 + 0.75 * acylindric**2) / rad_gyr_sqr**2
        return{
        "Rg^2": rad_gyr_sqr,
        "shape": [aspheric,
                  acylindric,
                  rel_shape_anis],
        "eva0": (w[order[0]], v[:, order[0]]),
        "eva1": (w[order[1]], v[:, order[1]]),
        "eva2": (w[order[2]], v[:, order[2]])}

    #
    # momentofinertiamatrix
    #

    def moment_of_inertia_matrix(self, p_type=None):
        """
        Returns the 3x3 moment of inertia matrix for particles of a given type.

        Parameters
        ----------
        p_type : :obj:`int`
            A particle :attr:`~espressomd.particle_data.ParticleHandle.type`

        Returns
        -------
        :obj:`ndarray`
            3x3 moment of inertia matrix.

        """

        cdef double[9] MofImatrix

        if p_type is None:
            raise ValueError(
                "The p_type keyword argument must be provided (particle type)")
        check_type_or_throw_except(p_type, 1, int, "p_type has to be an int")
        if (p_type < 0 or p_type >= analyze.max_seen_particle_type):
            raise ValueError("Particle type", p_type, "does not exist!")

        analyze.momentofinertiamatrix(
            analyze.partCfg(), p_type, MofImatrix)

        MofImatrix_np = np.empty((9))
        for i in range(9):
            MofImatrix_np[i] = MofImatrix[i]

        return MofImatrix_np.reshape((3, 3))

    #
    # Vkappa
    #

    _Vkappa = {
        "Vk1": 0.0,
        "Vk2": 0.0,
        "avk": 0.0
    }

    def v_kappa(self, mode=None, Vk1=None, Vk2=None, avk=None):
        """
        .. todo:: Looks to be incomplete

        Calculates the compressibility thought volume fluctuations.

        Parameters
        ----------
        mode : :obj:`str`, \{'read', 'set' or 'reset'\}
            Mode.
        Vk1 : :obj:`float`
            Volume.
        Vk2 : :obj:`float`
            Volume squared.
        avk : :obj:`float`
            Number of averages.

        """

        check_type_or_throw_except(mode, 1, str, "mode has to be a string")

        if (mode == "reset"):
            self._Vkappa["Vk1"] = 0.0
            self._Vkappa["Vk2"] = 0.0
            self._Vkappa["avk"] = 0.0
        elif (mode == "read"):
            return self._Vkappa
        elif (mode == "set"):
            check_type_or_throw_except(Vk1, 1, float, "Vk1 has to be a float")
            self._Vkappa["Vk1"] = Vk1
            check_type_or_throw_except(Vk2, 1, float, "Vk2 has to be a float")
            self._Vkappa["Vk2"] = Vk2
            check_type_or_throw_except(avk, 1, float, "avk has to be a float")
            self._Vkappa["avk"] = avk
            if (self._Vkappa["avk"] <= 0.0):
                result = self._Vkappa["Vk1"] = self._Vkappa[
                    "Vk2"] = self._Vkappa["avk"] = 0.0
                raise Exception(
                    "ERROR: # of averages <avk> has to be positive! Resetting values.")
            else:
                result = self._Vkappa["Vk2"] / self._Vkappa["avk"] - \
                    (self._Vkappa["Vk1"] / self._Vkappa["avk"])**2
            return result
        else:
            raise Exception("ERROR: Unknown mode.")
