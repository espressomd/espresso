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
include "myconfig.pxi"
from . cimport analyze
from libcpp.vector cimport vector  # import std::vector as vector
from .interactions cimport BONDED_IA_NONE
from .interactions cimport bonded_ia_params
import numpy as np
cimport numpy as np
from .globals import Globals

from collections import OrderedDict
from .system import System
from .utils import array_locked, is_valid_type
from .utils cimport Vector3i, Vector3d, Vector9d
from .utils cimport handle_errors, check_type_or_throw_except
from .utils cimport create_nparray_from_double_array
from .particle_data cimport get_n_part


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
        assert get_n_part(), "No particles to append!"
        if get_n_configs() > 0:
            assert analyze.get_n_part_conf() == get_n_part(), \
                "All configurations stored must have the same length"

        analyze.analyze_append(analyze.partCfg())

    def min_dist(self, p1='default', p2='default'):
        """Minimal distance between two sets of particle types.

        Parameters
        ----------
        p1, p2 : lists of :obj:`int`
            Particle :attr:`~espressomd.particle_data.ParticleHandle.type` in
            both sets. If both are set to ``'default'``, the minimum distance
            of all pairs is returned.

        """

        if p1 == 'default' and p2 == 'default':
            return analyze.mindist(analyze.partCfg(), [], [])
        elif p1 == 'default' or p2 == 'default':
            raise ValueError("Both p1 and p2 have to be specified")
        else:
            for i in range(len(p1)):
                if not is_valid_type(p1[i], int):
                    raise TypeError(
                        "Particle types in p1 have to be of type int, got: " + repr(p1[i]))

            for i in range(len(p2)):
                if not is_valid_type(p2[i], int):
                    raise TypeError(
                        "Particle types in p2 have to be of type int, got: " + repr(p2[i]))

            return analyze.mindist(analyze.partCfg(), p1, p2)

    #
    # Analyze Linear Momentum
    #

    def linear_momentum(self, include_particles=True,
                        include_lbfluid=True):
        """
        Calculates the system's linear momentum.

        Parameters
        ----------
        include_particles : :obj:`bool`, optional
            whether to include the particles contribution to the linear
            momentum.
        include_lbfluid : :obj:`bool`, optional
            whether to include the lattice-Boltzmann fluid contribution
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
        Calculate the system's center of mass.

        Note that virtual sites are not included, as they do not have a
        meaningful mass.

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
        if p_type < 0 or p_type >= analyze.max_seen_particle_type:
            raise ValueError("Particle type {} does not exist!".format(p_type))

        return analyze.centerofmass(analyze.partCfg(), p_type)

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
            If given, ``r_catch`` is the distance to the respective plane.

        Returns
        -------
        array of :obj:`int`
            The neighbouring particle ids.

        """

        cdef Vector3i planedims
        cdef Vector3d c_pos

        check_type_or_throw_except(
            pos, 3, float, "pos=(float,float,float) must be passed to nbhood")
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
            raise ValueError(
                'Invalid argument for specifying plane, must be xy, xz, or yz plane')

        for i in range(3):
            c_pos[i] = pos[i]

        return analyze.nbhood(analyze.partCfg(), c_pos, r_catch, planedims)

    def pressure(self, v_comp=False):
        """Calculate the instantaneous pressure (in parallel). This is only
        sensible in an isotropic system which is homogeneous (on average)! Do
        not use this in an anisotropic or inhomogeneous system. In order to
        obtain the pressure, the ensemble average needs to be calculated.

        Returns
        -------
        :obj:`dict`
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

        # Update in ESPResSo core if necessary
        if analyze.total_pressure.init_status != 1 + v_comp:
            analyze.update_pressure(v_comp)

        # Individual components of the pressure

        # Total pressure
        cdef int i
        total = 0
        for i in range(analyze.total_pressure.data.size()):
            total += analyze.total_pressure.data[i]

        p["total"] = total

        # kinetic
        p["kinetic"] = analyze.total_pressure.data[0]

        # Bonded
        cdef double total_bonded
        total_bonded = 0
        for i in range(bonded_ia_params.size()):
            if bonded_ia_params[i].type != BONDED_IA_NONE:
                p["bonded", i] = analyze.obsstat_bonded( & analyze.total_pressure, i)[0]
                total_bonded += analyze.obsstat_bonded( & analyze.total_pressure, i)[0]
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
                p["non_bonded", i, j] = analyze.obsstat_nonbonded( & analyze.total_pressure, i, j)[0]
                total_non_bonded += analyze.obsstat_nonbonded( & analyze.total_pressure, i, j)[0]
                total_intra += analyze.obsstat_nonbonded_intra( & analyze.total_pressure_non_bonded, i, j)[0]
                p["non_bonded_intra", i, j] = analyze.obsstat_nonbonded_intra( & analyze.total_pressure_non_bonded, i, j)[0]
                p["non_bonded_inter", i, j] = analyze.obsstat_nonbonded_inter( & analyze.total_pressure_non_bonded, i, j)[0]
                total_inter += analyze.obsstat_nonbonded_inter( & analyze.total_pressure_non_bonded, i, j)[0]
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
        """Calculate the instantaneous stress tensor (in parallel). This is
        sensible in an anisotropic system. Still it assumes that the system is
        homogeneous since the volume averaged stress tensor is used. Do not use
        this stress tensor in an (on average) inhomogeneous system. If the
        system is (on average inhomogeneous) then use a local stress tensor.
        In order to obtain the stress tensor, the ensemble average needs to be
        calculated.

        Returns
        -------
        :obj:`dict`
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

        # Update in ESPResSo core if necessary
        if analyze.total_p_tensor.init_status != 1 + v_comp:
            analyze.update_pressure(v_comp)

        # Individual components of the pressure

        # Total pressure
        cdef int i
        total = np.zeros(9)
        for i in range(9):
            for k in range(analyze.total_p_tensor.data.size() // 9):
                total[i] += analyze.total_p_tensor.data[9 * k + i]

        p["total"] = total.reshape((3, 3))

        # kinetic
        p["kinetic"] = create_nparray_from_double_array(
            analyze.total_p_tensor.data.data(), 9)
        p["kinetic"] = p["kinetic"].reshape((3, 3))

        # Bonded
        total_bonded = np.zeros((3, 3))
        for i in range(bonded_ia_params.size()):
            if bonded_ia_params[i].type != BONDED_IA_NONE:
                p["bonded", i] = np.reshape(create_nparray_from_double_array(
                    analyze.obsstat_bonded( & analyze.total_p_tensor, i), 9),
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
            A dictionary with keys ``total``, ``kinetic``, ``bonded``, ``nonbonded``,
            ``coulomb``, ``external_fields``.


        Examples
        --------
        >>> energy = system.analysis.energy()
        >>> print(energy["total"])
        >>> print(energy["kinetic"])
        >>> print(energy["bonded"])
        >>> print(energy["non_bonded"])
        >>> print(energy["external_fields"])

        """
        #  assert len(self._system.part), "no particles in the system"

        e = OrderedDict()

        if analyze.total_energy.init_status == 0:
            analyze.init_energies( & analyze.total_energy)
            analyze.master_energy_calc()
            handle_errors("calc_long_range_energies failed")

        # Individual components of the pressure

        # Total energy
        cdef int i
        total = analyze.total_energy.data[0]  # kinetic energy
        total += calculate_current_potential_energy_of_system()

        e["total"] = total
        e["external_fields"] = analyze.total_energy.external_fields[0]

        # Kinetic energy
        e["kinetic"] = analyze.total_energy.data[0]

        # Non-bonded
        cdef double total_bonded
        total_bonded = 0
        for i in range(bonded_ia_params.size()):
            if bonded_ia_params[i].type != BONDED_IA_NONE:
                e["bonded", i] = analyze.obsstat_bonded( & analyze.total_energy, i)[0]
                total_bonded += analyze.obsstat_bonded( & analyze.total_energy, i)[0]
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
                e["non_bonded", i, j] = analyze.obsstat_nonbonded( & analyze.total_energy, i, j)[0]
                if i <= j:
                    total_non_bonded += analyze.obsstat_nonbonded( & analyze.total_energy, i, j)[0]
        #       total_intra +=analyze.obsstat_nonbonded_intra(&analyze.total_energy_non_bonded, i, j)[0]
        #       e["non_bonded_intra",i,j] =analyze.obsstat_nonbonded_intra(&analyze.total_energy_non_bonded, i, j)[0]
        #       e["nonBondedInter",i,j] =analyze.obsstat_nonbonded_inter(&analyze.total_energy_non_bonded, i, j)[0]
        #       total_inter+= analyze.obsstat_nonbonded_inter(&analyze.total_energy_non_bonded, i, j)[0]
        # e["nonBondedIntra"]=total_intra
        # e["nonBondedInter"]=total_inter
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
        Calculate the mean end-to-end distance of chains and its
        standard deviation, as well as mean square end-to-end distance of
        chains and its standard deviation.

        This requires that a set of chains of equal length which start
        with the particle number ``chain_start`` and are consecutively
        numbered, the last particle in that topology having id number
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
        re = analyze.calc_re(
            chain_start,
            number_of_chains,
            chain_length)
        return np.array([re[0], re[1], re[2], re[3]])

    def calc_rg(self, chain_start=None, number_of_chains=None,
                chain_length=None):
        """
        Calculate the mean radius of gyration of chains and its standard
        deviation, as well as the mean square radius of gyration of chains
        and its standard deviation.

        This requires that a set of chains of equal length which start
        with the particle number ``chain_start`` and are consecutively
        numbered, the last particle in that topology having id number
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
            Where [0] is the mean radius of gyration of the chains and [1] its
            standard deviation, [2] the mean square radius of gyration and [3]
            its standard deviation.

        """
        self.check_topology(chain_start, number_of_chains, chain_length)
        rg = analyze.calc_rg(
            chain_start,
            number_of_chains,
            chain_length)
        return np.array([rg[0], rg[1], rg[2], rg[3]])

    def calc_rh(self, chain_start=None, number_of_chains=None,
                chain_length=None):
        """
        Calculate the hydrodynamic mean radius of chains and its standard
        deviation.

        This requires that a set of chains of equal length which start
        with the particle number ``chain_start`` and are consecutively
        numbered, the last particle in that topology having id number
        ``chain_start + number_of_chains * chain_length - 1``.

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
        rh = analyze.calc_rh(
            chain_start,
            number_of_chains,
            chain_length)
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
            if not self._system.part.exists(i):
                raise ValueError('particle with id {0} does not exist\n'
                                 'cannot perform analysis on the range '
                                 'chain_start={1}, number_of_chains={2}, chain_length={3}\n'
                                 'please provide a contiguous range of particle ids'.format(
                                     i, chain_start, number_of_chains, chain_length))

    #
    # Structure factor
    #

    def structure_factor(self, sf_types=None, sf_order=None):
        """
        Calculate the structure factor for given types.  Returns the
        spherically averaged structure factor of particles specified in
        ``sf_types``.  The structure factor is calculated for all possible wave
        vectors q up to ``sf_order``. Do not choose parameter ``sf_order`` too
        large because the number of calculations grows as ``sf_order`` to the
        third power.

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

        sf = analyze.calc_structurefactor(
            analyze.partCfg(), sf_types, sf_order)

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
        rdf_type : :obj:`str`, \{'rdf', '<rdf>'\}
            Type of analysis.
        type_list_a : lists of :obj:`int`
            Left :attr:`~espressomd.particle_data.ParticleHandle.type` of the rdf.
        type_list_b : lists of :obj:`int`, optional
            Right :attr:`~espressomd.particle_data.ParticleHandle.type` of the rdf.
        r_min : :obj:`float`
            Minimal distance to consider.
        r_max : :obj:`float`
            Maximal distance to consider.
        r_bins : :obj:`int`
            Number of bins.
        n_conf : :obj:`int`, optional
            If ``rdf_type`` is ``'<rdf>'`` this determines
            the number of stored configs that are used (if
            ``None``, all configurations are used).

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
        if (type_list_b is not None) and (
                not hasattr(type_list_b, '__iter__')):
            raise ValueError("type_list_b has to be a list!")
        if type_list_b is None:
            type_list_b = type_list_a

        if rdf_type != 'rdf':
            if get_n_configs() == 0:
                raise ValueError("No configurations founds!\n",
                                 "Use `analyze.append()` to save configurations,",
                                 "or `analyze.rdf('rdf')` to only look at current RDF!""")
            if n_conf is None:
                n_conf = get_n_configs()

        if r_max is None:
            r_max = min(Globals().box_l) / 2

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
            raise ValueError("Unknown rdf_type value {!r}".format(rdf_type))

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
        log_flag : :obj:`bool`
            When set to ``False``, the bins are linearly equidistant; when set
            to ``True``, the bins are logarithmically equidistant.
        int_flag : :obj:`bool`
            When set to ``True``, the result is an integrated distribution.

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
            r_max = min(Globals().box_l) / 2

        assert r_min >= 0.0, "r_min was chosen too small!"
        assert not log_flag or r_min != 0.0, "r_min cannot include zero"
        assert r_max > r_min, "r_max has to be greater than r_min!"
        assert r_bins >= 1, "r_bins has to be greater than zero!"

        cdef double low
        cdef vector[double] distribution
        distribution.resize(r_bins)

        analyze.calc_part_distribution(
            analyze.partCfg(), type_list_a, type_list_b,
            r_min, r_max, r_bins, < bint > log_flag, & low, distribution.data())

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
        """
        Calculates the system's angular momentum with respect to the origin.

        Note that virtual sites are not included, as they do not have a meaningful mass.

        Parameters
        ----------
        p_type : :obj:`int`
            Particle :attr:`~espressomd.particle_data.ParticleHandle.type` for
            which to calculate the center of mass.

        Returns
        -------
        (3,) :obj:`ndarray` of :obj:`float`
           The center of mass of the system.

        """
        check_type_or_throw_except(
            p_type, 1, int, "p_type has to be an int")

        cdef int p1 = p_type
        cdef Vector3d res = analyze.angularmomentum(analyze.partCfg(), p1)

        return np.array([res[0], res[1], res[2]])

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
        :obj:`dict`
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
        for ptype in p_type:
            check_type_or_throw_except(
                ptype, 1, int, "particle type has to be an int")
            if ptype < 0 or ptype >= analyze.max_seen_particle_type:
                raise ValueError(
                    "Particle type {} does not exist!".format(ptype))
        selection = self._system.part.select(lambda p: (p.type in p_type))
        cm = np.mean(selection.pos, axis=0)
        mat = np.zeros(shape=(3, 3))
        for i, j in np.ndindex((3, 3)):
            mat[i, j] = np.mean(((selection.pos)[:, i] - cm[i]) * (
                (selection.pos)[:, j] - cm[j]))
        w, v = np.linalg.eig(mat)
        # return eigenvalue/vector tuples in order of increasing eigenvalues
        order = np.argsort(np.abs(w))[::-1]
        rad_gyr_sqr = mat[0, 0] + mat[1, 1] + mat[2, 2]
        aspheric = w[order[0]] - 0.5 * (w[order[1]] + w[order[2]])
        acylindric = w[order[1]] - w[order[2]]
        rel_shape_anis = (aspheric**2 + 0.75 * acylindric**2) / rad_gyr_sqr**2
        return {
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
        if p_type < 0 or p_type >= analyze.max_seen_particle_type:
            raise ValueError("Particle type {} does not exist!".format(p_type))

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

        Calculate the compressibility through volume fluctuations.

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

        if mode == "reset":
            self._Vkappa["Vk1"] = 0.0
            self._Vkappa["Vk2"] = 0.0
            self._Vkappa["avk"] = 0.0
        elif mode == "read":
            return self._Vkappa
        elif mode == "set":
            check_type_or_throw_except(Vk1, 1, float, "Vk1 has to be a float")
            self._Vkappa["Vk1"] = Vk1
            check_type_or_throw_except(Vk2, 1, float, "Vk2 has to be a float")
            self._Vkappa["Vk2"] = Vk2
            check_type_or_throw_except(avk, 1, float, "avk has to be a float")
            self._Vkappa["avk"] = avk
            if self._Vkappa["avk"] <= 0.0:
                result = self._Vkappa["Vk1"] = self._Vkappa[
                    "Vk2"] = self._Vkappa["avk"] = 0.0
                raise ValueError(
                    "# of averages <avk> has to be positive! Resetting values.")
            else:
                result = self._Vkappa["Vk2"] / self._Vkappa["avk"] - \
                    (self._Vkappa["Vk1"] / self._Vkappa["avk"])**2
            return result
        else:
            raise ValueError("Unknown mode {!r}".format(mode))
