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
import numpy as np
cimport numpy as np
from .grid cimport box_geo

from .system import System
from .utils import array_locked, is_valid_type, handle_errors
from .utils cimport Vector3i, Vector3d, Vector9d
from .utils cimport make_Vector3d
from .utils cimport make_array_locked
from .utils cimport check_type_or_throw_except
from .utils cimport create_nparray_from_double_array


def autocorrelation(time_series):
    """
    Calculate the unnormalized autocorrelation function :math:`R_{XX}`
    of an observable :math:`X` measured at time :math:`t` with constant
    time step for lag times :math:`\\tau` such that:

    :math:`\\displaystyle R_{XX}(\\tau) = \\frac{1}{N - \\tau} \\sum_{t=0}^{N - \\tau} X_{t} \\cdot X_{t + \\tau}`

    This is a scipy implementation of the algorithm described
    in :cite:`debuyl18a` that uses FFT. This implementation:

    - doesn't subtract the mean of the time series before calculating the ACF,
    - doesn't normalize the ACF by the variance of the time series,
    - assumes the time series is aperiodic (i.e. uses zero-padding).

    Parameters
    ----------
    time_series : (N,) or (N, M) array_like of :obj:`float`
        Time series to correlate. For multi-dimensional data, M is the
        number of columns.

    Returns
    -------
    (N,) array_like of :obj:`float`
        The time series autocorrelation function.
    """
    import scipy.signal

    def acf_1d(signal, n_with_padding, n):
        acf = scipy.signal.correlate(signal, signal, mode="full", method="fft")
        acf = acf[-n_with_padding:][:n] / (n - np.arange(n))
        return acf

    if not isinstance(time_series, np.ndarray):
        time_series = np.array(time_series)

    n = time_series.shape[0]
    if n == 0:
        return np.array([], dtype=float)

    n_with_padding = 2**int(np.ceil(np.log2(n)) + 1)
    signal_padded = np.zeros(n_with_padding)

    if time_series.ndim == 1:
        signal_padded[:n] = time_series
        return acf_1d(signal_padded, n_with_padding, n)
    elif time_series.ndim == 2:
        signal_acf = np.zeros(n)
        for i in range(time_series.shape[1]):
            signal_padded[:n] = time_series[:, i]
            signal_acf += acf_1d(signal_padded, n_with_padding, n)
        return signal_acf
    else:
        raise ValueError(f"Only 1-dimensional and 2-dimensional time series "
                         f"are supported, got shape {time_series.shape}")


class Analysis:

    def __init__(self, system):
        if not isinstance(system, System):
            raise TypeError("An instance of System is required as argument")
        self._system = system

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
                        f"Particle types in p1 have to be of type int, got: {repr(p1[i])}")

            for i in range(len(p2)):
                if not is_valid_type(p2[i], int):
                    raise TypeError(
                        f"Particle types in p2 have to be of type int, got: {repr(p2[i])}")

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
            raise ValueError(f"Particle type {p_type} does not exist!")

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

        check_type_or_throw_except(pos, 3, float, "pos must be 3 floats")
        check_type_or_throw_except(
            r_catch, 1, float, "r_catch must be a float")

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

        return analyze.nbhood(
            analyze.partCfg(), make_Vector3d(pos), r_catch, planedims)

    def pressure(self):
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
            * ``"non_bonded"``: total non-bonded pressure
            * ``"non_bonded", <type_i>, <type_j>``: non-bonded pressure which arises from the interactions between type_i and type_j
            * ``"non_bonded_intra", <type_i>, <type_j>``: non-bonded pressure between short ranged forces between type i and j and with the same mol_id
            * ``"non_bonded_inter", <type_i>, <type_j>``: non-bonded pressure between short ranged forces between type i and j and different mol_ids
            * ``"coulomb"``: Coulomb pressure, how it is calculated depends on the method. It is equivalent to 1/3 of the trace of the Coulomb pressure tensor.
              For how the pressure tensor is calculated, see :ref:`Pressure Tensor`. The averaged value in an isotropic NVT simulation is equivalent to the average of
              :math:`E^{coulomb}/(3V)`, see :cite:`brown95a`.
            * ``"coulomb", <i>``: Coulomb pressure from particle pairs (``i=0``), electrostatics solvers (``i=1``)
            * ``"dipolar"``: not implemented
            * ``"virtual_sites"``: Pressure contribution from virtual sites
            * ``"external_fields"``: external fields contribution

        """

        obs = analyze.get_scalar_pressure()
        handle_errors("calculate_pressure() failed")
        return obs

    def pressure_tensor(self):
        """Calculate the instantaneous pressure_tensor (in parallel). This is
        sensible in an anisotropic system. Still it assumes that the system is
        homogeneous since the volume averaged pressure_tensor is used. Do not use
        this pressure_tensor in an (on average) inhomogeneous system. If the
        system is (on average inhomogeneous) then use a local pressure_tensor.
        In order to obtain the pressure_tensor, the ensemble average needs to be
        calculated.

        Returns
        -------
        :obj:`dict`
            A dictionary with the following keys:

            * ``"total"``: total pressure tensor
            * ``"kinetic"``: kinetic pressure tensor
            * ``"bonded"``: total bonded pressure tensor
            * ``"bonded", <bond_type>``: bonded pressure tensor which arises from the given bond_type
            * ``"non_bonded"``: total non-bonded pressure tensor
            * ``"non_bonded", <type_i>, <type_j>``: non-bonded pressure tensor which arises from the interactions between type_i and type_j
            * ``"non_bonded_intra", <type_i>, <type_j>``: non-bonded pressure tensor between short ranged forces between type i and j and with the same mol_id
            * ``"non_bonded_inter", <type_i>, <type_j>``: non-bonded pressure tensor between short ranged forces between type i and j and different mol_ids
            * ``"coulomb"``: Maxwell pressure tensor, how it is calculated depends on the method
            * ``"coulomb", <i>``: Maxwell pressure tensor from particle pairs (``i=0``), electrostatics solvers (``i=1``)
            * ``"dipolar"``: not implemented
            * ``"virtual_sites"``: pressure tensor contribution from virtual sites
            * ``"external_fields"``: external fields contribution

        """

        obs = analyze.get_pressure_tensor()
        handle_errors("calculate_pressure() failed")
        return obs

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
            A dictionary with the following keys:

            * ``"total"``: total energy
            * ``"kinetic"``: linear and rotational kinetic energy
            * ``"bonded"``: total bonded energy
            * ``"bonded", <bond_type>``: bonded energy which arises from the given bond_type
            * ``"non_bonded"``: total non-bonded energy
            * ``"non_bonded", <type_i>, <type_j>``: non-bonded energy which arises from the interactions between type_i and type_j
            * ``"non_bonded_intra", <type_i>, <type_j>``: non-bonded energy between short ranged forces between type i and j and with the same mol_id
            * ``"non_bonded_inter", <type_i>, <type_j>``: non-bonded energy between short ranged forces between type i and j and different mol_ids
            * ``"coulomb"``: Coulomb energy, how it is calculated depends on the method
            * ``"coulomb", <i>``: Coulomb energy from particle pairs (``i=0``), electrostatics solvers (``i=1``)
            * ``"dipolar"``: dipolar energy
            * ``"dipolar", <i>``: dipolar energy from particle pairs and magnetic field constraints (``i=0``), magnetostatics solvers (``i=1``)
            * ``"external_fields"``: external fields contribution


        Examples
        --------
        >>> energy = system.analysis.energy()
        >>> print(energy["total"])
        >>> print(energy["kinetic"])
        >>> print(energy["bonded"])
        >>> print(energy["non_bonded"])
        >>> print(energy["external_fields"])

        """

        obs = analyze.get_energy()
        handle_errors("calculate_energy() failed")
        return obs

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
                raise ValueError(f'particle with id {i} does not exist\n'
                                 f'cannot perform analysis on the range '
                                 f'chain_start={chain_start}, number_of_chains='
                                 f'{number_of_chains}, chain_length={chain_length}\n'
                                 f'please provide a contiguous range of particle ids')

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

        if sf_types is None or not hasattr(sf_types, '__iter__'):
            raise ValueError("sf_types has to be a list!")
        check_type_or_throw_except(
            sf_order, 1, int, "sf_order has to be an int!")

        cdef vector[double] wavevectors
        cdef vector[double] intensities
        analyze.calc_structurefactor(
            analyze.partCfg(), sf_types, sf_order, wavevectors, intensities)

        return np.vstack([wavevectors, intensities])

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
            box_l = make_array_locked(< Vector3d > box_geo.length())
            r_max = min(box_l) / 2

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

        cdef Vector3d res = analyze.angularmomentum(analyze.partCfg(), p_type)

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
                raise ValueError(f"Particle type {ptype} does not exist!")
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
            raise ValueError(f"Particle type {p_type} does not exist!")

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
            raise ValueError(f"Unknown mode {mode!r}")
