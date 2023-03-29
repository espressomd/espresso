#
# Copyright (C) 2013-2022 The ESPResSo project
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
import numpy as np

from . import utils
from .code_features import assert_features, has_features
from .script_interface import script_interface_register, ScriptInterfaceHelper


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


@script_interface_register
class _ObservableStat(ScriptInterfaceHelper):
    _so_name = "ScriptInterface::Analysis::ObservableStat"
    _so_creation_policy = "GLOBAL"

    def _generate_summary(self, obj, dim, calc_sp):
        """
        Compute derived quantities and reshape pressure tensors as 3x3 matrices.
        """

        def zero():
            if dim == 1 or calc_sp:
                return 0.
            return np.zeros(9, dtype=float)

        def reduction(obj, key):
            total = zero()
            for k in obj.keys():
                if isinstance(k, tuple) and k[0] == key:
                    total += obj[k]
            obj[key] = total

        out = {}
        out["bonded"] = zero()
        out["non_bonded_intra"] = zero()
        out["non_bonded_inter"] = zero()
        for k, v in obj.items():
            if "," not in k:
                out[k] = v
            else:
                k = k.split(",")
                k = (k[0], *map(int, k[1:]))
                out[k] = v
                if k[0] == "bonded":
                    out["bonded"] += v
                elif k[0].startswith("non_bonded_"):
                    if k[0] == "non_bonded_intra":
                        out["non_bonded_intra"] += v
                    else:
                        out["non_bonded_inter"] += v
                    k = ("non_bonded", *k[1:])
                    if k not in out:
                        out[k] = zero()
                    out[k] += v

        out["non_bonded"] = out["non_bonded_intra"] + out["non_bonded_inter"]
        if has_features("ELECTROSTATICS"):
            reduction(out, "coulomb")
        if has_features("DIPOLES"):
            reduction(out, "dipolar")
        if has_features("VIRTUAL_SITES"):
            reduction(out, "virtual_sites")

        if dim == 1 or calc_sp:
            return out
        return {k: np.reshape(v, (3, 3)) for k, v in out.items()}


@script_interface_register
class Analysis(ScriptInterfaceHelper):
    """
    Methods
    -------
    linear_momentum()
        Calculates the system's linear momentum.

        Parameters
        ----------
        include_particles : :obj:`bool`, optional
            Whether to include the particles contribution to the linear
            momentum. True by default.
        include_lbfluid : :obj:`bool`, optional
            whether to include the lattice-Boltzmann fluid contribution
            to the linear momentum. True by default.

        Returns
        -------
        (3,) array_like of :obj:`float`
            The linear momentum of the system.

    center_of_mass()
        Calculate the center of mass of particles of a given type.

        Note that virtual sites are not included, as they do not have a
        meaningful mass.

        Parameters
        ----------
        p_type : :obj:`int`
            Particle :attr:`~espressomd.particle_data.ParticleHandle.type`
            for which to calculate the center of mass.

        Returns
        -------
        (3,) array_like of :obj:`float`
            The center of mass of the particles.

    nbhood()
        Get all particles in a defined neighborhood.

        Parameters
        ----------
        pos : array of :obj:`float`
            Reference position for the neighborhood.
        r_catch : :obj:`float`
            Radius of the region.

        Returns
        -------
        (N,) array_like of :obj:`int`
            The neighbouring particle ids.

    particle_neighbor_pids()
        Get a list of all short-range neighbors for each particle.

        Returns
        -------
        :obj: `dict`
            A dictionary where each item is a pair of a particle id and
            its respective neighboring particle ids.

    calc_re()
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

    calc_rg()
        Calculate the mean radius of gyration of chains and its standard
        deviation, as well as the mean square radius of gyration of chains
        and its standard deviation.

        This requires that a set of chains of equal length which start
        with the particle number ``chain_start`` and are consecutively
        numbered, the last particle in that topology having id number
        ``chain_start + number_of_chains * chain_length - 1``.

        The radius of gyration is the radius of a sphere which would have
        the same moment of inertia as a polymer, and is defined as

        .. math::

           R_{\\mathrm G}^2 = \\frac{1}{N} \\sum\\limits_{i=1}^{N} \\left(\\vec r_i - \\vec r_{\\mathrm{cm}}\\right)^2\\,,

        where :math:`\\vec r_i` are position vectors of individual particles
        constituting the polymer and :math:`\\vec r_{\\mathrm{cm}}` is the
        position of its center of mass. The sum runs over all :math:`N`
        particles comprising the polymer. For more information see any
        polymer science book, e.g. :cite:`rubinstein03a`.

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

    calc_rh()
        Calculate the mean hydrodynamic radius of chains and its standard
        deviation.

        This requires that a set of chains of equal length which start
        with the particle number ``chain_start`` and are consecutively
        numbered, the last particle in that topology having id number
        ``chain_start + number_of_chains * chain_length - 1``.

        The following formula is used for the calculation:

        .. math::

           \\frac{1}{R_{\\mathrm H}} = \\frac{2}{N(N-1)} \\sum\\limits_{i=1}^{N} \\sum\\limits_{j<i}^{N} \\frac{1}{|\\vec r_i - \\vec r_j|}\\,,

        This formula is only valid under certain assumptions. For more
        information, see chapter 4 and equation 4.102 in :cite:`doi86a`.
        Note that the hydrodynamic radius is sometimes defined in a similar
        fashion but with a denominator of :math:`N^2` instead of :math:`N(N-1)`
        in the prefactor. Both versions are equivalent in the
        :math:`N\\rightarrow \\infty` limit but give numerically different
        values for finite polymers.

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

    angular_momentum()
        Calculate the system's angular momentum with respect to the origin.

        Note that virtual sites are not included, as they do not have a
        meaningful mass.

        Parameters
        ----------
        p_type : :obj:`int`
            Particle :attr:`~espressomd.particle_data.ParticleHandle.type` for
            which to calculate the center of mass, or ``-1`` for all particles.

        Returns
        -------
        (3,) array_like of :obj:`float`
           The center of mass of the system.

    structure_factor()
        Calculate the structure factor for given types.  Returns the
        spherically averaged structure factor of particles specified in
        ``sf_types``.  The structure factor is calculated for all possible wave
        vectors q up to ``sf_order``. Do not choose parameter ``sf_order`` too
        large because the number of calculations grows as ``sf_order`` to the
        third power.

        Parameters
        ----------
        sf_types : list of :obj:`int`
            Particle :attr:`~espressomd.particle_data.ParticleHandle.type`
            which should be considered.
        sf_order : :obj:`int`
            Specifies the maximum wavevector.

        Returns
        -------
        :obj:`ndarray`
            Where [0] contains the wave vectors q
            and [1] contains the structure factors S(q)

    distribution()
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
            Minimum distance. Default is 0.
        r_max : :obj:`float`
            Maximum distance. By default, it is half the box size.
            A value larger than half the box size is allowed for systems
            with :ref:`open boundary conditions <Open boundaries>`.
        r_bins : :obj:`int`
            Number of bins. Default is 100.
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
    _so_name = "ScriptInterface::Analysis::Analysis"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = (
        "linear_momentum",
        "center_of_mass",
        "nbhood",
        "particle_neighbor_pids",
        "calc_re",
        "calc_rg",
        "calc_rh",
        "angular_momentum",
        "structure_factor",
        "distribution")

    def min_dist(self, p1='default', p2='default'):
        """
        Minimal distance between two sets of particle types.

        Parameters
        ----------
        p1, p2 : lists of :obj:`int`
            Particle :attr:`~espressomd.particle_data.ParticleHandle.type` in
            both sets. If both are set to ``'default'``, the minimum distance
            of all pairs is returned.

        Returns
        -------
        :obj:`float`
            Minimal distance.

        """

        if p1 == 'default' and p2 == 'default':
            p1 = []
            p2 = []
        elif p1 == 'default' or p2 == 'default':
            raise ValueError("Both p1 and p2 have to be specified")
        return self.call_method("min_dist", p_types1=p1, p_types2=p2)

    def pressure(self):
        """
        Calculate the instantaneous scalar pressure in parallel. This is only
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
            * ``"bonded", <bond_id>``: bonded pressure from the bond
              identified by ``bond_id``
            * ``"non_bonded"``: total non-bonded pressure
            * ``"non_bonded", <type_i>, <type_j>``: non-bonded pressure which
              arises from the interactions between ``type_i`` and ``type_j``
            * ``"non_bonded_intra", <type_i>, <type_j>``: non-bonded pressure
              from short-range forces between ``type_i`` and ``type_j``
              with the same ``mol_id``
            * ``"non_bonded_inter", <type_i>, <type_j>``: non-bonded pressure
              from short-range forces between ``type_i`` and ``type_j``
              with different ``mol_id``
            * ``"coulomb"``: Coulomb pressure, how it is calculated depends on
              the method. It is equivalent to 1/3 of the trace of the Coulomb
              pressure tensor. For how the pressure tensor is calculated,
              see :ref:`Pressure Tensor`. The averaged value in an isotropic
              NVT simulation is equivalent to the average of
              :math:`E^{\\mathrm{coulomb}}/(3V)`, see :cite:`brown95a`.
            * ``"coulomb", <i>``: Coulomb pressure from particle pairs
              (``i=0``), electrostatics solvers (``i=1``)
            * ``"dipolar"``: not implemented
            * ``"virtual_sites"``: Pressure contribution from virtual sites
            * ``"external_fields"``: external fields contribution

        """
        obs_stat = _ObservableStat()
        observable = obs_stat.call_method("calculate_scalar_pressure")
        utils.handle_errors("calculate_pressure() failed")
        return obs_stat._generate_summary(observable, 9, True)

    def pressure_tensor(self):
        """
        Calculate the instantaneous pressure tensor in parallel. This is
        sensible in an anisotropic system. Still it assumes that the system is
        homogeneous since the volume-averaged pressure tensor is used. Do not use
        this pressure tensor in an (on average) inhomogeneous system. If the
        system is (on average inhomogeneous) then use a local pressure tensor.
        In order to obtain the pressure tensor, the ensemble average needs to be
        calculated.

        Returns
        -------
        :obj:`dict`
            A dictionary with the following keys:

            * ``"total"``: total pressure tensor
            * ``"kinetic"``: kinetic pressure tensor
            * ``"bonded"``: total bonded pressure tensor
            * ``"bonded", <bond_id>``: bonded pressure tensor from the bond
              identified by ``bond_id``
            * ``"non_bonded"``: total non-bonded pressure tensor
            * ``"non_bonded", <type_i>, <type_j>``: non-bonded pressure tensor
              from short-range forces between ``type_i`` and ``type_j``
            * ``"non_bonded_intra", <type_i>, <type_j>``: non-bonded pressure
              tensor from short-range forces between ``type_i`` and ``type_j``
              with the same ``mol_id``
            * ``"non_bonded_inter", <type_i>, <type_j>``: non-bonded pressure
              tensor from short-range forces between ``type_i`` and ``type_j``
              with different ``mol_id``
            * ``"coulomb"``: Maxwell pressure tensor, how it is calculated
              depends on the method
            * ``"coulomb", <i>``: Maxwell pressure tensor from particle pairs
              (``i=0``), electrostatics solvers (``i=1``)
            * ``"dipolar"``: not implemented
            * ``"virtual_sites"``: pressure tensor contribution from virtual sites
            * ``"external_fields"``: external fields contribution

        """
        obs_stat = _ObservableStat()
        observable = obs_stat.call_method("calculate_pressure_tensor")
        utils.handle_errors("calculate_pressure() failed")
        return obs_stat._generate_summary(observable, 9, False)

    def energy(self):
        """
        Calculate the system energy in parallel.

        Returns
        -------
        :obj:`dict`
            A dictionary with the following keys:

            * ``"total"``: total energy
            * ``"kinetic"``: linear and rotational kinetic energy
            * ``"bonded"``: total bonded energy
            * ``"bonded", <bond_id>``: bonded energy from the bond
              identified by ``bond_id``
            * ``"non_bonded"``: total non-bonded energy
            * ``"non_bonded", <type_i>, <type_j>``: non-bonded energy
              from short-range interactions between ``type_i`` and ``type_j``
            * ``"non_bonded_intra", <type_i>, <type_j>``: non-bonded energy
              from short-range interactions between ``type_i`` and ``type_j``
              with the same ``mol_id``
            * ``"non_bonded_inter", <type_i>, <type_j>``: non-bonded energy
              from short-range interactions between ``type_i`` and ``type_j``
              with different ``mol_id``
            * ``"coulomb"``: Coulomb energy, how it is calculated depends
              on the method
            * ``"coulomb", <i>``: Coulomb energy from particle pairs
              (``i=0``), electrostatics solvers (``i=1``)
            * ``"dipolar"``: dipolar energy
            * ``"dipolar", <i>``: dipolar energy from particle pairs and
              magnetic field constraints (``i=0``), magnetostatics solvers
              (``i=1``)
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
        obs_stat = _ObservableStat()
        observable = obs_stat.call_method("calculate_energy")
        utils.handle_errors("calculate_energy() failed")
        return obs_stat._generate_summary(observable, 1, False)

    def particle_energy(self, particle):
        """
        Calculate the non-bonded energy of a single given particle.

        Parameters
        ----------
        particle : :class:`~espressomd.particle_data.ParticleHandle`

        Returns
        -------
        :obj: `float`
            Non-bonded energy of that particle

        """
        return self.call_method("particle_energy", pid=particle.id)

    def dpd_stress(self):
        assert_features("DPD")
        return np.reshape(self.call_method("dpd_stress"), (3, 3))

    def gyration_tensor(self, p_type=None):
        """
        Analyze the gyration tensor of particles of a given type.

        Parameters
        ----------
        p_type : list of :obj:`int`
            A particle :attr:`~espressomd.particle_data.ParticleHandle.type`,
            or list of particle types to be considered.

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
        vec = self.call_method("gyration_tensor", p_types=p_type)
        mat = np.reshape(vec, (3, 3))
        w, v = np.linalg.eig(mat)
        # return eigenvalue/vector tuples in order of increasing eigenvalues
        order = np.argsort(np.abs(w))[::-1]
        rad_gyr_sqr = np.trace(mat)
        aspheric = w[order[0]] - 0.5 * (w[order[1]] + w[order[2]])
        acylindric = w[order[1]] - w[order[2]]
        rel_shape_anis = (aspheric**2 + 0.75 * acylindric**2) / rad_gyr_sqr**2
        return {
            "Rg^2": rad_gyr_sqr,
            "shape": [aspheric, acylindric, rel_shape_anis],
            "eva0": (w[order[0]], v[:, order[0]]),
            "eva1": (w[order[1]], v[:, order[1]]),
            "eva2": (w[order[2]], v[:, order[2]])}

    def moment_of_inertia_matrix(self, p_type=None):
        """
        Returns the 3x3 moment of inertia matrix for particles of a given type.

        Parameters
        ----------
        p_type : :obj:`int`
            A particle :attr:`~espressomd.particle_data.ParticleHandle.type`.

        Returns
        -------
        (3,3) array_like of :obj:`float`
            Moment of inertia matrix.

        """
        vec = self.call_method("moment_of_inertia_matrix", p_type=p_type)
        return np.reshape(vec, (3, 3))
