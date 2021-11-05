# Copyright (C) 2010-2019 The ESPResSo project
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
include "myconfig.pxi"
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr
from cython.operator cimport dereference as deref
import numpy as np
import warnings


cdef class ReactionAlgorithm:
    """

    This class provides the base class for Reaction Algorithms like
    the Reaction Ensemble algorithm and the constant pH method.
    Initialize the reaction algorithm by setting the
    standard pressure, temperature, and the exclusion radius.

    Note: When creating particles the velocities of the new particles are set
    according the Maxwell-Boltzmann distribution. In this step the mass of the
    new particle is assumed to equal 1.


    Parameters
    ----------
    kT : :obj:`float`
        Thermal energy of the system in simulation units
    exclusion_radius : :obj:`float`
        Minimal distance from any particle, within which new particle will not
        be inserted. This is useful to avoid integrator failures if particles
        are too close and there is a diverging repulsive interaction, or to
        prevent two oppositely charged particles from being placed on top of
        each other. The Boltzmann factor :math:`\\exp(-\\beta E)` gives these
        configurations a small contribution to the partition function,
        therefore they can be neglected.
    seed : :obj:`int`
        Initial counter value (or seed) of the Mersenne Twister RNG.
    """
    cdef object _params
    cdef CReactionAlgorithm * RE

    def _valid_keys(self):
        return "kT", "exclusion_radius", "seed"

    def _required_keys(self):
        return "kT", "exclusion_radius", "seed"

    def _set_params_in_es_core(self):
        deref(self.RE).kT = self._params["kT"]
        # setting a volume is a side effect, sets the default volume of the
        # reaction ensemble as the volume of the cuboid simulation box. this
        # volume can be altered by the command "reaction ensemble volume
        # <volume>" if one wants to simulate e.g. in a system with constraint
        # (e.g. cuboid box with cylinder constraint, so that the particles are
        # only contained in the volume of the cylinder)
        if deref(self.RE).volume < 0:
            deref(self.RE).set_cuboid_reaction_ensemble_volume()
        deref(self.RE).exclusion_radius = self._params["exclusion_radius"]

    def remove_constraint(self):
        """
        Remove any previously defined constraint.
        Requires setting the volume using :meth:`set_volume`.

        """
        deref(self.RE).remove_constraint()

    def set_cylindrical_constraint_in_z_direction(self, center_x, center_y,
                                                  radius_of_cylinder):
        """
        Constrain the reaction moves within a cylinder defined by its axis
        passing through centres (:math:`x` and :math:`y`) and the radius.
        Requires setting the volume using :meth:`set_volume`.

        Parameters
        ----------
        center_x : :obj:`float`
            x coordinate of center of the cylinder.
        center_y : :obj:`float`
            y coordinate of center of the cylinder.
        radius_of_cylinder : :obj:`float`
            radius of the cylinder

        """
        deref(self.RE).set_cyl_constraint(
            center_x, center_y, radius_of_cylinder)

    def set_wall_constraints_in_z_direction(self, slab_start_z, slab_end_z):
        """
        Restrict the sampling area to a slab in z-direction. Requires setting
        the volume using :meth:`set_volume`. This constraint is necessary when
        working with :ref:`Electrostatic Layer Correction (ELC)`.

        Parameters
        ----------
        slab_start_z : :obj:`float`
            z coordinate of the bottom wall.
        slab_end_z : :obj:`float`
            z coordinate of the top wall.

        Examples
        --------

        >>> import espressomd
        >>> import espressomd.shapes
        >>> import espressomd.electrostatics
        >>> import espressomd.reaction_ensemble
        >>> import numpy as np
        >>> # setup a charged system
        >>> box_l = 20
        >>> elc_gap = 10
        >>> system = espressomd.System(box_l=[box_l, box_l, box_l + elc_gap])
        >>> system.time_step = 0.001
        >>> system.cell_system.skin = 0.4
        >>> types = {"HA": 0, "A-": 1, "H+": 2, "wall": 3}
        >>> charges = {types["HA"]: 0, types["A-"]: -1, types["H+"]: +1}
        >>> for i in range(10):
        ...     system.part.add(pos=np.random.random(3) * box_l, type=types["A-"], q=charges[types["A-"]])
        ...     system.part.add(pos=np.random.random(3) * box_l, type=types["H+"], q=charges[types["H+"]])
        >>> for particle_type in charges.keys():
        ...     system.non_bonded_inter[particle_type, types["wall"]].wca.set_params(epsilon=1.0, sigma=1.0)
        >>> # add ELC actor
        >>> p3m = espressomd.electrostatics.P3M(prefactor=1.0, accuracy=1e-2)
        >>> elc = espressomd.electrostatics.ELC(p3m_actor=p3m, maxPWerror=1.0, gap_size=elc_gap)
        >>> system.actors.add(elc)
        >>> # add constant pH method
        >>> RE = espressomd.reaction_ensemble.ConstantpHEnsemble(kT=1, exclusion_radius=1, seed=77)
        >>> RE.constant_pH = 2
        >>> RE.add_reaction(gamma=0.0088, reactant_types=[types["HA"]],
        ...                 product_types=[types["A-"], types["H+"]],
        ...                 default_charges=charges)
        >>> # add walls for the ELC gap
        >>> RE.set_wall_constraints_in_z_direction(0, box_l)
        >>> RE.set_volume(box_l**3)
        >>> system.constraints.add(shape=espressomd.shapes.Wall(dist=0, normal=[0, 0, 1]),
        ...                        particle_type=types["wall"])
        >>> system.constraints.add(shape=espressomd.shapes.Wall(dist=-box_l, normal=[0, 0, -1]),
        ...                        particle_type=types["wall"])


        """
        deref(self.RE).set_slab_constraint(slab_start_z, slab_end_z)

    def get_wall_constraints_in_z_direction(self):
        """
        Returns the restrictions of the sampling area in z-direction.

        """
        v = deref(self.RE).get_slab_constraint_parameters()
        return [v[0], v[1]]

    def set_volume(self, volume):
        """
        Set the volume to be used in the acceptance probability of the reaction
        ensemble. This can be useful when using constraints, if the relevant
        volume is different from the box volume. If not used the default volume
        which is used, is the box volume.

        """
        deref(self.RE).volume = volume

    def get_volume(self):
        """
        Get the volume to be used in the acceptance probability of the reaction
        ensemble.

        """
        return deref(self.RE).volume

    def get_acceptance_rate_configurational_moves(self):
        """
        Returns the acceptance rate for the configuration moves.

        """
        return deref(self.RE).get_acceptance_rate_configurational_moves()

    def get_acceptance_rate_reaction(self, reaction_id):
        """
        Returns the acceptance rate for the given reaction.

        """
        return deref(self.RE).reactions[reaction_id].get_acceptance_rate()

    def set_non_interacting_type(self, non_interacting_type):
        """
        Sets the particle type for non-interacting particles.
        Default value: 100.
        This is used to temporarily hide particles during a reaction trial
        move, if they are to be deleted after the move is accepted. Please
        change this value if you intend to use the type 100 for some other
        particle types with interactions, or if you need improved performance,
        as the default value of 100 causes some overhead.
        Please also note that particles
        in the current implementation of the Reaction Ensemble are only
        hidden with respect to Lennard-Jones and Coulomb interactions. Hiding
        of other interactions, for example a magnetic, needs to be implemented
        in the code.
        """
        deref(self.RE).non_interacting_type = non_interacting_type

    def get_non_interacting_type(self):
        """
        Returns the type which is used for hiding particles.

        """
        return deref(self.RE).non_interacting_type

    def add_reaction(self, *args, **kwargs):
        """
        Sets up a reaction in the forward and backward direction.

        Parameters
        ----------
        gamma : :obj:`float`
            Equilibrium constant of the reaction, :math:`\\gamma` (see the User
            guide, section 6.6 for the definition and further details).
        reactant_types : list of :obj:`int`
            List of particle types of reactants in the reaction.
        reactant_coefficients : list of :obj:`int`
            List of stoichiometric coefficients of the reactants in the same
            order as the list of their types.
        product_types : list of :obj:`int`
            List of particle types of products in the reaction.
        product_coefficients : list of :obj:`int`
            List of stoichiometric coefficients of products of the reaction in
            the same order as the list of their types
        default_charges : :obj:`dict`
            A dictionary of default charges for types that occur in the provided reaction.
        check_for_electroneutrality : :obj:`bool`
            Check for electroneutrality of the given reaction if ``True``.

        """
        self._params["check_for_electroneutrality"] = True
        for k in self._required_keys_add():
            if k not in kwargs:
                raise ValueError(
                    f"At least the following keys have to be given as keyword "
                    f"arguments: {self._required_keys()}, got {kwargs}")
            self._params[k] = kwargs[k]

        for k in self._valid_keys_add():
            if k in kwargs:
                self._params[k] = kwargs[k]
        self._check_lengths_of_arrays()
        self._validate_params_default_charge()
        self._set_params_in_es_core_add()

    def _valid_keys_add(self):
        return ("gamma", "reactant_types", "reactant_coefficients",
                "product_types", "product_coefficients", "default_charges",
                "check_for_electroneutrality")

    def _required_keys_add(self):
        return ["gamma", "reactant_types", "reactant_coefficients",
                "product_types", "product_coefficients", "default_charges"]

    def _check_lengths_of_arrays(self):
        if(len(self._params["reactant_types"]) != len(self._params["reactant_coefficients"])):
            raise ValueError(
                "Reactants: Number of types and coefficients have to be equal")
        if(len(self._params["product_types"]) != len(self._params["product_coefficients"])):
            raise ValueError(
                "Products: Number of types and coefficients have to be equal")

    def _set_params_in_es_core_add(self):
        cdef vector[int] reactant_types
        cdef vector[int] reactant_coefficients
        cdef vector[int] product_types
        cdef vector[int] product_coefficients
        for value in self._params["reactant_types"]:
            reactant_types.push_back(value)
        for value in self._params["reactant_coefficients"]:
            reactant_coefficients.push_back(value)
        for value in self._params["product_types"]:
            product_types.push_back(value)
        for value in self._params["product_coefficients"]:
            product_coefficients.push_back(value)

        # forward reaction
        deref(self.RE).add_reaction(
            self._params["gamma"], reactant_types, reactant_coefficients, product_types, product_coefficients)
        # backward reaction
        deref(self.RE).add_reaction(
            1.0 / self._params["gamma"], product_types, product_coefficients, reactant_types, reactant_coefficients)

        for key, value in self._params["default_charges"].items():
            deref(self.RE).charges_of_types[int(key)] = value
        deref(self.RE).check_reaction_method()

    def _validate_params_default_charge(self):
        if not isinstance(self._params["default_charges"], dict):
            raise ValueError(
                "No dictionary for relation between types and default charges provided.")
        # check electroneutrality of the provided reaction
        if self._params["check_for_electroneutrality"]:
            charges = np.array(list(self._params["default_charges"].values()))
            if np.count_nonzero(charges) == 0:
                # all particles have zero charge
                # no need to check electroneutrality
                return
            total_charge_change = 0.0
            for i in range(len(self._params["reactant_coefficients"])):
                type_here = self._params["reactant_types"][i]
                total_charge_change -= self._params["reactant_coefficients"][
                    i] * self._params["default_charges"][type_here]
            for j in range(len(self._params["product_coefficients"])):
                type_here = self._params["product_types"][j]
                total_charge_change += self._params["product_coefficients"][
                    j] * self._params["default_charges"][type_here]
            min_abs_nonzero_charge = np.min(
                np.abs(charges[np.nonzero(charges)[0]]))
            if abs(total_charge_change) / min_abs_nonzero_charge > 1e-10:
                raise ValueError("Reaction system is not charge neutral")

    def reaction(self, reaction_steps=1):
        """
        Performs randomly selected reactions.

        Parameters
        ----------
        reaction_steps : :obj:`int`, optional
            The number of reactions to be performed at once, defaults to 1.

        """
        deref(self.RE).do_reaction(int(reaction_steps))

    def displacement_mc_move_for_particles_of_type(self, type_mc,
                                                   particle_number_to_be_changed=1):
        """
        Performs a displacement Monte Carlo move for particles of given type.
        New positions of the displaced particles are chosen from the whole box
        with a uniform probability distribution. If there are multiple types,
        that are being moved in a simulation, they should be moved in a random
        order to avoid artefacts.

        Parameters
        ----------
        type_mc : :obj:`int`
            particle type which should be moved

        """

        deref(self.RE).do_global_mc_move_for_particles_of_type(
            type_mc, particle_number_to_be_changed)

    def get_status(self):
        """
        Returns the status of the reaction ensemble in a dictionary containing
        the used reactions, the used kT and the used exclusion radius.

        """
        deref(self.RE).check_reaction_method()
        reactions = []
        for single_reaction_i in range(deref(self.RE).reactions.size()):
            reactant_types = []
            for i in range(
                    deref(self.RE).reactions[single_reaction_i].reactant_types.size()):
                reactant_types.append(
                    deref(self.RE).reactions[single_reaction_i].reactant_types[i])
            reactant_coefficients = []
            for i in range(
                    deref(self.RE).reactions[single_reaction_i].reactant_types.size()):
                reactant_coefficients.append(
                    deref(self.RE).reactions[single_reaction_i].reactant_coefficients[i])

            product_types = []
            for i in range(
                    deref(self.RE).reactions[single_reaction_i].product_types.size()):
                product_types.append(
                    deref(self.RE).reactions[single_reaction_i].product_types[i])
            product_coefficients = []
            for i in range(
                    deref(self.RE).reactions[single_reaction_i].product_types.size()):
                product_coefficients.append(
                    deref(self.RE).reactions[single_reaction_i].product_coefficients[i])
            reaction = {"reactant_coefficients": reactant_coefficients,
                        "reactant_types": reactant_types,
                        "product_types": product_types,
                        "product_coefficients": product_coefficients,
                        "reactant_types": reactant_types,
                        "gamma": deref(self.RE).reactions[single_reaction_i].gamma}
            reactions.append(reaction)

        return {"reactions": reactions, "kT": deref(
            self.RE).kT, "exclusion_radius": deref(self.RE).exclusion_radius}

    def delete_particle(self, p_id):
        """
        Deletes the particle of the given p_id and makes sure that the particle
        range has no holes. This function has some restrictions, as e.g. bonds
        are not deleted. Therefore only apply this function to simple ions.

        """
        deref(self.RE).delete_particle(p_id)

    def change_reaction_constant(self, reaction_id, gamma):
        """
        Changes the reaction constant of a given reaction
        (for the forward and backward reaction).
        The ``reaction_id`` which is assigned to a reaction
        depends on the order in which :meth:`add_reaction` was called.
        The 0th reaction has ``reaction_id=0``, the next added
        reaction needs to be addressed with ``reaction_id=1``, etc.

        Parameters
        ----------
        reaction_id : :obj:`int`
            reaction_id
        gamma : :obj:`float`
            new reaction constant

        """
        reaction_id = int(reaction_id)
        if(reaction_id > deref(self.RE).reactions.size() / 2 - 1 or reaction_id < 0):
            raise ValueError(
                "You provided an invalid reaction_id, please provide a valid reaction_id")
        # for the forward reaction
        deref(self.RE).reactions[2 * reaction_id].gamma = gamma
        # for the backward reaction
        deref(self.RE).reactions[2 * reaction_id + 1].gamma = 1.0 / gamma

    def delete_reaction(self, reaction_id):
        """
        Delete a reaction from the set of used reactions
        (the forward and backward reaction).
        The ``reaction_id`` which is assigned to a reaction
        depends on the order in which :meth:`add_reaction` was called.
        The 0th reaction has ``reaction_id=0``, the next added
        reaction needs to be addressed with ``reaction_id=1``, etc.
        After the deletion of a reaction subsequent reactions
        take the ``reaction_id`` of the deleted reaction.

        Parameters
        ----------
        reaction_id : :obj:`int`
            reaction_id

        """
        reaction_id = int(reaction_id)
        if(reaction_id > deref(self.RE).reactions.size() / 2 - 1 or reaction_id < 0):
            raise ValueError(
                "You provided an invalid reaction_id, please provide a valid reaction_id")
        deref(self.RE).delete_reaction(2 * reaction_id + 1)
        deref(self.RE).delete_reaction(2 * reaction_id)

cdef class ReactionEnsemble(ReactionAlgorithm):
    """
    This class implements the Reaction Ensemble.
    """

    cdef unique_ptr[CReactionEnsemble] REptr

    def __init__(self, *args, **kwargs):
        self._params = {"kT": 1,
                        "exclusion_radius": 0}
        for k in self._required_keys():
            if k not in kwargs:
                raise ValueError(
                    f"At least the following keys have to be given as keyword "
                    f"arguments: {self._required_keys()}, got {kwargs}")
            self._params[k] = kwargs[k]

        self.REptr.reset(new CReactionEnsemble(int(self._params["seed"])))
        self.RE = <CReactionAlgorithm * > self.REptr.get()

        for k in kwargs:
            if k in self._valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError(f"{k} is not a valid key")

        self._set_params_in_es_core()

cdef class ConstantpHEnsemble(ReactionAlgorithm):
    """
    This class implements the constant pH Ensemble.

    When adding an acid-base reaction, the acid and base particle types
    are always assumed to be at index 0 of the lists passed to arguments
    ``reactant_types`` and ``product_types``.

    """
    cdef unique_ptr[CConstantpHEnsemble] constpHptr

    def __init__(self, *args, **kwargs):
        self._params = {"kT": 1, "exclusion_radius": 0}
        for k in self._required_keys():
            if k not in kwargs:
                raise ValueError(
                    f"At least the following keys have to be given as keyword "
                    f"arguments: {self._required_keys()}, got {kwargs}")
            self._params[k] = kwargs[k]

        self.constpHptr.reset(new CConstantpHEnsemble(int(self._params["seed"])))
        self.RE = <CReactionAlgorithm * > self.constpHptr.get()

        for k in kwargs:
            if k in self._valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError(f"{k} is not a valid key")

        self._set_params_in_es_core()

    def add_reaction(self, *args, **kwargs):
        warn_msg = (
            "arguments 'reactant_coefficients' and 'product_coefficients' "
            "are deprecated and are no longer necessary for the constant pH "
            "ensemble. They are kept for backward compatibility but might "
            "be deleted in future versions.")
        err_msg = ("All product and reactant coefficients must equal one in "
                   "the constant pH method as implemented in ESPResSo.")
        warn_user = False

        if "reactant_coefficients" in kwargs:
            if kwargs["reactant_coefficients"][0] != 1:
                raise ValueError(err_msg)
            else:
                warn_user = True
        else:
            kwargs["reactant_coefficients"] = [1]

        if "product_coefficients" in kwargs:
            if kwargs["product_coefficients"][0] != 1 or kwargs["product_coefficients"][1] != 1:
                raise ValueError(err_msg)
            else:
                warn_user = True
        else:
            kwargs["product_coefficients"] = [1, 1]

        if warn_user:
            warnings.warn(warn_msg, FutureWarning)

        if(len(kwargs["product_types"]) != 2 or len(kwargs["reactant_types"]) != 1):
            raise ValueError(
                "The constant pH method is only implemented for reactions "
                "with two product types and one adduct type.")

        super().add_reaction(*args, **kwargs)

    property constant_pH:
        """
        Sets the input pH for the constant pH ensemble method.

        """

        def __set__(self, double pH):
            """
            Sets the pH that the method assumes for the implicit pH bath.

            """

            deref(self.constpHptr).m_constant_pH = pH

cdef class WidomInsertion(ReactionAlgorithm):
    """
    This class implements the Widom insertion method in the canonical ensemble
    for homogeneous systems, where the excess chemical potential is not
    depending on the location.

    """

    cdef unique_ptr[CWidomInsertion] WidomInsertionPtr

    def _required_keys(self):
        return ("kT", "seed")

    def _valid_keys(self):
        return ("kT", "seed")

    def _valid_keys_add(self):
        return ("reactant_types", "reactant_coefficients", "product_types",
                "product_coefficients", "default_charges",
                "check_for_electroneutrality")

    def _required_keys_add(self):
        return ("reactant_types", "reactant_coefficients",
                "product_types", "product_coefficients", "default_charges")

    def __init__(self, *args, **kwargs):
        self._params = {"kT": 1}
        for k in self._required_keys():
            if k not in kwargs:
                raise ValueError(
                    f"At least the following keys have to be given as keyword "
                    f"arguments: {self._required_keys()}, got {kwargs}")
            self._params[k] = kwargs[k]
        self._params["exclusion_radius"] = 0.0  # not used by this method
        self._params["gamma"] = 1.0  # not used by this method

        self.WidomInsertionPtr.reset(new CWidomInsertion(int(self._params["seed"])))
        self.RE = <CReactionAlgorithm * > self.WidomInsertionPtr.get()
        for k in kwargs:
            if k in self._valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError(f"{k} is not a valid key")

        self._set_params_in_es_core()

    def calculate_particle_insertion_potential_energy(self, reaction_id=0):
        """
        Measures the potential energy when particles are inserted in the
        system following the reaction  provided ``reaction_id``. Please
        define the insertion moves first by calling the method
        :meth:`~ReactionAlgorithm.add_reaction` (with only product types
        specified).

        Note that although this function does not provide directly
        the chemical potential, it can be used to calculate it.
        For an example of such an application please check
        :file:`/samples/widom_insertion.py`.
        """
        # make inverse widom scheme (deletion of particles) inaccessible.
        # The deletion reactions are the odd reaction_ids
        if(reaction_id < 0 or reaction_id > (deref(self.WidomInsertionPtr).reactions.size() + 1) / 2):
            raise ValueError("This reaction is not present")
        return deref(self.WidomInsertionPtr).calculate_particle_insertion_potential_energy(
            deref(self.WidomInsertionPtr).reactions[int(2 * reaction_id)])

    def calculate_excess_chemical_potential(
            self, particle_insertion_potential_energy_samples, N_blocks=16):
        """
        Given a set of samples of the particle insertion potential energy,
        calculates the excess chemical potential and its statistical error.

        Parameters
        ----------
        particle_insertion_potential_energy_samples : array_like of :obj:`float`
            Samples of the particle insertion potential energy.
        N_blocks : :obj:`int`, optional
            Number of bins for binning analysis.

        Returns
        -------
        mean : :obj:`float`
            Mean excess chemical potential.
        error : :obj:`float`
            Standard error of the mean.
        """

        def do_block_analysis(samples, N_blocks=16):
            """
            Performs a binning analysis of samples.
            Divides the samples in ``N_blocks`` equispaced blocks
            and returns the mean and its uncertainty
            """
            size_block = int(len(samples) / N_blocks)
            block_list = []
            for block in range(N_blocks):
                block_list.append(
                    np.mean(samples[block * size_block:(block + 1) * size_block]))

            sample_mean = np.mean(block_list)
            sample_std = np.std(block_list, ddof=1)
            sample_uncertainty = sample_std / np.sqrt(N_blocks)

            return sample_mean, sample_uncertainty

        gamma_samples = np.exp(-1.0 * np.array(
            particle_insertion_potential_energy_samples) / self._params["kT"])

        gamma_mean, gamma_std = do_block_analysis(
            samples=gamma_samples, N_blocks=N_blocks)

        mu_ex_mean = -1 * np.log(gamma_mean) * self._params["kT"]

        # full propagation of error
        mu_ex_Delta = 0.5 * self._params["kT"] * abs(-np.log(gamma_mean + gamma_std) -
                                                     (-np.log(gamma_mean - gamma_std)))

        return mu_ex_mean, mu_ex_Delta
