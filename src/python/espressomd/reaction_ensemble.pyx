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


class WangLandauHasConverged(Exception):
    pass


cdef class ReactionAlgorithm:
    """

    This class provides the base class for Reaction Algorithms like the Reaction
    Ensemble algorithm, the Wang-Landau Reaction Ensemble algorithm and the
    constant pH method. Initialize the reaction algorithm by setting the
    standard pressure, temperature, and the exclusion radius.

    Note: When creating particles the velocities of the new particles are set
    according the Maxwell-Boltzmann distribution. In this step the mass of the
    new particle is assumed to equal 1.


    Parameters
    ----------
    temperature : :obj:`float`
        The temperature at which the reaction is performed.
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
        return "temperature", "exclusion_radius", "seed"

    def _required_keys(self):
        return "temperature", "exclusion_radius", "seed"

    def _set_params_in_es_core(self):
        deref(self.RE).temperature = self._params["temperature"]
        # setting a volume is a side effect, sets the default volume of the
        # reaction ensemble as the volume of the cuboid simulation box. this
        # volume can be altered by the command "reaction ensemble volume
        # <volume>" if one wants to simulate e.g. in a system with constraint
        # (e.g. cuboid box with cylinder constraint, so that the particles are
        # only contained in the volume of the cylinder)
        if(deref(self.RE).volume < 0):
            deref(self.RE).set_cuboid_reaction_ensemble_volume()
        deref(self.RE).exclusion_radius = self._params["exclusion_radius"]

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
        deref(self.RE).cyl_x = center_x
        deref(self.RE).cyl_y = center_y
        deref(self.RE).cyl_radius = radius_of_cylinder
        deref(self.RE).box_is_cylindric_around_z_axis = True

    def set_wall_constraints_in_z_direction(self, slab_start_z, slab_end_z):
        """
        Restrict the sampling area to a slab in z-direction. Requires setting
        the volume using :meth:`set_volume`.

        """
        deref(self.RE).slab_start_z = slab_start_z
        deref(self.RE).slab_end_z = slab_end_z
        deref(self.RE).box_has_wall_constraints = True

    def get_wall_constraints_in_z_direction(self):
        """
        Returns the restrictions of the sampling area in z-direction.

        """
        return deref(self.RE).slab_start_z, deref(self.RE).slab_end_z

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
        particle types with interactions. Please also note that particles
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
                raise ValueError("At least the following keys have to be given as keyword arguments: " +
                                 self._required_keys_add().__str__() + " got " + kwargs.__str__())
            self._params[k] = kwargs[k]

        for k in self._valid_keys_add():
            try:
                self._params[k] = kwargs[k]
            except BaseException:
                pass
        self._check_lengths_of_arrays()
        self._validate_params_default_charge()
        self._set_params_in_es_core_add()

    def _valid_keys_add(self):
        return "gamma", "reactant_types", "reactant_coefficients", "product_types", "product_coefficients", "default_charges", "check_for_electroneutrality"

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
        for i in range(len(self._params["reactant_types"])):
            reactant_types.push_back(self._params["reactant_types"][i])
        cdef vector[int] reactant_coefficients
        for i in range(len(self._params["reactant_coefficients"])):
            reactant_coefficients.push_back(
                self._params["reactant_coefficients"][i])
        cdef vector[int] product_types
        for i in range(len(self._params["product_types"])):
            product_types.push_back(self._params["product_types"][i])
        cdef vector[int] product_coefficients
        for i in range(len(self._params["product_coefficients"])):
            product_coefficients.push_back(
                self._params["product_coefficients"][i])
        deref(self.RE).add_reaction(
            self._params["gamma"], reactant_types, reactant_coefficients, product_types, product_coefficients)
        deref(self.RE).add_reaction(
            1.0 / self._params["gamma"], product_types, product_coefficients, reactant_types, reactant_coefficients)

        for key in self._params["default_charges"]:  # the keys are the types
            deref(self.RE).charges_of_types[
                int(key)] = self._params["default_charges"][key]
        deref(self.RE).check_reaction_ensemble()

    def _validate_params_default_charge(self):
        if(isinstance(self._params["default_charges"], dict) == False):
            raise ValueError(
                "No dictionary for relation between types and default charges provided.")
        # check electroneutrality of the provided reaction
        if(self._params["check_for_electroneutrality"]):
            charges = np.array(list(self._params["default_charges"].values()))
            if(np.count_nonzero(charges) == 0):
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

        use_wang_landau = False
        deref(self.RE).do_global_mc_move_for_particles_of_type(
            type_mc, particle_number_to_be_changed, use_wang_landau)

    def get_status(self):
        """
        Returns the status of the reaction ensemble in a dictionary containing
        the used reactions, the used temperature and the used exclusion radius.

        """
        deref(self.RE).check_reaction_ensemble()
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

        return {"reactions": reactions, "temperature": deref(
            self.RE).temperature, "exclusion_radius": deref(self.RE).exclusion_radius}

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
        self._params = {"temperature": 1,
                        "exclusion_radius": 0}
        for k in self._required_keys():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self._required_keys().__str__() + " got " + kwargs.__str__())
            self._params[k] = kwargs[k]

        self.REptr.reset(new CReactionEnsemble(int(self._params["seed"])))
        self.RE = <CReactionAlgorithm * > self.REptr.get()

        for k in kwargs:
            if k in self._valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a valid key" % k)

        self._set_params_in_es_core()

cdef class ConstantpHEnsemble(ReactionAlgorithm):
    cdef unique_ptr[CConstantpHEnsemble] constpHptr

    def __init__(self, *args, **kwargs):
        self._params = {"temperature": 1,
                        "exclusion_radius": 0}
        for k in self._required_keys():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self._required_keys().__str__() + " got " + kwargs.__str__())
            self._params[k] = kwargs[k]

        self.constpHptr.reset(new CConstantpHEnsemble(int(self._params["seed"])))
        self.RE = <CReactionAlgorithm * > self.constpHptr.get()

        for k in kwargs:
            if k in self._valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a valid key" % k)

        self._set_params_in_es_core()

    def add_reaction(self, *args, **kwargs):
        if(len(kwargs["product_types"]) != 2 or len(kwargs["reactant_types"]) != 1):
            raise ValueError(
                "The constant pH method is only implemented for reactions with two product types and one adduct type.")
        if(kwargs["reactant_coefficients"][0] != 1 or kwargs["product_coefficients"][0] != 1 or kwargs["product_coefficients"][1] != 1):
            raise ValueError(
                "All product and reactant coefficients must equal one in the constant pH method as implemented in ESPResSo.")
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

cdef class WangLandauReactionEnsemble(ReactionAlgorithm):
    """
    This Class implements the Wang-Landau Reaction Ensemble.
    """

    cdef unique_ptr[CWangLandauReactionEnsemble] WLRptr

    def __init__(self, *args, **kwargs):
        self._params = {"temperature": 1,
                        "exclusion_radius": 0}
        for k in self._required_keys():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self._required_keys().__str__() + " got " + kwargs.__str__())
            self._params[k] = kwargs[k]
        for k in kwargs:
            if k in self._valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a valid key" % k)

        self.WLRptr.reset(new CWangLandauReactionEnsemble(int(self._params["seed"])))
        self.RE = <CReactionAlgorithm * > self.WLRptr.get()

        self._set_params_in_es_core()

    def reaction(self, reaction_steps=1):
        """
        Performs reaction_steps reactions. Sets the number of reaction steps
        which are performed at once. Do not use too many reaction steps
        consecutively without having conformation-changing steps in-between
        (especially important for the Wang-Landau reaction ensemble). Providing
        a number for the parameter reaction steps reduces the need for the
        interpreter to be called between consecutive reactions.

        """
        status_wang_landau = deref(
            self.WLRptr).do_reaction(int(reaction_steps))
        if(status_wang_landau < 0):
            raise WangLandauHasConverged(
                "The Wang-Landau algorithm has converged.")

    def add_collective_variable_degree_of_association(self, *args, **kwargs):
        """
        Adds the degree of association as a collective variable (reaction coordinate) for the Wang-Landau Reaction Ensemble.
        Several collective variables can be set simultaneously.

        Parameters
        ----------
        associated_type : :obj:`int`
            Particle type of the associated state of the reacting species.
        min : :obj:`float`
            Minimum value of the collective variable.
        max : :obj:`float`
            Maximum value of the collective variable.
        corresponding_acid_types : list of :obj:`int`
            List of the types of the version of the species.

        """
        for k in kwargs:
            if k in self._valid_keys_add_collective_variable_degree_of_association():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a valid key" % k)

        for k in self._required_keys_add_collective_variable_degree_of_association():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self._required_keys_add_collective_variable_degree_of_association().__str__() + " got " + kwargs.__str__())
            self._params[k] = kwargs[k]

        cdef vector[int] _corresponding_acid_types
        for i in range(len(self._params["corresponding_acid_types"])):
            _corresponding_acid_types.push_back(
                self._params["corresponding_acid_types"][i])
        deref(self.WLRptr).add_new_CV_degree_of_association(
            self._params["associated_type"], self._params["min"], self._params["max"], _corresponding_acid_types)

    def _valid_keys_add_collective_variable_degree_of_association(self):
        return "associated_type", "min", "max", "corresponding_acid_types"

    def _required_keys_add_collective_variable_degree_of_association(self):
        return "associated_type", "min", "max", "corresponding_acid_types"

    def add_collective_variable_potential_energy(self, *args, **kwargs):
        """
        Adds the potential energy as a collective variable (reaction coordinate) for the Wang-Landau Reaction Ensemble.
        Several collective variables can be set simultaneously.

        Parameters
        ----------
        filename : :obj:`str`
            Filename of the energy boundary file which provides the
            potential energy boundaries (min E_pot, max E_pot) tabulated
            for all degrees of association. Make sure to only list the
            degrees of association which are used by the degree of
            association collective variable within this file. The energy
            boundary file can be created in a preliminary energy run. By
            the help of the functions
            :meth:`update_maximum_and_minimum_energies_at_current_state`
            and :meth:`write_out_preliminary_energy_run_results`. This
            file has to be obtained before being able to run a
            simulation with the energy as collective variable.
        delta : :obj:`float`
            Provides the discretization of the potential energy range. Only
            for small enough delta the results of the energy reweighted
            averages are correct. If delta is chosen too big there are
            discretization errors in the numerical integration which occurs
            during the energy reweighting process.

        """
        for k in kwargs:
            if k in self._valid_keys_add_collective_variable_potential_energy():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a valid key" % k)

            for k in self._required_keys_add_collective_variable_potential_energy():
                if k not in kwargs:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self._required_keys_add_collective_variable_degree_of_association().__str__() + " got " + kwargs.__str__())
                self._params[k] = kwargs[k]
        filname_potential_energy_boundaries_file = self._params[
            "filename"].encode("utf-8")
        deref(self.WLRptr).add_new_CV_potential_energy(
            filname_potential_energy_boundaries_file, self._params["delta"])

    def _valid_keys_add_collective_variable_potential_energy(self):
        return "filename", "delta"

    def _required_keys_add_collective_variable_potential_energy(self):
        return "filename", "delta"

    def set_wang_landau_parameters(self, *args, **kwargs):
        """
        Sets the final Wang-Landau parameter.

        Parameters
        ----------
        final_wang_landau_parameter : :obj:`float`
            Sets the final Wang-Landau parameter, which is the Wang-Landau
            parameter after which the simulation should stop.
        full_path_to_output_filename : :obj:`str`
            Sets the path to the output file of the Wang-Landau algorithm which
            contains the Wang-Landau potential
        do_not_sample_reaction_partition_function : :obj:`bool`
            Avoids sampling the Reaction ensemble partition function in the
            Wang-Landau algorithm. Therefore this option makes all degrees of
            association equally probable. This option may be used in the
            sweeping mode of the reaction ensemble, since the reaction ensemble
            partition function can be later added analytically.

        """
        for k in kwargs:
            if k in self._valid_keys_set_wang_landau_parameters():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a valid key" % k)

        deref(self.WLRptr).final_wang_landau_parameter = self._params[
            "final_wang_landau_parameter"]
        deref(self.WLRptr).output_filename = self._params[
            "full_path_to_output_filename"].encode("utf-8")
        deref(self.WLRptr).do_not_sample_reaction_partition_function = self._params[
            "do_not_sample_reaction_partition_function"]

    def _valid_keys_set_wang_landau_parameters(self):
        return "final_wang_landau_parameter", "full_path_to_output_filename", "do_not_sample_reaction_partition_function"

    def load_wang_landau_checkpoint(self):
        """
        Loads the dumped Wang-Landau potential file.

        """
        checkpoint_name = "checkpoint".encode("utf-8")
        deref(self.WLRptr).load_wang_landau_checkpoint(checkpoint_name)

    def write_wang_landau_checkpoint(self):
        """
        Dumps the Wang-Landau potential to a checkpoint file. Can be used to
        checkpoint the Wang-Landau histogram, potential, parameter and the
        number of executed trial moves.

        """
        checkpoint_name = "checkpoint".encode("utf-8")
        deref(self.WLRptr).write_wang_landau_checkpoint(checkpoint_name)

    def update_maximum_and_minimum_energies_at_current_state(self):
        """
        Records the minimum and maximum potential energy as a function of the
        degree of association in a preliminary Wang-Landau reaction ensemble
        simulation where the acceptance probability includes the factor
        :math:`\\exp(-\\beta \\Delta E_{pot})`. The minimal and maximal
        potential energies which occur in the system are needed for the energy
        reweighting simulations where the factor :math:`\\exp(-\\beta \\Delta E_{pot})`
        is not included in the acceptance probability in
        order to avoid choosing the wrong potential energy boundaries.

        """
        self.WLRptr.get(
        ).update_maximum_and_minimum_energies_at_current_state()

    def write_out_preliminary_energy_run_results(self):
        """
        This writes out the minimum and maximum potential energy as a function
        of the degree of association to a file. It requires that previously
        :meth:`update_maximum_and_minimum_energies_at_current_state` was used.

        """
        filename = "preliminary_energy_run_results".encode("utf-8")
        deref(self.WLRptr).write_out_preliminary_energy_run_results(filename)

    def write_wang_landau_results_to_file(self, filename):
        """
        This writes out the Wang-Landau potential as a function of the used
        collective variables.

        """
        deref(self.WLRptr).write_wang_landau_results_to_file(
            filename.encode("utf-8"))

    def displacement_mc_move_for_particles_of_type(self, type_mc,
                                                   particle_number_to_be_changed=1):
        """
        Performs an MC (Monte Carlo) move for ``particle_number_to_be_changed``
        particle of type ``type_mc``. Positions for the particles are drawn
        uniformly and randomly within the box. The command takes into account
        the Wang-Landau terms in the acceptance probability.
        If there are multiple types, that need to be moved, make sure to move
        them in a random order to avoid artefacts. For the Wang-Landau algorithm
        in the case of energy reweighting you would also need to move the
        monomers of the polymer with special moves for the MC part. Those
        polymer configuration-changing moves need to be implemented in the
        case of using Wang-Landau with energy reweighting and a polymer in the
        system. Polymer configuration-changing moves had been implemented
        before but were removed from ESPResSo.

        """
        use_wang_landau = True
        deref(self.WLRptr).do_global_mc_move_for_particles_of_type(
            type_mc, particle_number_to_be_changed, use_wang_landau)


cdef class WidomInsertion(ReactionAlgorithm):
    """
    This class implements the Widom insertion method in the canonical ensemble
    for homogeneous systems, where the excess chemical potential is not
    depending on the location.

    """

    cdef unique_ptr[CWidomInsertion] WidomInsertionPtr

    def _required_keys(self):
        return "temperature", "seed"

    def _valid_keys(self):
        return "temperature", "seed"

    def _valid_keys_add(self):
        return "reactant_types", "reactant_coefficients", "product_types", "product_coefficients", "default_charges", "check_for_electroneutrality"

    def _required_keys_add(self):
        return ["reactant_types", "reactant_coefficients",
                "product_types", "product_coefficients", "default_charges"]

    def __init__(self, *args, **kwargs):
        self._params = {"temperature": 1}
        for k in self._required_keys():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self._required_keys().__str__() + " got " + kwargs.__str__())
            self._params[k] = kwargs[k]
        self._params[
            "exclusion_radius"] = 0.0  # this is not used by the widom insertion method
        self._params[
            "gamma"] = 1.0  # this is not used by the widom insertion method

        self.WidomInsertionPtr.reset(new CWidomInsertion(int(self._params["seed"])))
        self.RE = <CReactionAlgorithm * > self.WidomInsertionPtr.get()
        for k in kwargs:
            if k in self._valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a valid key" % k)

        self._set_params_in_es_core()

    def measure_excess_chemical_potential(self, reaction_id=0):
        """
        Measures the excess chemical potential in a homogeneous system for
        the provided ``reaction_id``. Please define the insertion moves
        first by calling the method :meth:`~ReactionAlgorithm.add_reaction`
        (with only product types specified).
        Returns the excess chemical potential and the standard error for
        the excess chemical potential. The error estimate assumes that
        your samples are uncorrelated.

        """
        if(reaction_id < 0 or reaction_id > (deref(self.WidomInsertionPtr).reactions.size() + 1) / 2):  # make inverse widom scheme (deletion of particles) inaccessible
            raise ValueError("This reaction is not present")
        return deref(self.WidomInsertionPtr).measure_excess_chemical_potential(
            int(2 * reaction_id))  # make inverse widom scheme (deletion of particles) inaccessible. The deletion reactions are the odd reaction_ids
