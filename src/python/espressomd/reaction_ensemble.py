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

import numpy as np
import warnings
from .script_interface import ScriptInterfaceHelper, script_interface_register
from . import utils


@script_interface_register
class SingleReaction(ScriptInterfaceHelper):
    _so_name = "ReactionMethods::SingleReaction"
    _so_creation_policy = "LOCAL"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if not 'sip' in kwargs:
            utils.check_valid_keys(self.valid_keys(), kwargs.keys())

    def valid_keys(self):
        return ("reactant_types", "reactant_coefficients",
                "product_types", "product_coefficients", "gamma")

    def required_keys(self):
        return ("reactant_types", "reactant_coefficients",
                "product_types", "product_coefficients", "gamma")

    def make_backward_reaction(self):
        return SingleReaction(
            gamma=1. / self.gamma, reactant_types=self.product_types,
            reactant_coefficients=self.product_coefficients,
            product_types=self.reactant_types,
            product_coefficients=self.reactant_coefficients)


@script_interface_register
class ReactionAlgorithm(ScriptInterfaceHelper):
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

    Methods
    -------
    remove_constraint()
        Remove any previously defined constraint.
        Requires setting the volume using :meth:`set_volume`.

    set_cylindrical_constraint_in_z_direction()
        Constrain the reaction moves within a cylinder aligned on the z-axis.
        Requires setting the volume using :meth:`set_volume`.

        Parameters
        ----------
        center_x : :obj:`float`
            x coordinate of center of the cylinder.
        center_y : :obj:`float`
            y coordinate of center of the cylinder.
        radius_of_cylinder : :obj:`float`
            radius of the cylinder

    set_wall_constraints_in_z_direction()
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

    get_wall_constraints_in_z_direction()
        Returns the restrictions of the sampling area in z-direction.

    set_volume()
        Set the volume to be used in the acceptance probability of the reaction
        ensemble. This can be useful when using constraints, if the relevant
        volume is different from the box volume. If not used the default volume
        which is used, is the box volume.

        Parameters
        ----------
        volume : :obj:`float`
            Volume of the system in simulation units

    get_volume()
        Get the volume to be used in the acceptance probability of the reaction
        ensemble.

    get_acceptance_rate_configurational_moves():
        Returns the acceptance rate for the configuration moves.

    get_acceptance_rate_reaction()
        Returns the acceptance rate for the given reaction.

        Parameters
        ----------
        reaction_id : :obj:`int`
            Reaction id

    set_non_interacting_type()
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

        Parameters
        ----------
        type : :obj:`int`
            Particle type for the hidden particles

    get_non_interacting_type()
        Returns the type which is used for hiding particle

    reaction()
        Performs randomly selected reactions.

        Parameters
        ----------
        reaction_steps : :obj:`int`, optional
            The number of reactions to be performed at once, defaults to 1.

    displacement_mc_move_for_particles_of_type()
        Performs a displacement Monte Carlo move for particles of given type.
        New positions of the displaced particles are chosen from the whole box
        with a uniform probability distribution. If there are multiple types,
        that are being moved in a simulation, they should be moved in a random
        order to avoid artefacts.

        Parameters
        ----------
        type_mc : :obj:`int`
            Particle type which should be moved
        particle_number_to_be_changed : :obj:`int`
            Number of particles to move, defaults to 1.

    delete_particle()
        Deletes the particle of the given p_id and makes sure that the particle
        range has no holes. This function has some restrictions, as e.g. bonds
        are not deleted. Therefore only apply this function to simple ions.

        Parameters
        ----------
        p_id : :obj:`int`
            Id of the particle to be deleted.

    change_reaction_constant()
        Changes the reaction constant of a given reaction
        (for both the forward and backward reactions).
        The ``reaction_id`` which is assigned to a reaction
        depends on the order in which :meth:`add_reaction` was called.
        The 0th reaction has ``reaction_id=0``, the next added
        reaction needs to be addressed with ``reaction_id=1``, etc.

        Parameters
        ----------
        reaction_id : :obj:`int`
            Reaction id
        gamma : :obj:`float`
            New reaction constant

    delete_reaction()
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
            Reaction id

    """

    _so_name = "ReactionMethods::ReactionAlgorithm"
    _so_creation_policy = "LOCAL"
    _so_bind_methods = ("remove_constraint",
                        "get_wall_constraints_in_z_direction",
                        "set_wall_constraints_in_z_direction",
                        "set_cylindrical_constraint_in_z_direction",
                        "set_volume",
                        "get_volume",
                        "get_acceptance_rate_reaction",
                        "set_non_interacting_type",
                        "get_non_interacting_type",
                        "reaction",
                        "displacement_mc_move_for_particles_of_type",
                        "check_reaction_method",
                        "change_reaction_constant",
                        "delete_reaction",
                        "delete_particle",
                        )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if not 'sip' in kwargs:
            utils.check_valid_keys(self.valid_keys(), kwargs.keys())

    def valid_keys(self):
        return {"kT", "exclusion_radius", "seed"}

    def required_keys(self):
        return {"kT", "exclusion_radius", "seed"}

    def add_reaction(self, **kwargs):
        """
        Sets up a reaction in the forward and backward direction.

        Parameters
        ----------
        gamma : :obj:`float`
            Equilibrium constant of the reaction in simulation units, :math:`\\gamma` (see the User
            guide, section 16.7.1. for the definition and further details).
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
            A dictionary of default charges for types that occur
            in the provided reaction.
        check_for_electroneutrality : :obj:`bool`
            Check for electroneutrality of the given reaction.
            Default is ``True``.

        """
        default_charges = kwargs.pop("default_charges")
        neutrality_check = kwargs.pop("check_for_electroneutrality", True)
        forward_reaction = SingleReaction(**kwargs)
        backward_reaction = forward_reaction.make_backward_reaction()
        if neutrality_check:
            self._check_charge_neutrality(
                type2charge=default_charges,
                reaction=forward_reaction)

        self.call_method("add_reaction", reaction=forward_reaction)
        self.call_method("add_reaction", reaction=backward_reaction)

        for ptype, charge in default_charges.items():
            self.call_method("set_charge_of_type", type=ptype, charge=charge)
        self.call_method("check_reaction_method")

    def _check_charge_neutrality(self, type2charge, reaction):
        if not isinstance(type2charge, dict):
            raise ValueError(
                "No dictionary for relation between types and default charges provided.")
        charges = np.array(list(type2charge.values()))
        if np.count_nonzero(charges) == 0:
            # all particles have zero charge
            # no need to check electroneutrality
            return
        # calculate net change of electrical charge for the reaction
        net_charge_change = 0.0
        for coef, ptype in zip(
                reaction.reactant_coefficients, reaction.reactant_types):
            net_charge_change -= coef * type2charge[ptype]
        for coef, ptype in zip(
                reaction.product_coefficients, reaction.product_types):
            net_charge_change += coef * type2charge[ptype]
        min_abs_nonzero_charge = np.min(
            np.abs(charges[np.nonzero(charges)[0]]))
        if abs(net_charge_change) / min_abs_nonzero_charge > 1e-10:
            raise ValueError("Reaction system is not charge neutral")

    def get_status(self):
        """
        Returns the status of the reaction ensemble in a dictionary containing
        the used reactions, the used kT and the used exclusion radius.

        """

        self.call_method("check_reaction_method")
        reactions_list = []

        for core_reaction in self.reactions:
            reaction = {"reactant_coefficients": core_reaction.reactant_coefficients,
                        "reactant_types": core_reaction.reactant_types,
                        "product_types": core_reaction.product_types,
                        "product_coefficients": core_reaction.product_coefficients,
                        "reactant_types": core_reaction.reactant_types,
                        "gamma": core_reaction.gamma}
            reactions_list.append(reaction)

        return {"reactions": reactions_list, "kT": self.kT, 
                "exclusion_radius": self.exclusion_radius}


@script_interface_register
class ReactionEnsemble(ReactionAlgorithm):
    """
    This class implements the Reaction Ensemble.
    """

    _so_name = "ReactionMethods::ReactionEnsemble"
    _so_creation_policy = "LOCAL"


@script_interface_register
class ConstantpHEnsemble(ReactionAlgorithm):
    """
    This class implements the constant pH Ensemble.

    When adding an acid-base reaction, the acid and base particle types
    are always assumed to be at index 0 of the lists passed to arguments
    ``reactant_types`` and ``product_types``.

    """
    _so_name = "ReactionMethods::ConstantpHEnsemble"
    _so_creation_policy = "LOCAL"

    def valid_keys(self):
        return {"kT", "exclusion_radius", "seed", "constant_pH"}

    def required_keys(self):
        return {"kT", "exclusion_radius", "seed", "constant_pH"}

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


@script_interface_register
class WidomInsertion(ReactionAlgorithm):
    """
    This class implements the Widom insertion method in the canonical ensemble
    for homogeneous systems, where the excess chemical potential is not
    depending on the location.

    """

    _so_name = "ReactionMethods::WidomInsertion"
    _so_creation_policy = "LOCAL"

    def required_keys(self):
        return {"kT", "seed"}

    def valid_keys(self):
        return {"kT", "seed"}

    def add_reaction(self, **kwargs):
        kwargs['gamma'] = 1.
        super().add_reaction(**kwargs)

    def calculate_particle_insertion_potential_energy(self, **kwargs):
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

        return self.call_method(
            "calculate_particle_insertion_potential_energy", **kwargs)

    def calculate_excess_chemical_potential(
            self, **kwargs):
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

        def do_block_analysis(samples, N_blocks):
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

        kT = self.kT

        gamma_samples = np.exp(-1.0 * np.array(
            kwargs["particle_insertion_potential_energy_samples"]) / kT)

        gamma_mean, gamma_std = do_block_analysis(
            samples=gamma_samples, N_blocks=kwargs.get("N_blocks", 16))

        mu_ex_mean = -kT * np.log(gamma_mean)

        # full propagation of error
        mu_ex_Delta = 0.5 * kT * abs(-np.log(gamma_mean + gamma_std) -
                                     (-np.log(gamma_mean - gamma_std)))

        return mu_ex_mean, mu_ex_Delta
