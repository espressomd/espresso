# Copyright (C) 2010-2022 The ESPResSo project
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
import math
import sys
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .code_features import has_features
from . import utils


@script_interface_register
class SingleReaction(ScriptInterfaceHelper):
    _so_name = "ReactionMethods::SingleReaction"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("get_acceptance_rate",)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if not 'sip' in kwargs:
            utils.check_valid_keys(self.valid_keys(), kwargs.keys())
        self.reaction_types = tuple(
            sorted(set(tuple(self.reactant_types) + tuple(self.product_types))))

    def valid_keys(self):
        return self.required_keys()

    def required_keys(self):
        return {"reactant_types", "reactant_coefficients", "gamma",
                "product_types", "product_coefficients"}

    def make_backward_reaction(self):
        return SingleReaction(
            gamma=1. / self.gamma, reactant_types=self.product_types,
            reactant_coefficients=self.product_coefficients,
            product_types=self.reactant_types,
            product_coefficients=self.reactant_coefficients)


@script_interface_register
class ExclusionRadius(ScriptInterfaceHelper):
    """
    Neighbor search algorithm that detects when a particle enters the exclusion
    zone of another particle. The exclusion radii are particle type-dependent.

    During the neighbor search, the following cases can arise:

    * the central particle per-type exclusion radius is zero: return ``False``
    * the neighbor particle per-type exclusion radius is zero: return ``False``
    * the central and neighbor particles per-type exclusion radii are non-zero:
      return ``True`` if the inter-particle distance is smaller than the sum of
      their respective exclusion radii, ``False`` otherwise
    * either the central particle type or the neighbor particle type is not in
      ``exclusion_radius_per_type``: return ``True`` if the inter-particle
      distance is smaller than ``exclusion_range``, ``False`` otherwise

    Attributes
    ----------
    exclusion_radius_per_type : :obj:`dict`, optional
         Mapping of particle types to exclusion radii.
    exclusion_range : :obj:`float`
        Minimal distance from any particle whose type
        is not in ``exclusion_radius_per_type``.
    search_algorithm : :obj:`str`
        Pair search algorithm. Default is ``"order_n"``, which evaluates the
        distance between the queried particle and all other particles in the
        system, and scales with O(N). For MPI-parallel simulations, the
        ``"parallel"`` method is faster. The ``"parallel"`` method is not
        recommended for simulations on 1 MPI rank, since it comes with the
        overhead of a ghost particle update.

    Methods
    -------
    check_exclusion_range()
        Check the neighborhood of a central particle and detect if any neighbor
        is too close.

        Parameters
        -----------
        pid : :obj:`int`
            Particle id.
        ptype : :obj:`int`, optional
            Particle type. If not provided, it will be read from the particle
            and communicated to all MPI ranks.

        Returns
        -------
        :obj:`bool` :
            Whether the particle is within the exclusion radius
            of another particle.

    """
    _so_name = "ReactionMethods::ExclusionRadius"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("check_exclusion_range",)


class ReactionAlgorithm(ScriptInterfaceHelper):
    """

    This class provides the base class for Reaction Algorithms like
    the Reaction Ensemble algorithm and the constant pH method.
    Initialize the reaction algorithm by setting the
    standard pressure, temperature, and the exclusion range.
    The exclusion range mechanism is explained in more detail
    in :class:`~espressomd.reaction_methods.ExclusionRadius`.

    Note: When creating particles the velocities of the new particles are set
    according the Maxwell-Boltzmann distribution. In this step the mass of the
    new particle is assumed to equal 1.


    Parameters
    ----------
    kT : :obj:`float`
        Thermal energy of the system in simulation units
    exclusion_range : :obj:`float`
        Minimal distance from any particle, within which new particles will not
        be inserted.
    seed : :obj:`int`
        Initial counter value (or seed) of the Mersenne Twister RNG.
    exclusion_radius_per_type : :obj:`dict`, optional
         Mapping of particle types to exclusion radii.
    search_algorithm : :obj:`str`
        Pair search algorithm. Default is ``"order_n"``, which evaluates the
        distance between the inserted particle and all other particles in the
        system, which scales with O(N). For MPI-parallel simulations, the
        ``"parallel"`` method is faster. The ``"parallel"`` method is not
        recommended for simulations on 1 MPI rank, since it comes with the
        overhead of a ghost particle update.

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
        >>> import espressomd.reaction_methods
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
        >>> elc = espressomd.electrostatics.ELC(actor=p3m, maxPWerror=1.0, gap_size=elc_gap)
        >>> system.actors.add(elc)
        >>> # add constant pH method
        >>> RE = espressomd.reaction_methods.ConstantpHEnsemble(kT=1, exclusion_range=1, seed=77)
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

    get_acceptance_rate_configurational_moves()
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

    displacement_mc_move_for_particles_of_type()
        Performs displacement Monte Carlo moves for particles of a given type.
        New positions of the displaced particles are chosen from the whole box
        with a uniform probability distribution and new velocities are
        sampled from the Maxwell-Boltzmann distribution.

        The sequence of moves is only accepted if each individual move in
        the sequence was accepted. Particles are sampled without replacement.
        Therefore, calling this method once for 10 particles is not equivalent
        to calling this method 10 times for 1 particle.

        Parameters
        ----------
        type_mc : :obj:`int`
            Particle type which should be moved
        particle_number_to_be_changed : :obj:`int`
            Number of particles to move, defaults to 1.
            Particles are selected without replacement.

        Returns
        -------
        :obj:`bool`
            Whether all moves were accepted.

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
            Identifier of the reaction to modify.
            Will be multiplied by 2 internally!
        gamma : :obj:`float`
            New reaction constant for the forward reaction.

    """

    _so_name = "ReactionMethods::ReactionAlgorithm"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("remove_constraint",
                        "get_wall_constraints_in_z_direction",
                        "set_wall_constraints_in_z_direction",
                        "set_cylindrical_constraint_in_z_direction",
                        "set_volume",
                        "get_volume",
                        "get_acceptance_rate_reaction",
                        "set_non_interacting_type",
                        "get_non_interacting_type",
                        "displacement_mc_move_for_particles_of_type",
                        "change_reaction_constant",
                        "delete_particle",
                        )

    def __init__(self, **kwargs):
        if self._so_name == ReactionAlgorithm._so_name:
            raise RuntimeError(
                f"Base class '{self.__class__.__name__}' cannot be instantiated")
        if 'exclusion_radius' in kwargs:
            raise KeyError(
                'the keyword `exclusion_radius` is obsolete. Currently, the equivalent keyword is `exclusion_range`')
        super().__init__(exclusion=ExclusionRadius(**kwargs), **kwargs)
        if not 'sip' in kwargs:
            utils.check_valid_keys(self.valid_keys(), kwargs.keys())
        self._rebuild_reaction_cache()

    def valid_keys(self):
        return {"kT", "exclusion_range", "seed",
                "exclusion_radius_per_type", "search_algorithm"}

    def required_keys(self):
        return {"kT", "exclusion_range", "seed"}

    def add_reaction(self, **kwargs):
        """
        Sets up a reaction in the forward and backward direction.

        Parameters
        ----------
        gamma : :obj:`float`
            Equilibrium constant :math:`\\Gamma` of the reaction in simulation
            units (see section :ref:`Reaction Ensemble` for its definition).
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
        if not isinstance(default_charges, dict):
            raise TypeError("Argument 'default_charges' needs to be a dict")
        forward_reaction = SingleReaction(**kwargs)
        backward_reaction = forward_reaction.make_backward_reaction()
        if neutrality_check:
            self._check_charge_neutrality(
                type2charge=default_charges,
                reaction=forward_reaction)

        old_default_charges = self.default_charges
        for ptype, charge in default_charges.items():
            if ptype in old_default_charges:
                if abs(old_default_charges[ptype] - charge) > 1e-10:
                    raise ValueError(
                        f"Cannot override charge on particle type {ptype} (from {old_default_charges[ptype]} to {charge})")
        for ptype, charge in default_charges.items():
            if ptype not in old_default_charges:
                self.call_method(
                    "set_charge_of_type",
                    type=ptype,
                    charge=charge)

        self.call_method("add_reaction", reaction=forward_reaction)
        self.call_method("add_reaction", reaction=backward_reaction)
        self.check_reaction_method()
        self._rebuild_reaction_cache()

    def delete_reaction(self, **kwargs):
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
            Reaction id

        """
        self.call_method("delete_reaction", **kwargs)
        self._rebuild_reaction_cache()

    def check_reaction_method(self):
        if len(self.reactions) == 0:
            raise RuntimeError("Reaction system not initialized")

        # charges of all reactive types need to be known
        if has_features("ELECTROSTATICS"):
            for reaction in self.reactions:
                for p_type in reaction.reaction_types:
                    if p_type not in self.default_charges:
                        raise RuntimeError(
                            f"Forgot to assign charge to type {p_type}")

    def _check_charge_neutrality(self, type2charge, reaction):
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

    def _rebuild_reaction_cache(self):
        """
        Maintain a local cache of the Python representation of the
        script interface reactions. This helps reducing the overhead
        when calculating the acceptance probability.
        """
        self._reactions_cache = [x for x in self.reactions]

    def reaction(self, steps):
        """
        Performs randomly selected reactions.

        Parameters
        ----------
        steps : :obj:`int`, optional
            The number of reactions to be performed at once, defaults to 1.

        """
        self.call_method("setup_bookkeeping_of_empty_pids")
        E_pot = self.call_method("potential_energy")
        for _ in range(steps):
            reaction_id = self.call_method("get_random_reaction_index")
            E_pot = self.generic_oneway_reaction(reaction_id, E_pot)

    def calculate_acceptance_probability(self, reaction_id, E_pot_diff):
        """
        Calculate the acceptance probability of a Monte Carlo move.

        Parameters
        ----------
        reaction_id : :obj:`int`
            Identifier of the reaction that was carried out in the move.
        E_pot_diff : :obj:`float`
            The potential energy difference for the move.

        Returns
        -------
        :obj:`float`
            The acceptance probability.

        """
        factorial_expr = self.call_method("calculate_factorial_expression")
        reaction = self._reactions_cache[reaction_id]
        return self.get_volume()**reaction.nu_bar * reaction.gamma * \
            factorial_expr * math.exp(-E_pot_diff / self.kT)

    def generic_oneway_reaction(self, reaction_id, E_pot_old):
        """
        Carry out a generic one-way chemical reaction of the type
        `A+B+...+G +... --> K+...X + Z +...`

        You need to use `2A --> B` instead of `A+A --> B` since
        in the latter you assume distinctness of the particles, however both
        ways to describe the reaction are equivalent in the thermodynamic limit
        (large particle numbers). Furthermore, the order of the reactant and
        product types matters since particles will be replaced in that order!
        If there are less reactants than products, new product particles are
        created randomly in the box. Reactants get their type and charge
        changed to the corresponding type and charge of the products.
        If there are more reactants than products, excess reactant particles
        are deleted.

        Parameters
        ----------
        reaction_id : :obj:`int`
            Identifier of the reaction to attempt.
        E_pot_old : :obj:`float`
            The current potential energy.

        Returns
        -------
        E_pot_new : :obj:`float`
            The potential energy after the move.

        """
        try:
            E_pot_new = self.call_method(
                "create_new_trial_state", reaction_id=reaction_id)
            if E_pot_new is None:
                return E_pot_old
            E_pot_diff = E_pot_new - E_pot_old
            bf = self.calculate_acceptance_probability(reaction_id, E_pot_diff)
            return self.call_method("make_reaction_mc_move_attempt",
                                    reaction_id=reaction_id, bf=bf,
                                    E_pot_new=E_pot_new, E_pot_old=E_pot_old)
        except BaseException as err:
            tb = sys.exc_info()[2]
            raise RuntimeError(
                "An exception was raised by a chemical reaction; the particle "
                "state tracking is no longer guaranteed to be correct! -- "
                f"{err}").with_traceback(tb)

    def _check_reaction_index(self, reaction_index):
        if reaction_index < 0 or reaction_index >= len(self.reactions):
            raise IndexError(f"No reaction with id {reaction_index}")

    def get_status(self):
        """
        Returns the status of the reaction ensemble in a dictionary containing
        the used reactions, the used kT and the used exclusion radius.

        """

        self.check_reaction_method()
        property_keys = {"reactant_coefficients", "reactant_types",
                         "product_coefficients", "product_types", "gamma"}
        reactions_list = [{key: getattr(reaction, key) for key in property_keys}
                          for reaction in self.reactions]

        return {"reactions": reactions_list, "kT": self.kT,
                "exclusion_range": self.exclusion_range,
                "exclusion_radius_per_type": self.exclusion_radius_per_type}


@script_interface_register
class ReactionEnsemble(ReactionAlgorithm):
    """
    This class implements the Reaction Ensemble.
    """

    _so_name = "ReactionMethods::ReactionEnsemble"


@script_interface_register
class ConstantpHEnsemble(ReactionAlgorithm):
    """
    This class implements the constant pH Ensemble.

    When adding an acid-base reaction, the acid and base particle types
    are always assumed to be at index 0 of the lists passed to arguments
    ``reactant_types`` and ``product_types``.

    Attributes
    ----------
    constant_pH : :obj:`float`
        Constant pH value.

    """
    _so_name = "ReactionMethods::ConstantpHEnsemble"

    def valid_keys(self):
        return {"kT", "exclusion_range", "seed",
                "constant_pH", "exclusion_radius_per_type", "search_algorithm"}

    def required_keys(self):
        return {"kT", "exclusion_range", "seed", "constant_pH"}

    def calculate_acceptance_probability(self, reaction_id, E_pot_diff):
        """
        Calculate the acceptance probability of a Monte Carlo move.

        Parameters
        ----------
        reaction_id : :obj:`int`
            Identifier of the reaction that was carried out in the move.
        E_pot_diff : :obj:`float`
            The potential energy difference for the move.

        Returns
        -------
        :obj:`float`
            The acceptance probability.

        """
        factorial_expr = self.call_method("calculate_factorial_expression")
        reaction = self._reactions_cache[reaction_id]
        ln_bf = E_pot_diff - reaction.nu_bar * self.kT * math.log(10.) * (
            self.constant_pH + reaction.nu_bar * math.log10(reaction.gamma))
        return factorial_expr * math.exp(-ln_bf / self.kT)

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
        system following the reaction provided in ``reaction_id``. Please
        define the insertion moves by calling the method
        :meth:`~ReactionAlgorithm.add_reaction` (with only product types
        specified).

        Note that although this function does not provide directly
        the chemical potential, it can be used to calculate it.
        For an example of such an application please check
        :file:`/samples/widom_insertion.py`.

        Parameters
        ----------
        reaction_id : :obj:`int`
            Reaction identifier. Will be multiplied by 2 internally to
            skip reverse reactions, i.e. deletion reactions!

        Returns
        -------
        :obj:`float`
            The particle insertion potential energy.

        """
        return self.call_method(
            "calculate_particle_insertion_potential_energy", **kwargs)

    def calculate_excess_chemical_potential(self, **kwargs):
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
