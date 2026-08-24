#
# Copyright (C) 2010-2026 The ESPResSo project
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
import math
import sys
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .code_features import has_features
from . import utils


class SingleReaction:
    def __init__(self, **kwargs):
        utils.check_required_keys(self.required_keys(), kwargs.keys())
        utils.check_valid_keys(self.valid_keys(), kwargs.keys())
        self.reactant_types = kwargs["reactant_types"]
        self.reactant_coefficients = kwargs["reactant_coefficients"]
        self.product_types = kwargs["product_types"]
        self.product_coefficients = kwargs["product_coefficients"]
        self.gamma = kwargs["gamma"]
        if len(self.reactant_types) != len(self.reactant_coefficients):
            raise ValueError(
                "reactants: number of types and coefficients have to match")
        if len(self.product_types) != len(self.product_coefficients):
            raise ValueError(
                "products: number of types and coefficients have to match")
        if self.gamma <= 0.:
            raise ValueError("gamma needs to be a strictly positive value")
        self.accepted_moves = 0
        self.trial_moves = 0
        self.accumulator_potential_energy_difference_exponential = []
        self.nu_bar = sum(self.product_coefficients) - \
            sum(self.reactant_coefficients)

    def get_acceptance_rate(self):
        return self.accepted_moves / self.trial_moves

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

    Parameters
    ----------
    exclusion_radius_per_type : :obj:`dict`, optional
         Mapping of particle types to exclusion radii.
    exclusion_range : :obj:`float`
        Minimal distance from any particle whose type
        is not in ``exclusion_radius_per_type``.
    search_algorithm : :obj:`str`, optional
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

    def __init__(self, **kwargs):
        utils.check_required_keys(self.required_keys(), kwargs.keys())
        utils.check_valid_keys(self.valid_keys(), kwargs.keys())
        super().__init__(**kwargs)

    @staticmethod
    def required_keys():
        return {"exclusion_range"}

    @staticmethod
    def valid_keys():
        return {"exclusion_range",
                "exclusion_radius_per_type", "search_algorithm"}


class ReactionAlgorithm:
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
    system : :obj:`espressomd.system.System`
        Which system to modify.
    kT : :obj:`float`
        Thermal energy of the system in simulation units
    exclusion_range : :obj:`float`
        Minimal distance from any particle, within which new particles will not
        be inserted.
    seed : :obj:`int`
        Initial counter value (or seed) of the philox RNG.
    exclusion_radius_per_type : :obj:`dict`, optional
         Mapping of particle types to exclusion radii.
    search_algorithm : :obj:`str`
        Pair search algorithm. Default is ``"order_n"``, which evaluates the
        distance between the inserted particle and all other particles in the
        system, which scales with O(N). For MPI-parallel simulations, the
        ``"parallel"`` method is faster. The ``"parallel"`` method is not
        recommended for simulations on 1 MPI rank, since it comes with the
        overhead of a ghost particle update.

    """

    @script_interface_register
    class _ReactionAlgorithmHelper(ScriptInterfaceHelper):
        """
        Stateless class encapsulating helper functions.
        """
        _so_name = "ReactionMethods::ReactionAlgorithm"
        _so_creation_policy = "GLOBAL"

    @script_interface_register
    class _ParticleModifier(ScriptInterfaceHelper):
        """
        Stateful class encapsulating helper functions.
        """
        _so_name = "Particles::ParticleModifier"
        _so_checkpointable = False
        _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        if type(self) is ReactionAlgorithm:
            raise RuntimeError(
                f"Base class '{self.__class__.__name__}' cannot be instantiated")
        utils.check_required_keys(self.required_keys(), kwargs.keys())
        utils.check_valid_keys(self.valid_keys(), kwargs.keys())
        self.system = kwargs.pop("system")
        particle_modifier = self._ParticleModifier(
            id=-1, __cell_structure=self.system.cell_system)
        self._helper = self._ReactionAlgorithmHelper(
            system=self.system, particle_modifier=particle_modifier)
        self.kT = kwargs["kT"]
        if self.kT < 0.:
            raise ValueError("Invalid value for 'kT'")
        self.rng = np.random.Generator(np.random.Philox(kwargs["seed"]))
        self.exclusion_range_touched = False
        if "exclusion_range" not in kwargs:
            kwargs["exclusion_range"] = 0.
        self.exclusion = ExclusionRadius(
            **{k: kwargs[k] for k in ExclusionRadius.valid_keys() if k in kwargs})
        self.constraint_type = "none"
        self.params_boundaries = {}
        self.m_accepted_configurational_MC_moves = 0
        self.m_tried_configurational_MC_moves = 0
        self.non_interacting_type = 100
        self.reactions = []
        self.default_charges = {}
        self._empty_p_ids_smaller_than_max_seen_particle = []
        self._initialize_particle_changes()
        self._particle_numbers = {}
        self._analysis = self.system.analysis
        self._system_part = self.system.part

    def valid_keys(self):
        raise NotImplementedError("Derived classes must implement this method")

    def required_keys(self):
        raise NotImplementedError("Derived classes must implement this method")

    @property
    def exclusion_range(self):
        return self.exclusion.exclusion_range

    @exclusion_range.setter
    def exclusion_range(self, value):
        self.exclusion.exclusion_range = value

    @property
    def exclusion_radius_per_type(self):
        return self.exclusion.exclusion_radius_per_type

    @exclusion_radius_per_type.setter
    def exclusion_radius_per_type(self, value):
        self.exclusion.exclusion_radius_per_type = value

    @property
    def search_algorithm(self):
        return self.exclusion.search_algorithm

    @search_algorithm.setter
    def search_algorithm(self, value):
        self.exclusion.search_algorithm = value

    @classmethod
    def calculate_factorial_expression(cls, reaction, particle_numbers):
        raise NotImplementedError("Derived classes must implement this method")

    def _check_exclusion_range(self, pid, ptype=None):
        if ptype is None:
            result = self.exclusion.check_exclusion_range(pid=pid)
        else:
            result = self.exclusion.check_exclusion_range(pid=pid, ptype=ptype)
        self.exclusion_range_touched |= result

    def _check_exclusion_range_any(self, pids, ptype):
        self.exclusion_range_touched |= self.exclusion.call_method(
            b"check_exclusion_range_any", pids=pids, ptype=ptype)

    def get_random_positions_in_box(self, n):
        """
        Return random positions in the simulation box.
        If any constraint is active, the random positions will be sampled inside
        the constraint volume.
        """
        box_l = np.copy(self.system.box_l)
        if self.constraint_type == "none":
            return self.rng.random(size=(n, 3)) * box_l
        if self.constraint_type == "slab":
            box_l[2] = self.params_boundaries["slab_end_z"] - \
                self.params_boundaries["slab_start_z"]
            positions = self.rng.random(size=(n, 3)) * box_l
            positions[:, 2] += self.params_boundaries["slab_start_z"]
            return positions
        if self.constraint_type == "cylinder":
            # see https://mathworld.wolfram.com/DiskPointPicking.html
            # for uniform disk point picking in cylindrical coordinates
            radii = self.params_boundaries["radius"] * \
                np.sqrt(self.rng.uniform(size=(n,)))
            phi = 2 * np.pi * self.rng.uniform(size=(n,))
            z = box_l[2] * self.rng.uniform(size=(n,))
            return np.vstack((
                self.params_boundaries["center_x"] + radii * np.cos(phi),
                self.params_boundaries["center_y"] + radii * np.sin(phi),
                z)).T
        raise NotImplementedError(
            f"Constraint type {self.constraint_type} is not implemented")

    def get_random_pids(self, ptype, size):
        pids = self.system.call_method(b"get_pids_of_type", ptype=ptype)
        if size == 1:
            return (pids[self.rng.integers(len(pids))],)
        return self.rng.choice(pids, size=size, replace=False)

    def remove_constraint(self):
        """
        Remove any previously defined constraint.
        """
        self.constraint_type = "none"
        self.params_boundaries = {}

    def set_cylindrical_constraint_in_z_direction(
            self, center_x, center_y, radius):
        """
        Constrain the reaction moves within a cylinder aligned with the z-axis.

        Parameters
        ----------
        center_x : :obj:`float`
            x coordinate of center of the cylinder.
        center_y : :obj:`float`
            y coordinate of center of the cylinder.
        radius : :obj:`float`
            radius of the cylinder.

        """
        if center_x < 0. or center_x > self.system.box_l[0]:
            raise ValueError(f"center_x is outside the box")
        if center_y < 0. or center_y > self.system.box_l[1]:
            raise ValueError(f"center_y is outside the box")
        if radius < 0.:
            raise ValueError(f"radius is invalid")

        self.constraint_type = "cylinder"
        self.params_boundaries = {"radius": radius,
                                  "center_x": center_x,
                                  "center_y": center_y}

    def set_wall_constraints_in_z_direction(self, slab_start_z, slab_end_z):
        """
        Restrict the sampling area to a slab in z-direction.
        This constraint is necessary when
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
        >>> RE = espressomd.reaction_methods.ConstantpHEnsemble(kT=1., exclusion_range=1., seed=77, system=system)
        >>> RE.constant_pH = 2
        >>> RE.add_reaction(gamma=0.0088, reactant_types=[types["HA"]],
        ...                 product_types=[types["A-"], types["H+"]],
        ...                 default_charges=charges)
        >>> # add walls for the ELC gap
        >>> RE.set_wall_constraints_in_z_direction(0, box_l)
        >>> system.constraints.add(shape=espressomd.shapes.Wall(dist=0, normal=[0, 0, 1]),
        ...                        particle_type=types["wall"])
        >>> system.constraints.add(shape=espressomd.shapes.Wall(dist=-box_l, normal=[0, 0, -1]),
        ...                        particle_type=types["wall"])

        """
        if slab_start_z < 0. or slab_start_z > self.system.box_l[2]:
            raise ValueError("slab_start_z is outside the box")
        if slab_end_z < 0. or slab_end_z > self.system.box_l[2]:
            raise ValueError("slab_end_z is outside the box")
        if slab_end_z < slab_start_z:
            raise ValueError("slab_end_z must be >= slab_start_z")

        self.constraint_type = "slab"
        self.params_boundaries = {"slab_start_z": slab_start_z,
                                  "slab_end_z": slab_end_z}

    def get_wall_constraints_in_z_direction(self):
        """
        Get the restrictions of the sampling area in z-direction.
        """
        if self.constraint_type != "slab":
            raise RuntimeError("no slab constraint is currently active")
        return [self.params_boundaries["slab_start_z"],
                self.params_boundaries["slab_end_z"]]

    def get_volume(self):
        """
        Get the volume to be used in the acceptance probability of the reaction
        ensemble.
        """
        if self.constraint_type == "slab":
            box_l = np.copy(self.system.box_l)
            box_l[2] = self.params_boundaries["slab_end_z"] - \
                self.params_boundaries["slab_start_z"]
            return float(np.prod(box_l))
        if self.constraint_type == "cylinder":
            radius = self.params_boundaries["radius"]
            height = float(self.system.box_l[2])
            return np.pi * radius**2 * height
        return self.system.volume()

    def get_acceptance_rate_configurational_moves(self):
        """
        Return the acceptance rate for the configuration moves.
        """
        return self.m_accepted_configurational_MC_moves / \
            self.m_tried_configurational_MC_moves

    def get_acceptance_rate_reaction(self, reaction_id):
        """
        Return the acceptance rate for the given reaction.

        Parameters
        ----------
        reaction_id : :obj:`int`
            Identifier of the reaction to modify.
            Will *not* be multiplied by 2 internally!
        """
        if reaction_id < 0 or reaction_id >= len(self.reactions):
            raise IndexError(f"No reaction with id {reaction_id}")
        return self.reactions[reaction_id].get_acceptance_rate()

    def set_non_interacting_type(self, type):
        """
        Set the particle type for non-interacting particles.
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
        """
        if type < 0:
            raise ValueError(f"Invalid type: {type}")
        self.non_interacting_type = type

    def get_non_interacting_type(self):
        """
        Return the type which is used for hiding particles.
        """
        return self.non_interacting_type

    def _displacement_mc_move(self, ptype, n_particles):
        # draw particle ids at random without replacement
        p_id = -1
        drawn_pids = [p_id]
        for _ in range(n_particles):
            # draw a new particle id
            while p_id in drawn_pids:
                p_id = self.get_random_pids(ptype, 1)[0]
            drawn_pids.append(p_id)
            # write new position and new velocity
            p = self._system_part.by_id(p_id)
            self._particle_changes["changed"].append(
                {"pid": p_id, "pos": p.pos, "v": p.v})
            new_pos = self.get_random_positions_in_box(1)[0]
            new_vel = self.rng.normal(size=3) * math.sqrt(self.kT / p.mass)
            p.update({"pos": new_pos, "v": new_vel})
            self._check_exclusion_range(p_id, ptype)
            if self.exclusion_range_touched:
                break

    def displacement_mc_move_for_particles_of_type(
            self, type_mc, particle_number_to_be_changed=1):
        """
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
        """
        if type_mc < 0:
            raise ValueError("Parameter 'type_mc' must be >= 0")
        if particle_number_to_be_changed < 0:
            raise ValueError(
                "Parameter 'particle_number_to_be_changed' must be >= 0")
        if particle_number_to_be_changed == 0:
            # reject
            return False
        self.m_tried_configurational_MC_moves += 1
        self.exclusion_range_touched = False

        n_particles_of_type = self.system.number_of_particles(type=type_mc)
        if particle_number_to_be_changed > n_particles_of_type:
            # reject
            return False

        E_pot_old = self._analysis.potential_energy()
        self._displacement_mc_move(type_mc, particle_number_to_be_changed)
        E_pot_new = float("inf")
        if not self.exclusion_range_touched:
            E_pot_new = self._analysis.potential_energy()
        exp_min = -708.4  # for IEEE-compatible double
        exponent = -(E_pot_new - E_pot_old) / self.kT
        exponential = 0. if (exponent < exp_min) else math.exp(exponent)

        # Metropolis algorithm since proposal density is symmetric
        bf = min(1., exponential)

        # // correct for enhanced proposal of small radii by using the
        # // Metropolis-Hastings algorithm for asymmetric proposal densities
        # double old_radius =
        #     std::sqrt(std::pow(particle_positions[0][0]-cyl_x,2) +
        #               std::pow(particle_positions[0][1]-cyl_y,2));
        # double new_radius =
        #     std::sqrt(std::pow(new_pos[0]-cyl_x,2)+std::pow(new_pos[1]-cyl_y,2));
        # auto const bf = std::min(1.0,
        #     exp(-beta*(E_pot_new-E_pot_old))*new_radius/old_radius);

        # Metropolis-Hastings algorithm for asymmetric proposal density
        if self.rng.random(1) < bf:
            # accept
            self.m_accepted_configurational_MC_moves += 1
            self._initialize_particle_changes()
            return True

        # reject: restore original particle properties
        self._restore_system()
        return False

    def get_reaction_index(self, reaction_id):
        """
        Check reaction id is within the reaction container bounds.
        Since each reaction has a corresponding backward reaction,
        the total number of reactions is doubled. Return the
        corresponding index for chosen reaction.
        """
        index = 2 * reaction_id
        if index < 0 or index >= len(self.reactions):
            raise IndexError(f"No reaction with id {reaction_id}")
        return index

    def delete_particle(self, p_id):
        """
        Deletes the particle of the given p_id and makes sure that the particle
        range has no holes. This function has some restrictions, as e.g. bonds
        are not deleted. Therefore only apply this function to simple ions.

        Parameters
        ----------
        p_id : :obj:`int`
            Id of the particle to be deleted.
        """
        if p_id < 0:
            raise ValueError(f"Invalid particle id: {p_id}")
        self._free_particle_id(p_id, precheck=True)
        self._helper.call_method(b"delete_particle", pid=p_id)

    def change_reaction_constant(self, reaction_id, gamma):
        """
        Changes the reaction constant of a given reaction
        (for both the forward and backward reactions).
        The ``reaction_id`` which is assigned to a reaction depends on the order
        in which :meth:`~espressomd.reaction_methods.ReactionAlgorithm.add_reaction`
        was called.
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
        if gamma <= 0.:
            raise ValueError("gamma needs to be a strictly positive value")
        index = self.get_reaction_index(reaction_id)
        self.reactions[index + 0].gamma = gamma
        self.reactions[index + 1].gamma = 1. / gamma

    def __reduce__(self):
        raise NotImplementedError(
            "Reaction methods do not support checkpointing")

    def _initialize_particle_changes(self):
        self._particle_changes = {"created": [], "changed": [], "hidden": []}

    def add_reaction(self, **kwargs):
        """
        Set up a reaction in the forward and backward directions.

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
        self.default_charges.update(default_charges)
        self.reactions.append(forward_reaction)
        self.reactions.append(backward_reaction)
        self.check_reaction_method()

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
            Reaction id
        """
        index = self.get_reaction_index(reaction_id)
        del self.reactions[index + 1]
        del self.reactions[index + 0]

    @classmethod
    def _factorial_Ni0_by_factorial_Ni0_plus_nu_i(cls, nu_i, N_i0):
        value = 1.
        if nu_i > 0:
            value /= math.factorial(N_i0 + nu_i) // math.factorial(N_i0)
        elif nu_i < 0:
            value *= math.factorial(N_i0) // math.factorial(N_i0 + nu_i)
        return value

    def make_reaction_attempt(self, reaction):
        """
        Carry out a chemical reaction and save the old system state.
        """
        minimum_number_of_types = min(len(reaction.reactant_types),
                                      len(reaction.product_types))
        maximum_number_of_types = max(len(reaction.reactant_types),
                                      len(reaction.product_types))

        for index in range(minimum_number_of_types):
            r_type = reaction.reactant_types[index]
            p_type = reaction.product_types[index]
            r_charge = self.default_charges[r_type]
            p_charge = self.default_charges[p_type]

            # change reactant particles to product particles
            size = min(reaction.reactant_coefficients[index],
                       reaction.product_coefficients[index])
            if self._particle_numbers:
                self._particle_numbers[r_type] -= size
                self._particle_numbers[p_type] += size
            pids = self.get_random_pids(r_type, size)
            self._helper.call_method(
                b"batch_update", pids=pids, properties={"type": p_type, "q": p_charge})
            for random_pid in pids:
                self._particle_changes["changed"].append(
                    {"pid": random_pid, "type": r_type, "q": r_charge})

            # measure stoichiometric excess
            delta_n = reaction.product_coefficients[index] - \
                reaction.reactant_coefficients[index]

            if delta_n > 0:
                # create product particles
                pids = self._create_particles(delta_n, p_type)
                self._check_exclusion_range_any(pids, p_type)
                for pid in pids:
                    self._particle_changes["created"].append(
                        {"pid": pid, "type": p_type, "q": p_charge})
            elif delta_n < 0:
                # hide reactant particles
                pids = self.get_random_pids(r_type, -delta_n)
                self._check_exclusion_range_any(pids, r_type)
                self._hide_particles(pids, r_type)
                for random_pid in pids:
                    self._particle_changes["hidden"].append(
                        {"pid": random_pid, "type": r_type, "q": r_charge})

        # create/hide particles with non-corresponding replacement types
        for index in range(minimum_number_of_types, maximum_number_of_types):
            if len(reaction.product_types) < len(reaction.reactant_types):
                r_type = reaction.reactant_types[index]
                r_charge = self.default_charges[r_type]
                size = reaction.reactant_coefficients[index]
                # hide superfluous reactant particles
                pids = self.get_random_pids(r_type, size)
                self._check_exclusion_range_any(pids, r_type)
                self._hide_particles(pids, r_type)
                for random_pid in pids:
                    self._particle_changes["hidden"].append(
                        {"pid": random_pid, "type": r_type, "q": r_charge})
            else:
                p_type = reaction.product_types[index]
                p_charge = self.default_charges[p_type]
                # create additional product particles
                delta_n = reaction.product_coefficients[index]
                pids = self._create_particles(delta_n, p_type)
                self._check_exclusion_range_any(pids, p_type)
                for pid in pids:
                    self._particle_changes["created"].append(
                        {"pid": pid, "type": p_type, "q": p_charge})

    def all_reactant_particles_exist(self, reaction):
        for r_type in reaction.reactant_types:
            r_index = reaction.reactant_types.index(r_type)
            r_coef = reaction.reactant_coefficients[r_index]
            if self._particle_numbers[r_type] < r_coef:
                return False
        return True

    def count_number_of_particles_per_type(self):
        types = [self.non_interacting_type]
        for reaction in self.reactions:
            types += reaction.reactant_types + reaction.product_types
        types = list(set(types))
        numbers = self._helper.call_method(
            b"count_number_of_particles_per_type", types=types)
        return dict(zip(types, numbers))

    def _free_particle_id(self, p_id, precheck=False):
        old_max_seen_id = self.system.call_method(
            b"reaction_get_maximal_particle_id")
        if p_id == old_max_seen_id:
            self._empty_p_ids_smaller_than_max_seen_particle = [
                x for x in self._empty_p_ids_smaller_than_max_seen_particle if x < old_max_seen_id]
        elif p_id <= old_max_seen_id:
            self._empty_p_ids_smaller_than_max_seen_particle.append(p_id)
        elif precheck:
            raise RuntimeError(
                "Particle id is greater than the max seen particle id")

    def _delete_created_particles(self):
        pids = []
        for particle_info in self._particle_changes["created"]:
            pids.append(particle_info["pid"])
            if self._particle_numbers:
                self._particle_numbers[particle_info["type"]] -= 1
            self._free_particle_id(particle_info["pid"])
        self._helper.call_method(b"delete_particles", pids=pids)

    def _delete_hidden_particles(self):
        pids = []
        for particle_info in self._particle_changes["hidden"]:
            pids.append(particle_info["pid"])
            if self._particle_numbers:
                self._particle_numbers[self.non_interacting_type] -= 1
            self._free_particle_id(particle_info["pid"])
        self._helper.call_method(b"delete_particles", pids=pids)

    def _restore_system(self):
        # restore properties of changed and hidden particles
        for particle_info in self._particle_changes["changed"] + \
                self._particle_changes["hidden"]:
            pid = particle_info.pop("pid")
            ptype = self._helper.call_method(
                b"single_update", pid=pid, properties=particle_info)
            if self._particle_numbers:
                self._particle_numbers[ptype] -= 1
                self._particle_numbers[particle_info["type"]] += 1
        # destroy created particles
        self._delete_created_particles()
        self._initialize_particle_changes()

    def _hide_particle(self, pid):
        ptype = self._helper.call_method(
            b"single_update", pid=pid,
            properties={"type": self.non_interacting_type, "q": 0.})
        if self._particle_numbers:
            self._particle_numbers[ptype] -= 1
            self._particle_numbers[self.non_interacting_type] += 1

    def _hide_particles(self, pids, ptype):
        self._helper.call_method(
            b"batch_update", pids=pids,
            properties={"type": self.non_interacting_type, "q": 0.})
        if self._particle_numbers:
            self._particle_numbers[ptype] -= len(pids)
            self._particle_numbers[self.non_interacting_type] += len(pids)

    def _create_particles(self, size, ptype):
        pids = []
        highest_particle_id = self._system_part.highest_particle_id
        for _ in range(size):
            if len(self._empty_p_ids_smaller_than_max_seen_particle) == 0:
                pid = highest_particle_id + 1
                highest_particle_id = pid
            else:
                pid = min(self._empty_p_ids_smaller_than_max_seen_particle)
                self._empty_p_ids_smaller_than_max_seen_particle.remove(pid)
                highest_particle_id = max(highest_particle_id, pid)
            pids.append(pid)
        new_pos = self.get_random_positions_in_box(size)
        new_v = self.rng.normal(size=(size, 3)) * math.sqrt(self.kT)
        for i in range(size):
            self._system_part.add(
                id=pids[i], type=ptype, q=self.default_charges[ptype],
                pos=new_pos[i], v=new_v[i])
        if self._particle_numbers:
            self._particle_numbers[ptype] += size
        return pids

    def _setup_bookkeeping_of_empty_pids(self):
        particle_ids = self._system_part.all().id
        available_pids = self._find_missing_pids(pids_list=particle_ids)
        self._empty_p_ids_smaller_than_max_seen_particle = available_pids

    def _find_missing_pids(self, pids_list):
        """
        Finds the missing particles ids in `pids_list`.
        NOTE: ``pids_list`` must be a sorted list [0,1,3,5,7..]
        """
        return [i for x, y in zip(pids_list, pids_list[1:])
                for i in range(x + 1, y) if y - x > 1]

    def check_reaction_method(self):
        if len(self.reactions) == 0:
            raise RuntimeError("Reaction system not initialized")

        # charges of all reactive types need to be known
        if has_features("ELECTROSTATICS"):
            for reaction in self.reactions:
                for p_type in reaction.reactant_types:
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

    def _setup_cache(self):
        self._setup_bookkeeping_of_empty_pids()
        self._particle_numbers = self.count_number_of_particles_per_type()

    def reaction(self, steps=1):
        """
        Perform reaction steps. Chemical reactions are selected at random.

        Parameters
        ----------
        steps : :obj:`int`, optional
            The number of reactions to be performed at once, defaults to 1.

        """
        self._setup_cache()
        E_pot = self._analysis.potential_energy()
        n_reactions = len(self.reactions)
        for i in self.rng.choice(n_reactions, size=steps, replace=True):
            E_pot = self.generic_oneway_reaction(self.reactions[i], E_pot)

    def calculate_log_acceptance_probability(
            self, reaction, E_pot_diff, old_particle_numbers):
        """
        Calculate the logarithmic acceptance probability of a Monte Carlo move.

        Parameters
        ----------
        reaction : :class:`SingleReaction`
            The reaction that was carried out in the move.
        E_pot_diff : :obj:`float`
            The potential energy difference for the move.
        old_particle_numbers : :obj:`dict`
            The particle numbers before the move.

        Returns
        -------
        :obj:`float`
            The acceptance probability.

        """
        raise NotImplementedError("Derived classes must implement this method")

    def generic_oneway_reaction(self, reaction, E_pot_old):
        """
        Carry out a generic one-way chemical reaction of the type
        ``A + B + ... --> C + D + ...`` and return the new potential energy.

        You need to use ``2A --> B`` instead of ``A+A --> B`` since
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
        reaction : :obj:`SingleReaction`
            The reaction to attempt.
        E_pot_old : :obj:`float`
            The current potential energy.

        Returns
        -------
        E_pot_new : :obj:`float`
            The potential energy in the new configuration if the trial move
            was accepted, otherwise the original potential energy.
        """
        try:
            reaction.trial_moves += 1
            self.exclusion_range_touched = False
            if not self.all_reactant_particles_exist(reaction):
                return E_pot_old

            types = reaction.reactant_types + reaction.product_types
            old_particle_numbers = {
                k: v for k, v in self._particle_numbers.items() if k in types}
            self.make_reaction_attempt(reaction)

            if self.exclusion_range_touched:
                # reject trial move
                self._restore_system()
                self.exclusion_range_touched = False
                return E_pot_old

            E_pot_new = self._analysis.potential_energy()
            E_pot_diff = E_pot_new - E_pot_old
            ln_bf = self.calculate_log_acceptance_probability(
                reaction, E_pot_diff, old_particle_numbers)
            reaction.accumulator_potential_energy_difference_exponential.append(
                math.exp(-E_pot_diff / self.kT))

            if -self.rng.standard_exponential() >= ln_bf:
                # reject trial move
                self._restore_system()
                return E_pot_old

            # accept trial move
            self._delete_hidden_particles()
            reaction.accepted_moves += 1
            self._initialize_particle_changes()
            return E_pot_new
        except BaseException as err:
            tb = sys.exc_info()[2]
            raise RuntimeError(
                "An exception was raised by a chemical reaction; the particle "
                "state tracking is no longer guaranteed to be correct! -- "
                f"{err}").with_traceback(tb)

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

    @staticmethod
    def _ln_factorial_Ni0_div_factorial_Ni0_nu_i(N_i0, nu_i):
        """
        Calculate :math:`\\frac{N_i^0!}{(N_i^0+\\nu_{i}\\xi)!}`
        """
        if nu_i == 0:
            return 0.
        if N_i0 + nu_i < 0:
            return -float("inf")
        return math.lgamma(N_i0 + 1) - math.lgamma(N_i0 + nu_i + 1)


class ReactionEnsemble(ReactionAlgorithm):
    """
    This class implements the Reaction Ensemble.
    """

    def valid_keys(self):
        return {"system", "kT", "exclusion_range", "seed",
                "exclusion_radius_per_type", "search_algorithm"}

    def required_keys(self):
        return {"system", "kT", "exclusion_range", "seed"}

    def _setup_cache(self):
        super()._setup_cache()
        self.volume = self.get_volume()

    def calculate_log_acceptance_probability(
            self, reaction, E_pot_diff, old_particle_numbers):
        __doc__ = ReactionAlgorithm.__doc__  # pylint: disable=unused-variable
        ln_factorial = self.calculate_factorial_expression(
            reaction, old_particle_numbers)
        ln_bf = -E_pot_diff / self.kT + reaction.nu_bar * \
            math.log(self.volume) + math.log(reaction.gamma)
        return ln_factorial + ln_bf

    @classmethod
    def calculate_factorial_expression(cls, reaction, particle_numbers):
        """
        Calculate the logarithm of the product of factorial expressions which
        occur in the reaction ensemble acceptance probability :cite:`smith94c`.

        :math:`\\log\\left(\\prod_{i=1}\\frac{N_i^0!}{(N_i^0+\\nu_{i}\\xi)!}\\right)`
        """
        value = 0.
        # factorial contribution of reactants
        for i in range(len(reaction.reactant_types)):
            nu_i = -reaction.reactant_coefficients[i]
            N_i0 = particle_numbers[reaction.reactant_types[i]]
            value += cls._ln_factorial_Ni0_div_factorial_Ni0_nu_i(N_i0, nu_i)
        # factorial contribution of products
        for i in range(len(reaction.product_types)):
            nu_i = reaction.product_coefficients[i]
            N_i0 = particle_numbers[reaction.product_types[i]]
            value += cls._ln_factorial_Ni0_div_factorial_Ni0_nu_i(N_i0, nu_i)
        return value


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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.constant_pH = kwargs["constant_pH"]

    def valid_keys(self):
        return {"system", "kT", "exclusion_range", "seed",
                "constant_pH", "exclusion_radius_per_type", "search_algorithm"}

    def required_keys(self):
        return {"system", "kT", "exclusion_range", "seed", "constant_pH"}

    def calculate_log_acceptance_probability(
            self, reaction, E_pot_diff, old_particle_numbers):
        __doc__ = ReactionAlgorithm.__doc__  # pylint: disable=unused-variable
        ln_factorial_expr = self.calculate_factorial_expression(
            reaction, old_particle_numbers)
        ln_bf = E_pot_diff - reaction.nu_bar * self.kT * math.log(10.) * (
            self.constant_pH + reaction.nu_bar * math.log10(reaction.gamma))

        return ln_factorial_expr - ln_bf / self.kT

    def add_reaction(self, **kwargs):
        kwargs["reactant_coefficients"] = [1]
        kwargs["product_coefficients"] = [1, 1]
        super().add_reaction(**kwargs)

    @classmethod
    def calculate_factorial_expression(cls, reaction, particle_numbers):
        """
        Calculate the logarithm of the product of factorial expressions which
        occur in the constant pH method with symmetric proposal probability
        :cite:`landsgesell17b`. Only the acid and conjugated base are involved
        in the product.

        :math:`\\log\\left(\\prod_{i=1}\\frac{N_i^0!}{(N_i^0+\\nu_{i}\\xi)!}\\right)`
        """
        value = 0.
        # factorial contribution of reactants
        nu_i = -reaction.reactant_coefficients[0]
        N_i0 = particle_numbers[reaction.reactant_types[0]]
        value += cls._ln_factorial_Ni0_div_factorial_Ni0_nu_i(N_i0, nu_i)
        # factorial contribution of products
        nu_i = reaction.product_coefficients[0]
        N_i0 = particle_numbers[reaction.product_types[0]]
        value += cls._ln_factorial_Ni0_div_factorial_Ni0_nu_i(N_i0, nu_i)
        return value


class WidomInsertion(ReactionAlgorithm):
    """
    This class implements the Widom insertion method in the canonical ensemble
    for homogeneous systems, where the excess chemical potential is not
    dependent on the particle position.

    """

    @property
    def exclusion_range(self):
        return self.exclusion.exclusion_range

    @exclusion_range.setter
    def exclusion_range(self, value):
        raise RuntimeError("No search algorithm for WidomInsertion")

    @property
    def exclusion_radius_per_type(self):
        return self.exclusion.exclusion_radius_per_type

    @exclusion_radius_per_type.setter
    def exclusion_radius_per_type(self, value):
        raise RuntimeError("No search algorithm for WidomInsertion")

    @property
    def search_algorithm(self):
        return None

    @search_algorithm.setter
    def search_algorithm(self, value):
        raise RuntimeError("No search algorithm for WidomInsertion")

    def required_keys(self):
        return {"system", "kT", "seed"}

    def valid_keys(self):
        return {"system", "kT", "seed"}

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
        self._setup_cache()
        index = self.get_reaction_index(kwargs.pop("reaction_id"))
        reaction = self.reactions[index]
        if not self.all_reactant_particles_exist(reaction):
            raise RuntimeError("Trying to remove some non-existing particles "
                               "from the system via the inverse Widom scheme.")
        self._setup_bookkeeping_of_empty_pids()
        E_pot_old = self._analysis.potential_energy()
        self.make_reaction_attempt(reaction)
        E_pot_new = self._analysis.potential_energy()
        self._restore_system()
        return E_pot_new - E_pot_old

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
