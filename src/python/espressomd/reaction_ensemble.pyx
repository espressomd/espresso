include "myconfig.pxi"
from libcpp.vector cimport vector
from libcpp.string cimport string
import numpy as np


class WangLandauHasConverged(Exception):
    pass


cdef class ReactionAlgorithm(object):
    """

    This class provides the base class for Reaction Algorithms like the Reaction Ensemble algorithm, the Wang-Landau
    Reaction Ensemble algorithm and the constant pH method. Initialize the
    reaction algorithm by setting the standard pressure, temperature, and the
    exclusion radius.

    Note: When creating particles the velocities of the new particles are set
    according the Maxwell-Boltzmann distribution. In this step the mass of the
    new particle is assumed to equal 1.


    Parameters
    ----------
    temperature : :obj:`float`
                  The temperature at which the reaction is performed.
    exclusion_radius : :obj:`float`
                       Minimal distance from any particle, within which new
                       particle will not be inserted. This is useful to avoid
                       integrator failures if particles are too close and there
                       is a diverging repulsive interaction, or to prevent two
                       oppositely charged particles from being placed on top of
                       each other.  The Boltzmann factor :math:`\exp(-\\beta
                       E)` gives these configurations a small contribution to
                       the partition function, therefore they can be neglected.
    """
    cdef object _params
    cdef CReactionAlgorithm * RE

    def _valid_keys(self):
        return "temperature", "exclusion_radius"

    def _required_keys(self):
        return "temperature", "exclusion_radius"

    def _set_params_in_es_core(self):
        self.RE.temperature = self._params[
            "temperature"]
        # setting a volume is a side effect, sets the default volume of the
        # reaction ensemble as the volume of the cuboid simulation box. this
        # volume can be altered by the command "reaction ensemble volume
        # <volume>" if one wants to simulate e.g. in a system with constraint
        # (e.g. cuboid box with cylinder constraint, so that the particles are
        # only contained in the volume of the cylinder)
        if(self.RE.volume < 0):
            self.RE.set_cuboid_reaction_ensemble_volume()
        self.RE.exclusion_radius = self._params[
            "exclusion_radius"]

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
        self.RE.cyl_x = center_x
        self.RE.cyl_y = center_y
        self.RE.cyl_radius = radius_of_cylinder
        self.RE.box_is_cylindric_around_z_axis = True

    def set_wall_constraints_in_z_direction(self, slab_start_z, slab_end_z):
        """
        Restrict the sampling area to a slab in z-direction. Requires setting
        the volume using :meth:`set_volume`.

        """
        self.RE.slab_start_z = slab_start_z
        self.RE.slab_end_z = slab_end_z
        self.RE.box_has_wall_constraints = True

    def get_wall_constraints_in_z_direction(self):
        """
        Returns the restrictions of the sampling area in z-direction.

        """
        return self.RE.slab_start_z, self.RE.slab_end_z

    def set_volume(self, volume):
        """
        Set the volume to be used in the acceptance probability of the reaction
        ensemble. This can be useful when using constraints, if the relevant
        volume is different from the box volume. If not used the default volume
        which is used, is the box volume.

        """
        self.RE.volume = volume

    def get_volume(self):
        """
        Get the volume to be used in the acceptance probability of the reaction
        ensemble.

        """
        return self.RE.volume

    def get_acceptance_rate_configurational_moves(self):
        """
        Returns the acceptance rate for the configuration moves.

        """
        return self.RE.get_acceptance_rate_configurational_moves()

    def get_acceptance_rate_reaction(self, reaction_id):
        """
        Returns the acceptance rate for the given reaction.

        """
        return self.RE.reactions[reaction_id].get_acceptance_rate()

    def set_non_interacting_type(self, non_interacting_type):
        """
        Sets the particle type for non-interacting particles.
        Default value: 100.
        This is used to temporarily hide
        particles during a reaction trial move, if they are to be deleted after
        the move is accepted.
        Please change this value if you intend to use the type 100 for some other
        particle types with interactions. Please also note that particles in the current
        implementation of the Reaction Ensemble are only hidden with respect to
        Lennard-Jones and Coulomb interactions. Hiding of other interactions,
        for example a magnetic, needs to be implemented in the code.
        """
        self.RE.non_interacting_type = non_interacting_type

    def get_non_interacting_type(self):
        """
        Returns the type which is used for hiding particles.

        """
        return self.RE.non_interacting_type

    def add_reaction(self, *args, **kwargs):
        """
        Sets up a reaction in the forward and backward direction.

        Parameters
        ----------
        gamma : :obj:`float`
                               Equilibrium constant of the reaction, :math:`\gamma`
                               (see the User guide, section 6.6 for the definition and further details).
        reactant_types : list of :obj:`int`
                                List of particle types of reactants in the reaction.
        reactant_coefficients : list
                                List of stoichiometric coefficients of the
                                reactants in the same order as the list of
                                their types.
        product_types : list
                               List of particle types of products in the reaction.
        product_coefficients : list
                               List of stoichiometric coefficients of
                               products of the reaction in the same order as
                               the list of their types
        default_charges : dictionary
                        A dictionary of default charges for types that occur in the provided reaction.

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
            except:
                pass
        self._check_lengths_of_arrays()
        self._validate_params_default_charge()
        self._set_params_in_es_core_add()

    def _valid_keys_add(self):
        return "gamma", "reactant_types", "reactant_coefficients", "product_types", "product_coefficients", "default_charges", "check_for_electroneutrality"

    def _required_keys_add(self):
        return ["gamma", "reactant_types", "reactant_coefficients", "product_types", "product_coefficients", "default_charges"]

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
        self.RE.add_reaction(
            self._params["gamma"], reactant_types, reactant_coefficients, product_types, product_coefficients)
        self.RE.add_reaction(
            1.0 / self._params["gamma"], product_types, product_coefficients, reactant_types, reactant_coefficients)

        for key in self._params["default_charges"]:  # the keys are the types
            self.RE.charges_of_types[
                int(key)] = self._params["default_charges"][key]
        self.RE.check_reaction_ensemble()

    def _validate_params_default_charge(self):
        if(isinstance(self._params["default_charges"], dict) == False):
            raise ValueError(
                "No dictionary for relation between types and default charges provided.")
        #check electroneutrality of the provided reaction
        if(self._params["check_for_electroneutrality"]):
            charges = np.array(list(self._params["default_charges"].values()))
            if(np.count_nonzero(charges) == 0):
                # all partices have zero charge
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
        self.RE.do_reaction(int(reaction_steps))

    def displacement_mc_move_for_particles_of_type(self, type_mc,
                                                   particle_number_to_be_changed=1):
        """
        Performs a displacement Monte Carlo move for particles of given type. New positions
        of the displaced particles are chosen from the whole box with a uniform probability distribution.
        If there are multiple types, that are being moved in a simulation, they should be moved in a
        random order to avoid artefacts.

        Parameters
        ----------
        type_mc : :obj:`int`
            particle type which should be moved

        """

        use_wang_landau = False
        self.RE.do_global_mc_move_for_particles_of_type(
            type_mc, particle_number_to_be_changed, use_wang_landau)

    def get_status(self):
        """
        Returns the status of the reaction ensemble in a dictionary containing
        the used reactions, the used temperature and the used exclusion radius.

        """
        self.RE.check_reaction_ensemble()
        reactions = []
        for single_reaction_i in range(self.RE.reactions.size()):
            reactant_types = []
            for i in range(self.RE.reactions[single_reaction_i].reactant_types.size()):
                reactant_types.append(
                    self.RE.reactions[single_reaction_i].reactant_types[i])
            reactant_coefficients = []
            for i in range(self.RE.reactions[single_reaction_i].reactant_types.size()):
                reactant_coefficients.append(
                    self.RE.reactions[single_reaction_i].reactant_coefficients[i])

            product_types = []
            for i in range(self.RE.reactions[single_reaction_i].product_types.size()):
                product_types.append(
                    self.RE.reactions[single_reaction_i].product_types[i])
            product_coefficients = []
            for i in range(self.RE.reactions[single_reaction_i].product_types.size()):
                product_coefficients.append(
                    self.RE.reactions[single_reaction_i].product_coefficients[i])
            reaction = {"reactant_coefficients": reactant_coefficients, "reactant_types": reactant_types, "product_types": product_types, "product_coefficients":
                        product_coefficients, "reactant_types": reactant_types, "gamma": self.RE.reactions[single_reaction_i].gamma}
            reactions.append(reaction)

        return {"reactions": reactions, "temperature": self.RE.temperature, "exclusion_radius": self.RE.exclusion_radius}

    def delete_particle(self, p_id):
        """
        Deletes the particle of the given p_id and makes sure that the particle
        range has no holes. This function has some restrictions, as e.g. bonds
        are not deleted. Therefore only apply this function to simple ions.

        """
        self.RE.delete_particle(p_id)


cdef class ReactionEnsemble(ReactionAlgorithm):
    """
    This class implements the Reaction Ensemble.
    """

    cdef CReactionEnsemble * REptr

    def __init__(self, *args, **kwargs):
        self.RE = <CReactionAlgorithm * > new CReactionEnsemble()
        self.REptr = <CReactionEnsemble * > self.RE
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
                raise KeyError("%s is not a vaild key" % k)

        self._set_params_in_es_core()

    def __dealloc__(self):
        del(self.REptr)

cdef class ConstantpHEnsemble(ReactionAlgorithm):
    cdef CConstantpHEnsemble * constpHptr

    def __init__(self, *args, **kwargs):
        self.RE = <CReactionAlgorithm * > new CConstantpHEnsemble()
        self.constpHptr = <CConstantpHEnsemble * > self.RE
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
                raise KeyError("%s is not a vaild key" % k)

        self._set_params_in_es_core()

    def __dealloc__(self):
        del(self.constpHptr)

    def add_reaction(self, *args, **kwargs):
        if(len(kwargs["product_types"]) != 2 or len(kwargs["reactant_types"]) != 1):
            raise ValueError(
                "The constant pH method is only implemented for reactionw with two product types and one educt type.")
        if(kwargs["reactant_coefficients"][0] != 1 or kwargs["product_coefficients"][0] != 1 or kwargs["product_coefficients"][1] != 1):
            raise ValueError(
                "All product and reactant coefficients must equal one in the constant pH method as implemented in Espresso.")
        super(ConstantpHEnsemble, self).add_reaction(*args, **kwargs)

    property constant_pH:
        """
        Sets the input pH for the constant pH ensemble method.

        """

        def __set__(self, double pH):
            """
            Sets the pH that the method assumes for the implicit pH bath.

            """

            self.constpHptr.m_constant_pH = pH

cdef class WangLandauReactionEnsemble(ReactionAlgorithm):
    """
    This Class implements the Wang-Landau Reaction Ensemble.
    """

    cdef CWangLandauReactionEnsemble * WLRptr

    def __init__(self, *args, **kwargs):
        self.RE = <CReactionAlgorithm * > new CWangLandauReactionEnsemble()
        self.WLRptr = <CWangLandauReactionEnsemble * > self.RE
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
                raise KeyError("%s is not a vaild key" % k)

        self._set_params_in_es_core()

    def __dealloc__(self):
        del(self.WLRptr)

    def reaction(self, reaction_steps=1):
        """
        Performs reaction_steps reactions. Sets the number of reaction steps which are
        performed at once. Do not use too many reaction steps
        steps consecutively without having conformation
        changing steps in between (especially important for the Wang-Landau reaction ensemble). Providing a number for the parameter reaction steps reduces the need for the interpreter to be
        called between consecutive reactions.

        """
        status_wang_landau = self.WLRptr.do_reaction(int(reaction_steps))
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
        corresponding_acid_types : list
                                   List of the types of the version of the
                                   species.

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
        self.WLRptr.add_new_CV_degree_of_association(
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
        self.WLRptr.add_new_CV_potential_energy(
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
                                      Sets the final Wang-Landau parameter, which is the Wang-Landau parameter after which the simulation should stop.).
        full_path_to_output_filename : :obj:`str`
                                       Sets the path to the output file of the
                                       Wang-Landau algorithm which contains the
                                       Wang-Landau potential
        do_not_sample_reaction_partition_function : :obj:`bool`
                                                    Avoids sampling the
                                                    Reaction ensemble partition
                                                    function in the Wang-Landau
                                                    algorithm. Therefore this
                                                    option makes all degrees of
                                                    association equally
                                                    probable. This option may
                                                    be used in the sweeping
                                                    mode of the reaction
                                                    ensemble, since the
                                                    reaction ensemble partition
                                                    function can be later added
                                                    analytically.

        """
        for k in kwargs:
            if k in self._valid_keys_set_wang_landau_parameters():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a valid key" % k)

        self.WLRptr.final_wang_landau_parameter = self._params[
            "final_wang_landau_parameter"]
        self.WLRptr.output_filename = self._params[
            "full_path_to_output_filename"].encode("utf-8")
        self.WLRptr.do_not_sample_reaction_partition_function = self._params[
            "do_not_sample_reaction_partition_function"]

    def _valid_keys_set_wang_landau_parameters(self):
        return "final_wang_landau_parameter", "full_path_to_output_filename", "do_not_sample_reaction_partition_function"

    def load_wang_landau_checkpoint(self):
        """
        Loads the dumped Wang-Landau potential file.

        """
        checkpoint_name = "checkpoint".encode("utf-8")
        self.WLRptr.load_wang_landau_checkpoint(checkpoint_name)

    def write_wang_landau_checkpoint(self):
        """
        Dumps the Wang-Landau potential to a checkpoint file. Can be used to
        checkpoint the Wang-Landau histogram, potential, parameter and the
        number of executed trial moves.

        """
        checkpoint_name = "checkpoint".encode("utf-8")
        self.WLRptr.write_wang_landau_checkpoint(checkpoint_name)

    def update_maximum_and_minimum_energies_at_current_state(self):
        """
        Records the minimum and maximum potential energy as a function of the
        degree of association in a preliminary Wang-Landau reaction ensemble
        simulation where the acceptance probability includes the factor
        :math:`\exp(-\\beta \\Delta E_{pot})`. The minimal and maximal
        potential energies which occur in the system are needed for the energy
        reweighting simulations where the factor :math:`\exp(-\\beta \\Delta E_{pot})`
        is not included in the acceptance probability in
        order to avoid choosing the wrong potential energy boundaries.

        """
        self.WLRptr.update_maximum_and_minimum_energies_at_current_state()

    def write_out_preliminary_energy_run_results(self):
        """
        This writes out the minimum and maximum potential energy as a function
        of the degree of association to a file. It requires that previously
        :meth:`update_maximum_and_minimum_energies_at_current_state` was used.

        """
        filename = "preliminary_energy_run_results".encode("utf-8")
        self.WLRptr.write_out_preliminary_energy_run_results(filename)

    def write_wang_landau_results_to_file(self, filename):
        """
        This writes out the Wang-Landau potential as a function of the used
        collective variables.

        """
        self.WLRptr.write_wang_landau_results_to_file(filename.encode("utf-8"))

    def displacement_mc_move_for_particles_of_type(self, type_mc,
                                                   particle_number_to_be_changed=1):
        """
        Performs an MC (Monte Carlo) move for particle_number_to_be_changed particle of type type_mc. Positions for the particles are drawn uniformly random within the box. The command takes into account the Wang-Landau terms in the acceptance probability.
        If there are multiple types, that
        need to be moved, make sure to move them in a random order to avoid
        artefacts. For the Wang-Landau algorithm in the case of energy reweighting you would also need to move the monomers of the polymer with special moves for the MC part. Those polymer configuration changing moves need to be implemented in the case of using Wang-Landau with energy reweighting and a polymer in the system. Polymer configuration changing moves had been implemented before but were removed from espresso.

        """
        use_wang_landau = True
        self.WLRptr.do_global_mc_move_for_particles_of_type(
            type_mc, particle_number_to_be_changed, use_wang_landau)


cdef class WidomInsertion(ReactionAlgorithm):
    """
    This class implements the Widom Insertion Method for homogeneous systems, where the excess chemical potential is not depending on the location.

    """

    cdef CWidomInsertion * WidomInsertionPtr

    def __init__(self, *args, **kwargs):
        self.WidomInsertionPtr = new CWidomInsertion()
        self.RE = <CReactionAlgorithm * > self.WidomInsertionPtr
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
                raise KeyError("%s is not a vaild key" % k)

        self._set_params_in_es_core()

    def __dealloc__(self):
        del(self.WidomInsertionPtr)

    def measure_excess_chemical_potential(self, reaction_id=0):
        """
        Measures the excess chemical potential in a homogeneous system. Returns the excess chemical potential and the standard error for the excess chemical potential. It assumes that your samples are uncorrelated in estimating the standard error.

        """
        return self.WidomInsertionPtr.measure_excess_chemical_potential(int(reaction_id))
