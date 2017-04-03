include "myconfig.pxi"
from libc.stdlib cimport malloc, realloc, calloc

IF REACTION_ENSEMBLE:
    es_error=1
    class Wang_Landau_has_converged(Exception):
        pass

    cdef _set_params_in_es_core_add_reverse(equilibrium_constant, int* reactant_types, len_reactant_types, int* reactant_coefficients, int* product_types, len_product_types, int* product_coefficients, nu_bar):
        #add reverse reaction based on last reaction that was added            
        cdef single_reaction* new_reaction =<single_reaction *>malloc(sizeof(single_reaction))     
        new_reaction.equilibrium_constant=equilibrium_constant
        new_reaction.reactant_types=reactant_types
        new_reaction.len_reactant_types=len_reactant_types
        new_reaction.reactant_coefficients=reactant_coefficients
        new_reaction.product_types=product_types
        new_reaction.len_product_types=len_product_types
        new_reaction.product_coefficients=product_coefficients
        new_reaction.nu_bar=nu_bar
        #if everything is fine:
        current_reaction_system.reactions=<single_reaction**> realloc(current_reaction_system.reactions,sizeof(single_reaction*)*(current_reaction_system.nr_single_reactions+1)); #enlarge current_reaction_system
        current_reaction_system.reactions[current_reaction_system.nr_single_reactions]=new_reaction;
        current_reaction_system.nr_single_reactions+=1;

    cdef class ReactionEnsemble:

        def __init__(self,*args,**kwargs):
            """
            initialize the reaction ensemble by setting the standard pressure, temperature, and the exclusion radius. 

            Parameters
            ----------
            standard_pressure : float
                                the pressure in simulation units where the reactions should occur. This is an input parameter of the reaction ensemble
            temperature : float
                          the temperature at which the reaction is performed
            exclusion_radius : float
                               Exclusion radius is the minimal distance that a new particle must have towards another particle when the new particle is inserted. This is valid if there is some repulsive potential in the system, that brings the energy to (approximately) infinity if particles are too close and therefore :math:`\exp(-\\beta E)` gives these configurations aproximately zero contribution in the partition function. The exclusion radius needs to be set in order to avoid oppositley charged particles to be set too close to each other or in order to avoid too steep gradients from the short ranged interaction potential when using the Reaction ensemble together with a MD scheme.
            
            """
            self._params=self._default_params()

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

        def _default_params(self):
            return {"standard_pressure":0,
                    "temperature":1,
                    "exclusion_radius":0}

        def _valid_keys(self):
            return "standard_pressure", "temperature", "exclusion_radius"

        def _required_keys(self):
            return "standard_pressure", "temperature", "exclusion_radius"

        def _get_params_from_es_core(self):
            params = {"temperature":current_reaction_system.temperature_reaction_ensemble,
                      "standard_pressure":current_reaction_system.standard_pressure_in_simulation_units,
                      "exclusion_radius":current_reaction_system.exclusion_radius}
            return params

        def _set_params_in_es_core(self):
            current_reaction_system.temperature_reaction_ensemble = self._params["temperature"]
            #setting a volume is a side effect, sets the default volume of the reaction ensemble as the volume of the cuboid simulation box. this volume can be altered by the command "reaction ensemble volume <volume>" if one wants to simulate e.g. in a system with constraint (e.g. cuboid box with cylinder constraint, so that the particles are only contained in the volume of the cylinder)
            if(current_reaction_system.volume<0):
                set_cuboid_reaction_ensemble_volume()
            current_reaction_system.standard_pressure_in_simulation_units = self._params["standard_pressure"]
            current_reaction_system.exclusion_radius = self._params["exclusion_radius"]

        def set_cylindrical_constraint_in_z_direction(center_x, center_y, radius_of_cylinder):
            """
            constrain the reaction moves within a cylinder defined by its axis
            passing through centres (:math:`x` and :math:`y`) and the radius.
            Requires setting the volume using :meth:`set_volume`.

            Parameters
            ----------
            center_x : float
                       x coordinate of center of the cylinder.
            center_y : float
                       y coordinate of center of the cylinder.
            radius_of_cylinder : float
                       radius of the cylinder

            """
            current_reaction_system.cyl_x=center_x
            current_reaction_system.cyl_y=center_y
            current_reaction_system.cyl_radius=radius_of_cylinder
            current_reaction_system.box_is_cylindric_around_z_axis=True
            
        def set_wall_constraints_in_z_direction(slab_start_z,slab_end_z):
            """
            restrict the sampling area to a slab in z-direction. Requires setting the volume using :meth:`set_volume`.
            
            """
            current_reaction_system.slab_start_z=slab_start_z
            current_reaction_system.slab_end_z=slab_end_z
         
        def set_volume(volume):
            """
            set the volume to be used in the acceptance probability of the reaction ensemble. This can be useful when using constraints, if the relevant volume is different from the box volume. If not used the default volume which is used, is the box volume.

            """
            current_reaction_system.volume=volume

        def print_acceptance_rate_configurational_moves(self):
            """
            prints the acceptance rate for the configuration changing moves
            """
            print "acceptance rate for configurational moves is", (1.0*accepted_configurational_MC_moves)/tried_configurational_MC_moves

        def set_non_interacting_type(non_interacting_type):
            """
            sets a type which is assumed to be non interacting in order to hide particles temporarily during a reaction trial move if they are to be deleted. The default value for this non_interacting type is 100. Please change this value if you intend to use a type 100 which has interactions. Please also note that particles in the current implementation of the Reaction Ensemble are only hidden with respect to Lennard-Jones interactions and Coulomb interactions. If there is for example a magnetic interaction hiding for this needs to be implemented in the code.
            """
            current_reaction_system.non_interacting_type=non_interacting_type

        def add(self,*args,**kwargs):
            """
            sets up a reaction in the forward and backward direction.
            
            Parameters
            ----------
            equilibrium_constant : float
                                   dimensionless (thermodynamic) equilibrium constant of the reaction, see Eq.[Keq].
            reactant_types : list 
                          a list of types of reactants in the reaction
            reactant_coefficients : list 
                                 a list of stoichiometric coefficients of the reactants in the same order as the list of their types
            product_types : list
                            a list of product types of the reaction.
            product_coefficients : list
                                   a list of stoichiometric coefficients of products of the reaction in the same order as the list of their types

            """            
            for k in self._required_keys_add():
                if k not in kwargs:
                    raise ValueError("At least the following keys have to be given as keyword arguments: " + self._required_keys_add().__str__() + " got " + kwargs.__str__())
                self._params[k] = kwargs[k]
            self._set_params_in_es_core_add()

        def _valid_keys_add(self):
            return "equilibrium_constant", "reactant_types", "reactant_coefficients", "product_types", "product_coefficients"

        def _required_keys_add(self):
            return ["equilibrium_constant", "reactant_types", "reactant_coefficients", "product_types", "product_coefficients"]

        def _set_params_in_es_core_add(self):
            
            cdef single_reaction* new_reaction =<single_reaction *>malloc(sizeof(single_reaction))
            #initialize values of reactant/ products in reaction to reasonable value
            new_reaction.len_reactant_types=0
            new_reaction.len_product_types=0
            
            new_reaction.equilibrium_constant=self._params["equilibrium_constant"]
            cdef int *reactant_types = <int *>malloc(len(self._params["reactant_types"]) * sizeof(int))
            for i in range(len(self._params["reactant_types"])):
                reactant_types[i]=self._params["reactant_types"][i]
            new_reaction.reactant_types=reactant_types
            new_reaction.len_reactant_types=len(self._params["reactant_types"]);
            
            
            cdef int *reactant_coefficients = <int *>malloc(len(self._params["reactant_coefficients"]) * sizeof(int))
            for i in range(len(self._params["reactant_coefficients"])):
                reactant_coefficients[i]=self._params["reactant_coefficients"][i]
            new_reaction.reactant_coefficients=reactant_coefficients
            
            cdef int *product_types = <int *>malloc(len(self._params["product_types"]) * sizeof(int))
            for i in range(len(self._params["product_types"])):
                product_types[i]=self._params["product_types"][i]
            new_reaction.product_types=product_types
            new_reaction.len_product_types=len(self._params["product_types"]);
            
            
            cdef int *product_coefficients = <int *>malloc(len(self._params["product_coefficients"]) * sizeof(int))
            for i in range(len(self._params["product_coefficients"])):
                product_coefficients[i]=self._params["product_coefficients"][i]
            new_reaction.product_coefficients=product_coefficients
            
            new_reaction.nu_bar=calculate_nu_bar(new_reaction.reactant_coefficients, new_reaction. len_reactant_types,  new_reaction.product_coefficients, new_reaction.len_product_types);
            
            #if everything is fine:
            current_reaction_system.reactions=<single_reaction**> realloc(current_reaction_system.reactions,sizeof(single_reaction*)*(current_reaction_system.nr_single_reactions+1)); #enlarge current_reaction_system
            current_reaction_system.reactions[current_reaction_system.nr_single_reactions]=new_reaction;
            current_reaction_system.nr_single_reactions+=1;
            
            #assign different types an index in a growing list that starts at and is incremented by 1 for each new type
            status=update_type_index(new_reaction.reactant_types, new_reaction.len_reactant_types, new_reaction.product_types, new_reaction.len_product_types);
            if(status==es_error):
                raise Exception("could not initialize gc particle list for types")
            
            #add reverse reaction
            _set_params_in_es_core_add_reverse(1.0/new_reaction.equilibrium_constant, new_reaction.product_types, new_reaction.len_product_types, new_reaction.product_coefficients, new_reaction.reactant_types, new_reaction.len_reactant_types, new_reaction.reactant_coefficients, -1*new_reaction.nu_bar )

        def default_charges(self,*args,**kwargs):
            """
            sets the charges of the particle types that are created. Note that it has to be called for each type that occurs in the reaction system individually.
            """
            for k in kwargs:
                if k in self._valid_keys_default_charge():
                    self._params[k] = kwargs[k]
                else:
                    raise KeyError("%s is not a vaild key" % k)
    
            self._validate_params_default_charge()
            
            for key in self._params["dictionary"]:
                    if(find_index_of_type(int(key))>=0):
                        current_reaction_system.charges_of_types[find_index_of_type(int(key))] = self._params["dictionary"][key]
            
        def _valid_keys_default_charge(self):
            return "dictionary"
        
        def _validate_params_default_charge(self):
            if(isinstance(self._params["dictionary"],dict)==False):
                raise ValueError("No dictionary for relation between types and default charges provided.")
        
        def reaction(self):
            """
            performs one randomly selected reaction of the provided reaction system
            """
            do_reaction()
        
        def do_global_mc_move_for_one_particle_of_type(type_mc):
            """
            performs a global mc move for one particle of type type_mc.
            If there are multiple types, that need to be moved, make sure to move them in a random order to avoid artefacts.
            """
            do_global_mc_move_for_particles_of_type(type_mc, -10, -10, 1, False)

        def print_status(self):
            """
            prints the status of the reaction ensemble, e.g. the used reactions
            """
            if(current_reaction_system.nr_single_reactions == 0):
                print("Reaction System is not initialized")
            else:
                print("Reaction System is the following:")
                print("reactant Types")
                for single_reaction_i in range(current_reaction_system.nr_single_reactions):
                    print "Reaction Nr.", single_reaction_i
                    print "reactant Types"
                    for i in range(current_reaction_system.reactions[single_reaction_i].len_reactant_types):
                        print current_reaction_system.reactions[single_reaction_i].reactant_types[i]

                    print "reactant coefficients"
                    for i in range(current_reaction_system.reactions[single_reaction_i].len_reactant_types):
                        print current_reaction_system.reactions[single_reaction_i].reactant_coefficients[i]
                    
                    print "Product Types"
                    for i in range(current_reaction_system.reactions[single_reaction_i].len_product_types):
                        print current_reaction_system.reactions[single_reaction_i].product_types[i]

                    print "Product coefficients"
                    for i in range(current_reaction_system.reactions[single_reaction_i].len_product_types):
                        print current_reaction_system.reactions[single_reaction_i].product_coefficients[i]
                    print "Equilibrium constant"
                    print current_reaction_system.reactions[single_reaction_i].equilibrium_constant
                    
            print "Reaction ensemble temperature:", current_reaction_system.temperature_reaction_ensemble
            print "Exclusion radius:", current_reaction_system.exclusion_radius
            if(check_reaction_ensemble()==es_error):
                raise ValueError("")
        
        def free(self):
            """
            Frees the reaction ensemble data structures in the core.
            """
            free_reaction_ensemble()

        #//////////////////////////Wang-Landau algorithm
        def add_collective_variable_degree_of_association(self,*args,**kwargs):
            """
            adds a reaction coordinate of the type degree of association
            
            Parameters
            ----------
            associated_type : int
                              type of the associated version of the species
                         
            min : float
                  minimum value of the collective variable
            max : float
                  maximum value of the collective variable
            corresponding_acid_types : list
                  list of the types of the version of the species
            
            """
            for k in kwargs:
                if k in self._valid_keys_add_collective_variable_degree_of_association():
                    self._params[k]=kwargs[k]
                else: KeyError("%s is not a valid key" %k)
            cdef collective_variable* new_collective_variable=<collective_variable*> calloc(1,sizeof(collective_variable))

            for k in self._required_keys_add_collective_variable_degree_of_association():
                if k not in kwargs:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self._required_keys_add_collective_variable_degree_of_association().__str__() + " got " + kwargs.__str__())
                self._params[k] = kwargs[k]

            new_collective_variable.associated_type=self._params["associated_type"]
            new_collective_variable.CV_minimum=self._params["min"]
            new_collective_variable.CV_maximum=self._params["max"]

            cdef int *corresponding_acid_types = <int *>malloc(len(self._params["corresponding_acid_types"]) * sizeof(int))
            for i in range(len(self._params["corresponding_acid_types"])):
                corresponding_acid_types[i]=self._params["corresponding_acid_types"][i]
                new_collective_variable.corresponding_acid_types=corresponding_acid_types
            new_collective_variable.nr_corresponding_acid_types=len(self._params["corresponding_acid_types"])


            current_wang_landau_system.collective_variables=<collective_variable**> realloc(current_wang_landau_system.collective_variables,sizeof(collective_variable*)*(current_wang_landau_system.nr_collective_variables+1))
            current_wang_landau_system.collective_variables[current_wang_landau_system.nr_collective_variables]=new_collective_variable
            current_wang_landau_system.nr_collective_variables+=1
            
            initialize_wang_landau()

        def _valid_keys_add_collective_variable_degree_of_association(self):
            return "associated_type", "min", "max", "corresponding_acid_types"

        def _required_keys_add_collective_variable_degree_of_association(self):
            return "associated_type", "min", "max", "corresponding_acid_types"

        def add_collective_variable_potential_energy(self,*args,**kwargs):
            """
            adds a reaction coordinate of the type potential energy
            
            Parameters
            ----------
            filename : str
                       filename of the energy boundary file which provides the potential energy boundaries (min E_pot, max E_pot) tabulated for all degrees of association. Make sure to only list the degrees of association which are used by the degree of association collective variable within this file. The energy boundary file can be created in a preliminary energy run. By the help of the functions :meth:`update_maximum_and_minimum_energies_at_current_state` and :meth:`write_out_preliminary_energy_run_results`. This file has to be obtained before being able to run a simulation with the energy as collective variable.
            delta : float
                    provides the discretization of the potential energy range. Only for small enough delta the results of the energy reweighted averages are correct. If delta is chosen too big there are discretization errors in the numerical integration which occurs during the energy reweighting process.
            """
            for k in kwargs:
                if k in self._valid_keys_add_collective_variable_potential_energy():
                    self._params[k]=kwargs[k]
                else: KeyError("%s is not a valid key" %k)
                
                for k in self._required_keys_add_collective_variable_potential_energy():
                    if k not in kwargs:
                        raise ValueError(
                            "At least the following keys have to be given as keyword arguments: " + self._required_keys_add_collective_variable_degree_of_association().__str__() + " got " + kwargs.__str__())
                    self._params[k] = kwargs[k]
            cdef collective_variable* new_collective_variable=<collective_variable*> calloc(1,sizeof(collective_variable)*(current_wang_landau_system.nr_collective_variables+1))

            new_collective_variable.energy_boundaries_filename=self._params["filename"]
            new_collective_variable.delta_CV=self._params["delta"]
            
            current_wang_landau_system.collective_variables=<collective_variable**> realloc(current_wang_landau_system.collective_variables,sizeof(collective_variable*)*(current_wang_landau_system.nr_collective_variables+1));
            current_wang_landau_system.collective_variables[current_wang_landau_system.nr_collective_variables]=new_collective_variable
            current_wang_landau_system.nr_collective_variables+=1;
            
            initialize_wang_landau()

        def _valid_keys_add_collective_variable_potential_energy(self):
            return "filename","delta"

        def _required_keys_add_collective_variable_potential_energy(self):
            return "filename","delta"

        def set_wang_landau_parameters(self,*args,**kwargs):
            """
            sets the final Wang-Landau parameter
            
            Parameters
            ----------
            final_wang_landau_parameter : float
                                          sets the final Wang-Landau parameter, which is the Wang-Landau parameter after which the simulation should stop.).
            wang_landau_steps : float
                                sets the number of Wang-Landau steps which are performed at once. Do not use too many Wang-Landau steps consequetively without having conformation changing steps in between. Number of Wang-Landau steps performed at once. This is for performance. It reduces the need for the interpreter to be called.
            full_path_to_output_filename : string
                                           sets the path to the output file of the Wang-Landau algorithm which contains the Wang-Landau potential
            do_not_sample_reaction_partition_function : bool
                                                        avoids sampling the Reaction ensemble partition function in the Wang-Landau algorithm. Therefore this option makes all degrees of association equally probable. This option may be used in the sweeping mode of the reaction ensemble, since the reaction ensemble partition function can be later added analytically.
            use_hybrid_monte_carlo : bool
                                     this is an experimental implementation only and per default it is turned off! Check the implementation again before using HMC here. Make sure not to use an MD thermostat in the case of using the Wang-Landau algorithm with Hybrid-Monte-Carlo moves. Wang-Landau moves with the Hybrid-Monte-Carlo moves are interesting for polymer systems since they avoid trapping in the energy reweighting case. However it is stressed here again that the implementation is experimental only. Sets whether the conformation changing Monte-Carlo moves should use a hybrid Monte Carlo scheme (use MD to propose new configurations and accept these proposed configurations with a probability proportional to :math:`\exp(-\\beta \\Delta E_\\text{pot})`).            
            """
            for k in kwargs:
                if k in self._valid_keys_set_wang_landau_parameters():
                    self._params[k]=kwargs[k]
                else: KeyError("%s is not a valid key" %k)
            
            current_wang_landau_system.final_wang_landau_parameter=self._params["final_wang_landau_parameter"]
            current_wang_landau_system.wang_landau_steps=self._params["wang_landau_steps"]
            current_wang_landau_system.output_filename=self._params["full_path_to_output_filename"]
            current_wang_landau_system.do_not_sample_reaction_partition_function=self._params["do_not_sample_reaction_partition_function"]
            current_wang_landau_system.use_hybrid_monte_carlo=self._params["use_hybrid_monte_carlo"]

        def _valid_keys_set_wang_landau_parameters(self):
            return "final_wang_landau_parameter", "wang_landau_steps", "full_path_to_output_filename", "do_not_sample_reaction_partition_function", "use_hybrid_monte_carlo"
            
        def load_wang_landau_checkpoint(self):
            """
            Loads the dumped wang landau potential file
            """
            load_wang_landau_checkpoint("checkpoint")
        def write_wang_landau_checkpoint(self):
            """
            Dumps the wang landau potential to a checkpoint file. Can be used to checkpoint the Wang-Landau histogram, potential, parameter and the number of executed trial moves
            """
            write_wang_landau_checkpoint("checkpoint")
            
        def update_maximum_and_minimum_energies_at_current_state(self):
            """
            records the minimum and maximum potential energy as a function of the degree of association in a preliminary Wang-Landau reaction ensemble simulation where the acceptance probability includes the factor :math:`\exp(-\\beta \\Delta E_{pot})`. The minimal and maximal potential energys which occur in the system are needed for the energy reweighting simulations wehere the factor :math:`\exp(-\\beta \\Delta E_{pot})` is not included in the acceptance probability in order to avoid choosing the wrong potential energy boundaries.
            """
            update_maximum_and_minimum_energies_at_current_state()
        
        def write_out_preliminary_energy_run_results(self):
            """
            this writes out the minimum and maximum potential energy as a function of the degree of association to a file. It requires that previously :meth:`update_maximum_and_minimum_energies_at_current_state` was used.
            """
            write_out_preliminary_energy_run_results("preliminary_energy_run_results")
            
        def do_reaction_wang_landau(self):
            """
            performs a reaction in the Wang-Landau reaction ensemble.
            """
            status_wang_landau=do_reaction_wang_landau()
            if(status_wang_landau<0):
                    raise Wang_Landau_has_converged("The Wang-Landau algorithm has converged.")

        def do_global_mc_move_for_one_particle_of_type_wang_landau(self,type_mc):
            """
            performs a global mc move for one particle of type type_mc (depending on the energy reweighting scheme)
            If there are multiple types, that need to be moved, make sure to move them in a random order to avoid artefacts.
            """
            do_global_mc_move_for_particles_of_type(type_mc, current_wang_landau_system.polymer_start_id,current_wang_landau_system.polymer_end_id, 1, True)

        def wang_landau_free(self):
            """
            frees the Wang-Landau data structures in the core
            
            """
            free_wang_landau()
        
        ##specify information for configuration changing monte carlo move
        property polymer_start_id:
            """
            Optional: since you might not want to change the configuration of your polymer, e.g. if you are trying to simulate a rigid conformation. Sets the start id of the polymer, optional. Should be set when you have a non fixed polymer and want it to be moved by MC trail moves in order to sample its configuration space.. MC moves for free particles and polymer particles may be very different.
            """
            def __set__(self, int start_id):
                current_wang_landau_system.polymer_start_id=start_id
            def __get__(self):
                        return current_wang_landau_system.polymer_start_id
        property polymer_end_id:
            """
            Optional: since you might not want to change the configuration of your polymer, e.g. if you are trying to simulate a rigid conformation. Sets the end id of the polymer, optional. Should be set when you have a non fixed polymer and want it to be moved by MC trail moves in order to sample its configuration space. MC moves for free particles and polymer particles may be very different.
            """
            def __set__(self, int end_id):
                current_wang_landau_system.polymer_end_id=end_id
            def __get__(self):
                return current_wang_landau_system.polymer_end_id
            
        property fix_polymer_monomers:
            """
            fixes the polymer monomers in the Monte Carlo moves
            """
            def __set__(self, bool fix_polymer):
                current_wang_landau_system.fix_polymer=fix_polymer
            def __get__(self):
                return current_wang_landau_system.fix_polymer

        
        #//////////////////////////constant pH ensemble
        def set_pH_core(self,pH):
                """
                sets the pH that the method assumes for the implicit pH bath
                """
                set_pH(pH)      
        def do_reaction_constant_pH(self):
            """
            performs a reaction according to the constant pH method
            """
            do_reaction_constant_pH()
        

