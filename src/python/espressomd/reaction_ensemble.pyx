include "myconfig.pxi"
from libc.stdlib cimport malloc, realloc, calloc

IF REACTION_ENSEMBLE:
    class Wang_Landau_has_converged(Exception):
        pass

    cdef class ReactionEnsemble:

        def __init__(self,*args,**kwargs):
            self._params=self.default_params()

            for k in self.required_keys():
                if k not in kwargs:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__() + " got " + kwargs.__str__())
                self._params[k] = kwargs[k]

            for k in kwargs:
                if k in self.valid_keys():
                    self._params[k] = kwargs[k]
                else:
                    raise KeyError("%s is not a vaild key" % k)
                
            self.validate_params()
            self._set_params_in_es_core()

        def default_params(self):
            return {"standard_pressure":0,
                    "temperature":1,
                    "exclusion_radius":0}

        def validate_params(self):
            return -1       

        def validate_params_add(self):
            return -1

        def valid_keys(self):
            return "standard_pressure", "temperature", "exclusion_radius"

        def required_keys(self):
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
            current_reaction_system.cyl_x=center_x
            current_reaction_system.cyl_y=center_y
            current_reaction_system.cyl_radius=radius_of_cylinder
            current_reaction_system.box_is_cylindric_around_z_axis=True
            
        def set_wall_constraints_in_z_direction(slab_start_z,slab_end_z):
            current_reaction_system.slab_start_z=slab_start_z
            current_reaction_system.slab_end_z=slab_end_z
         
        def print_acceptance_rate_configurational_moves(self):
                print "acceptance rate for configurational moves is", (1.0*accepted_configurational_MC_moves)/tried_configurational_MC_moves



        def add(self,*args,**kwargs):
             for k in self.required_keys_add():
                if k not in kwargs:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self.required_keys_add().__str__() + " got " + kwargs.__str__())
                self._params[k] = kwargs[k]
             self._set_params_in_es_core_add()

        def valid_keys_add(self):
            return "equilibrium_constant", "educt_types", "educt_coefficients", "product_types", "product_coefficients"

        def required_keys_add(self):
            return ["equilibrium_constant", "educt_types", "educt_coefficients", "product_types", "product_coefficients"]

        def default_params_add(self):
            
            return {}

        def _set_params_in_es_core_add(self):
            
            cdef single_reaction* new_reaction =<single_reaction *>malloc(sizeof(single_reaction))
            
            new_reaction.equilibrium_constant=self._params["equilibrium_constant"]
            cdef int *educt_types = <int *>malloc(len(self._params["educt_types"]) * sizeof(int))
            for i in range(len(self._params["educt_types"])):
                educt_types[i]=self._params["educt_types"][i]
            new_reaction.educt_types=educt_types
            new_reaction.len_educt_types=len(self._params["educt_types"]);
            
            
            cdef int *educt_coefficients = <int *>malloc(len(self._params["educt_coefficients"]) * sizeof(int))
            for i in range(len(self._params["educt_coefficients"])):
                educt_coefficients[i]=self._params["educt_coefficients"][i]
            new_reaction.educt_coefficients=educt_coefficients
            
            cdef int *product_types = <int *>malloc(len(self._params["product_types"]) * sizeof(int))
            for i in range(len(self._params["product_types"])):
                product_types[i]=self._params["product_types"][i]
            new_reaction.product_types=product_types
            new_reaction.len_product_types=len(self._params["product_types"]);
            
            
            cdef int *product_coefficients = <int *>malloc(len(self._params["product_coefficients"]) * sizeof(int))
            for i in range(len(self._params["product_coefficients"])):
                product_coefficients[i]=self._params["product_coefficients"][i]
            new_reaction.product_coefficients=product_coefficients
            
            new_reaction.nu_bar=calculate_nu_bar(new_reaction.educt_coefficients, new_reaction. len_educt_types,  new_reaction.product_coefficients, new_reaction.len_product_types);
            
            #if everything is fine:
            current_reaction_system.reactions=<single_reaction**> realloc(current_reaction_system.reactions,sizeof(single_reaction*)*(current_reaction_system.nr_single_reactions+1)); #enlarge current_reaction_system
            current_reaction_system.reactions[current_reaction_system.nr_single_reactions]=new_reaction;
            current_reaction_system.nr_single_reactions+=1;
            
            #assign different types an index in a growing list that starts at and is incremented by 1 for each new type
            update_type_index(new_reaction.educt_types, new_reaction.len_educt_types, new_reaction.product_types, new_reaction.len_product_types);
            
            #add reverse reaction
            self._set_params_in_es_core_add_reverse(1.0/new_reaction.equilibrium_constant, new_reaction.product_types, new_reaction.len_product_types, new_reaction.product_coefficients, new_reaction.educt_types, new_reaction.len_educt_types, new_reaction.educt_coefficients, -1*new_reaction.nu_bar )

        cdef _set_params_in_es_core_add_reverse(self, equilibrium_constant, int* educt_types, len_educt_types, int* educt_coefficients, int* product_types, len_product_types, int* product_coefficients, nu_bar):
            #add reverse reaction based on last reaction that was added            
            cdef single_reaction* new_reaction =<single_reaction *>malloc(sizeof(single_reaction))     
            new_reaction.equilibrium_constant=equilibrium_constant
            new_reaction.educt_types=educt_types
            new_reaction.len_educt_types=len_educt_types
            new_reaction.educt_coefficients=educt_coefficients
            new_reaction.product_types=product_types
            new_reaction.len_product_types=len_product_types
            new_reaction.product_coefficients=product_coefficients
            new_reaction.nu_bar=nu_bar
            #if everything is fine:
            current_reaction_system.reactions=<single_reaction**> realloc(current_reaction_system.reactions,sizeof(single_reaction*)*(current_reaction_system.nr_single_reactions+1)); #enlarge current_reaction_system
            current_reaction_system.reactions[current_reaction_system.nr_single_reactions]=new_reaction;
            current_reaction_system.nr_single_reactions+=1;
            
        def default_charges(self,*args,**kwargs):
                for k in kwargs:
                    if k in self.valid_keys_default_charge():
                        self._params[k] = kwargs[k]
                    else:
                        raise KeyError("%s is not a vaild key" % k)
        
                self.validate_params_default_charge()
                
                for key in self._params["dictionary"]:
                        current_reaction_system.charges_of_types[find_index_of_type(int(key))] = self._params["dictionary"][key]
                
        def valid_keys_default_charge(self):
            return "dictionary"
        
        def validate_params_default_charge(self):
            if(isinstance(self._params["dictionary"],dict)==False):
                raise ValueError("No dictionary for relation between types and default charges provided.")
        
        def reaction(self):
            do_reaction()

        def print_status(self):
            if(current_reaction_system.nr_single_reactions == 0):
                print("Reaction System is not initialized")
            else:
                print("Reaction System is the following:")
                print("Educt Types")
                for single_reaction_i in range(current_reaction_system.nr_single_reactions):
                    print "Reaction Nr.", single_reaction_i
                    print "Educt Types"
                    for i in range(current_reaction_system.reactions[single_reaction_i].len_educt_types):
                        print current_reaction_system.reactions[single_reaction_i].educt_types[i]

                    print "Educt coefficients"
                    for i in range(current_reaction_system.reactions[single_reaction_i].len_educt_types):
                        print current_reaction_system.reactions[single_reaction_i].educt_coefficients[i]
                    
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
            check_reaction_ensemble()
            
        #//////////////////////////Wang-Landau algorithm
        def add_collective_variable_degree_of_association(self,*args,**kwargs):
            for k in kwargs:
                if k in self.valid_keys_add_collective_variable_degree_of_association():
                    self._params[k]=kwargs[k]
                else: KeyError("%s is not a valid key" %k)
            cdef collective_variable* new_collective_variable=<collective_variable*> calloc(1,sizeof(collective_variable))

            for k in self.required_keys_add_collective_variable_degree_of_association():
                if k not in kwargs:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self.required_keys_add_collective_variable_degree_of_association().__str__() + " got " + kwargs.__str__())
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

        def valid_keys_add_collective_variable_degree_of_association(self):
            return "associated_type", "min", "max", "corresponding_acid_types"

        def required_keys_add_collective_variable_degree_of_association(self):
            return "associated_type", "min", "max", "corresponding_acid_types"

        def add_collective_variable_potential_energy(self,*args,**kwargs):
            for k in kwargs:
                if k in self.valid_keys_add_collective_variable_potential_energy():
                    self._params[k]=kwargs[k]
                else: KeyError("%s is not a valid key" %k)
                
                for k in self.required_keys_add_collective_variable_potential_energy():
                    if k not in kwargs:
                        raise ValueError(
                            "At least the following keys have to be given as keyword arguments: " + self.required_keys_add_collective_variable_degree_of_association().__str__() + " got " + kwargs.__str__())
                    self._params[k] = kwargs[k]
            cdef collective_variable* new_collective_variable=<collective_variable*> calloc(1,sizeof(collective_variable)*(current_wang_landau_system.nr_collective_variables+1))

            new_collective_variable.energy_boundaries_filename=self._params["filename"]
            new_collective_variable.delta_CV=self._params["delta"]
            
            current_wang_landau_system.collective_variables=<collective_variable**> realloc(current_wang_landau_system.collective_variables,sizeof(collective_variable*)*(current_wang_landau_system.nr_collective_variables+1));
            current_wang_landau_system.collective_variables[current_wang_landau_system.nr_collective_variables]=new_collective_variable
            current_wang_landau_system.nr_collective_variables+=1;
            
            initialize_wang_landau()

        def valid_keys_add_collective_variable_potential_energy(self):
            return "filename","delta"

        def required_keys_add_collective_variable_potential_energy(self):
            return "filename","delta"

        def set_wang_landau_parameters(self,*args,**kwargs):
            for k in kwargs:
                if k in self.valid_keys_set_wang_landau_parameters():
                    self._params[k]=kwargs[k]
                else: KeyError("%s is not a valid key" %k)
            
            current_wang_landau_system.final_wang_landau_parameter=self._params["final_wang_landau_parameter"]
            current_wang_landau_system.wang_landau_steps=self._params["wang_landau_steps"]
            current_wang_landau_system.output_filename=self._params["full_path_to_output_filename"]
            current_wang_landau_system.do_not_sample_reaction_partition_function=self._params["do_not_sample_reaction_partition_function"]
            current_wang_landau_system.use_hybrid_monte_carlo=self._params["use_hybrid_monte_carlo"]

        def valid_keys_set_wang_landau_parameters(self):
            return "final_wang_landau_parameter", "wang_landau_steps", "full_path_to_output_filename", "do_not_sample_reaction_partition_function", "use_hybrid_monte_carlo"
            
        def load_wang_landau_checkpoint(self):
            load_wang_landau_checkpoint("checkpoint")
        def write_wang_landau_checkpoint(self):
            write_wang_landau_checkpoint("checkpoint")
            
        def update_maximum_and_minimum_energies_at_current_state(self):
            update_maximum_and_minimum_energies_at_current_state()
        
        def write_out_preliminary_energy_run_results(self):
            write_out_preliminary_energy_run_results("preliminary_energy_run_results")
            
        
        ##specify information for configuration changing monte carlo move
        property counter_ion_type:
            def __set__(self, int c_type):
                current_wang_landau_system.counter_ion_type=c_type
            def __get__(self):
                return current_wang_landau_system.counter_ion_type
        property polymer_start_id:
            def __set__(self, int start_id):
                current_wang_landau_system.polymer_start_id=start_id
            def __get__(self):
                        return current_wang_landau_system.polymer_start_id
        property polymer_end_id:
            def __set__(self, int end_id):
                current_wang_landau_system.polymer_end_id=end_id
            def __get__(self):
                return current_wang_landau_system.polymer_end_id
            
        property fix_polymer_monomers:
            def __set__(self, bool fix_polymer):
                current_wang_landau_system.fix_polymer=fix_polymer
            def __get__(self):
                return current_wang_landau_system.fix_polymer

        def do_reaction_wang_landau(self):
            status_wang_landau=do_reaction_wang_landau()
            if(status_wang_landau<0):
                    raise Wang_Landau_has_converged("The Wang-Landau algorithm has converged.")

        #//////////////////////////constant pH ensemble
        def set_pH_core(self,pH):
                set_pH(pH)      
        def do_reaction_constant_pH(self):
            do_reaction_constant_pH()
        

