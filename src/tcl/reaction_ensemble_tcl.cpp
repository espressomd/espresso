#include "parser.hpp"

#ifdef REACTION_ENSEMBLE
#include "reaction_ensemble.hpp"
#include "particle_data.hpp" //for particle creation, modification

//function declartion

int tclcommand_reaction_ensemble_print_status(Tcl_Interp *interp){
	char buffer[3000];
	if(current_reaction_system.nr_single_reactions == 0){
		sprintf(buffer,"Reaction System is not initialized\n");
		Tcl_AppendResult(interp, buffer, (char *)NULL);
	}else{
		sprintf(buffer,"Reaction System is the following:\n");
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		for(int single_reaction_i=0;single_reaction_i<current_reaction_system.nr_single_reactions;single_reaction_i++){
			sprintf(buffer, "#Reaction %d# \n", single_reaction_i);
			Tcl_AppendResult(interp, buffer, "educt types:\n", (char *)NULL);
			//educt types and coefficients
			for(int i=0;i<current_reaction_system.reactions[single_reaction_i]->len_educt_types;i++){
				sprintf(buffer, "%d " , current_reaction_system.reactions[single_reaction_i]->educt_types[i]);
				Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
			}
			sprintf(buffer, "\neduct coefficients: \n");
			Tcl_AppendResult(interp, buffer, (char *)NULL);
			for(int i=0;i<current_reaction_system.reactions[single_reaction_i]->len_educt_types;i++){
				sprintf(buffer, "%d " , current_reaction_system.reactions[single_reaction_i]->educt_coefficients[i]);
				Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
			}
				
    			Tcl_AppendResult(interp, "\n", (char *)NULL);
			//product types and coefficients
			Tcl_AppendResult(interp, "product types:\n", (char *)NULL);
			for(int i=0;i<current_reaction_system.reactions[single_reaction_i]->len_product_types;i++){
				sprintf(buffer, "%d " , current_reaction_system.reactions[single_reaction_i]->product_types[i]);
				Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
			}
			sprintf(buffer, "\nproduct coefficients: \n");
			Tcl_AppendResult(interp, buffer, (char *)NULL);
			for(int i=0;i<current_reaction_system.reactions[single_reaction_i]->len_product_types;i++){
				sprintf(buffer, "%d " , current_reaction_system.reactions[single_reaction_i]->product_coefficients[i]);
				Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
			}				
			//reaction_constant
			sprintf(buffer, "\nequilibrium constant: %f " , current_reaction_system.reactions[single_reaction_i]->equilibrium_constant);
			Tcl_AppendResult(interp, buffer, "\n", (char *)NULL);
		}

		//temperature
		sprintf(buffer, "\nreaction ensemble temperature: %f " , current_reaction_system.temperature_reaction_ensemble);
		Tcl_AppendResult(interp, buffer, "\n", (char *)NULL);
		//exclusion radius
		sprintf(buffer, "\nexclusion radius: %f " , current_reaction_system.exclusion_radius);
		Tcl_AppendResult(interp, buffer, "\n", (char *)NULL);
		//volume
		sprintf(buffer, "\nvolume: %f " , current_reaction_system.volume);
		Tcl_AppendResult(interp, buffer, "\n", (char *)NULL);

		if(check_reaction_ensemble()==ES_ERROR)
                        return TCL_ERROR;

	}
	return TCL_OK;
}

int tclcommand_add_reaction(Tcl_Interp *interp, int argc, char **argv){
	single_reaction* new_reaction =(single_reaction *) malloc(sizeof(single_reaction));
	//initialize values of educt/ products in reaction to reasonable value 
	new_reaction->len_educt_types=0;
	new_reaction->len_product_types=0;
	
	argc -= 1; argv += 1;
	if(ARG1_IS_S("equilibrium_constant")){
		ARG_IS_D(2, new_reaction->equilibrium_constant);
		argc-=2; argv+=2;
	}
	if(ARG1_IS_S("educt_types")){
		argc-=1; argv+=1;
		int next_added_type;
		int* educt_types=NULL;
		int educt_type_counter=0;
		while(ARG_IS_I(1,next_added_type)){
			educt_types=(int*) realloc(educt_types,sizeof(int)*(educt_type_counter+1));
			educt_types[educt_type_counter]=next_added_type;
			educt_type_counter+=1;
			//check for terminus of string
			if(argc<3) {
				break;
			}
			argc-=1; argv+=1;
		}
		new_reaction->len_educt_types=educt_type_counter;
		new_reaction->educt_types=educt_types;
	}
	if(ARG1_IS_S("educt_coefficients")){
		argc-=1; argv+=1;
		int next_added_type_coeff;
		int* educt_coefficients=NULL;
		int educt_type_counter=0;
		while(ARG_IS_I(1,next_added_type_coeff)&& argc>1){
			educt_coefficients=(int*) realloc(educt_coefficients,sizeof(int)*(educt_type_counter+1));
			educt_coefficients[educt_type_counter]=next_added_type_coeff;
			educt_type_counter+=1;
			//check for terminus of string
			if(argc<3) {
				break;
			}
			argc-=1; argv+=1;
		}
		new_reaction->educt_coefficients=educt_coefficients;
	}
	if(ARG1_IS_S("product_types")){
		argc-=1; argv+=1;
		int next_added_type;
		int* product_types=NULL;
		int product_type_counter=0;
		while(ARG_IS_I(1,next_added_type)){
			product_types=(int*) realloc(product_types,sizeof(int)*(product_type_counter+1));
			product_types[product_type_counter]=next_added_type;
			product_type_counter+=1;
			//check for terminus of string
			if(argc<3) {
				break;
			}
			argc-=1; argv+=1;
		}
		new_reaction->len_product_types=product_type_counter;
		new_reaction->product_types=product_types;
	}
	if(ARG1_IS_S("product_coefficients")){
		argc-=1; argv+=1;
		int next_added_type_coeff;
		int* product_coefficients=NULL;
		int product_type_counter=0;
		while(ARG_IS_I(1,next_added_type_coeff)&& argc>1){
			product_coefficients=(int*) realloc(product_coefficients,sizeof(int)*(product_type_counter+1));
			product_coefficients[product_type_counter]=next_added_type_coeff;
			product_type_counter+=1;
			new_reaction->product_coefficients=product_coefficients;
			//check for terminus of string
			if(argc<3) {
				break;
			}
			argc-=1; argv+=1;
		}
	
	}
	
	//if everything is fine:
	new_reaction->nu_bar=calculate_nu_bar(new_reaction->educt_coefficients, new_reaction-> len_educt_types,  new_reaction->product_coefficients, new_reaction->len_product_types);
	current_reaction_system.reactions=(single_reaction**) realloc(current_reaction_system.reactions,sizeof(single_reaction*)*(current_reaction_system.nr_single_reactions+1)); //enlarge current_reaction_system
	current_reaction_system.reactions[current_reaction_system.nr_single_reactions]=new_reaction;
	current_reaction_system.nr_single_reactions+=1;
	//assign different types an index in a growing list that starts at and is incremented by 1 for each new type
	update_type_index(new_reaction->educt_types, new_reaction->len_educt_types, new_reaction->product_types, new_reaction->len_product_types);
	return TCL_OK;
}

int tclcommand_add_reaction_coordinate(Tcl_Interp *interp, int argc, char **argv) {
	argc-=1;
	argv+=1;
	collective_variable* new_collective_variable=(collective_variable*) calloc(1,sizeof(collective_variable));
	if(ARG1_IS_S("degree_of_association")){
		argc-=1;
		argv+=1;
		if(ARG1_IS_S("associated_type")){
			argc-=1;
			argv+=1;
			ARG_IS_I(1,new_collective_variable->associated_type);
			argc-=1;
			argv+=1;
		}
		if(ARG1_IS_S("min")){
			argc-=1;
			argv+=1;
			ARG_IS_D(1,new_collective_variable->CV_minimum);
			argc-=1;
			argv+=1;
		}
		if(ARG1_IS_S("max")){
			argc-=1; argv+=1;
			ARG_IS_D(1,new_collective_variable->CV_maximum);
		}

		argc-=1; argv+=1;
		if(ARG1_IS_S("corresponding_acid_types")){
			argc-=1; argv+=1;
			int next_added_acid_type;
			int* corresponding_acid_types=NULL;
			int corresponding_type_counter=0;
			while(ARG_IS_I(1,next_added_acid_type)){
				corresponding_acid_types=(int*) realloc(corresponding_acid_types,sizeof(int)*(corresponding_type_counter+1));
				corresponding_acid_types[corresponding_type_counter]=next_added_acid_type;
				corresponding_type_counter+=1;
				//check for terminus of string
				if(argc<3) {
					break;
				}
				argc-=1; argv+=1;
			}
			new_collective_variable->corresponding_acid_types=corresponding_acid_types;
			new_collective_variable->nr_corresponding_acid_types=corresponding_type_counter;
		}
	}

	if(ARG1_IS_S("energy")){
		//needs to be called after all other collective variables are known
		argc-=1; argv+=1;
		if(ARG1_IS_S("filename")){ //full path to file which saves the energies in the format nbar_i \t nbar_j \t ... \t energy_min \t energy_max
			argc-=1; argv+=1;
			char* energy_boundaries_filename =strdup(argv[1]);
			new_collective_variable->energy_boundaries_filename=energy_boundaries_filename;
		}else{
			return TCL_ERROR;
		}
		argc-=1; argv+=1;
		if(ARG1_IS_S("delta")){
			argc-=1; argv+=1;
			ARG_IS_D(1,new_collective_variable->delta_CV);
		}else{
			return TCL_ERROR;
		}
	}

	current_wang_landau_system.collective_variables=(collective_variable**) realloc(current_wang_landau_system.collective_variables,sizeof(collective_variable*)*(current_wang_landau_system.nr_collective_variables+1));
	current_wang_landau_system.collective_variables[current_wang_landau_system.nr_collective_variables]=new_collective_variable;
	current_wang_landau_system.nr_collective_variables+=1;
	int return_code =TCL_OK;
	if (initialize_wang_landau()==ES_ERROR)
		return_code=TCL_ERROR;
	return TCL_OK;

}

int tclcommand_reaction_ensemble(ClientData data, Tcl_Interp *interp, int argc, char **argv){
	int err = TCL_OK;
	bool provided_unknown_command=true;
	if(argc == 1){
		provided_unknown_command=false;
		err= tclcommand_reaction_ensemble_print_status(interp);
	} else {
		if(ARG1_IS_S("do")) {
			provided_unknown_command=false;
			do_reaction();
		} else { //for performance reasons skip the other checks if do is found as argument
			if( ARG1_IS_S("add_reaction") ) {
				provided_unknown_command=false;
				tclcommand_add_reaction(interp, argc, argv);
		
			}
			if( ARG1_IS_S("temperature_reaction_ensemble") ){
				provided_unknown_command=false;
				argc-=1; argv+=1;
				ARG_IS_D(1,current_reaction_system.temperature_reaction_ensemble);
				
				//setting a volume is a side effect, sets the default volume of the reaction ensemble as the volume of the cuboid simulation box. this volume can be altered by the command "reaction ensemble volume <volume>" if one wants to simulate e.g. in a system with constraint (e.g. cuboid box with cylinder constraint, so that the particles are only contained in the volume of the cylinder)
				if(current_reaction_system.volume<0)
					set_cuboid_reaction_ensemble_volume();
			}
			if( ARG1_IS_S("exclusion_radius") ){
				provided_unknown_command=false;
				argc-=1; argv+=1;
				ARG_IS_D(1,current_reaction_system.exclusion_radius);
			}
			
			if(ARG1_IS_S("set_default_charge_of_type")) {
				provided_unknown_command=false;
				//needs to be called for each type individually after reaction was added
				int type;
				ARG_IS_I(2,type);
				double charge;
				ARG_IS_D(3,charge);
				int index_of_type=find_index_of_type(type);
				if(index_of_type>=0){
					current_reaction_system.charges_of_types[index_of_type]=charge;
				}
				
			}
		
			if( ARG1_IS_S("free_memory")) {
				provided_unknown_command=false;
				free_reaction_ensemble();
			}
			
			if (ARG1_IS_S("monte_carlo_move_for_type")) {
				provided_unknown_command=false;
				int mc_type;
				ARG_IS_I(2,mc_type);
				
				int particle_number_of_type;
				number_of_particles_with_type(mc_type, &(particle_number_of_type));
				for(int i=0; i<particle_number_of_type; i++) {
					do_global_mc_move_for_one_particle_of_type(mc_type,-10,-10);
				}
			}
		
			if( ARG1_IS_S("standard_pressure_in_simulation_units")) {
				provided_unknown_command=false;
				ARG_IS_D(2,current_reaction_system.standard_pressure_in_simulation_units);
			}
			if( ARG1_IS_S("volume")) {
				provided_unknown_command=false;
				ARG_IS_D(2,current_reaction_system.volume);
			}
			if( ARG1_IS_S("cylindrical_constraint_in_z_direction")) {
				provided_unknown_command=false;
				argc-=1; argv+=1;
				if(ARG1_IS_S("center")){
					argc-=1; argv+=1;
					ARG_IS_D(1,current_reaction_system.cyl_x);
					argc-=1; argv+=1;
					ARG_IS_D(1,current_reaction_system.cyl_y);
				}
				argc-=1; argv+=1;
				if(ARG1_IS_S("radius")){
					argc-=1; argv+=1;
					ARG_IS_D(1,current_reaction_system.cyl_radius);
				}
				current_reaction_system.box_is_cylindric_around_z_axis=true;
			}
			if( ARG1_IS_S("wall_constraints_in_z_direction")) {
				provided_unknown_command=false;
				argc-=1; argv+=1;
				if(ARG1_IS_S("start")){
					argc-=1; argv+=1;
					ARG_IS_D(1,current_reaction_system.slab_start_z);
				}
				argc-=1; argv+=1;
				if(ARG1_IS_S("end")){
					argc-=1; argv+=1;
					ARG_IS_D(1,current_reaction_system.slab_end_z);
				}
				current_reaction_system.box_has_wall_constraints=true;
			}			
			
			if(ARG1_IS_S("print_acceptance_rate")){
				provided_unknown_command=false;
				char buffer[3000];
				sprintf(buffer, "acceptance rate for configurational moves is %f", (1.0*accepted_configurational_MC_moves)/tried_configurational_MC_moves);
				Tcl_AppendResult(interp, buffer, "\n", (char *)NULL);
			}
		}
		///////////////////////////////////////////// Wang-Landau algorithm
		if (ARG1_IS_S("wang_landau")){ // for performance reasons skip other tests
			argc-=1;
			argv+=1;
			if(ARG1_IS_S("do")) {
				provided_unknown_command=false;
				int status_wang_landau=do_reaction_wang_landau();
				if (status_wang_landau <0)
					exit(0); // Wang-Landau simulation has converged, end simulation
			} else { //for performance reasons skip the other checks if do is found as argument
				//do other checks here

				if(ARG1_IS_S("add")){
					provided_unknown_command=false;
					tclcommand_add_reaction_coordinate(interp, argc, argv);
				}

				if(ARG1_IS_S("final_wang_landau_parameter")) {
					provided_unknown_command=false;
					argc-=1; argv+=1;
					ARG_IS_D(1,current_wang_landau_system.final_wang_landau_parameter);
				}
				
				if(ARG1_IS_S("wang_landau_steps")) {
					provided_unknown_command=false;
					argc-=1; argv+=1;
					ARG_IS_I(1,current_wang_landau_system.wang_landau_steps);
				}		
	
				if(ARG1_IS_S("full_path_to_output_filename")){
					provided_unknown_command=false;
					argc-=1; argv+=1;
					char* output_filename =strdup(argv[1]);
					current_wang_landau_system.output_filename=output_filename;
				}

				if(ARG1_IS_S("free")) {
					provided_unknown_command=false;
					//needs to be called after all observables are added
					free_wang_landau();
				}
				
				if(ARG1_IS_S("update_maximum_and_minimum_energies_at_current_state")){
					provided_unknown_command=false;
					update_maximum_and_minimum_energies_at_current_state();
				}
				if(ARG1_IS_S("write_out_preliminary_energy_run_results")){
					provided_unknown_command=false;
					argc-=1; argv+=1;
					write_out_preliminary_energy_run_results(argv[1]);
				}
				if(ARG1_IS_S("do_not_sample_reaction_partition_function")){
					provided_unknown_command=false;
		current_wang_landau_system.do_not_sample_reaction_partition_function=true;
				}
				if(ARG1_IS_S("counter_ion_type")){
					provided_unknown_command=false;
					argc-=1; argv+=1;
					ARG_IS_I(1,current_wang_landau_system.counter_ion_type);
				}
				if(ARG1_IS_S("polymer_start_id")){
					provided_unknown_command=false;
					argc-=1; argv+=1;
					ARG_IS_I(1,current_wang_landau_system.polymer_start_id);
					if(current_wang_landau_system.polymer_start_id<0){
						printf("negative start id for polymer");
						return TCL_ERROR;
					}
					argc-=1; argv+=1;
				}
				if(ARG1_IS_S("polymer_end_id")){
					provided_unknown_command=false;
					argc-=1; argv+=1;
					ARG_IS_I(1,current_wang_landau_system.polymer_end_id);
					if(current_wang_landau_system.polymer_end_id<0){
						printf("negative end id for polymer");
						return TCL_ERROR;
					}
					argc-=1; argv+=1;			
					
				}
				if(ARG1_IS_S("fix_polymer_monomers")){
					provided_unknown_command=false;
					current_wang_landau_system.fix_polymer=true;
				}
				if(ARG1_IS_S("use_hybrid_monte_carlo")){
					provided_unknown_command=false;
					current_wang_landau_system.use_hybrid_monte_carlo=true;
					printf("Make sure to use a bigger time step when using hybrid monte carlo steps!\n");
				}
				if(ARG1_IS_S("write_wang_landau_checkpoint")){
					provided_unknown_command=false;					
					argc-=1; argv+=1;
					char* identifier =strdup(argv[1]);
					write_wang_landau_checkpoint(identifier);
				}
				if(ARG1_IS_S("load_wang_landau_checkpoint")){
					provided_unknown_command=false;
					argc-=1; argv+=1;
					char* identifier =strdup(argv[1]);
					load_wang_landau_checkpoint(identifier);
				}

				if(ARG1_IS_S("monte_carlo_move_for_type")) {
					provided_unknown_command=false;
					int mc_type;
					ARG_IS_I(2,mc_type);
					int particle_number_of_type;
					number_of_particles_with_type(mc_type, &(particle_number_of_type));
					for(int i=0; i<std::min(30,particle_number_of_type); i++) {
						do_global_mc_move_for_one_particle_of_type_wang_landau(mc_type,current_wang_landau_system.polymer_start_id,current_wang_landau_system.polymer_end_id);
					}
					
					if(current_wang_landau_system.use_hybrid_monte_carlo==true){
						do_HMC_move_wang_landau();
					}
					
				}

				
				
			}
		}
		
		
		
		///////////////////////////////////////////// constant_pH
		if (ARG1_IS_S("constant_pH")){
			argc-=1; argv+=1;
			if(ARG1_IS_S("do")) {
				provided_unknown_command=false;
				do_reaction_constant_pH();
			} else {
				if(ARG1_IS_S("pH")){
					provided_unknown_command=false;				
					argc-=1; argv+=1;
					ARG_IS_D(1,constant_pH);
				}
			}			
		
		}
		
	}

	if (provided_unknown_command==true){
		printf("Provided unknown parameter to reaction ensemble.\n");
		return TCL_ERROR;	
	}
		
	return gather_runtime_errors(interp,err);
}

#endif
