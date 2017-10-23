/** reaction ensemble method according to smith94x for the reaction ensemble at constant volume and temperature, for the reaction ensemble at constant pressure additionally employ a barostat!
*NOTE: a chemical reaction consists of a forward and backward reaction. Here both reactions have to be defined seperately.
*The extent of the reaction is here chosen to be +1.
*If the reaction trial move for a dissociation of HA is accepted then there is one more dissociated ion pair H+ and A-
*/

/** @file */ 

#include "reaction_ensemble.hpp"
#include "random.hpp" //for random numbers
#include "energy.hpp"	//for energies
#include "external_potential.hpp" //for energies
#include "global.hpp" //for access to global variables
#include "particle_data.hpp" //for particle creation, modification
#include "statistics.hpp" //for distto
#include "integrate.hpp" //for time_step
#include <stdio.h> //for getline()
#include <iostream> //for std::cout
#include <fstream> //for std::ifstream, std::ofstream for input output into files
#include "utils.hpp" // for PI and random vectors
#include "partCfg_global.hpp"

namespace ReactionEnsemble{

ReactionEnsemble::ReactionEnsemble(){}

ReactionEnsemble::~ReactionEnsemble(){
    this->free_reaction_ensemble();
    if(this->m_current_wang_landau_system.len_histogram!=0)
        this->free_wang_landau();
}

/**
* Performs a randomly selected reaction in the reaction ensemble
*/
int ReactionEnsemble::do_reaction(int reaction_steps){
    for (int i = 0; i< reaction_steps; i++){
        int reaction_id=i_random(m_current_reaction_system.nr_single_reactions);
    	generic_oneway_reaction(reaction_id, reaction_ensemble_mode);
    }
	return 0;
}

/**
* Adds a reaction to the reaction system
*/
void ReactionEnsemble::add_reaction(double equilibrium_constant, std::vector<int> _reactant_types, std::vector<int> _reactant_coefficients, std::vector<int> _product_types, std::vector<int> _product_coefficients){
            single_reaction* new_reaction =(single_reaction *)malloc(sizeof(single_reaction));
            int len_reactant_types=_reactant_types.size();
            int len_product_types=_product_types.size();
            new_reaction->len_reactant_types=len_reactant_types;
            new_reaction->len_product_types=len_product_types;
            
            new_reaction->equilibrium_constant=equilibrium_constant;
            int *reactant_types = (int *)malloc(len_reactant_types * sizeof(int));
            for(int i=0; i< len_reactant_types; i++)
                reactant_types[i]=_reactant_types[i];
            new_reaction->reactant_types=reactant_types;
            new_reaction->len_reactant_types=len_reactant_types;
            
            
            int* reactant_coefficients = (int*) malloc(len_reactant_types * sizeof(int));
            for(int i=0; i< len_reactant_types; i++)
                reactant_coefficients[i]=_reactant_coefficients[i];
            new_reaction->reactant_coefficients=reactant_coefficients;
            
            int* product_types = (int *)malloc(len_product_types * sizeof(int));
            for(int i=0; i< len_product_types; i++)
                product_types[i]=_product_types[i];
            new_reaction->product_types=product_types;
            new_reaction->len_product_types=len_product_types;
            
            
            int* product_coefficients = (int *)malloc(len_product_types * sizeof(int));
            for(int i=0; i< len_product_types; i++)
                product_coefficients[i]=_product_coefficients[i];
            new_reaction->product_coefficients=product_coefficients;
            
            new_reaction->nu_bar=calculate_nu_bar(new_reaction->reactant_coefficients, new_reaction->len_reactant_types,  new_reaction->product_coefficients, new_reaction->len_product_types);
            
            //if everything is fine:
            m_current_reaction_system.reactions=(single_reaction**) realloc(m_current_reaction_system.reactions,sizeof(single_reaction*)*(m_current_reaction_system.nr_single_reactions+1)); //enlarge RE.m_current_reaction_system
            m_current_reaction_system.reactions[m_current_reaction_system.nr_single_reactions]=new_reaction;
            m_current_reaction_system.nr_single_reactions+=1;
            
            //assign different types an index in a growing list that starts at and is incremented by 1 for each new type
            int status=update_type_index(new_reaction->reactant_types, new_reaction->len_reactant_types, new_reaction->product_types, new_reaction->len_product_types);
            if(status==ES_ERROR)
                 throw std::runtime_error("could not initialize gc particle list for types\n");

}

/**
* Checks whether all necessary variables for the reaction ensemble have been set.
*/
int ReactionEnsemble::check_reaction_ensemble(){
	/**checks the reaction_ensemble struct for valid parameters */
	int check_is_successfull =ES_OK;
	if(m_current_reaction_system.standard_pressure_in_simulation_units<0 and not (std::abs(m_constant_pH-(-10)) > std::numeric_limits<double>::epsilon() ) ){
	    throw std::runtime_error("Please initialize your reaction ensemble standard pressure before calling initialize.\n");
		check_is_successfull=ES_ERROR;
	}
	if(m_current_reaction_system.temperature_reaction_ensemble<0){
	    throw std::runtime_error("Temperatures cannot be negative. Please provide a temperature (in k_B T) to the simulation. Normally it should be 1.0. This will be used directly to calculate beta:=1/(k_B T) which occurs in the exp(-beta*E)\n");
		check_is_successfull=ES_ERROR;
	}
	#ifdef ELECTROSTATICS
	//check for the existence of default charges for all types that take part in the reactions
	for (int i=0; i<m_current_reaction_system.nr_different_types; i++) {
		if(m_current_reaction_system.charges_of_types[m_current_reaction_system.type_index[i]]==m_invalid_charge) {
		    std::string message = std::string("Forgot to assign charge to type") +std::to_string(m_current_reaction_system.charges_of_types[m_current_reaction_system.type_index[i]]);
		    throw std::runtime_error(message);
			check_is_successfull=ES_ERROR;
		}
	}
	#endif
	return check_is_successfull;
}

/**
* Frees the data structures that were allocated for the reaction ensemble
*/
int ReactionEnsemble::free_reaction_ensemble(){
	/**needs to be called at the end of the simulation*/
	for(int single_reaction_i=0;single_reaction_i<m_current_reaction_system.nr_single_reactions;single_reaction_i++){
		//free reactant types and coefficients
		free(m_current_reaction_system.reactions[single_reaction_i]->reactant_types);
		free(m_current_reaction_system.reactions[single_reaction_i]->reactant_coefficients);	
		//free product types and coefficients
		free(m_current_reaction_system.reactions[single_reaction_i]->product_types);
		free(m_current_reaction_system.reactions[single_reaction_i]->product_coefficients);
		free(m_current_reaction_system.reactions[single_reaction_i]);	
	}
	free(m_current_reaction_system.reactions);
	free(m_current_reaction_system.type_index);
	free(m_current_reaction_system.charges_of_types);

	return 0;
}


//boring helper functions
/**
* Automatically sets the volume which is used by the reaction ensemble to the volume of a cuboid box
*/
void ReactionEnsemble::set_cuboid_reaction_ensemble_volume(){
	if (m_current_reaction_system.volume <0)
		m_current_reaction_system.volume=box_l[0]*box_l[1]*box_l[2];
}

/**
* Calculates the factorial expression which occurs in the reaction ensemble acceptance probability
*/
float ReactionEnsemble::factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i) {
	float value=1;
	if (nu_i == 0) {
		value=1.0;
	} else {
		value=1.0;
		if (nu_i > 0) {
			for(int i=1;i <= nu_i;i++) {
				value=value*1.0/(Ni0+i);
			}
		} else {
			int abs_nu_i =-1.0*nu_i;
			for(int i=0;i < abs_nu_i;i++) {
				value= value*(Ni0-i);
			}
		}
	}
	return value;
}

/**
* Checks wether all particles exist for the provided reaction.
*/
bool ReactionEnsemble::all_reactant_particles_exist(int reaction_id) {
	bool enough_particles=true;
	for(int i=0;i<m_current_reaction_system.reactions[reaction_id]->len_reactant_types;i++){
		int current_number;
		number_of_particles_with_type(m_current_reaction_system.reactions[reaction_id]->reactant_types[i], &current_number);
		if(current_number<m_current_reaction_system.reactions[reaction_id]->reactant_coefficients[i]){
			enough_particles=false;
			break;
		}
	}
	return enough_particles;
}

/**
* Calculates the current potential energy of the system and returns it. The function arguments are only dummy arguments needed to provide the same function signature as the function calculate_degree_of_association()
*/
double calculate_current_potential_energy_of_system_wrap(int unimportant_int, void* unimportant_wang_landau_system){
	return calculate_current_potential_energy_of_system();
}

/**
* Stores the particle property of a random particle of the provided type into the provided vector
*/
void ReactionEnsemble::append_particle_property_of_random_particle(int type, std::vector<stored_particle_property>& list_of_particles){
	int p_id ;
	find_particle_type(type, &p_id);
	stored_particle_property property_of_part={p_id,\
						   m_current_reaction_system.charges_of_types[find_index_of_type(type)],\
						   type
						  };
	list_of_particles.push_back(property_of_part);
}

/**
*Performs a trial reaction move
*/
void ReactionEnsemble::make_reaction_attempt(single_reaction* current_reaction, std::vector<stored_particle_property>& changed_particles_properties, std::vector<int>& p_ids_created_particles, std::vector<stored_particle_property>& hidden_particles_properties){
	const int number_of_saved_properties=3;//save p_id, charge and type of the reactant particle, only thing we need to hide the particle and recover it
	//create or hide particles of types with corresponding types in reaction
	for(int i=0;i<std::min(current_reaction->len_product_types,current_reaction->len_reactant_types);i++){
		//change std::min(reactant_coefficients(i),product_coefficients(i)) many particles of reactant_types(i) to product_types(i)
		for(int j=0;j<std::min(current_reaction->product_coefficients[i],current_reaction->reactant_coefficients[i]);j++){
			append_particle_property_of_random_particle(current_reaction->reactant_types[i], changed_particles_properties);
			replace(changed_particles_properties.back().p_id,current_reaction->product_types[i]);
		}
		//create product_coefficients(i)-reactant_coefficients(i) many product particles iff product_coefficients(i)-reactant_coefficients(i)>0,
		//iff product_coefficients(i)-reactant_coefficients(i)<0, hide this number of reactant particles
		if ( current_reaction->product_coefficients[i]-current_reaction->reactant_coefficients[i] >0) {
			for(int j=0; j< current_reaction->product_coefficients[i]-current_reaction->reactant_coefficients[i] ;j++) {
				int p_id=create_particle(current_reaction->product_types[i]);
				p_ids_created_particles.push_back(p_id);
			}
		} else if (current_reaction->reactant_coefficients[i]-current_reaction->product_coefficients[i] >0) {
			for(int j=0;j<current_reaction->reactant_coefficients[i]-current_reaction->product_coefficients[i];j++) {
				append_particle_property_of_random_particle(current_reaction->reactant_types[i], hidden_particles_properties);			
				hide_particle(hidden_particles_properties.back().p_id,current_reaction->reactant_types[i]);
			}
		}

	}
	//create or hide particles of types with noncorresponding replacement types
	for(int i=std::min(current_reaction->len_product_types,current_reaction->len_reactant_types);i< std::max(current_reaction->len_product_types,current_reaction->len_reactant_types);i++ ) {
		if(current_reaction->len_product_types<current_reaction->len_reactant_types){
			//hide superfluous reactant_types particles
			for(int j=0;j<current_reaction->reactant_coefficients[i];j++){
				append_particle_property_of_random_particle(current_reaction->reactant_types[i], hidden_particles_properties);
				hide_particle(hidden_particles_properties.back().p_id,current_reaction->reactant_types[i]);
			}
		} else {
			//create additional product_types particles
			for(int j=0;j<current_reaction->product_coefficients[i];j++){
				int p_id= create_particle(current_reaction->product_types[i]);
				p_ids_created_particles.push_back(p_id);
			}
		}
	}	
	
}

/**
* Calculates the whole product of factorial expressions which occur in the reaction ensemble acceptance probability
*/
double ReactionEnsemble::calculate_factorial_expression(single_reaction* current_reaction, int* old_particle_numbers){
	double factorial_expr=1.0;	
	//factorial contribution of reactants
	for(int i=0;i<current_reaction->len_reactant_types;i++) {
		int nu_i=-1*current_reaction->reactant_coefficients[i];
		int N_i0= old_particle_numbers[find_index_of_type(current_reaction->reactant_types[i])];
		factorial_expr=factorial_expr*factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N_i0,nu_i); //zeta = 1 (see smith paper) since we only perform one reaction at one call of the function
	}
	//factorial contribution of products
	for(int i=0;i<current_reaction->len_product_types;i++) {
		int nu_i=current_reaction->product_coefficients[i];
		int N_i0= old_particle_numbers[find_index_of_type(current_reaction->product_types[i])];
		factorial_expr=factorial_expr*factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N_i0,nu_i); //zeta = 1 (see smith paper) since we only perform one reaction at one call of the function
	}
	return factorial_expr;
}

/**
* Restores the previosly stored particle properties. This funtion is invoked when a reaction attempt is rejected.
*/
void ReactionEnsemble::restore_properties(std::vector<stored_particle_property> property_list ,const int number_of_saved_properties){
	//this function restores all properties of all particles provided in the property list, the format of the property list is (p_id,charge,type) repeated for each particle that occurs in that list
	for(int i=0;i<property_list.size();i++) {
		double charge= property_list[i].charge;
		int type=(int) property_list[i].type;
		#ifdef ELECTROSTATICS
		//set charge
		set_particle_q(property_list[i].p_id, charge);
		#endif
		//set type
		set_particle_type(property_list[i].p_id, type);
	}
}

/**
* Calculates the expression in the acceptance probability in the reaction ensemble
*/
double ReactionEnsemble::calculate_boltzmann_factor_reaction_ensemble(single_reaction* current_reaction, double E_pot_old, double E_pot_new, std::vector<int>& old_particle_numbers){
	/**calculate the acceptance probability in the reaction ensemble */
	const double volume = m_current_reaction_system.volume;
	const double factorial_expr=calculate_factorial_expression(current_reaction, old_particle_numbers.data());

	const double beta =1.0/m_current_reaction_system.temperature_reaction_ensemble;
	const double standard_pressure_in_simulation_units=m_current_reaction_system.standard_pressure_in_simulation_units;
	//calculate boltzmann factor
	return std::pow(volume*beta*standard_pressure_in_simulation_units, current_reaction->nu_bar) * current_reaction->equilibrium_constant * factorial_expr * exp(-beta * (E_pot_new - E_pot_old));
}

/**
*generic one way reaction
*A+B+...+G +... --> K+...X + Z +...
*you need to use 2A --> B instead of A+A --> B since in the last case you assume distinctness of the particles, however both ways to describe the reaction are equivalent in the thermodynamic limit
*further it is crucial for the function in which order you provide the reactant and product types since particles will be replaced correspondingly! If there are less reactants than products, new product particles are created randomly in the box. Matching particles simply change the types. If there are more reactants than products, old reactant particles are deleted.
 */
int ReactionEnsemble::generic_oneway_reaction(int reaction_id, int reaction_modus){

	single_reaction* current_reaction=m_current_reaction_system.reactions[reaction_id];
	//Wang-Landau begin
	int old_state_index;
	if(reaction_modus==reaction_ensemble_wang_landau_mode){
		old_state_index=get_flattened_index_wang_landau_of_current_state();
		if(old_state_index>=0){
			if(m_current_wang_landau_system.histogram[old_state_index]>=0)
				m_current_wang_landau_system.monte_carlo_trial_moves+=1;
		}		
	}
	//Wang-Landau end
	
	if (all_reactant_particles_exist(reaction_id) ==false ) {
		//makes sure, no incomplete reaction is performed -> only need to consider rollback of complete reactions
		
		//Wang-Landau begin
		//increase the wang landau potential and histogram at the current nbar (this case covers the cases nbar=0 or nbar=1)
		if(reaction_modus==reaction_ensemble_wang_landau_mode)
			update_wang_landau_potential_and_histogram(old_state_index);
		//Wang-Landau end
		
		return 0;
	}
	
	//calculate potential energy
	const double E_pot_old=calculate_current_potential_energy_of_system_wrap(0, NULL); //only consider potential energy since we assume that the kinetic part drops out in the process of calculating ensemble averages (kinetic part may be seperated and crossed out)
	
	//find reacting molecules in reactants and save their properties for later recreation if step is not accepted
	//do reaction
	//save old particle_numbers
	std::vector<int> old_particle_numbers(m_current_reaction_system.nr_different_types);
	if(reaction_modus!=constant_pH_mode){
		for(int type_index=0;type_index<m_current_reaction_system.nr_different_types;type_index++)
			number_of_particles_with_type(m_current_reaction_system.type_index[type_index], &(old_particle_numbers[type_index])); // here could be optimized by not going over all types but only the types that occur in the reaction
	}

	std::vector<int> p_ids_created_particles;
	std::vector<stored_particle_property> hidden_particles_properties;
	std::vector<stored_particle_property> changed_particles_properties;
	const int number_of_saved_properties=3; //save p_id, charge and type of the reactant particle, only thing we need to hide the particle and recover it
	make_reaction_attempt(current_reaction, changed_particles_properties, p_ids_created_particles, hidden_particles_properties);
	
	const double E_pot_new=calculate_current_potential_energy_of_system_wrap(0, NULL);


	//Wang-Landau begin
	//save new_state_index
	int new_state_index;
	if(reaction_modus==reaction_ensemble_wang_landau_mode)
		new_state_index=get_flattened_index_wang_landau_of_current_state();
	double bf;
	if(reaction_modus==reaction_ensemble_mode){
		bf=calculate_boltzmann_factor_reaction_ensemble(current_reaction, E_pot_old, E_pot_new, old_particle_numbers);
	}else if (reaction_modus==reaction_ensemble_wang_landau_mode)
		bf=calculate_boltzmann_factor_reaction_ensemble_wang_landau(current_reaction, E_pot_old, E_pot_new, old_particle_numbers, old_state_index, new_state_index, false);
	else if (reaction_modus==constant_pH_mode)
		bf=calculate_boltzmann_factor_consant_pH(current_reaction, E_pot_old, E_pot_new);
	else
		throw std::runtime_error("Reaction mode is unknown");
	int reaction_is_accepted=false;
	//Wang-Landau begin
	int accepted_state;
	//Wang-Landau end
	if ( d_random() < bf ) {
		//accept
		if(reaction_modus==reaction_ensemble_wang_landau_mode)
			accepted_state=new_state_index;
		
		//delete hidden reactant_particles (remark: dont delete changed particles)
		//extract ids of to be deleted particles and sort them. needed since delete_particle changes particle p_ids. start deletion from the largest p_id onwards
		int len_hidden_particles_properties=hidden_particles_properties.size();
		int to_be_deleted_hidden_ids[len_hidden_particles_properties];
		for(int i=0;i<len_hidden_particles_properties;i++) {
			int p_id = (int) hidden_particles_properties[i].p_id;
			to_be_deleted_hidden_ids[i]=p_id;
		}
		std::sort(to_be_deleted_hidden_ids,to_be_deleted_hidden_ids+len_hidden_particles_properties,std::greater<int>());
		
		for(int i=0;i<len_hidden_particles_properties;i++)
			delete_particle(to_be_deleted_hidden_ids[i]); //delete particle
		reaction_is_accepted= true;
	} else {
		//reject
		if(reaction_modus==reaction_ensemble_wang_landau_mode)
			accepted_state=old_state_index;
		//reverse reaction
		//1) delete created product particles
		std::sort(p_ids_created_particles.begin(),p_ids_created_particles.end(),std::greater<int>());// needed since delete_particle changes particle p_ids. start deletion from the largest p_id onwards
		for(int i=0;i<p_ids_created_particles.size();i++){
			delete_particle(p_ids_created_particles[i]);
		}
		//2)restore previously hidden reactant particles
		restore_properties(hidden_particles_properties, number_of_saved_properties);
		//2)restore previously changed reactant particles
		restore_properties(changed_particles_properties, number_of_saved_properties);
		reaction_is_accepted= false;
	}
	if(reaction_modus==reaction_ensemble_wang_landau_mode)
		update_wang_landau_potential_and_histogram(accepted_state);

	return reaction_is_accepted;
}

/**
* Calculates the change in particle numbers for the given reaction
*/
int ReactionEnsemble::calculate_nu_bar(int* reactant_coefficients, int len_reactant_types,  int* product_coefficients, int len_product_types){
	//should only be used at when defining a new reaction
	int nu_bar =0;
	for(int i=0;i<len_reactant_types;i++){
		nu_bar-=reactant_coefficients[i];
	}
	for(int i=0;i<len_product_types;i++){
		nu_bar+=product_coefficients[i];
	}
	return nu_bar;
}

/**
* Adds types to an index. the index is later used inside of the Reaction ensemble algorithm to determine e.g. the charge which is associated to a type.
*/
int ReactionEnsemble::add_types_to_index(int* type_list, int len_type_list, int status_gc_init){
	int status_gc_init_temp=0;
	for (int i =0; i<len_type_list;i++){
		bool type_i_is_known=is_in_list(type_list[i],m_current_reaction_system.type_index,m_current_reaction_system.nr_different_types);
		if (type_i_is_known==false){
			m_current_reaction_system.type_index=(int*) realloc(m_current_reaction_system.type_index, sizeof(int)*(m_current_reaction_system.nr_different_types+1));
			m_current_reaction_system.type_index[m_current_reaction_system.nr_different_types]=type_list[i];
			m_current_reaction_system.nr_different_types+=1;
			status_gc_init_temp=init_type_array(type_list[i]); //make types known in espresso
			status_gc_init=status_gc_init || status_gc_init_temp;
		}
	}
	return status_gc_init;
}

/**
* Adds types to an index and copes with the case that no index was present before. the index is later used inside of the Reaction ensemble algorithm to determine e.g. the charge which is associated to a type.
*/
int ReactionEnsemble::update_type_index(int* reactant_types, int len_reactant_types, int* product_types, int len_product_types){
	//should only be used when defining a new reaction
	int status_gc_init=0;
	if(m_current_reaction_system.type_index==NULL){
		m_current_reaction_system.type_index=(int*) calloc(1,sizeof(int));
		if(len_reactant_types>0)
			m_current_reaction_system.type_index[0]=reactant_types[0];
		else
			m_current_reaction_system.type_index[0]=product_types[0];
		m_current_reaction_system.nr_different_types=1;
		status_gc_init=init_type_array(m_current_reaction_system.type_index[0]); //make types known in espresso
	}
	status_gc_init=add_types_to_index(reactant_types, len_reactant_types,status_gc_init);
	status_gc_init=add_types_to_index(product_types, len_product_types,status_gc_init);
	
	//increase m_current_reaction_system.charges_of_types length
	m_current_reaction_system.charges_of_types =(double*) realloc(m_current_reaction_system.charges_of_types,sizeof(double)*m_current_reaction_system.nr_different_types);
	m_current_reaction_system.charges_of_types[m_current_reaction_system.nr_different_types-1]=m_invalid_charge;
	return status_gc_init;
}

/**
* Finds the index of a given type.
*/
int ReactionEnsemble::find_index_of_type(int type){
	int index =-100; //initialize to invalid index
	for(int i=0; i<m_current_reaction_system.nr_different_types;i++){
		if(type==m_current_reaction_system.type_index[i]){
			index=i;
			break;
		}
	}
	if(index<0){
		throw std::runtime_error("Invalid Index");
		
	}
	return index;
}

/**
* Replaces a particle with the given particle id to be of a certain type. This especially means that the particle type and the particle charge are changed.
*/
int ReactionEnsemble::replace(int p_id, int desired_type){
	int err_code_type=set_particle_type(p_id, desired_type);
	int err_code_q;
	#ifdef ELECTROSTATICS
	err_code_q=set_particle_q(p_id, (double) m_current_reaction_system.charges_of_types[find_index_of_type(desired_type)]);
	#endif
	return (err_code_q bitor err_code_type);
}

/**
* Hides a particle from short ranged interactions and from the electrostatic interaction. Additional hiding from interactions would need to be implemented here.
*/
int ReactionEnsemble::hide_particle(int p_id, int previous_type){
	/**
    *remove_charge and put type to a non existing one --> no interactions anymore it is as if the particle was non existing (currently only type-based interactions are swithced off, as well as the electrostatic interaction)  
    *hide_particle() does not break bonds for simple reactions. as long as there are no reactions like 2A -->B where one of the reacting A particles occurs in the polymer (think of bond breakages if the monomer in the polymer gets deleted in the reaction). This constraint is not of fundamental reason, but there would be a need for a rule for such "collision" reactions (a reaction like the one above).
    */
	#ifdef ELECTROSTATICS
	//set charge
	set_particle_q(p_id, 0.0);
	#endif
	//set type
	int err_code_type=set_particle_type(p_id, m_current_reaction_system.non_interacting_type);
	return err_code_type;
}

/**
* Deletes the particle with the given p_id. This method is intended to only delete unbonded particles since bonds are coupled to ids. The methods is aware of holes in the particle id range and fills them if a particle gets created.
*/

int ReactionEnsemble::delete_particle(int p_id) {
	/**deletes the particle with the provided id and stores if it created a hole at that position in the particle id range */
	if (p_id == max_seen_particle) {
		// last particle, just delete
		remove_particle(p_id);
		// remove all saved empty p_ids which are greater than the max_seen_particle this is needed in order to avoid the creation of holes
        for (auto p_id_iter = m_empty_p_ids_smaller_than_max_seen_particle.begin(); p_id_iter != m_empty_p_ids_smaller_than_max_seen_particle.end(); ){
            if( (*p_id_iter) >=max_seen_particle)
                p_id_iter=m_empty_p_ids_smaller_than_max_seen_particle.erase(p_id_iter); //update iterator after container was modified
            else
                ++p_id_iter;
        }
	} else if (p_id <= max_seen_particle){
        remove_particle(p_id);
        m_empty_p_ids_smaller_than_max_seen_particle.push_back(p_id);
	} else {
	    throw std::runtime_error("Particle id is greater than the max seen particle id");
	}
	return 0;
}

/**
* Writes a random position inside the central box into the provided array.
*/
void ReactionEnsemble::get_random_position_in_box (double* out_pos) {
	if(m_current_reaction_system.box_is_cylindric_around_z_axis==true) {
		//see http://mathworld.wolfram.com/DiskPointPicking.html
		double random_radius=m_current_reaction_system.cyl_radius*std::sqrt(d_random()); //for uniform disk point picking in cylinder
		double phi=2.0*PI*d_random();
		out_pos[0]=random_radius*cos(phi);
		out_pos[1]=random_radius*sin(phi);
		while (std::pow(out_pos[0],2)+std::pow(out_pos[1],2)<=std::pow(m_current_reaction_system.exclusion_radius,2)){
			random_radius=m_current_reaction_system.cyl_radius*std::sqrt(d_random());
			out_pos[0]=random_radius*cos(phi);
			out_pos[1]=random_radius*sin(phi);		
		}
		out_pos[0]+=m_current_reaction_system.cyl_x;
		out_pos[1]+=m_current_reaction_system.cyl_y;
		out_pos[2]=box_l[2]*d_random();
	} else if (m_current_reaction_system.box_has_wall_constraints==true) {
		out_pos[0]=box_l[0]*d_random();
		out_pos[1]=box_l[1]*d_random();	
		out_pos[2]=m_current_reaction_system.slab_start_z+(m_current_reaction_system.slab_end_z-m_current_reaction_system.slab_start_z)*d_random();
	}else{
		//cubic case
		out_pos[0]=box_l[0]*d_random();
		out_pos[1]=box_l[1]*d_random();
		out_pos[2]=box_l[2]*d_random();
	
	}
} 

/**
* Writes a random position inside the central box into the provided array. Additionally it proposes points with a small radii more often than a uniform random probability density would do it.
*/
void ReactionEnsemble::get_random_position_in_box_enhanced_proposal_of_small_radii (double* out_pos) {
	double random_radius=m_current_reaction_system.cyl_radius*d_random(); //for enhanced proposal of small radii, needs correction within metropolis hasting algorithm, proposal density is p(x,y)=1/(2*pi*cyl_radius*r(x,y)), that means small radii are proposed more often
	double phi=2.0*PI*d_random();
	out_pos[0]=random_radius*cos(phi);
	out_pos[1]=random_radius*sin(phi);
	while (std::pow(out_pos[0],2)+std::pow(out_pos[1],2)<=std::pow(m_current_reaction_system.exclusion_radius,2) or std::pow(out_pos[0],2)+std::pow(out_pos[1],2) > std::pow(m_current_reaction_system.cyl_radius,2)){
		random_radius=m_current_reaction_system.cyl_radius*d_random();
		out_pos[0]=random_radius*cos(phi);
		out_pos[1]=random_radius*sin(phi);		
	}
	out_pos[0]+=m_current_reaction_system.cyl_x;
	out_pos[1]+=m_current_reaction_system.cyl_y;
	out_pos[2]=box_l[2]*d_random();
}

/**
* Creates a particle at the end of the observed particle id range.
*/
int ReactionEnsemble::create_particle(int desired_type){
    int p_id;
    if(m_empty_p_ids_smaller_than_max_seen_particle.size()>0){
        auto p_id_iter= std::min_element(std::begin(m_empty_p_ids_smaller_than_max_seen_particle), std::end(m_empty_p_ids_smaller_than_max_seen_particle));
        p_id=*p_id_iter;
        m_empty_p_ids_smaller_than_max_seen_particle.erase(p_id_iter);
    }else{
        p_id=max_seen_particle+1;
    }
	double pos_vec[3];
	
	//create random velocity vector according to Maxwell Boltzmann distribution for components
	double vel[3];
	//we usse mass=1 for all particles, think about adapting this
	vel[0]=std::pow(2*PI*m_current_reaction_system.temperature_reaction_ensemble,-3.0/2.0)*gaussian_random()*time_step;//scale for internal use in espresso
	vel[1]=std::pow(2*PI*m_current_reaction_system.temperature_reaction_ensemble,-3.0/2.0)*gaussian_random()*time_step;//scale for internal use in espresso
	vel[2]=std::pow(2*PI*m_current_reaction_system.temperature_reaction_ensemble,-3.0/2.0)*gaussian_random()*time_step;//scale for internal use in espresso
	double charge= (double) m_current_reaction_system.charges_of_types[find_index_of_type(desired_type)];
	bool particle_inserted_too_close_to_another_one=true;
	int max_insert_tries=1000;
	int insert_tries=0;
	double min_dist=m_current_reaction_system.exclusion_radius; //setting of a minimal distance is allowed to avoid overlapping configurations if there is a repulsive potential. States with very high energies have a probability of almost zero and therefore do not contribute to ensemble averages.
	if(min_dist!=0){
		while(particle_inserted_too_close_to_another_one && insert_tries<max_insert_tries) {
			get_random_position_in_box(pos_vec);
			place_particle(p_id,pos_vec);
			//set type
			set_particle_type(p_id, desired_type);
			#ifdef ELECTROSTATICS
			//set charge
			set_particle_q(p_id, charge);
			#endif
			//set velocities
			set_particle_v(p_id,vel);
			double d_min=distto(partCfg(), pos_vec,p_id); //TODO also catch constraints with an IFDEF CONSTRAINTS here, but only interesting, when doing MD/ HMC because then the system might explode easily here due to high forces
			insert_tries+=1;
			if(d_min>m_current_reaction_system.exclusion_radius)
				particle_inserted_too_close_to_another_one=false;
		}
	}else{
		get_random_position_in_box(pos_vec);
		place_particle(p_id,pos_vec);
		//set type
		set_particle_type(p_id, desired_type);	
		//set velocities
		set_particle_v(p_id,vel);
		#ifdef ELECTROSTATICS
		//set charge
		set_particle_q(p_id, charge);
		#endif
	}
	if(insert_tries>max_insert_tries){
		throw std::runtime_error("No particle inserted");
		return -1;	
	}
	return p_id;
}

/**
* Checks wether an integer is in an array of integers.
*/
bool ReactionEnsemble::is_in_list(int value, int* list, int len_list){
	for(int i=0;i<len_list;i++){
		if(list[i]==value)
			return true;
	}
	return false;
}

//the following 2 functions are directly taken from ABHmath.tcl
/**
* Calculates the normed vector of a given vector
*/
std::vector<double> vecnorm(std::vector<double> vec, double desired_length){
	for(int i=0;i<vec.size();i++){
		vec[i]=vec[i]/utils::veclen(vec)*desired_length;	
	}
	return vec;
}

/**
* Calculates a uniformly distributed vector on a sphere of given radius.
*/
std::vector<double> vec_random(double desired_length){
	/**returns a random vector of length len
	*(uniform distribution on a sphere)
	*This is done by chosing 3 uniformly distributed random numbers [-1,1]
	*If the length of the resulting vector is <= 1.0 the vector is taken and normalized
	*to the desired length, otherwise the procedure is repeated until succes.
	*On average the procedure needs 5.739 random numbers per vector.
	*(This is probably not the most efficient way, but it works!)
	*/
	std::vector<double> vec;
	while(1){
		for(int i=0;i<3;i++){
			vec.push_back(2*d_random()-1.0);
		}
		if (utils::veclen(vec)<=1)
			break;
	}
	vecnorm(vec,desired_length);
	return vec;
}

/**
* Adds a random vector of given length to the provided array named vector.
*/
void ReactionEnsemble::add_random_vector(double* vector, int len_vector, double length_of_displacement){
	//adds a vector which is uniformly distributed on a sphere
	std::vector<double> random_direction_vector = vec_random(length_of_displacement);
	for(int i=0;i<len_vector;i++){
		vector[i]+=random_direction_vector[i];
	}
}

/**
* Performs a global mc move for a particle of the provided type.
*/
bool ReactionEnsemble::do_global_mc_move_for_particles_of_type(int type, int start_id_polymer, int end_id_polymer, int particle_number_of_type, bool use_wang_landau){
	
	m_tried_configurational_MC_moves+=1;
	bool got_accepted=false;

	int old_state_index;
	if(use_wang_landau==true){
		old_state_index=get_flattened_index_wang_landau_of_current_state();
		if(old_state_index>=0){
			if(m_current_wang_landau_system.histogram[old_state_index]>=0)
				m_current_wang_landau_system.monte_carlo_trial_moves+=1;
		}
	}

	int p_id;
	number_of_particles_with_type(type, &(particle_number_of_type));
	if(particle_number_of_type==0){
		//reject
		if(m_current_wang_landau_system.do_energy_reweighting==true&& use_wang_landau== true){
			update_wang_landau_potential_and_histogram(old_state_index);
		}
		return got_accepted;
	}


	const double E_pot_old=calculate_current_potential_energy_of_system_wrap(0, NULL);

	double particle_positions[3*particle_number_of_type];
	int changed_particle_counter=0;
	int p_id_s_changed_particles[particle_number_of_type];

	//save old_position
	double temp_pos[3];
	get_random_position_in_box(temp_pos);
	while(changed_particle_counter<particle_number_of_type){
		if(changed_particle_counter==0){
			find_particle_type(type, &p_id);
		}else{
			//determine a p_id you have not touched yet
			while(is_in_list(p_id,p_id_s_changed_particles,changed_particle_counter) or changed_particle_counter==0){
				find_particle_type(type, &p_id); //check wether you already touched this p_id
			}
		}

		auto part = get_particle_data(p_id);
		double ppos[3];
		memmove(ppos, part->r.p, 3*sizeof(double));

		particle_positions[3*changed_particle_counter]=ppos[0];
		particle_positions[3*changed_particle_counter+1]=ppos[1];
		particle_positions[3*changed_particle_counter+2]=ppos[2];
		place_particle(p_id,temp_pos);
		p_id_s_changed_particles[changed_particle_counter]=p_id;
		changed_particle_counter+=1;
	}
	
	//propose new positions
	changed_particle_counter=0;
	int max_tries=100*particle_number_of_type;//important for very dense systems
	int attempts=0;
	double new_pos[3];
	while(changed_particle_counter<particle_number_of_type){
		p_id=p_id_s_changed_particles[changed_particle_counter];
		bool particle_inserted_too_close_to_another_one=true;
		while(particle_inserted_too_close_to_another_one==true&& attempts<max_tries){
			//change particle position
			get_random_position_in_box(new_pos);
//			get_random_position_in_box_enhanced_proposal_of_small_radii(new_pos); //enhanced proposal of small radii
			place_particle(p_id,new_pos);
			double d_min=distto(partCfg(), new_pos,p_id);
			if(d_min>m_current_reaction_system.exclusion_radius){
				particle_inserted_too_close_to_another_one=false;
			}
			attempts+=1;
		}
		changed_particle_counter+=1;
	}
		if(attempts==max_tries){
		//reversing
		//create particles again at the positions they were
		for(int i=0;i<particle_number_of_type;i++)
			place_particle(p_id_s_changed_particles[i],&particle_positions[3*i]);
	}
	
	
	//change polymer conformation if start and end id are provided
	double old_pos_polymer_particle[3*(end_id_polymer-start_id_polymer+1)];
	if(start_id_polymer>=0 && end_id_polymer >=0 ){
		
		for(int i=start_id_polymer;i<=end_id_polymer;i++){
			auto part = get_particle_data(i);
			//move particle to new position nearby
			const double length_of_displacement=0.05;
			add_random_vector(part->r.p, 3, length_of_displacement);
			place_particle(i,part->r.p);
		}
		
	}
	
	const double E_pot_new=calculate_current_potential_energy_of_system_wrap(0, NULL);
	double beta =1.0/m_current_reaction_system.temperature_reaction_ensemble;
	
	int new_state_index;
	double bf=1.0;
	if(use_wang_landau==true){
		new_state_index=get_flattened_index_wang_landau_of_current_state();
		std::vector<int> dummy_old_particle_numbers;
		bf=calculate_boltzmann_factor_reaction_ensemble_wang_landau(NULL, E_pot_old, E_pot_new, dummy_old_particle_numbers, old_state_index, new_state_index, true);
	}else{
		bf=std::min(1.0, bf*exp(-beta*(E_pot_new-E_pot_old))); //Metropolis Algorithm since proposal density is symmetric
	}
	
	
//	//correct for enhanced proposal of small radii by using the metropolis hastings algorithm for asymmetric proposal densities.
//	double old_radius=std::sqrt(std::pow(particle_positions[0]-m_current_reaction_system.cyl_x,2)+std::pow(particle_positions[1]-m_current_reaction_system.cyl_y,2));
//	double new_radius=std::sqrt(std::pow(new_pos[0]-m_current_reaction_system.cyl_x,2)+std::pow(new_pos[1]-m_current_reaction_system.cyl_y,2));
//	bf=std::min(1.0, bf*exp(-beta*(E_pot_new-E_pot_old))*new_radius/old_radius); //Metropolis-Hastings Algorithm for asymmetric proposal density

	if(d_random()<bf){
		//accept
		m_accepted_configurational_MC_moves+=1;
		got_accepted=true;
		if(use_wang_landau==true && m_current_wang_landau_system.do_energy_reweighting==true){
			//modify wang_landau histogram and potential
			update_wang_landau_potential_and_histogram(new_state_index);		
		}
	}else{
		//reject
		//modify wang_landau histogram and potential
		if(use_wang_landau==true && m_current_wang_landau_system.do_energy_reweighting==true)
			update_wang_landau_potential_and_histogram(old_state_index);

		//create particles again at the positions they were
		for(int i=0;i<particle_number_of_type;i++)
			place_particle(p_id_s_changed_particles[i],&particle_positions[3*i]);
		//restore polymer particle again at original position
		if(start_id_polymer>=0 && end_id_polymer >=0 ){
			//place_particle(random_polymer_particle_id, old_pos_polymer_particle);
			for(int i=start_id_polymer;i<=end_id_polymer;i++)
				place_particle(i,&old_pos_polymer_particle[3*i]);
		}
	}
	return got_accepted;
}

///////////////////////////////////////////// Wang-Landau algorithm

/**
* Adds a new collective variable (CV) of the type degree of association to the Wang-Landau sampling
*/
void ReactionEnsemble::add_new_CV_degree_of_association(int associated_type, double CV_minimum, double CV_maximum, std::vector<int> _corresponding_acid_types){
            collective_variable* new_collective_variable=(collective_variable*) calloc(1,sizeof(collective_variable));
            new_collective_variable->associated_type=associated_type;
            new_collective_variable->CV_minimum=CV_minimum;
            new_collective_variable->CV_maximum=CV_maximum;

            int* corresponding_acid_types = (int *) malloc(_corresponding_acid_types.size() * sizeof(int));
            for(int i=0; i<_corresponding_acid_types.size();i++)
                corresponding_acid_types[i]=_corresponding_acid_types[i];
            new_collective_variable->corresponding_acid_types=corresponding_acid_types;
            new_collective_variable->nr_corresponding_acid_types=_corresponding_acid_types.size();


            m_current_wang_landau_system.collective_variables=(collective_variable**) realloc(m_current_wang_landau_system.collective_variables,sizeof(collective_variable*)*(m_current_wang_landau_system.nr_collective_variables+1));
            m_current_wang_landau_system.collective_variables[m_current_wang_landau_system.nr_collective_variables]=new_collective_variable;
            m_current_wang_landau_system.nr_collective_variables+=1;
            initialize_wang_landau();
}

/**
* Adds a new collective variable (CV) of the type potential energy to the Wang-Landau sampling
*/
void ReactionEnsemble::add_new_CV_potential_energy(std::string filename, double delta_CV){
            collective_variable* new_collective_variable=(collective_variable*) calloc(1,sizeof(collective_variable));
            new_collective_variable->energy_boundaries_filename=strdup(filename.c_str());
            new_collective_variable->delta_CV=delta_CV;
            m_current_wang_landau_system.collective_variables=(collective_variable**) realloc(m_current_wang_landau_system.collective_variables,sizeof(collective_variable*)*(m_current_wang_landau_system.nr_collective_variables+1));
            m_current_wang_landau_system.collective_variables[m_current_wang_landau_system.nr_collective_variables]=new_collective_variable;
            m_current_wang_landau_system.nr_collective_variables+=1;
            initialize_wang_landau();
}

/**
* Returns the flattened index of the multidimensional Wang-Landau histogram
*/
int ReactionEnsemble::get_flattened_index_wang_landau(double* current_state, double* collective_variables_minimum_values, double* collective_variables_maximum_values, double* delta_collective_variables_values, int nr_collective_variables){
	int index=-10; //negative number is not allowed as index and therefore indicates error
	int individual_indices[nr_collective_variables]; //pre result
	memset(individual_indices, -1, sizeof(individual_indices)); //initialize individual_indices to -1
	int* nr_subindices_of_collective_variable =m_current_wang_landau_system.nr_subindices_of_collective_variable;

	//check for the current state to be an allowed state in the [range collective_variables_minimum_values:collective_variables_maximum_values], else return a negative index
	for(int collective_variable_i=0;collective_variable_i<nr_collective_variables;collective_variable_i++){
		if(current_state[collective_variable_i]>collective_variables_maximum_values[collective_variable_i]+delta_collective_variables_values[collective_variable_i]*0.98 || current_state[collective_variable_i]<collective_variables_minimum_values[collective_variable_i]-delta_collective_variables_values[collective_variable_i]*0.01)
			return -10;
	}

	for(int collective_variable_i=0;collective_variable_i<nr_collective_variables;collective_variable_i++){
		if(collective_variable_i==m_current_wang_landau_system.nr_collective_variables-1 && m_current_wang_landau_system.do_energy_reweighting==true)	//for energy collective variable (simple truncating conversion desired)
			individual_indices[collective_variable_i]=(int) ((current_state[collective_variable_i]-collective_variables_minimum_values[collective_variable_i])/delta_collective_variables_values[collective_variable_i]);
		else	//for degree of association collective variables (rounding conversion desired)
			individual_indices[collective_variable_i]=(int) ((current_state[collective_variable_i]-collective_variables_minimum_values[collective_variable_i])/delta_collective_variables_values[collective_variable_i]+0.5);
		if(individual_indices[collective_variable_i]<0 or individual_indices[collective_variable_i]>=nr_subindices_of_collective_variable[collective_variable_i]) // sanity check
			return -10;
	}
	
	//get flattened index from individual_indices
	index=0; //this is already part of the algorithm to find the correct index
	for(int collective_variable_i=0;collective_variable_i<nr_collective_variables;collective_variable_i++){
		int factor=1;
		for(int j=collective_variable_i+1;j<nr_collective_variables;j++){
			factor*=nr_subindices_of_collective_variable[j];
		}
		index+=factor*individual_indices[collective_variable_i];
		
	}	
	return index;
}

/**
* Returns the flattened index of the multidimensional Wang-Landau histogram for the current state of the simulation.
*/
int ReactionEnsemble::get_flattened_index_wang_landau_of_current_state(){
	int nr_collective_variables=m_current_wang_landau_system.nr_collective_variables;
	//get current state
	double current_state[nr_collective_variables];
	for(int CV_i=0;CV_i<nr_collective_variables;CV_i++){
		current_state[CV_i]=(m_current_wang_landau_system.collective_variables[CV_i]->determine_current_state_in_collective_variable_with_index)(CV_i,&m_current_wang_landau_system);	
	}

	//get collective_variables_minimum_values
	double collective_variables_minimum_values[nr_collective_variables];
	for(int CV_i=0;CV_i<nr_collective_variables;CV_i++){
		collective_variables_minimum_values[CV_i]=m_current_wang_landau_system.collective_variables[CV_i]->CV_minimum;	
	}
	//get collective_variables_maximum_values
	double collective_variables_maximum_values[nr_collective_variables];
	for(int CV_i=0;CV_i<nr_collective_variables;CV_i++){
		collective_variables_maximum_values[CV_i]=m_current_wang_landau_system.collective_variables[CV_i]->CV_maximum;	
	}
	//get delta_collective_variables_values
	double delta_collective_variables_values[nr_collective_variables];
	for(int CV_i=0;CV_i<nr_collective_variables;CV_i++){
		delta_collective_variables_values[CV_i]=m_current_wang_landau_system.collective_variables[CV_i]->delta_CV;	
	}
	int index=get_flattened_index_wang_landau(current_state, collective_variables_minimum_values, collective_variables_maximum_values, delta_collective_variables_values, nr_collective_variables);
	return index;
}

/**
* Returns the minimum value of the collective variable on a delta_CV spaced grid which starts at 0
*/
double get_minimum_CV_value_on_delta_CV_spaced_grid(double min_CV_value, double delta_CV) {
	//assume grid has it s origin at 0
	double minimum_CV_value_on_delta_CV_spaced_grid=floor(min_CV_value/delta_CV)*delta_CV;
	return minimum_CV_value_on_delta_CV_spaced_grid;
};


/**
* Calculates the smallest difference in the degree of association which can be observed when changing the degree of association by one single reaction.
*/
double ReactionEnsemble::calculate_delta_degree_of_association(int index_of_current_collective_variable){
	//calculate Delta in the degree of association so that EVERY reaction step is driven.
	collective_variable* current_collective_variable=m_current_wang_landau_system.collective_variables[index_of_current_collective_variable];
	int total_number_of_corresponding_acid=0;
	for(int corresponding_type_i=0; corresponding_type_i<current_collective_variable->nr_corresponding_acid_types;corresponding_type_i++){
		int num_of_current_type;
		number_of_particles_with_type(current_collective_variable->corresponding_acid_types[corresponding_type_i],&num_of_current_type);
		total_number_of_corresponding_acid+=num_of_current_type;
	}
	double delta=1.0/total_number_of_corresponding_acid;
	//now modify the minimum value of the CV to lie on the grid
	current_collective_variable->CV_minimum=get_minimum_CV_value_on_delta_CV_spaced_grid(current_collective_variable->CV_minimum,delta);
	return delta ;
}

/**
* Initializes the Wang-Landau histogram.
*/
int* ReactionEnsemble::initialize_histogram(){
	int needed_bins=1;
	for(int CV_i=0;CV_i<m_current_wang_landau_system.nr_collective_variables;CV_i++){
		collective_variable* current_collective_variable=m_current_wang_landau_system.collective_variables[CV_i];
		needed_bins*=int((current_collective_variable->CV_maximum-current_collective_variable->CV_minimum)/current_collective_variable->delta_CV)+1; // plus 1 needed for degrees of association related part of histogram (think of only one acid particle)
	}
	int* histogram =(int*) calloc(1,sizeof(int)*needed_bins); //calloc initializes everything to zero
	m_current_wang_landau_system.len_histogram=needed_bins;
	return histogram;
}

/**
* Initializes the Wang-Landau potential.
*/
double* ReactionEnsemble::initialize_wang_landau_potential(){
	int needed_bins=1;
	for(int CV_i=0;CV_i<m_current_wang_landau_system.nr_collective_variables;CV_i++){
		collective_variable* current_collective_variable=m_current_wang_landau_system.collective_variables[CV_i];
		needed_bins*=int((current_collective_variable->CV_maximum-current_collective_variable->CV_minimum)/current_collective_variable->delta_CV)+1; // plus 1 needed for degrees of association related part of histogram (think of only one acid particle) 
	}
	double* wang_landau_potential =(double*) calloc(1,sizeof(double)*needed_bins); //calloc initializes everything to zero
	return wang_landau_potential;
}

/**
* Returns the degree of association for a given index of the collective variable. This is needed since you may use multiple degrees of association as collective variable for the Wang-Landau algorithm.
*/
double calculate_degree_of_association(int index_of_current_collective_variable, void* _m_current_wang_landau_system){
    wang_landau_system* current_wang_landau_system=(wang_landau_system*) _m_current_wang_landau_system;
	collective_variable* current_collective_variable=current_wang_landau_system->collective_variables[index_of_current_collective_variable];
	int total_number_of_corresponding_acid=0;
	for(int corresponding_type_i=0; corresponding_type_i<current_collective_variable->nr_corresponding_acid_types;corresponding_type_i++){
		int num_of_current_type;
		number_of_particles_with_type(current_collective_variable->corresponding_acid_types[corresponding_type_i],&num_of_current_type);
		total_number_of_corresponding_acid+=num_of_current_type;
	}
	if(total_number_of_corresponding_acid==0){
	    throw std::runtime_error("Have you forgotten to specify all corresponding acid types? Total particle number of corresponding acid type is zero\n");
	}
	int num_of_associated_acid;
	number_of_particles_with_type(current_collective_variable->associated_type,&num_of_associated_acid);
	double degree_of_association=(double) num_of_associated_acid/total_number_of_corresponding_acid; //cast to double because otherwise any fractional part is lost
	return degree_of_association;
}

/**
* Finds the minimum non negative value in the provided double array and returns this value.
*/
double find_minimum_non_negative_value(double* list, int len){
	double minimum =list[0];
	for (int i=0;i<len;i++){
		if(minimum<0)
			minimum=list[i];//think of negative histogram values that indicate not allowed energies in the case of an energy observable
		if(list[i]<minimum && list[i]>=0)
			minimum=list[i];	
	}
	return minimum;
}

/**
* Finds the minimum in a double array and returns it.
*/
double find_minimum(double* list, int len){
	double minimum =list[0];
	for (int i=0;i<len;i++){
		if(list[i]<minimum)
			minimum=list[i];	
	}
	return minimum;
}

/**
* Finds the maximum in a double array and returns it.
*/
double find_maximum(double* list, int len){
	double maximum =list[0];
	for (int i=0;i<len;i++){
		if(list[i]>maximum)
			maximum=list[i];	
	}
	return maximum;
}

/**
* Initializes the current Wang-Landau system.
*/
int ReactionEnsemble::initialize_wang_landau(){
	if(m_current_wang_landau_system.nr_subindices_of_collective_variable!=NULL){
		//initialize_wang_landau() has been called before, free everything that was allocated
		free(m_current_wang_landau_system.nr_subindices_of_collective_variable);
	}
	
	//initialize deltas for collective variables which are of the type of a degree of association
	int energy_collective_variable_index=-10;
	double* min_boundaries_energies=NULL;
	double* max_boundaries_energies=NULL;
	for(int collective_variable_i=0; collective_variable_i<m_current_wang_landau_system.nr_collective_variables;collective_variable_i++){
		collective_variable* current_collective_variable=m_current_wang_landau_system.collective_variables[collective_variable_i];
		if(current_collective_variable->corresponding_acid_types!=NULL){
			//found a collective variable which is of the type of a degree_of_association
			current_collective_variable->delta_CV=calculate_delta_degree_of_association(collective_variable_i);
		}
		
		int flattened_index_previous_run=0; //len_histogram of energy preparation run
		if(current_collective_variable->energy_boundaries_filename!=NULL){
			//found a collective variable which is not of the type of an energy
			m_current_wang_landau_system.do_energy_reweighting=true;
			energy_collective_variable_index=collective_variable_i;
			//load energy boundaries from file
			FILE* pFile;
			pFile = fopen(current_collective_variable->energy_boundaries_filename,"r");
			if (pFile==NULL){
			    throw std::runtime_error("ERROR: energy boundaries file for the specific system could not be read.\n");
				// Note that you cannot change the other collective variables in the pre-production run and the production run
				return ES_ERROR;
			}
			//save minimum and maximum energies as a function of the other collective variables under m_current_wang_landau_system.energ...
			char *line = NULL;
			size_t len = 0;
			ssize_t length_line;
			const char* delim="\t ";
			getline(&line, &len, pFile);//dummy call of getline to get rid of header line (first line in file)
			while ((length_line = getline(&line, &len, pFile)) != -1) {
				int counter_words_in_line=0;
				for(char* word=strtok(line,delim);word!=NULL;word=strtok(NULL,delim)){
					if(counter_words_in_line<m_current_wang_landau_system.nr_collective_variables-1){
						counter_words_in_line+=1;
						continue;
					}else if(counter_words_in_line==m_current_wang_landau_system.nr_collective_variables-1){
						double energy_boundary_minimum=atof(word);
						counter_words_in_line+=1;
						min_boundaries_energies=(double*) realloc(min_boundaries_energies,sizeof(double)*(flattened_index_previous_run+1));			
						min_boundaries_energies[flattened_index_previous_run]=energy_boundary_minimum;
					}else if(counter_words_in_line==m_current_wang_landau_system.nr_collective_variables){
						double energy_boundary_maximum=atof(word);
						max_boundaries_energies=(double*) realloc(max_boundaries_energies,sizeof(double)*(flattened_index_previous_run+1));				
						max_boundaries_energies[flattened_index_previous_run]=energy_boundary_maximum;
						counter_words_in_line+=1;							
					}
					
				}
				flattened_index_previous_run+=1;		
			}
			
			current_collective_variable->CV_minimum=find_minimum(min_boundaries_energies,flattened_index_previous_run);
			current_collective_variable->CV_maximum=find_maximum(max_boundaries_energies,flattened_index_previous_run);
		}
	}


	int* nr_subindices_of_collective_variable=(int*) malloc(m_current_wang_landau_system.nr_collective_variables*sizeof(int));
	for(int collective_variable_i=0;collective_variable_i<m_current_wang_landau_system.nr_collective_variables;collective_variable_i++){
		nr_subindices_of_collective_variable[collective_variable_i]=int((m_current_wang_landau_system.collective_variables[collective_variable_i]->CV_maximum-m_current_wang_landau_system.collective_variables[collective_variable_i]->CV_minimum)/m_current_wang_landau_system.collective_variables[collective_variable_i]->delta_CV)+1; //+1 for collecive variables which are of type degree of association
	}
	m_current_wang_landau_system.nr_subindices_of_collective_variable=nr_subindices_of_collective_variable;

	//construct (possibly higher dimensional) histogram over Gamma (the room which should be equally sampled when the wang-landau algorithm has converged)
	m_current_wang_landau_system.histogram=initialize_histogram();

	//construct (possibly higher dimensional) wang_landau potential over Gamma (the room which should be equally sampled when the wang-landau algorithm has converged)
	m_current_wang_landau_system.wang_landau_potential=initialize_wang_landau_potential();
	
	m_current_wang_landau_system.used_bins=m_current_wang_landau_system.len_histogram; //initialize for 1/t wang_landau algorithm
	
	if(energy_collective_variable_index>=0){
		//make values in histogram and wang landau potential negative if they are not allowed at the given degree of association, because the energy boundaries prohibit them

		int empty_bins_in_memory=0;

		for(int flattened_index=0;flattened_index<m_current_wang_landau_system.len_histogram;flattened_index++){
			//unravel index
			int unraveled_index[m_current_wang_landau_system.nr_collective_variables];
			unravel_index(nr_subindices_of_collective_variable,m_current_wang_landau_system.nr_collective_variables,flattened_index,unraveled_index);
			//use unraveled index
			double current_energy=unraveled_index[energy_collective_variable_index]*m_current_wang_landau_system.collective_variables[energy_collective_variable_index]->delta_CV+m_current_wang_landau_system.collective_variables[energy_collective_variable_index]->CV_minimum;
			if(current_energy>max_boundaries_energies[get_flattened_index_wang_landau_without_energy_collective_variable(flattened_index,energy_collective_variable_index)] || current_energy<min_boundaries_energies[get_flattened_index_wang_landau_without_energy_collective_variable(flattened_index,energy_collective_variable_index)]-m_current_wang_landau_system.collective_variables[energy_collective_variable_index]->delta_CV ){
				m_current_wang_landau_system.histogram[flattened_index]=m_current_wang_landau_system.int_fill_value;
				m_current_wang_landau_system.wang_landau_potential[flattened_index]=m_current_wang_landau_system.double_fill_value;
				empty_bins_in_memory+=1;	
			}
		}
		
		m_current_wang_landau_system.used_bins=m_current_wang_landau_system.len_histogram-empty_bins_in_memory;
		
	}

	free(max_boundaries_energies);
	free(min_boundaries_energies);

	//assign determine_current_state_in_this_collective_variable function pointers to correct function
	for(int collective_variable_i=0; collective_variable_i<m_current_wang_landau_system.nr_collective_variables;collective_variable_i++){
		collective_variable* current_collective_variable=m_current_wang_landau_system.collective_variables[collective_variable_i];
		if(current_collective_variable->corresponding_acid_types!=NULL){
			//found a collective variable which is not of the type of a degree_of_association association)	
			current_collective_variable->determine_current_state_in_collective_variable_with_index=&calculate_degree_of_association;
		}
		if(current_collective_variable->energy_boundaries_filename!=NULL){
			//found a collective variable which is not of the type of an energy
			current_collective_variable->determine_current_state_in_collective_variable_with_index=&calculate_current_potential_energy_of_system_wrap;
		}
		
	}
	
	return ES_OK;
}

/**
* Calculates the expression which occurs in the Wang-Landau acceptance probability.
*/
double ReactionEnsemble::calculate_boltzmann_factor_reaction_ensemble_wang_landau(single_reaction* current_reaction, double E_pot_old, double E_pot_new, std::vector<int>& old_particle_numbers, int old_state_index, int new_state_index, bool only_make_configuration_changing_move){
	/**determine the acceptance probabilities of the reaction move 
	* in wang landau reaction ensemble
	*/
	double beta =1.0/m_current_reaction_system.temperature_reaction_ensemble;
	double bf;
	if(m_current_wang_landau_system.do_not_sample_reaction_partition_function==true || only_make_configuration_changing_move==true){
		bf=1.0;
	}else{
		const double volume = m_current_reaction_system.volume;
		double factorial_expr=calculate_factorial_expression(current_reaction, old_particle_numbers.data());
		double standard_pressure_in_simulation_units=m_current_reaction_system.standard_pressure_in_simulation_units;
		bf=std::pow(volume*beta*standard_pressure_in_simulation_units, current_reaction->nu_bar) * current_reaction->equilibrium_constant * factorial_expr;
	}
	
	if(m_current_wang_landau_system.do_energy_reweighting==false){
		bf= bf * exp(-beta * (E_pot_new - E_pot_old));
	} else {
		//pass
	}
	//look wether the proposed state lies in the reaction coordinate space Gamma and add the Wang-Landau modification factor, this is a bit nasty due to the energy collective variable case (memory layout of storage array of the histogram and the wang_landau_potential values is "cuboid")
	if(old_state_index>=0 && new_state_index>=0){
		if(m_current_wang_landau_system.histogram[new_state_index]>=0 &&m_current_wang_landau_system.histogram[old_state_index]>=0 ){
			bf=std::min(1.0, bf*exp(m_current_wang_landau_system.wang_landau_potential[old_state_index]-m_current_wang_landau_system.wang_landau_potential[new_state_index])); //modify boltzmann factor according to wang-landau algorithm, according to grand canonical simulation paper "Density-of-states Monte Carlo method for simulation of fluids"
			//this makes the new state being accepted with the conditinal probability bf (bf is a transition probability = conditional probability from the old state to move to the new state)
		}else{
			if(m_current_wang_landau_system.histogram[new_state_index]>=0 &&m_current_wang_landau_system.histogram[old_state_index]<0 )
				bf=10;//this makes the reaction get accepted, since we found a state in Gamma
			else if (m_current_wang_landau_system.histogram[new_state_index]<0 &&m_current_wang_landau_system.histogram[old_state_index]<0)
				bf=10;//accept, in order to be able to sample new configs, which might lie in Gamma
			else if(m_current_wang_landau_system.histogram[new_state_index]<0 &&m_current_wang_landau_system.histogram[old_state_index]>=0)
				bf=-10;//this makes the reaction get rejected, since the new state is not in Gamma while the old sate was in Gamma
		}
	}else if(old_state_index<0 && new_state_index>=0){
		bf=10;	//this makes the reaction get accepted, since we found a state in Gamma
	}else if(old_state_index<0 && new_state_index<0){
		bf=10;	//accept, in order to be able to sample new configs, which might lie in Gamma
	}else if(old_state_index>=0 && new_state_index<0){
		bf=-10; //this makes the reaction get rejected, since the new state is not in Gamma while the old sate was in Gamma
	}
	return bf;	
}

void ReactionEnsemble::update_wang_landau_potential_and_histogram(int index_of_state_after_acceptance_or_rejection){
	/**increase the wang landau potential and histogram at the current nbar */
	if(index_of_state_after_acceptance_or_rejection>=0 ){
		if(m_current_wang_landau_system.histogram[index_of_state_after_acceptance_or_rejection]>=0){
			m_current_wang_landau_system.histogram[index_of_state_after_acceptance_or_rejection]+=1;
			m_current_wang_landau_system.wang_landau_potential[index_of_state_after_acceptance_or_rejection]+=m_current_wang_landau_system.wang_landau_parameter;
		}
	}

}


/** Performs a randomly selected reaction using the Wang-Landau algorithm.
*make sure to perform additional configuration changing steps, after the reaction step! like in Density-of-states Monte Carlo method for simulation of fluids Yan, De Pablo. this can be done with MD in the case of the no-energy-reweighting case, or with the functions do_global_mc_move_for_particles_of_type
*perform additional Monte-carlo moves to to sample configurational partition function
*according to "Density-of-states Monte Carlo method for simulation of fluids"
do as many steps as needed to get to a new conformation (compare Density-of-states Monte Carlo method for simulation of fluids Yan, De Pablo)*/
int ReactionEnsemble::do_reaction_wang_landau(){
	m_WL_tries+=m_current_wang_landau_system.wang_landau_steps;
	bool got_accepted=false;
	for(int step=0;step<m_current_wang_landau_system.wang_landau_steps;step++){
		int reaction_id=i_random(m_current_reaction_system.nr_single_reactions);
		got_accepted=generic_oneway_reaction(reaction_id, reaction_ensemble_wang_landau_mode);
		if(got_accepted==true){
			m_WL_accepted_moves+=1;
		}

		if(can_refine_wang_landau_one_over_t()&& m_WL_tries%10000==0){
			//check for convergence
			if(achieved_desired_number_of_refinements_one_over_t()==true){
				ReactionEnsemble::write_wang_landau_results_to_file(m_current_wang_landau_system.output_filename);
				return -10; //return negative value to indicate that the Wang-Landau algorithm has converged
			}
			refine_wang_landau_parameter_one_over_t();
		}
	}
	
	//shift wang landau potential minimum to zero
	if(m_WL_tries%(std::max(90000,9*m_current_wang_landau_system.wang_landau_steps))==0){
		//for numerical stability here we also subtract the minimum positive value of the wang_landau_potential from the wang_landau potential, allowed since only the difference in the wang_landau potential is of interest.
		double minimum_wang_landau_potential=find_minimum_non_negative_value(m_current_wang_landau_system.wang_landau_potential,m_current_wang_landau_system.len_histogram);
		for(int i=0;i<m_current_wang_landau_system.len_histogram;i++){
			if(m_current_wang_landau_system.wang_landau_potential[i]>=0)//check for wether we are in the valid range of the collective variable
				m_current_wang_landau_system.wang_landau_potential[i]-=minimum_wang_landau_potential;	
		}
		
		//write out preliminary wang-landau potential results
		write_wang_landau_results_to_file(m_current_wang_landau_system.output_filename);
	}
	return 0;	
};

/**
*Frees the Wang-Landau data structures.
*/
void ReactionEnsemble::free_wang_landau(){
	free(m_current_wang_landau_system.histogram);
	free(m_current_wang_landau_system.wang_landau_potential);
	for(int CV_i=0;CV_i<m_current_wang_landau_system.nr_collective_variables;CV_i++){
		collective_variable* current_collective_variable=m_current_wang_landau_system.collective_variables[CV_i];
		if(current_collective_variable->corresponding_acid_types!=NULL) { //check wether we have a collective variable which is of the type of a degree of association
			free(current_collective_variable->corresponding_acid_types);
		}
		if(current_collective_variable->energy_boundaries_filename!=NULL){//check wether we have a collective variable which is of the type of an energy
			free(current_collective_variable->energy_boundaries_filename);
		}
		free(current_collective_variable);
	}
	free(m_current_wang_landau_system.collective_variables);
	free(m_current_wang_landau_system.output_filename);
	free(m_current_wang_landau_system.nr_subindices_of_collective_variable);

	if(m_current_wang_landau_system.minimum_energies_at_flat_index!=NULL) //only present in energy preparation run
		free(m_current_wang_landau_system.minimum_energies_at_flat_index);
	if(m_current_wang_landau_system.maximum_energies_at_flat_index!=NULL)
		free(m_current_wang_landau_system.maximum_energies_at_flat_index);
}

//boring helper functions
/**
*Calculates the average of an integer array (used for the histogram of the Wang-Landau algorithm). It excludes values which are initialized to be negative. Those values indicate that the Wang-Landau algorithm should not sample those values. The values still occur in the list because we can only store "rectangular" value ranges.
*/
double ReactionEnsemble::average_int_list(int* int_number_list, int len_int_nr_list){
	double result=0.0;
	int counter_allowed_entries=0;
	for(int i=0;i<len_int_nr_list;i++){
		if(int_number_list[i]>=0){ //checks for validity of index i (think of energy collective variables, in a cubic memory layout there will be indices which are not allowed by the energy boundaries. These values will be initalized with a negative fill value)
			result+=(double) int_number_list[i];
			counter_allowed_entries+=1;
		}
	}
	return result/counter_allowed_entries;
}

/**
*finds the minimum in an integer array and returns it
*/
int ReactionEnsemble::find_minimum_in_int_list(int* list, int len){
	double minimum =list[0];
	for (int i=0;i<len;i++){
		if(minimum<0)
			minimum=list[i];//think of negative histogram values that indicate not allowed energies in the case of an energy observable
		if(list[i]<minimum && list[i]>=0)
			minimum=list[i];
	}
	return minimum;
}
/**
*Determines wether we can reduce the Wang-Landau parameter
*/
bool ReactionEnsemble::can_refine_wang_landau_one_over_t(){
	double minimum_required_value=0.80*average_int_list(m_current_wang_landau_system.histogram,m_current_wang_landau_system.len_histogram); // This is an additional constraint to sample configuration space better. Use flatness criterion according to 1/t algorithm as long as you are not in 1/t regime.
	if(m_current_wang_landau_system.do_energy_reweighting==true)
		minimum_required_value=20; //get faster in energy reweighting case

	if(find_minimum_in_int_list(m_current_wang_landau_system.histogram,m_current_wang_landau_system.len_histogram)>minimum_required_value || m_system_is_in_1_over_t_regime==true){
		return true;
	}else{
		return false;	
	}
}

/**
*Reset the Wang-Landau histogram.
*/
void ReactionEnsemble::reset_histogram(){
	printf("Histogram is flat. Refining. Previous Wang-Landau modification parameter was %f.\n",m_current_wang_landau_system.wang_landau_parameter);
	fflush(stdout);
	
	for(int i=0;i<m_current_wang_landau_system.len_histogram;i++){
		if(m_current_wang_landau_system.histogram[i]>=0){//checks for validity of index i (think of energy collective variables, in a cubic memory layout there will be indices which are not allowed by the energy boundaries. These values will be initalized with a negative fill value)
			m_current_wang_landau_system.histogram[i]=0;
		}	
	}
	
}

/**
*Refine the Wang-Landau parameter using the 1/t rule.
*/
void ReactionEnsemble::refine_wang_landau_parameter_one_over_t(){
	double monte_carlo_time = (double) m_current_wang_landau_system.monte_carlo_trial_moves/m_current_wang_landau_system.used_bins;
	if ( m_current_wang_landau_system.wang_landau_parameter/2.0 <=1.0/monte_carlo_time || m_system_is_in_1_over_t_regime==true){
		m_current_wang_landau_system.wang_landau_parameter= 1.0/monte_carlo_time;
		if(m_system_is_in_1_over_t_regime==false){		
			m_system_is_in_1_over_t_regime=true;
			printf("Refining: Wang-Landau parameter is now 1/t.\n");
		}
	} else {
		reset_histogram();
		m_current_wang_landau_system.wang_landau_parameter= m_current_wang_landau_system.wang_landau_parameter/2.0;
		
	}
}

/**
*Determine whether the desired number of refinements was achieved.
*/
bool ReactionEnsemble::achieved_desired_number_of_refinements_one_over_t() {
	if(m_current_wang_landau_system.wang_landau_parameter < m_current_wang_landau_system.final_wang_landau_parameter) {
		printf("Achieved desired number of refinements\n");
		return true;
	} else {
		return false;
	}

}

/**
*Returns the unraveled index of the provided flattened index (needed for writing the Wang-Landau results to file)
*/
void ReactionEnsemble::unravel_index(int* len_dims, int ndims, int flattened_index, int* unraveled_index_out){
	//idea taken from http://codinghighway.com/2014/02/22/c-multi-dimensional-arrays-part-2-flattened-to-unflattened-index/
	int mul[ndims];
	mul[ndims-1]=1;
	for (int j = ndims-2; j >= 0; j--)
		mul[j] = mul[j+1]*len_dims[j+1];
	for (int j = 0; j < ndims; j++)
		unraveled_index_out[j]=(flattened_index/mul[j])%len_dims[j];
}


/**
*Writes the Wang-Landau potential to file.
*/
void ReactionEnsemble::write_wang_landau_results_to_file(char* full_path_to_output_filename){

	FILE* pFile;
	pFile = fopen(full_path_to_output_filename,"w");
	if (pFile==NULL){
	    throw std::runtime_error("ERROR: Wang-Landau file could not be written\n");
	}else{
		int* nr_subindices_of_collective_variable =m_current_wang_landau_system.nr_subindices_of_collective_variable;
		for(int flattened_index=0;flattened_index<m_current_wang_landau_system.len_histogram;flattened_index++){
			//unravel index
			if(std::abs(m_current_wang_landau_system.wang_landau_potential[flattened_index]-m_current_wang_landau_system.double_fill_value)>1){ //only output data if they are not equal to m_current_reaction_system.double_fill_value. This if ensures that for the energy observable not allowed energies (energies in the interval [global_E_min, global_E_max]) in the multidimensional wang landau potential are printed out, since the range [E_min(nbar), E_max(nbar)] for each nbar may be a different one
				int unraveled_index[m_current_wang_landau_system.nr_collective_variables];
				unravel_index(nr_subindices_of_collective_variable,m_current_wang_landau_system.nr_collective_variables,flattened_index,unraveled_index);
				//use unraveled index
				for(int i=0;i<m_current_wang_landau_system.nr_collective_variables;i++){
					fprintf(pFile, "%f ",unraveled_index[i]*m_current_wang_landau_system.collective_variables[i]->delta_CV+m_current_wang_landau_system.collective_variables[i]->CV_minimum);
				}
				fprintf(pFile, "%f \n", m_current_wang_landau_system.wang_landau_potential[flattened_index]);
			}
		}
		fflush(pFile);
		fclose(pFile);
	}

}

/**
*Update the minimum and maximum observed energies using the current state. Needed for perliminary energy reweighting runs.
*/
int ReactionEnsemble::update_maximum_and_minimum_energies_at_current_state(){
	if(m_current_wang_landau_system.minimum_energies_at_flat_index==NULL || m_current_wang_landau_system.maximum_energies_at_flat_index==NULL){
		m_current_wang_landau_system.minimum_energies_at_flat_index=(double*) calloc(1,sizeof(double)*m_current_wang_landau_system.len_histogram);
		m_current_wang_landau_system.maximum_energies_at_flat_index=(double*) calloc(1,sizeof(double)*m_current_wang_landau_system.len_histogram);
		for (int i = 0; i < m_current_wang_landau_system.len_histogram; i++){
 	 		m_current_wang_landau_system.minimum_energies_at_flat_index[i] =m_current_wang_landau_system.double_fill_value;
 	 		m_current_wang_landau_system.maximum_energies_at_flat_index[i] =m_current_wang_landau_system.double_fill_value;
		}
	}
	
	const double E_pot_current=calculate_current_potential_energy_of_system_wrap(0, NULL);
	int index=get_flattened_index_wang_landau_of_current_state();

	//update stored energy values
	if( (( E_pot_current<m_current_wang_landau_system.minimum_energies_at_flat_index[index])|| std::abs(m_current_wang_landau_system.minimum_energies_at_flat_index[index] -m_current_wang_landau_system.double_fill_value)<std::numeric_limits<double>::epsilon()) ) {
		m_current_wang_landau_system.minimum_energies_at_flat_index[index]=E_pot_current;
	}
	if( ((E_pot_current>m_current_wang_landau_system.maximum_energies_at_flat_index[index]) || std::abs(m_current_wang_landau_system.maximum_energies_at_flat_index[index] -m_current_wang_landau_system.double_fill_value)<std::numeric_limits<double>::epsilon()) ) {
		m_current_wang_landau_system.maximum_energies_at_flat_index[index]= E_pot_current;
	}
	

	return 0;
}

/**
*Write out an energy boundary file using the energy boundaries observed in a preliminary energy reweighting run.
*/
void ReactionEnsemble::write_out_preliminary_energy_run_results (char* full_path_to_output_filename) {
	FILE* pFile;
	pFile = fopen(full_path_to_output_filename,"w");
	if(pFile==NULL){
	    throw std::runtime_error("ERROR: Wang-Landau file could not be written\n");
	}else{
		fprintf(pFile, "#nbar E_min E_max\n");
		int* nr_subindices_of_collective_variable =m_current_wang_landau_system.nr_subindices_of_collective_variable;

		for(int flattened_index=0;flattened_index<m_current_wang_landau_system.len_histogram;flattened_index++){
			//unravel index
			int unraveled_index[m_current_wang_landau_system.nr_collective_variables];
			unravel_index(nr_subindices_of_collective_variable,m_current_wang_landau_system.nr_collective_variables,flattened_index,unraveled_index);
			//use unraveled index
			for(int i=0;i<m_current_wang_landau_system.nr_collective_variables;i++){
				fprintf(pFile, "%f ",unraveled_index[i]*m_current_wang_landau_system.collective_variables[i]->delta_CV+m_current_wang_landau_system.collective_variables[i]->CV_minimum);
			}
			fprintf(pFile, "%f %f \n", m_current_wang_landau_system.minimum_energies_at_flat_index[flattened_index], m_current_wang_landau_system.maximum_energies_at_flat_index[flattened_index]);
		}
		fflush(pFile);
		fclose(pFile);
	}
}


/**
*Returns the flattened index of a given flattened index without the energy collective variable.
*/
int ReactionEnsemble::get_flattened_index_wang_landau_without_energy_collective_variable(int flattened_index_with_energy_collective_variable, int collective_variable_index_energy_observable){
	int* nr_subindices_of_collective_variable=m_current_wang_landau_system.nr_subindices_of_collective_variable;
	//unravel index
	int unraveled_index[m_current_wang_landau_system.nr_collective_variables];
	unravel_index(nr_subindices_of_collective_variable,m_current_wang_landau_system.nr_collective_variables,flattened_index_with_energy_collective_variable,unraveled_index);
	//use unraveled index
	const int nr_collective_variables=m_current_wang_landau_system.nr_collective_variables-1; //forget the last collective variable (the energy collective variable)
	double current_state[nr_collective_variables];
	for(int i=0;i<nr_collective_variables;i++){
		current_state[i]=unraveled_index[i]*m_current_wang_landau_system.collective_variables[i]->delta_CV+m_current_wang_landau_system.collective_variables[i]->CV_minimum;
	}
	
	//get collective_variables_minimum_values
	double collective_variables_minimum_values[nr_collective_variables];
	for(int CV_i=0;CV_i<nr_collective_variables;CV_i++){
		collective_variables_minimum_values[CV_i]=m_current_wang_landau_system.collective_variables[CV_i]->CV_minimum;	
	}
	//get collective_variables_maximum_values
	double collective_variables_maximum_values[nr_collective_variables];
	for(int CV_i=0;CV_i<nr_collective_variables;CV_i++){
		collective_variables_maximum_values[CV_i]=m_current_wang_landau_system.collective_variables[CV_i]->CV_maximum;	
	}
	//get delta_collective_variables_values
	double delta_collective_variables_values[nr_collective_variables];
	for(int CV_i=0;CV_i<nr_collective_variables;CV_i++){
		delta_collective_variables_values[CV_i]=m_current_wang_landau_system.collective_variables[CV_i]->delta_CV;	
	}
	int index=get_flattened_index_wang_landau(current_state, collective_variables_minimum_values, collective_variables_maximum_values, delta_collective_variables_values, nr_collective_variables);
	return index;
}

/** remove bins from the range of to be sampled values if they have not been sampled.
*use with caution otherwise you produce unpyhsical results, do only use when you know what you want to do. This can make wang landau converge on a reduced set Gamma. use this function e.g. in do_reaction_wang_landau() for the diprotonic acid
*compare "Wang-Landau sampling with self-adaptive range" by Troester and Dellago
*/
void ReactionEnsemble::remove_bins_that_have_not_been_sampled(){
	int removed_bins=0;
	double beta=1.0/m_current_reaction_system.temperature_reaction_ensemble;
	double largest_wang_landau_potential_at_given_particle_number=find_maximum(m_current_wang_landau_system.wang_landau_potential,m_current_wang_landau_system.len_histogram);
	for(int k=0;k<m_current_wang_landau_system.len_histogram;k++){
		if(m_current_wang_landau_system.wang_landau_potential[k]==0){
			removed_bins+=1;
			// criterion is derived from the canonical partition function and the ration of two summands for the same particle number
			m_current_wang_landau_system.histogram[k]=m_current_wang_landau_system.int_fill_value;
			m_current_wang_landau_system.wang_landau_potential[k]=m_current_wang_landau_system.double_fill_value;
		}
		
	}
	printf("Removed %d bins from the Wang-Landau spectrum\n",removed_bins);
	//update used bins
	m_current_wang_landau_system.used_bins-=removed_bins;
}

/**
*Writes the Wang-Landau parameter, the histogram and the potential to a file. You can restart a Wang-Landau simulation using this information. Additionally you should store the positions of the particles. Not storing them introduces small, small statistical errors.
*/
int ReactionEnsemble::write_wang_landau_checkpoint(char* identifier){
	std::ofstream outfile;

	//write current wang landau parameters (wang_landau_parameter, monte_carlo_trial_moves, flat_index_of_current_state)
	outfile.open(std::string("checkpoint_wang_landau_parameters_")+identifier);
	outfile << m_current_wang_landau_system.wang_landau_parameter << " " << m_current_wang_landau_system.monte_carlo_trial_moves << " " << get_flattened_index_wang_landau_of_current_state() << "\n" ;	
	outfile.close();

	//write histogram
	outfile.open(std::string("checkpoint_wang_landau_histogram_")+identifier);
	for(int i=0;i<m_current_wang_landau_system.len_histogram;i++){
		outfile << m_current_wang_landau_system.histogram[i] <<"\n" ;
	}
	outfile.close();
	//write wang landau potential
	outfile.open(std::string("checkpoint_wang_landau_potential_")+identifier);
	for(int i=0;i<m_current_wang_landau_system.len_histogram;i++){
		outfile << m_current_wang_landau_system.wang_landau_potential[i]<<"\n";	
	}
	outfile.close();
	return 0;
}


/**
*Loads the Wang-Landau checkpoint
*/
int ReactionEnsemble::load_wang_landau_checkpoint(char* identifier){
	std::ifstream infile;
	
	//restore wang landau parameters
	infile.open(std::string("checkpoint_wang_landau_parameters_")+identifier);
	if(infile.is_open()) {
		
		double wang_landau_parameter_entry;
		int wang_landau_monte_carlo_trial_moves_entry;
		int flat_index_of_state_at_checkpointing;
		int line=0;
		while (infile >> wang_landau_parameter_entry >> wang_landau_monte_carlo_trial_moves_entry >> flat_index_of_state_at_checkpointing)
		{
			m_current_wang_landau_system.wang_landau_parameter=wang_landau_parameter_entry;
			m_current_wang_landau_system.monte_carlo_trial_moves=wang_landau_monte_carlo_trial_moves_entry;
			line+=1;
		}
		infile.close();
	} else {
		std::cout << "Exception opening " << std::string("checkpoint_wang_landau_parameters_")+identifier  << "\n" << std::flush;
	}
	
	//restore histogram
	infile.open(std::string("checkpoint_wang_landau_histogram_")+identifier);
	if(infile.is_open()){ 
		int hist_entry;
		int line=0;
		while (infile >> hist_entry)
		{
			m_current_wang_landau_system.histogram[line]=hist_entry;
			line+=1;
		}
		infile.close();
	} else {
		std::cout << "Exception opening/ reading " << std::string("checkpoint_wang_landau_histogram_")+identifier << "\n" << std::flush;
	}
	
	//restore wang landau potential
	infile.open(std::string("checkpoint_wang_landau_potential_")+identifier);
	if(infile.is_open()) {	
		double wang_landau_potential_entry;
		int line=0;
		while (infile >> wang_landau_potential_entry)
		{
			m_current_wang_landau_system.wang_landau_potential[line]=wang_landau_potential_entry;
			line+=1;
		}
		infile.close();
	} else {
		std::cout << "Exception opening " << std::string("checkpoint_wang_landau_potential_")+identifier  << "\n" << std::flush;
	}
	
	//possible task: restore state in which the system was when the checkpoint was written. However as long as checkpointing and restoring the system form the checkpoint is rare this should not matter statistically.
	
	return 0;
}



int ReactionEnsemble::get_random_p_id(){
    int random_p_id = i_random(max_seen_particle); 
    while(is_in_list(random_p_id, m_empty_p_ids_smaller_than_max_seen_particle.data(), m_empty_p_ids_smaller_than_max_seen_particle.size()))
        random_p_id = i_random(max_seen_particle);
    return random_p_id;
}



/**
* Constant-pH Ensemble, for derivation see Reed and Reed 1992
* For the constant pH reactions you need to provide the deprotonation and afterwards the corresponding protonation reaction (in this order). If you want to deal with multiple reactions do it multiple times.
* Note that there is a difference in the usecase of the constant pH reactions and the above reaction ensemble. For the constant pH simulation directily the **apparent equilibrium constant which carries a unit** needs to be provided -- this is different from the reaction ensemble above, where the dimensionless reaction constant needs to be provided. Again: For the constant-pH algorithm not the dimensionless reaction constant needs to be provided here, but the apparent reaction constant.
*/

/**
*Performs a reaction in the constant pH ensemble
*/
int ReactionEnsemble::do_reaction_constant_pH(){
	//get a list of reactions where a randomly selected particle type occurs in the reactant list. the selection probability of the particle types has to be proportional to the number of occurances of the number of particles with this type
	
	//for optimizations this list could be determined during the initialization
	int* list_of_reaction_ids_with_given_reactant_type=NULL;
	int found_reactions_with_given_reactant_type=0;
	while(found_reactions_with_given_reactant_type==0) { // avoid selecting a (e.g. salt) particle which does not take part in a reaction
		int random_p_id =get_random_p_id(); // only used to determine which reaction is attempted.
		auto part = get_particle_data(random_p_id);

		if(!part)
			continue;

		int type_of_random_p_id = part->p.type;

		//construct list of reactions with the above reactant type
		for(int reaction_i=0;reaction_i<m_current_reaction_system.nr_single_reactions;reaction_i++){
			single_reaction* current_reaction=m_current_reaction_system.reactions[reaction_i];
			for(int reactant_i=0; reactant_i< 1; reactant_i++){ //reactant_i<1 since it is assumed in this place that the types A, and HA occur in the first place only. These are the types that should be switched, H+ should not be switched
				if(current_reaction->reactant_types[reactant_i]== type_of_random_p_id){
					found_reactions_with_given_reactant_type+=1;
					list_of_reaction_ids_with_given_reactant_type=(int*) realloc(list_of_reaction_ids_with_given_reactant_type, sizeof(int)*found_reactions_with_given_reactant_type);
					list_of_reaction_ids_with_given_reactant_type[found_reactions_with_given_reactant_type-1]=reaction_i;
					break;
				}
			}
		}
	}
	
	//randomly select a reaction to be performed
	int reaction_id=list_of_reaction_ids_with_given_reactant_type[i_random(found_reactions_with_given_reactant_type)];
	free(list_of_reaction_ids_with_given_reactant_type);
	generic_oneway_reaction(reaction_id, constant_pH_mode);
	return 0;
}

double ReactionEnsemble::calculate_boltzmann_factor_consant_pH(single_reaction* current_reaction, double E_pot_old, double E_pot_new){
	/**
    *Calculates the expression in the acceptance probability of the constant pH method.
    */
	double ln_bf;
	double pKa;
	const double beta =1.0/m_current_reaction_system.temperature_reaction_ensemble;
	if(current_reaction->nu_bar > 0){ //deprotonation of monomer
		pKa = -log10(current_reaction->equilibrium_constant);
		ln_bf= (E_pot_new - E_pot_old)- 1.0/beta*log(10)*(m_constant_pH-pKa) ;
	}else{ //protonation of monomer (yields neutral monomer)
		pKa = -(-log10(current_reaction->equilibrium_constant)); //additional minus, since in this case 1/Ka is stored in the equilibrium constant
		
		ln_bf= (E_pot_new - E_pot_old)+ 1.0/beta*log(10)*(m_constant_pH-pKa);
	}
	double bf=exp(-beta*ln_bf);
	return bf;
}

}
