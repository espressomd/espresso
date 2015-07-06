//this file is just a reference and contains all necessary information on how to use the reaction ensemble
//method according to smith94x
//as example a acid dissociation is implemented
//NOTE: a reaction here is one trial move to dissociate one acid molecule to its dissociated form 
//so if the reaction is accepted then there is one more dissociated ion pair H+ and A-
//NOTE: generic_oneway_reaction does not break bonds for simple reactions. as long as there are no reactions like 2A -->B where one of the reacting A particles occurs in the polymer
//NOTE: paricle types have to start at zero and have to increase by one for every type otherwise the function hide_particle cannot work correctly. Through adding 100 in the function hide_particle we ensure that the function works correctly if the particle types are monotonically increasing and if the largest particle type is smaller than 100.

#include "reaction_ensemble.hpp"
#include "random.hpp"
#include "energy.hpp"	//for energies
#include "external_potential.hpp" //for energies
#include "global.hpp" //for access to global variables
#include "particle_data.hpp" //for particle creation, modification
#include "statistics.hpp" //for mindist
#include "integrate.hpp" //for integrate (hack for updating particle lists)
#include <stdlib.h>  /* qsort() */

//For now the reaction ensemble is only implemented for the reaction VT ensemble. The reaction PT ensemble is also possible to implement.

reaction_system current_reaction_system={.nr_single_reactions=0, .reactions=NULL,.volume=0 , .type_index=NULL, .nr_different_types=0, .charges_of_types=NULL, .water_type=-100}; //initialize watertype to negative number, for checking wether it has been assigned


//declaration of borin helper functions
float factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i);
int calculate_nu_bar(int* educt_coefficients, int len_educt_types,  int* product_coefficients, int len_product_types); //should onlny be used at when defining a new reaction
bool all_educt_particles_exist(int reaction_id);
int generic_oneway_reaction(int reaction_id);
int find_index_of_type(int type);
int replace(int p_id, int desired_type);
double convert_conc_mol_to_vol(double len_sim, double len_real);
int create_particle(int desired_type);
int vec_random(double* vecrandom, double desired_length);
int hide_particle(int p_id, int previous_type);
int delete_particle (int p_id);
int intcmp(const void *aa, const void *bb);



int create_current_reaction_system_struct(){
	single_reaction** reactions_list=(single_reaction**) malloc(sizeof(single_reaction*));
	current_reaction_system.reactions=reactions_list;
	return 0;
}

int do_reaction(){
	int reaction_id=i_random(current_reaction_system.nr_single_reactions);
	generic_oneway_reaction(reaction_id);
	return 0;
}


int _type_is_in_list(int type, int* list, int len_list){
	bool is_in_list=false;
	for(int i=0;i<len_list;i++){
		if(list[i]==type){
			is_in_list=true;
			break;
		}
	}
	return is_in_list;
}

int initialize(){
	int* already_initialized_types=(int*)malloc(sizeof(int));
	int len_already_initialized_types=0;
	for(int single_reaction_i=0; single_reaction_i<current_reaction_system.nr_single_reactions;single_reaction_i++){
		for(int educt_types_i=0; educt_types_i<current_reaction_system.reactions[single_reaction_i]->len_educt_types;educt_types_i++){
			int type=current_reaction_system.reactions[single_reaction_i]->educt_types[educt_types_i];
			if(_type_is_in_list(type,already_initialized_types,len_already_initialized_types)==false)
				init_type_array(type);
		}
	}
	current_reaction_system.charges_of_types =(double*) malloc(sizeof(double)*current_reaction_system.nr_different_types);

	return 0;
}

int free_reaction_ensemble(){
		//needs to be called at the end of the simulation
		for(int single_reaction_i=0;single_reaction_i<current_reaction_system.nr_single_reactions;single_reaction_i++){
			//free educt types and coefficients
			free(current_reaction_system.reactions[single_reaction_i]->educt_types);
			free(current_reaction_system.reactions[single_reaction_i]->educt_coefficients);	
			//free product types and coefficients
			free(current_reaction_system.reactions[single_reaction_i]->product_types);
			free(current_reaction_system.reactions[single_reaction_i]->product_coefficients);	
		}
	free(current_reaction_system.reactions);
	free(current_reaction_system.type_index);
	free(current_reaction_system.charges_of_types);
	return 0;
}


//boring helper functions

float factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i) {
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


bool all_educt_particles_exist(int reaction_id) {
	bool enough_particles=true;
	for(int i=0;i<current_reaction_system.reactions[reaction_id]->len_educt_types;i++){
		int current_number;
		number_of_particles_with_type(current_reaction_system.reactions[reaction_id]->educt_types[i], &current_number);
		if(current_number<current_reaction_system.reactions[reaction_id]->educt_coefficients[i]){
			enough_particles=false;
			break;
		}
	}
	return enough_particles;
}

int generic_oneway_reaction(int reaction_id){
	float volume = current_reaction_system.volume;
	single_reaction* current_reaction=current_reaction_system.reactions[reaction_id];
	//type_H2O=current_reaction_system.water_type;
	
	//generic one way reaction
	//A+B+...+G +... --> K+...X + Z +...
	//you need to use 2A --> B instead of A+A --> B since in the last case you assume distinctness of the particles
	//further it is crucial for the function in which order you provide the educt and product types since particles will be replaced correspondingly!
	
	if (all_educt_particles_exist(reaction_id) ==false ) {
		//makes sure, no incomplete reaction is performed -> only need to consider rollback of complete reactions
		return 0;
	}
	
	//calculate potential energy
	if (total_energy.init_status == 0) {
		init_energies(&total_energy);
		master_energy_calc();
	}
  	int num_energies=total_energy.data.n;
	double kinetic_energy =total_energy.data.e[0];
	double sum_all_energies=0;
	for(int i=0;i<num_energies;i++){
		sum_all_energies+= total_energy.data.e[i];
	}
	for (int i = 0; i < n_external_potentials; i++) {
        	sum_all_energies += external_potentials[i].energy;
        }
	
	double E_pot_old=sum_all_energies-kinetic_energy; //only consider potential energy since we assume that the kinetic part drops out in the process of calculating ensemble averages (kinetic part may be seperated and crossed out)
	
	//find reacting molecules in educts and save their properties for later recreation if step is not accepted
	//do reaction
	//save old particle_numbers
	int* old_particle_numbers=(int*) malloc(sizeof(int) *current_reaction_system.nr_different_types);
	for(int type_index=0;type_index<current_reaction_system.nr_different_types;type_index++){
		number_of_particles_with_type(current_reaction_system.type_index[type_index], &(old_particle_numbers[type_index])); // here could be optimized by not going over all types but only the types that occur in the reaction
	}
	if(current_reaction_system.water_type>=0){
		//set number of water molecules to typical value 55.5 mol/l
		//see https://de.wikipedia.org/wiki/Eigenschaften_des_Wassers#Ionenprodukt
		int index_of_water_type=find_index_of_type(current_reaction_system.water_type);
		old_particle_numbers[index_of_water_type]=int(convert_conc_mol_to_vol(1,1) *55.5*current_reaction_system.volume); // TODO need for correct call of convert_conc_mol_to_vol, setup current_reaction_system.len_sim, current_reaction_system.len_real
	}
		
	int* p_ids_created_particles =NULL;
	int len_p_ids_created_particles=0;
	float* hidden_particles_properties=NULL;
	int len_hidden_particles_properties=0;
	float* changed_particles_properties=NULL;
	int len_changed_particles_properties=0;
	int number_of_saved_properties=3; //save p_id, charge and type of the educt particle, only thing we need to hide the particle and recover it
	
	//create or hide particles of types with corresponding types in reaction
	for(int i=0;i<min(current_reaction->len_product_types,current_reaction->len_educt_types);i++){
		//change min(educt_coefficients(i),product_coefficients(i)) many particles of educt_types(i) to product_types(i)
		for(int j=0;j<min(current_reaction->product_coefficients[i],current_reaction->educt_coefficients[i]);j++){
			int p_id ;
			find_particle_type(current_reaction->educt_types[i], &p_id);
			if(changed_particles_properties==NULL){
				changed_particles_properties=(float*) malloc(sizeof(float)*(len_changed_particles_properties+1)*number_of_saved_properties);			
			}else{
				realloc(changed_particles_properties,sizeof(float)*(len_changed_particles_properties+1)*number_of_saved_properties);
			}
			changed_particles_properties[len_changed_particles_properties*number_of_saved_properties]=(float) p_id;
			changed_particles_properties[len_changed_particles_properties*number_of_saved_properties+1]= (float) current_reaction_system.charges_of_types[find_index_of_type(current_reaction->educt_types[i])];
			changed_particles_properties[len_changed_particles_properties*number_of_saved_properties+2]=(float) current_reaction->educt_types[i];
			len_changed_particles_properties+=3;
			replace(p_id,current_reaction->product_types[i]);
		}
		//create product_coefficients(i)-educt_coefficients(i) many product particles iff product_coefficients(i)-educt_coefficients(i)>0,
		//iff product_coefficients(i)-educt_coefficients(i)<0, hide this number of educt particles
		if ( current_reaction->product_coefficients[i]-current_reaction->educt_coefficients[i] >0) {
			for(int j=0; j< current_reaction->product_coefficients[i]-current_reaction->educt_coefficients[i] ;j++) {
				int p_id=-1;
				while(p_id == -1) {
					p_id=create_particle(current_reaction->product_types[i]);
				}
				if(p_ids_created_particles==NULL){
					p_ids_created_particles=(int*) malloc(sizeof(int)*(len_p_ids_created_particles+1));
				}else{
					realloc(p_ids_created_particles,sizeof(int)*(len_p_ids_created_particles+1));
				}
				p_ids_created_particles[len_p_ids_created_particles]=p_id;
				len_p_ids_created_particles+=1;
			}
		} else if (current_reaction->educt_coefficients[i]-current_reaction->product_coefficients[i] >0) {
			for(int j=0;j<current_reaction->educt_coefficients[i]-current_reaction->product_coefficients[i];j++) {
				int p_id ;
				find_particle_type(current_reaction->educt_types[i], &p_id);
				if(hidden_particles_properties==NULL){
					hidden_particles_properties=(float*) malloc(sizeof(float)*(len_hidden_particles_properties+1)*number_of_saved_properties);			
				}else{
					realloc(hidden_particles_properties,sizeof(float)*(len_hidden_particles_properties+1)*number_of_saved_properties);
				}
				hidden_particles_properties[len_hidden_particles_properties*number_of_saved_properties]=(float) p_id;
				hidden_particles_properties[len_hidden_particles_properties*number_of_saved_properties+1]= (float) current_reaction_system.charges_of_types[find_index_of_type(current_reaction->educt_types[i])];
				hidden_particles_properties[len_hidden_particles_properties*number_of_saved_properties+2]=(float) current_reaction->educt_types[i];
				len_hidden_particles_properties+=3;
				hide_particle(p_id,current_reaction->educt_types[i]);
			}
		}

	}

	//create or hide particles of types with noncorresponding replacement types
	for(int i=min(current_reaction->len_product_types,current_reaction->len_educt_types);i< max(current_reaction->len_product_types,current_reaction->len_educt_types);i++ ) {
		if(current_reaction->len_product_types<current_reaction->len_educt_types){
			//hide superfluous educt_types particles
			for(int j=0;j<current_reaction->educt_coefficients[i];j++){
				int p_id ;
				find_particle_type(current_reaction->educt_types[i], &p_id);
				if(hidden_particles_properties==NULL){
					hidden_particles_properties=(float*) malloc(sizeof(float)*(len_hidden_particles_properties+1)*number_of_saved_properties);			
				}else{
					realloc(hidden_particles_properties,sizeof(float)*(len_hidden_particles_properties+1)*number_of_saved_properties);
				}
				hidden_particles_properties[len_hidden_particles_properties*number_of_saved_properties]=(float) p_id;
				hidden_particles_properties[len_hidden_particles_properties*number_of_saved_properties+1]= (float) current_reaction_system.charges_of_types[find_index_of_type(current_reaction->educt_types[i])];
				hidden_particles_properties[len_hidden_particles_properties*number_of_saved_properties+2]=(float) current_reaction->educt_types[i];
				len_hidden_particles_properties+=3;
				hide_particle(p_id,current_reaction->educt_types[i]);
			}
		} else {
			//create additional product_types particles
			for(int j=0;j<current_reaction->product_coefficients[i];j++){
				int p_id = -1;
				while (p_id == -1) {
					p_id= create_particle(current_reaction->product_types[i]);
				}
				if(p_ids_created_particles==NULL){
					p_ids_created_particles=(int*) malloc(sizeof(int)*(len_p_ids_created_particles+1));				
				}else{
					realloc(p_ids_created_particles,sizeof(int)*(len_p_ids_created_particles+1));
				}
				p_ids_created_particles[len_p_ids_created_particles]=p_id;
				len_p_ids_created_particles+=1;
			}
		}
	}
	
	for (int i=0;i<current_reaction_system.nr_different_types;i++)
		update_particle_array(current_reaction_system.type_index[i]);

	if (total_energy.init_status == 0) {
		init_energies(&total_energy);
		master_energy_calc();
	}	
	//calculate potential energy
	//master_energy_calc(); //not needed ??? only called once in tcl energy.cpp
  	num_energies=total_energy.data.n;
	kinetic_energy =total_energy.data.e[0];
	sum_all_energies=0;
	for(int i=0;i<num_energies;i++){
		sum_all_energies+= total_energy.data.e[i];
	}
	for (int i = 0; i < n_external_potentials; i++) {
        	sum_all_energies += external_potentials[i].energy;
        }
	
	double E_pot_new=sum_all_energies-kinetic_energy;
	double factorial_expr=1.0;
	//factorial contribution of educts
	for(int i=0;i<current_reaction->len_educt_types;i++) {
		int nu_i=-1*current_reaction->educt_coefficients[i];
		int N_i0= old_particle_numbers[find_index_of_type(current_reaction->educt_types[i])];
		factorial_expr=factorial_expr*factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N_i0,nu_i); //zeta = 1 (see smith paper) since we only perform one reaction at one call of the function
	}
	//factorial contribution of products
	for(int i=0;i<current_reaction->len_product_types;i++) {
		int nu_i=current_reaction->product_coefficients[i];
		int N_i0= old_particle_numbers[find_index_of_type(current_reaction->product_types[i])];
		factorial_expr=factorial_expr*factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N_i0,nu_i); //zeta = 1 (see smith paper) since we only perform one reaction at one call of the function
	}

	double beta =1.0;//XXX has to be obtained from global scope
	double standard_pressure_in_simulation_units=0.00108;// XXX has to be obtained from tcl script
	//calculate boltzmann factor
	double bf= pow(current_reaction_system.volume*beta*standard_pressure_in_simulation_units, current_reaction->nu_bar) * current_reaction->equilibrium_constant * factorial_expr * exp(-beta * (E_pot_new - E_pot_old));
	bf=1;
	if ( d_random() < bf ) {
		//accept
		//delete hidden educt_particles (remark: dont delete changed particles)
		int n_HA;
		int type=1;
		number_of_particles_with_type(type,&n_HA);
		if(n_HA>0){
		Particle part;
		int p_id_HA;
		find_particle_type(type,&p_id_HA);
		get_particle_data(p_id_HA,&part);
		}
		for(int i=0;i<len_hidden_particles_properties;i+=number_of_saved_properties) {
			int p_id = (int) hidden_particles_properties[i];
			delete_particle(p_id); //delete particle
		}
		for (int i=0;i<current_reaction_system.nr_different_types;i++)
			update_particle_array(current_reaction_system.type_index[i]);	
		return 1;
	} else {
		//reject
		//reverse reaction
		//1) delete created product particles
		qsort(p_ids_created_particles, len_p_ids_created_particles, sizeof(int), intcmp); // needed since delete_particle changes particle p_ids. start deletion from the largest p_id onwards
		for(int i=0;i<len_p_ids_created_particles;i++){
			delete_particle(p_ids_created_particles[i]);
		}
		for (int i=0;i<current_reaction_system.nr_different_types;i++)
			update_particle_array(current_reaction_system.type_index[i]);
		//2)restore previously hidden educt particles
		for(int i=0;i<len_hidden_particles_properties;i+=number_of_saved_properties) {
			int p_id = (int) hidden_particles_properties[i];
			double charge=(double) hidden_particles_properties[i+1];
			double type=(double) hidden_particles_properties[i+2];
			//set charge
			set_particle_q(p_id, charge);
			//set type
			set_particle_type(p_id, type);
		}
		//2)restore previously changed educt particles
		for(int i=0;i<len_changed_particles_properties;i+=number_of_saved_properties) {
			int p_id = (int) changed_particles_properties[i];
			double charge=(double) changed_particles_properties[i+1];
			double type=(double) changed_particles_properties[i+2];
			//set charge
			set_particle_q(p_id, charge);
			//set type
			set_particle_type(p_id, type);
		}
		for (int i=0;i<current_reaction_system.nr_different_types;i++)
			update_particle_array(current_reaction_system.type_index[i]);
		return 0;
	}


}

int convert_conc_mol_to_vol();
int convert_apparent_to_dimensionless_equilibrium_constant();

int calculate_nu_bar(int* educt_coefficients, int len_educt_types,  int* product_coefficients, int len_product_types){
	//should only be used at when defining a new reaction
	int nu_bar =0;
	for(int i=0;i<len_educt_types;i++){
		nu_bar-=educt_coefficients[i];
	}
	for(int i=0;i<len_product_types;i++){
		nu_bar-=product_coefficients[i];
	}
	return nu_bar;
}

int update_type_index(int* educt_types, int len_educt_types, int* product_types, int len_product_types){
	//should only be used at when defining a new reaction
	if(current_reaction_system.type_index==NULL){
		current_reaction_system.type_index=(int*) malloc(sizeof(int));
		current_reaction_system.type_index[0]=educt_types[0];
		current_reaction_system.nr_different_types=1;
	}
	for (int i =0; i<len_educt_types;i++){
		bool educt_type_i_is_known=false;
		for (int known_type_i=0;known_type_i<current_reaction_system.nr_different_types;known_type_i++){
			if(current_reaction_system.type_index[known_type_i]==educt_types[i]){
				educt_type_i_is_known=true;
				break;
			}
		}
		if (educt_type_i_is_known==false){
			current_reaction_system.type_index=(int*) realloc(current_reaction_system.type_index, sizeof(int)*(current_reaction_system.nr_different_types+1));
			current_reaction_system.type_index[current_reaction_system.nr_different_types]=educt_types[i];
			current_reaction_system.nr_different_types+=1;
		}
	}
	
	for (int i =0; i<len_product_types;i++){
		bool product_type_i_is_known=false;
		for (int known_type_i=0;known_type_i<current_reaction_system.nr_different_types;known_type_i++){
			if(current_reaction_system.type_index[known_type_i]==product_types[i]){
				product_type_i_is_known=true;
				break;
			}
		}
		if (product_type_i_is_known==false){
			current_reaction_system.type_index=(int*) realloc(current_reaction_system.type_index, sizeof(int)*(current_reaction_system.nr_different_types+1));
			current_reaction_system.type_index[current_reaction_system.nr_different_types]=product_types[i];
			current_reaction_system.nr_different_types+=1;
		}
	}
	
	return 0;
}

int find_index_of_type(int type){
	int index =-100; //initialize to invalid index
	for(int i=0; i<current_reaction_system.nr_different_types;i++){
		if(type==current_reaction_system.type_index[i]){
			index=i;
			break;
		}
	}
	return index;
}

double convert_conc_mol_to_vol(double len_sim, double len_real){
	//TODO needs to be implemented
	return 1;
}

int replace(int p_id, int desired_type){
	set_particle_type(p_id, desired_type);
	int err_code=set_particle_q(p_id, (double) current_reaction_system.charges_of_types[find_index_of_type(desired_type)]);
	return 0;
}

int create_particle(int desired_type){
	//remark might only work for cubic box
	int p_id=max_seen_particle+1;
	double pos_x=box_l[0]*d_random();
	double pos_y=box_l[1]*d_random();
	double pos_z=box_l[2]*d_random();

	//create random velocity vector
	double random_vel_vec[3];
	double pi=3.14159265359;
	double temperature =1; //XXX needs to be obtained from global
	double mean_abs_velocity=sqrt(8*temperature/pi);
	vec_random(random_vel_vec, mean_abs_velocity);
	
	double charge= (double) current_reaction_system.charges_of_types[find_index_of_type(desired_type)];
	double d_min=0;
	double d_old=-1;
	int max_insert_tries=1000;
	int insert_tries=0;
	double sig=1;//XXX needs to be obtained from global
	double min_dist=0.9*sig;
	int err_code=ES_PART_ERROR;
	if(min_dist!=0){
		while(d_min<min_dist && insert_tries<max_insert_tries) {
			double pos_vec[3]={pos_x,pos_y,pos_z};
			err_code=place_particle(p_id,pos_vec);
			//set type
			set_particle_type(p_id, desired_type);
			//set charge
			set_particle_q(p_id, charge);
			//set velocities
			set_particle_v(p_id,random_vel_vec);

			d_min=mindist(NULL,NULL); //TODO this seems to be buggy in gdb
			if(d_old==d_min)
				break; //accept, mininal distance caused by other particle
			insert_tries+=1;
			pos_x=box_l[0]*d_random();
			pos_y=box_l[1]*d_random();
			pos_z=box_l[2]*d_random();
			d_old=d_min;
		}
	}else{
		double pos_vec[3]={pos_x,pos_y,pos_z};
		err_code=place_particle(p_id,pos_vec);
		//set type
		set_particle_type(p_id, desired_type);	
		//set velocities
		set_particle_v(p_id,random_vel_vec);
		//set charge
		set_particle_q(p_id, charge);
	}
	
	if(insert_tries>max_insert_tries){
		printf("Error: Particle not inserted, ES_PART err_code %d\n", err_code);
		return -1;	
	}
	return p_id;
}


int hide_particle(int p_id, int previous_type){
	//remove_charge and put type to a non existing one --> no interactions anymore (not even bonds contribute to energy) it is as if the particle was non existing
	//set charge
	set_particle_q(p_id, 0);
	//set type
	int desired_type=previous_type+50;//+50 in order to assign types that are out of the "usual" range of types
	set_particle_type(p_id, desired_type);
	return 0;
}


//copies last particle to that position and deletes the then last one
int delete_particle (int p_id) {
  if (p_id == max_seen_particle) {
	// last particle, just delete
	remove_particle(p_id);
  } else {
	// otherwise, copy properties of last particle to particle with given pid, delete last particle
	// this avoids that the particle identities get excessive

	//read from last particle
	Particle last_particle;
	get_particle_data(max_seen_particle,&last_particle);
	//set pos
	double ppos[3];
	int img[3];
	memmove(ppos, last_particle.r.p, 3*sizeof(double));
	memmove(img, last_particle.l.i, 3*sizeof(int));
	unfold_position(ppos, img);

	//write to particle with p_id
	//set pos
	place_particle(p_id,ppos);
	//set velocities
	set_particle_v(p_id,last_particle.m.v);
	//set charge
	set_particle_q(p_id,last_particle.p.q);
	//set type
	set_particle_type(p_id,last_particle.p.type);


	p_id=max_seen_particle;
	remove_particle(p_id);
  }
return 0;
}

//the following 3 functions are directly taken from ABHmath.tcl
double vec_len(double* vec){
	double veclen=0;
	for(int i=0;i<3;i++){
		veclen+=pow(vec[i],2);	
	}
	return sqrt(veclen);
}

int vecnorm(double* vec, double desired_length){
	double veclen=vec_len(vec);
	for(int i=0;i<3;i++){
		vec[i]=vec[i]/veclen*desired_length;	
	}
	return 0;
}

int vec_random(double* vecrandom, double desired_length){
	//returns a random vector of length len
	//(uniform distribution on a sphere)
	//This is done by chosing 3 uniformly distributed random numbers [-1,1]
	//If the length of the resulting vector is <= 1.0 the vector is taken and normalized
	//to the desired length, otherwise the procedure is repeated until succes.
	//On average the procedure needs 5.739 random numbers per vector.
	//(This is probably not the most efficient way, but it works!)
	//Ask your favorit mathematician for a proof!
	while(1){
		for(int i=0;i<3;i++){
			vecrandom[i]=2*d_random()-1.0;
		}
		if (vec_len(vecrandom)<=1)
			break;
	}
	vecnorm(vecrandom,desired_length);
	return 0;
}

//for sorting integer lists, derived from http://rosettacode.org/wiki/Sort_an_integer_array#C, function is licensed under GPL
int intcmp(const void *aa, const void *bb){ 
    const int *a = (int*) aa, *b =(int*)  bb;
    return (*a > *b) ? -1 : (*a < *b);
}


///////////////////////////////////////////// Reaction Ensemble Wang-Landau algorithm

