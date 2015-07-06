#ifndef REACTION_ENSEMBLE_H
#define REACTION_ENSEMBLE_H

#include "utils.hpp"

typedef struct single_reaction{
	//strict input to the algorithm
	int* educt_types;
	int len_educt_types;
	int* educt_coefficients;
	int* product_types;
	int len_product_types;
	int* product_coefficients;
	double equilibrium_constant;
	//calculated values that are stored for performance reasons
	int nu_bar;
}  single_reaction;

typedef struct reaction_system {
	int nr_single_reactions;
	single_reaction** reactions;
	float volume;
	int* type_index;
	int nr_different_types; // is equal to length type_index
	double* charges_of_types;
	int water_type; //TODO needs to be used
} reaction_system;

extern reaction_system current_reaction_system;

int do_reaction();

int create_current_reaction_system_struct();

int free_reaction_ensemble();

int initialize();

int calculate_nu_bar(int* educt_coefficients, int len_educt_types,  int* product_coefficients, int len_product_types);

int update_type_index(int* educt_types, int len_educt_types , int* product_types, int len_product_types); //assign different types an index in a growing list that starts at 0 and is incremented by 1 for each new type. the entry in the index at place i is the "type_value". therefore the type of type "typevalue" has the index i; 

int generic_oneway_reaction(int reaction_id);

int find_index_of_type(int type);

#endif /* ifdef REACTION_H */
