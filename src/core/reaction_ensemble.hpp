#include "utils.hpp"
namespace ReactionEnsemble{

typedef struct single_reaction{
	//strict input to the algorithm
	int* reactant_types;
	int len_reactant_types;
	int* reactant_coefficients;
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
	int* type_index;
	int nr_different_types; // is equal to length type_index
	double* charges_of_types;
	double standard_pressure_in_simulation_units;
	double temperature_reaction_ensemble;
	double exclusion_radius; //this is used as a kind of hard sphere radius, if particles are closer than that it is assumed that their interaction energy gets approximately infinite => these configurations do not contribute to the partition function and ensemble averages.
	double volume;
	bool box_is_cylindric_around_z_axis;
	double cyl_radius;
	double cyl_x;
	double cyl_y;
	bool box_has_wall_constraints;
	double slab_start_z;
	double slab_end_z;
	int non_interacting_type;
} reaction_system;

extern reaction_system current_reaction_system;

extern int accepted_configurational_MC_moves;
extern int tried_configurational_MC_moves;

void set_cuboid_reaction_ensemble_volume();

int do_reaction();

int free_reaction_ensemble();

int check_reaction_ensemble();

int calculate_nu_bar(int* reactant_coefficients, int len_reactant_types,  int* product_coefficients, int len_product_types);

int update_type_index(int* reactant_types, int len_reactant_types , int* product_types, int len_product_types); //assign different types an index in a growing list that starts at 0 and is incremented by 1 for each new type. the entry in the index at place i is the "type_value". therefore the type of type "typevalue" has the index i; 

int generic_oneway_reaction(int reaction_id, int reaction_modus);

int find_index_of_type(int type);

struct invalid_index : public std::exception {
    const char* what () const throw ()
    {
	return "Type not indexed!";
    }
};

struct no_particle_inserted : public std::exception {
    const char* what () const throw ()
    {
	return "Error: Particle not inserted. System too dense?";
    }
};

struct reaction_mode_unknown : public std::exception {
    const char* what () const throw ()
    {
	return "Error: This reaction mode is unknown";
    }
};

bool do_global_mc_move_for_particles_of_type(int type, int start_id_polymer, int end_id_polymer, int particle_number_of_type, const bool use_wang_landau);

///////////////////////////////////////////// Wang-Landau algorithm

typedef struct collective_variable{
	double CV_minimum;
	double CV_maximum;
	double delta_CV;
	double (*determine_current_state_in_collective_variable_with_index) (int);	//declare a function pointer with name determine_current_state_in_this_collective_variable that has to point to a function that has 1 int input parameter (for the collective_variable_index) and returns a double
	//for collective variables of type degree of association
	int* corresponding_acid_types; // is NULL if the current collective variable has nothing to do with a degree of association
	int nr_corresponding_acid_types;
	int associated_type;
	//for collective variables of type energy
	char* energy_boundaries_filename;
} collective_variable;

typedef struct wang_landau_system {
	int* histogram;
	int len_histogram; //equals also len_wang_landau_potential and also len_Gamma
	double* wang_landau_potential; //equals the logarithm to basis e of the degeneracy of the states
	int nr_collective_variables;
	collective_variable** collective_variables;
	int* nr_subindices_of_collective_variable;
	double wang_landau_parameter; //equals the logarithm to basis e of the modification factor of the degeneracy of states when the state is visited
	double initial_wang_landau_parameter;
	
	int int_fill_value;
	double double_fill_value;
	
	int number_of_monte_carlo_moves_between_check_of_convergence;
	double final_wang_landau_parameter;
	int used_bins; //for 1/t algorithm
	int monte_carlo_trial_moves; //for 1/t algorithm

	int wang_landau_steps; //may be used for performance improvements, when you do not want to record other observables in the tcl script
	char* output_filename;
	
	double* minimum_energies_at_flat_index; //only present in energy preparation run
	double* maximum_energies_at_flat_index; //only present in energy preparation run
	
	bool do_energy_reweighting;
	int polymer_start_id;
	int polymer_end_id;
	bool fix_polymer;
	bool do_not_sample_reaction_partition_function;
	bool use_hybrid_monte_carlo;
} wang_landau_system;

extern wang_landau_system current_wang_landau_system;

int initialize_wang_landau(); //may first be called after all collective variables are added
int do_reaction_wang_landau();
void free_wang_landau();
int update_maximum_and_minimum_energies_at_current_state(); //use for preliminary energy reweighting runs
void write_out_preliminary_energy_run_results(char* filename);
bool do_HMC_move_wang_landau();

//checkpointing, only designed to reassign values of a previous simulation to a new simulation with the same initialization process
int write_wang_landau_checkpoint(char* identifier);
int load_wang_landau_checkpoint(char* identifier);








/////////////////////////////////////////////  Constant-pH Reactions
int do_reaction_constant_pH();
double extern constant_pH;
void set_pH(double pH);
}
