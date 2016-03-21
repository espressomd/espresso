cdef extern from "particle_data.hpp":
	int init_type_array(int type)
	int delete_particle_of_type(int type)
	int find_particle_type(int type, int *id)
	int number_of_particles_with_type(int type, int *number)
