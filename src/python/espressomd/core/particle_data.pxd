include "myconfig.pxi"

IF EXCLUSIONS:
    cdef extern from "particle_data.hpp":
        void auto_exclusions(int distance)

cdef extern from "particle_data.hpp":
    int max_seen_particle_type
    int init_type_map(int type) except +
    int get_random_p_id(int type) except +
    int number_of_particles_with_type(int type) except +
