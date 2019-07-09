include "myconfig.pxi"
from libcpp cimport bool

IF EXCLUSIONS:
    cdef extern from "particle_data.hpp":
        void auto_exclusions(int distance)

cdef extern from "particle_data.hpp":
    int max_seen_particle_type
    int init_type_map(int type) except +
    int get_random_p_id(int type) except +
    int number_of_particles_with_type(int type) except +
    extern int n_part
    extern bool swimming_particles_exist
