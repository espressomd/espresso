from __future__ import print_function, absolute_import
cdef extern from "particle_data.hpp":
    int init_type_array(int type)
    int find_particle_type(int type, int * id)
    int number_of_particles_with_type(int type, int * number)
    
    ctypedef struct TypeList:
        int max_entry
        int cur_size
        int *id_list
    
    cdef extern TypeList *type_array;
    
    
    ctypedef struct TypeOfIndex:
        int max_entry
        int *index
    cdef extern TypeOfIndex Type

    ctypedef struct IndexOfType:
        int max_entry
        int *type
    
    cdef extern IndexOfType Index
