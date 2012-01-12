cdef extern from "global.h":
    int FIELD_MAXNUMCELLS

cdef extern from "communication.h":
    void mpi_set_time_step(double time_step)
    void mpi_bcast_parameter(int p)

cdef extern from "integrate.h":
    double time_step
    extern int integ_switch

cdef extern from "verlet.h":
    double skin

cdef extern from "../src/domain_decomposition.h":
    ctypedef struct IA_Neighbor:
        pass
    ctypedef struct IA_Neighbor_List:
        pass
    ctypedef struct  DomainDecomposition:
          int cell_grid[3]
          double cell_size[3]
    
    extern DomainDecomposition dd
    extern int max_num_cells

  
cdef extern from "interaction_data.h":
    double dpd_gamma
    double dpd_r_cut
    extern double max_cut

cdef extern from "thermostat.h":
    double langevin_gamma

cdef extern from "grid.h":
    double box_l[3]
    double local_box_l[3]

