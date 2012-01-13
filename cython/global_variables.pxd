cdef extern from "global.h":
    int FIELD_MAXNUMCELLS
    int FIELD_MINNUMCELLS
    int FIELD_NODEGRID
    int FIELD_NPTISO_PISTON
    int FIELD_NPTISO_PDIFF
    int FIELD_PERIODIC
    int FIELD_SIMTIME

cdef extern from "communication.h":
    extern int n_nodes
    void mpi_set_time_step(double time_step)
    void mpi_bcast_parameter(int p)

cdef extern from "integrate.h":
    double time_step
    extern int integ_switch
    extern double sim_time
    extern double verlet_reuse


cdef extern from "verlet.h":
    double skin

cdef extern from "lattice.h":
    extern int lattice_switch

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
    extern int min_num_cells
    extern double max_skin
    int calc_processor_min_num_cells()

  
cdef extern from "particle_data.h":
    extern int  n_total_particles


cdef extern from "interaction_data.h":
    double dpd_gamma
    double dpd_r_cut
    extern double max_cut
    extern int max_seen_particle
    extern int n_particle_types
    extern double max_cut_nonbonded
    extern double max_cut_bonded


cdef extern from "thermostat.h":
    double langevin_gamma
    extern double nptiso_gamma0
    extern double nptiso_gammav
    extern double temperature 
    extern int thermo_switch     


cdef extern from "dpd.h":
    extern int dpd_wf
    extern double dpd_tgamma
    extern double dpd_tr_cut
    extern int dpd_twf


#FIXME this will only work with conditional compilation of address
#cdef extern from "adresso.h":
#    extern double adress_vars[7]


cdef extern from "cells.h":
    extern double max_range

cdef extern from "layered.h":
    extern int n_layers

cdef extern from "rattle.h":
    extern int n_rigidbonds


cdef extern from "tuning.h":
    extern int timing_samples

cdef extern from "imd.h":
    extern int transfer_rate


cdef extern from "grid.h":
    double box_l[3]
    double local_box_l[3]
    extern int node_grid[3]
    extern int periodic

cdef extern from "npt.h":
    ctypedef struct nptiso_struct:
        double p_ext
        double p_inst
        double p_inst_av
        double p_diff
        double piston
    extern nptiso_struct nptiso

