include "myconfig.pxi"

cimport global_variables
import numpy as np

cdef class GlobalsHandle:
    def __init__(self):
        pass
  
    property box_l:
        def __set__(self, _box_l):
            global box_l
            if len(_box_l) != 3:
                raise ValueError("Box length must be of length 3")
            for i in range(3):
                if _box_l[i] <= 0:
                    raise ValueError("Box length must be > 0 in all directions")
                box_l[i]=_box_l[i]
      
            mpi_bcast_parameter(0)
    
        def __get__(self):
            return np.array([box_l[0], box_l[1], box_l[2]])
  
    property cell_grid:
        def __get__(self):
            return np.array( [ dd.cell_grid[0], dd.cell_grid[1], dd.cell_grid[2]  ] )
  
    property cell_size:
        def __get__(self):
            return np.array( [ dd.cell_size[0], dd.cell_size[1], dd.cell_size[2] ] )
    
    property dpd_gamma:
        def __get__(self):
            return dpd_gamma;
  
    property dpd_r_cut:
        def __get__(self):
            return dpd_r_cut;

    property gamma:
        def __get__(self):
            return langevin_gamma;
  
    property integ_switch:
        def __get__(self):
            return integ_switch;
      
    property local_box_l:
        def __get__(self):
            return np.array([local_box_l[0], local_box_l[1], local_box_l[2]])
  
    property max_cut:
        def __get__(self):
            return max_cut;
      
    property max_num_cells:
        def __set__(self, int _max_num_cells):
            global max_num_cells
            if _max_num_cells < min_num_cells:
                raise ValueError("max_num_cells must be >= min_num_cells (currently "+str(min_num_cells)+")")
            max_num_cells=_max_num_cells
            mpi_bcast_parameter(FIELD_MAXNUMCELLS);
        def __get__(self):
            return max_num_cells;

    property min_num_cells:
        def __set__(self, int _min_num_cells): 
            global min_num_cells 
            min = calc_processor_min_num_cells()
            if _min_num_cells < min:
                raise ValueError("min_num_cells must be >= processor_min_num_cells (currently "+str(min)+")")
            if _min_num_cells > max_num_cells:
                raise ValueError("min_num_cells must be <= max_num_cells (currently "+str(max_num_cells)+")")
            min_num_cells=_min_num_cells
            mpi_bcast_parameter(FIELD_MINNUMCELLS);
        def __get__(self):
            return min_num_cells;
    
    property max_part:
        def __get__(self):
            return max_seen_particle;
  
    property max_range:
        def __get__(self):
            return max_range;
  
    property max_skin:
        def __get__(self):
            return max_skin;
  
    property n_layers:
        def __get__(self):
            return n_layers;
  
    property n_nodes:
        def __get__(self):
            return n_nodes;
  
    property n_part:
        def __get__(self):
            return n_part;
  
    property n_part_types:
        def __get__(self):
            return n_particle_types;
  
    property n_rigidbonds:
        def __get__(self):
            return n_rigidbonds;
  
    property node_grid:
        def __set__(self, _node_grid):
            global node_grid
            if len(_node_grid) != 3:
                raise ValueError("node_grid must be of length 3")
            for i in range(3):
                if _node_grid[i] <= 0:
                    raise ValueError("node_grid must be > 0 in all directions")
                node_grid[i]=_node_grid[i]
            if _node_grid[0]*_node_grid[1]*_node_grid[2] != n_nodes:
                raise ValueError("node_grid does not fit n_nodes ("+str(n_nodes)+")")
            for i in range(3):
                node_grid[i]=_node_grid[i]
            mpi_bcast_parameter(FIELD_NODEGRID)
        def __get__(self):
            return np.array([node_grid[0], node_grid[1], node_grid[2]])
  
    property nptiso_gamma0:
        def __get__(self):
            return nptiso_gamma0;
  
    property nptiso_gammav:
        def __get__(self):
            return nptiso_gammav;
  
    property npt_p_ext:
        def __get__(self):
            return nptiso.p_ext;
  
    property npt_p_inst:
        def __get__(self):
            return nptiso.p_inst;
  
    property npt_p_inst_av:
        def __get__(self):
            return nptiso.p_inst_av;
  
    property npt_piston:
        def __set__(self, _npt_piston):
            global npt_piston
            if _npt_piston < 0:
                raise ValueError("npt_piston must be > 0")
            nptiso.piston=_npt_piston
            mpi_bcast_parameter(FIELD_NPTISO_PISTON)
        def __get__(self):
            global npt_piston
            return nptiso.piston

    property npt_p_diff:
        def __set__(self, _npt_p_diff):
            global npt_p_diff
            nptiso.p_diff=_npt_p_diff
            mpi_bcast_parameter(FIELD_NPTISO_PDIFF)
        def __get__(self):
            global npt_p_diff
            return nptiso.p_diff

    property periodicity:
        def __set__(self, _periodic):
            global periodic
            if len(_periodic) != 3:
                raise ValueError("periodicity must be of length 3, got length "+str(len(_periodic)))
            periodicity=np.zeros(3);
            for i in range(3):
                if _periodic[i] != 1:
                    raise ValueError("Until we can handle conditional compilation, only periodicity [1,1,1] is supported in python interface")
            for i in range(3):
                periodicity[i]=_periodic[i];
            periodic=4*_periodic[2] + 2*_periodic[1] + _periodic[0]; 
            # first 3 bits of periodic determine the periodicity
            # until we can handle contitional compilatio, periodic=7 is the only value which makes sense
            mpi_bcast_parameter(FIELD_PERIODIC);
        def __get__(self):
            global periodic 
            periodicity=np.zeros(3);
            periodicity[0]=periodic%2; 
            periodicity[1]=int(periodic/2)%2; 
            periodicity[2]=int(periodic/4)%2; 
            return periodicity

    property skin:
        def __set__(self, double _skin):
            if _skin <= 0:
                raise ValueError("Skin must be >= 0")
            global skin
            skin=_skin
            mpi_bcast_parameter(28)
        def __get__(self): 
            global skin
            return skin

    property temperature:
        def __get__(self):
            global temperature;
            return temperature;
  
    property thermo_switch:
        def __get__(self):
            global thermo_switch;
            return thermo_switch;
  
    property time:
        def __set__(self, double _time):
            if _time <= 0:
                raise ValueError("Simulation time must be >= 0")
            global sim_time
            sim_time=_time
            mpi_bcast_parameter(FIELD_SIMTIME)
        def __get__(self): 
            global sim_time
            return sim_time

    property time_step:
        def __set__(self, double _time_step):
            IF LB:
                global lbpar
            IF LB_GPU:
                global lbpar_gpu
            if _time_step <= 0:
              raise ValueError("Time Step must be positive")
            IF LB:
                if lbpar.tau >= 0.0 and _time_step > lbpar.tau:
                  raise ValueError("Time Step must be > LB_time_step ("+str(lbpar.tau)+")")
            IF LB_GPU:
                if lbpar_gpu.tau >= 0.0 and _time_step > lbpar_gpu.tau:
                  raise ValueError("Time Step must be > LB_time_step ("+str(lbpar_gpu.tau)+")")
            mpi_set_time_step(_time_step)
        def __get__(self):
            global time_step
            return time_step
  
    property timings:
        def __set__(self, int _timings):
            global timing_samples;
            if _timings <= 0:
                timing_samples=0;
            else:
                timing_samples=_timings;
        def __get__(self):
            global timing_samples
            return timing_samples
  
    property transfer_rate:
        def __get__(self):
            global transfer_rate;
            return transfer_rate;
  
    property max_cut_nonbonded:
        def __get__(self):
            global max_cut_nonbonded;
            return max_cut_nonbonded;
  
    property verlet_reuse:
        def __get__(self):
            global verlet_reuse;
            return verlet_reuse;


    property lattice_switch:
        def __get__(self):
            global lattice_switch;
            return lattice_switch;
  
    property dpd_tgamma:
        def __get__(self):
            global dpd_tgamma;
            return dpd_tgamma;
  
    property dpd_tr_cut:
        def __get__(self):
            global dpd_tr_cut;
            return dpd_tr_cut;
  
    property dpd_twf:
        def __get__(self):
            global dpd_twf;
            return dpd_twf;
  
    property dpd_wf:
        def __get__(self):
            global dpd_wf;
            return dpd_wf;

    property adress_vars:
        def __get__(self):
            global adress_vars;
            return np.array( [ \
            adress_vars[0], \
            adress_vars[1], \
            adress_vars[2], \
            adress_vars[3], \
            adress_vars[4], \
            adress_vars[5], \
            adress_vars[6] \
            ])
  
    property max_cut_bonded:
        def __get__(self):
            global max_cut_bonded;
            return max_cut_bonded;
  
