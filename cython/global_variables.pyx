cimport global_variables
import numpy as np

cpdef int marcello_test(int x):
    return x**2

cdef class GlobalsHandle:
    def __init__(self):
        pass
  
    property time_step:
        def __set__(self, double _time_step):
            if _time_step <= 0:
              raise ValueError("Time Step must be positive")
            mpi_set_time_step(_time_step)
        def __get__(self):
            global time_step
            return time_step
  
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
	    #FIXME add check if max_num_cells > min_num_cells
            if _max_num_cells <= 0:
                raise ValueError("max_num_cells must be >= 0")
            max_num_cells=_max_num_cells
            mpi_bcast_parameter(FIELD_MAXNUMCELLS);
        def __get__(self):
            return max_num_cells;
      
  
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

