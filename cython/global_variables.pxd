cdef extern from "communication.h":
  void mpi_set_time_step(double time_step)
  void mpi_bcast_parameter(int p)

cdef extern from "integrate.h":
  double time_step

cdef extern from "verlet.h":
  double skin

cdef extern from "../src/domain_decomposition.h":
  ctypedef struct IA_Neighbor:
    pass
  ctypedef struct IA_Neighbor_List:
    pass
  ctypedef struct  DomainDecomposition:
      int cell_grid[3]
  
  extern DomainDecomposition dd

cdef extern from "grid.h":
  double box_l[3]

