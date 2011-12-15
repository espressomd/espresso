

cdef extern from "../../src/communication.h":
  void mpi_set_time_step(double time_step)
  void mpi_bcast_parameter(int p)

cdef extern from "../../src/integrate.h":
  double time_step

cdef extern from "../../src/verlet.h":
  double skin

cdef extern from "../../src/grid.h":
  double box_l[3]

