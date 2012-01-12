cdef extern from "communication.h":
  void mpi_set_time_step(double time_step)
  void mpi_bcast_parameter(int p)

cdef extern from "global.h":
  int FIELD_TEMPERATURE
