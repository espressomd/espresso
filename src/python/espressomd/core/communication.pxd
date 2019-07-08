include "myconfig.pxi"
from libcpp.vector cimport vector

from core.utils cimport Vector3d

cdef extern from "communication.hpp":
    void mpi_bcast_cell_structure(int cs)
    int n_nodes
    vector[int] mpi_resort_particles(int global_flag)
    void mpi_kill_particle_motion(int rotation)
    void mpi_kill_particle_forces(int torque)
    Vector3d mpi_system_CMS()
    Vector3d mpi_system_CMS_velocity()
    void mpi_galilei_transform()
    void mpi_init()
    void mpi_loop()
    int this_node

cdef extern from "communication.hpp" namespace "Random":
    void mpi_random_seed(int cnt, vector[int] & seed)


IF ELECTROSTATICS and P3M:
    cdef extern from "communication.hpp":
        int mpi_iccp3m_init()

IF ELECTROSTATICS:
    cdef extern from "communication.hpp":
        void mpi_bcast_coulomb_params()
