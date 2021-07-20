include "myconfig.pxi"

from libcpp cimport bool

from .actors cimport Actor
from .utils cimport Vector3i

cdef class EKinWalberla(Actor):
    pass

cdef class EKinRoutines:
    cdef Vector3i node

IF EK_WALBERLA:
    cdef extern from "grid_based_algorithms/ekin_walberla_instance.hpp":
        void mpi_init_ekin_walberla(double diffusion, double kT, double density, double tau) except +
        void mpi_destruct_ekin_walberla() except +

cdef extern from "grid_based_algorithms/ek_interface.hpp" namespace "EK":
    double get_diffusion() except +
    void set_diffusion(double diffusion) except +

    double get_kT() except +
    void set_kT(double kT) except +

    double get_tau() except +

    Vector3i get_shape() except +

    bool node_is_index_valid(const Vector3i & ind) except +

    double get_density(const Vector3i & ind) except +
    void set_density(const Vector3i & ind, double density) except +

    bool get_is_boundary(const Vector3i & ind) except +
