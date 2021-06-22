
from libcpp cimport bool

from .actors cimport Actor
from .utils cimport Vector3i

cdef class EKinWalberla(Actor):
    pass

cdef class EKinRoutines:
    cdef Vector3i node

cdef extern from "grid_based_algorithms/ekin_walberla_instance.hpp":
    void mpi_init_ekin_walberla(double diffusion, double kT, double density) except +
    void mpi_destruct_ekin_walberla() except +

cdef extern from "grid_based_algorithms/ekin_walberla_interface.hpp" namespace "walberla":
    double ek_get_diffusion() except +
    void ek_set_diffusion(double diffusion) except +

    double ek_get_kT() except +
    void ek_set_kT(double kT) except +

    Vector3i ek_get_shape() except +

    bool ek_node_is_index_valid(const Vector3i & ind) except +

    double ek_get_density(const Vector3i & ind) except +
    void ek_set_node_density(const Vector3i & ind, double density) except +
