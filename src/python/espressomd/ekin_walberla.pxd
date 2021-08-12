include "myconfig.pxi"

from libcpp cimport bool
from libcpp.string cimport string

from .actors cimport Actor
from .utils cimport Vector3i

cdef class EKinWalberla(Actor):
    pass

cdef class EKinRoutines:
    cdef Vector3i node

IF EK_WALBERLA:
    # TODO: figure out if this extern is even necessary
    cdef extern from "grid_based_algorithms/ek_interface.hpp":
        cdef enum OutputVTK:
            pass

    cdef extern from "grid_based_algorithms/ek_interface.hpp" namespace 'EKOutputVTK':
        cdef OutputVTK output_vtk_density 'EKOutputVTK::density'

    cdef extern from "grid_based_algorithms/ekin_walberla_instance.hpp":
        void mpi_init_ekin_walberla(double diffusion, double kT, double density, double tau) except +
        void mpi_destruct_ekin_walberla() except +

cdef extern from "grid_based_algorithms/ek_interface.hpp" namespace "EK":
    double get_diffusion(unsigned int id) except +
    void set_diffusion(unsigned int id, double diffusion) except +

    double get_kT(unsigned int id) except +
    void set_kT(unsigned int id, double kT) except +

    double get_tau() except +

    Vector3i get_shape(unsigned int id) except +

    bool node_is_index_valid(unsigned int id, const Vector3i & ind) except +

    double get_density(unsigned int id, const Vector3i & ind) except +
    void set_density(unsigned int id, const Vector3i & ind, double density) except +

    bool get_is_boundary(unsigned int id, const Vector3i & ind) except +

    void create_vtk(unsigned int id, unsigned delta_N, unsigned initial_count, unsigned flag_observables, const string identifier, const string base_folder, const string execution_folder) except +
    void switch_vtk(unsigned int id, const string vtk_uid, int status) except +
    void write_vtk(unsigned int id, const string vtk_uid) except +
