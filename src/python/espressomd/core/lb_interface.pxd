from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc cimport stdint

from core.utils cimport Vector3d
from core.utils cimport Vector3i
from core.utils cimport Vector6d
from core.utils cimport Vector19d

cdef extern from "grid_based_algorithms/lb_interface.hpp" namespace "ActiveLB":
    cdef ActiveLB NONE
    cdef ActiveLB CPU
    cdef ActiveLB GPU

cdef extern from "grid_based_algorithms/lb_interface.hpp":

    cdef enum ActiveLB:
        pass
    void lb_lbfluid_set_tau(double c_tau) except +
    double lb_lbfluid_get_tau() except +
    void lb_lbfluid_set_density(double c_dens) except +
    double lb_lbfluid_get_density() except +
    void lb_lbfluid_set_viscosity(double c_visc) except +
    double lb_lbfluid_get_viscosity() except +
    void lb_lbfluid_set_agrid(double c_agrid) except +
    double lb_lbfluid_get_agrid() except +
    void lb_lbfluid_set_gamma_odd(double c_gamma_odd) except +
    double lb_lbfluid_get_gamma_odd() except +
    void lb_lbfluid_set_gamma_even(double c_gamma_even) except +
    double lb_lbfluid_get_gamma_even() except +
    void lb_lbfluid_set_ext_force_density(const Vector3d forcedensity) except +
    const Vector3d lb_lbfluid_get_ext_force_density() except +
    void lb_lbfluid_set_bulk_viscosity(double c_bulk_visc) except +
    double lb_lbfluid_get_bulk_viscosity() except +
    void lb_lbfluid_print_vtk_velocity(string filename) except +
    void lb_lbfluid_print_vtk_velocity(string filename, vector[int] bb1, vector[int] bb2) except +
    void lb_lbfluid_print_vtk_boundary(string filename) except +
    void lb_lbfluid_print_velocity(string filename) except +
    void lb_lbfluid_print_boundary(string filename) except +
    void lb_lbfluid_save_checkpoint(string filename, int binary) except +
    void lb_lbfluid_load_checkpoint(string filename, int binary) except +
    void lb_lbfluid_set_lattice_switch(ActiveLB local_lattice_switch) except +
    ActiveLB lb_lbfluid_get_lattice_switch() except +
    Vector6d lb_lbfluid_get_stress() except +
    bool lb_lbnode_is_index_valid(Vector3i & ind) except +
    Vector3i lb_lbfluid_get_shape() except +
    const Vector3d lb_lbnode_get_velocity(const Vector3i & ind) except +
    void lb_lbnode_set_velocity(const Vector3i & ind, const Vector3d & u) except +
    double lb_lbnode_get_density(const Vector3i & ind) except +
    void lb_lbnode_set_density(const Vector3i & ind, double density) except +
    const Vector6d lb_lbnode_get_stress(const Vector3i & ind) except +
    const Vector6d lb_lbnode_get_stress_neq(const Vector3i & ind) except +
    const Vector19d lb_lbnode_get_pop(const Vector3i & ind) except +
    void lb_lbnode_set_pop(const Vector3i & ind, const Vector19d & populations) except +
    int lb_lbnode_get_boundary(const Vector3i & ind) except +
    stdint.uint64_t lb_lbfluid_get_rng_state() except +
    void lb_lbfluid_set_rng_state(stdint.uint64_t) except +
    void lb_lbfluid_set_kT(double) except +
    double lb_lbfluid_get_kT() except +
    double lb_lbfluid_get_lattice_speed() except +

