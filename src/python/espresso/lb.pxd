include "myconfig.pxi" 



IF LB_GPU==1:
  
  
  cdef extern from "../src/config.hpp":
    pass
  
  #cdef extern from "../src/lattice.hpp":
  #  int lattice_switch
  
  cdef extern from "lb.hpp":
    int lb_lbfluid_set_tau(double p_tau)
    int lb_lbfluid_set_density(double p_dens)
    int lb_lbfluid_get_density(double* p_dens)
    int lb_lbfluid_set_visc(double p_visc)
    int lb_lbfluid_get_visc(double* p_visc)
    int lb_lbfluid_set_agrid(double p_agrid)
    int lb_lbfluid_get_agrid(double* p_agrid)
    int lb_lbfluid_set_friction(double friction)
    int lb_lbfluid_get_friction(double* p_friction)
    int lb_lbfluid_set_gamma_odd(double p_gamma_odd)
    int lb_lbfluid_get_gamma_odd(double* p_gamma_odd)
    int lb_lbfluid_set_gamma_even(double p_gamma_even)
    int lb_lbfluid_get_gamma_even(double* p_gamma_even)
    int lb_lbfluid_set_ext_force(double p_fx, double p_fy, double p_fz)
    int lb_lbfluid_get_ext_force(double* p_fx, double* p_fy, double* p_fz)
    int lb_lbfluid_set_bulk_visc(double p_bulk_visc)
    int lb_lbfluid_get_bulk_visc(double* p_bulk_visc)
    int lb_lbfluid_print_vtk_velocity(char* filename)
    int lb_lbfluid_print_vtk_boundary(char* filename)
    int lb_lbfluid_print_velocity(char* filename)
    int lb_lbfluid_print_boundary(char* filename)
    int lb_lbfluid_save_checkpoint(char* filename, int binary)
    int lb_lbfluid_load_checkpoint(char* filename, int binary)
    void cython_lb_init(int switch)
