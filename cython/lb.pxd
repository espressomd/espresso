cdef extern from "config.h":
  pass

cdef extern from "lb.h":
  int lb_lbfluid_set_tau(double p_tau)
  int lb_lbfluid_set_density(double p_dens)
  int lb_lbfluid_get_density(double* p_dens)
  int lb_lbfluid_set_visc(double p_visc)
  int lb_lbfluid_get_visc(double* p_visc)
  int lb_lbfluid_set_agrid(double agrid)
  int lb_lbfluid_get_agrid(double* p_agrid)
  int lb_lbfluid_set_friction(double friction)
  int lb_lbfluid_get_friction(double* p_friction)
  void cython_lb_init(int switch)
