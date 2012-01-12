
cdef extern from "../myconfig.h":
  pass

cdef extern from "../src/lb.h":
  int lb_lbfluid_set_tau(double _tau)
  int lb_lbfluid_set_density(double _dens)
  int lb_lbfluid_get_density(double* _p_dens)
  int lb_lbfluid_set_visc(double _visc)
  int lb_lbfluid_get_visc(double* _p_visc)
  int lb_lbfluid_set_agrid(double _agrid)
  int lb_lbfluid_get_agrid(double* _p_agrid)
  int lb_lbfluid_set_friction(double _friction)
  int lb_lbfluid_get_friction(double* _p_friction)

