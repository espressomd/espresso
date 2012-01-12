#cimport global_variables
import numpy as np
cimport lb

cdef extern from "../myconfig.h":
  pass
cdef extern from "../src/lattice.h":
  int lattice_switch
#don t work within header lb.pxd
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
  void cython_lb_init(int switch)

cdef class LBparaHandle:
  cdef int switch
  def __init__(self, _dev):
    if _dev == "gpu" :
      switch=1
      cython_lb_init(switch)
    else: 
      switch=0
      cython_lb_init(switch)

  property tau:
    def __set__(self, double _tau):
      lb_lbfluid_set_tau(_tau)
    def __get__(self):
      raise Exception("tau can not be printed")

  property dens:
    def __set__(self, double _dens):
      lb_lbfluid_set_density(_dens)
    def __get__(self): 
      cdef double _p_dens
      return lb_lbfluid_get_density(&_p_dens)

  property visc:
    def __set__(self, _visc):
      lb_lbfluid_set_visc(_visc)
    def __get__(self):
      cdef double _p_visc
      return lb_lbfluid_get_visc(&_p_visc)

  property agrid:
    def __set__(self, double _agrid):
      lb_lbfluid_set_agrid(_agrid)
    def __get__(self):
      cdef double _p_agrid
      return lb_lbfluid_get_agrid(&_p_agrid)
      
  property friction:
    def __set__(self, double _friction):
      lb_lbfluid_set_friction(_friction)
    def __get__(self):
      cdef double _p_friction
      raise Exception("get friction c function not implemented in lb.c")
      #return lb_lbfluid_get_friction(&_p_friction)

class DeviceList:
  def __getitem__(self, _dev):
    return LBparaHandle(_dev)
    
  
