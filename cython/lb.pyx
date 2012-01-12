#cimport global_variables
import numpy as np
cimport lb

#cdef extern from "../src/lattice.h":
#  int lattice_switch
#don t work within header lb.pxd no idea why
cdef extern from "../src/lb.h":
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

cdef class LBparaHandle:
  cdef int switch
  cdef int checkpoint_binary
  cdef char* checkpoint_filename
  
  def __init__(self, _dev):
    if _dev == "gpu" :
      switch=1
      cython_lb_init(switch)
    else: 
      switch=0
      cython_lb_init(switch)

  property tau:
    def __set__(self, double p_tau):
      if lb_lbfluid_set_tau(p_tau):
        raise Exception("lb_lbfluid_set_tau error")
    def __get__(self):
      raise Exception("get friction c function not implemented in lb.c")

  property dens:
    def __set__(self, double p_dens):
      if lb_lbfluid_set_density(p_dens):
        raise Exception("lb_lbfluid_set_density error")
    def __get__(self): 
      cdef double _p_dens
      if lb_lbfluid_get_density(&_p_dens):
        raise Exception("lb_lbfluid_get_density error")
      return _p_dens

  property visc:
    def __set__(self, _visc):
      if lb_lbfluid_set_visc(_visc):
        raise Exception("lb_lbfluid_set_visc error")
    def __get__(self):
      cdef double _p_visc
      if lb_lbfluid_get_visc(&_p_visc):
        raise Exception("lb_lbfluid_get_visc error")
      return _p_visc
      
  property agrid:
    def __set__(self, double _agrid):
      if lb_lbfluid_set_agrid(_agrid):
        raise Exception("lb_lbfluid_set_agrid error")
    def __get__(self):
      cdef double _p_agrid
      if lb_lbfluid_get_agrid(&_p_agrid):
        raise Exception("lb_lbfluid_get_agrid error")
      return _p_agrid
      
  property friction:
    def __set__(self, double _friction):
      if lb_lbfluid_set_friction(_friction):
        raise Exception("lb_lbfluid_set_friction error")
    def __get__(self):
      cdef double _p_friction
      raise Exception("get friction c function not implemented in lb.c")
      #return lb_lbfluid_get_friction(&_p_friction)

  property gamma_odd:
    def __set__(self, double _gamma_odd):
      if lb_lbfluid_set_gamma_odd(_gamma_odd):
        raise Exception("lb_lbfluid_set_gamma_odd error")
    def __get__(self):
      cdef double _p_gamma_odd
      if lb_lbfluid_get_gamma_odd(&_p_gamma_odd):
        raise Exception("lb_lbfluid_get_gamma_odd error")
      return _p_gamma_odd
      
  property gamma_even:
    def __set__(self, double _gamma_even):
      if lb_lbfluid_set_gamma_even(_gamma_even):
        raise Exception("lb_lbfluid_set_gamma_even error")
    def __get__(self):
      cdef double _p_gamma_even
      if lb_lbfluid_get_gamma_even(&_p_gamma_even):
        raise Exception("lb_lbfluid_get_gamma_even error")
      return _p_gamma_even
      
  property ext_force:
    def __set__(self, _ext_force):
      if lb_lbfluid_set_ext_force(_ext_force[0], _ext_force[1], _ext_force[2]):
        raise Exception("lb_lbfluid_set_ext_force error")
    def __get__(self):
      cdef double _p_ext_force[3]
      if lb_lbfluid_get_ext_force(&_p_ext_force[0], &_p_ext_force[1], &_p_ext_force[2]):
        raise Exception("lb_lbfluid_get_ext_force error")
      return np.array([_p_ext_force[0], _p_ext_force[1], _p_ext_force[2]])
      
  property bulk_visc:
    def __set__(self, double _bulk_visc):
      if lb_lbfluid_set_bulk_visc(_bulk_visc):
        raise Exception("lb_lbfluid_set_bulk_visc error")
    def __get__(self):
      cdef double _p_bulk_visc
      if lb_lbfluid_get_bulk_visc(&_p_bulk_visc):
        raise Exception("lb_lbfluid_get_bulk_visc error")
      return _p_bulk_visc      
  
  property print_vtk_velocity:
    def __set__(self, char* _filename):
      if lb_lbfluid_print_vtk_velocity(_filename):
        raise Exception("lb_lbfluid_print_vtk_velocity error")
            
  property print_vtk_boundary:
    def __set__(self, char* _filename):
      if lb_lbfluid_print_vtk_boundary(_filename):
        raise Exception("lb_lbfluid_print_vtk_boundary error")   

  property print_velocity:
    def __set__(self, char* _filename):
      if lb_lbfluid_print_velocity(_filename):
        raise Exception("lb_lbfluid_print_vtk_velocity error")
            
  property print_boundary:
    def __set__(self, char* _filename):
      if lb_lbfluid_print_boundary(_filename):
        raise Exception("lb_lbfluid_print_vtk_boundary error")
        
  property checkpoint:
    def __set__(self, char* checkpoint_filename):
      self.checkpoint_filename=checkpoint_filename
      if lb_lbfluid_save_checkpoint(checkpoint_filename, self.checkpoint_binary):
        raise Exception("lb_lbfluid_save_checkpoint error")
    def __get__(self):
      if lb_lbfluid_load_checkpoint(self.checkpoint_filename, self.checkpoint_binary):
        raise Exception("lb_lbfluid_load_checkpoint error")  

  property checkpoint_style:
    def __set__(self, int _binary):
      self.checkpoint_binary=_binary
    def __get__(self):
      return self.checkpoint_binary
            
class DeviceList:
  def __getitem__(self, _dev):
    return LBparaHandle(_dev)
    
  
