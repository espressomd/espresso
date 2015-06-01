#
# Copyright (C) 2013,2014 The ESPResSo project
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#  
include "myconfig.pxi"
import numpy as np
from actors import Actor
#cimport cuda_init
#import cuda_init
#cimport lb
#cimport global_variables

class HydrodynamicInteraction(Actor):
  def _lb_init(self):
    raise Exception("Subclasses of HydrodynamicInteraction must define the _lb_init() method.")

IF LB_GPU or LB:
  class LB(HydrodynamicInteraction):

    #int switch
    #char* checkpoint_filename

    def __init__(self, _dev):
      if _dev == "gpu":
        switch=1
        #cython_lb_init(switch)
      else: 
        switch=0
        #cython_lb_init(switch)

    def validateParams(self):
      default_lb_params=self.defaultParams()

      if not (self._lb_params["agrid"] > 0.0):
        raise ValueError("Grid spacing should be a positive double")

    def validKeys(self):
      return "gpu", "agrid","dens","fric","ext_force","visc","tau"

    def requiredKeys(self):
      return ["agrid","visc", "tau"]

    def defaultParams(self):
      return {"gpu":0,\
              "agrid":-1,\
              "dens":[-1,-1],\
              "fric":[-1,-1],\
              "ext_force":[-1,-1,-1],\
              "visc":[-1,-1],\
              "bulk_visc":[-1,-1],\
              "tau":-1}

    def _setParamsInEsCore(self):
      #if lb_lbfluid_set_tau(self._lb_params["tau"]):
      #  raise Exception("lb_lbfluid_set_tau error")

      if python_lbfluid_set_density(self._lb_params["dens"]):
        raise Exception("lb_lbfluid_set_density error")

      #if lb_lbfluid_set_visc(self._lb_params["visc"]):
      #  raise Exception("lb_lbfluid_set_visc error")

      #if lb_lbfluid_set_agrid(self._lb_params["agrid"]):
      #  raise Exception("lb_lbfluid_set_agrid error")

      #if lb_lbfluid_set_friction(self._lb_params["friction"]):
      #  raise Exception("lb_lbfluid_set_friction error")
  
      #if lb_lbfluid_set_ext_force(1, self._params["ext_force",0], self._params["ext_force",1], self._params["ext_force",2]):
      #   raise Exception("lb_lbfluid_set_ext_force error")

      #if lb_lbfluid_set_bulk_visc(self._lb_params["bulk_visc"]):
      #  raise Exception("lb_lbfluid_set_bulk_visc error")
    
#    property print_vtk_velocity:
#      def __set__(self, char* _filename):
#        if lb_lbfluid_print_vtk_velocity(_filename):
#          raise Exception("lb_lbfluid_print_vtk_velocity error")
#              
#    property print_vtk_boundary:
#      def __set__(self, char* _filename):
#        if lb_lbfluid_print_vtk_boundary(_filename):
#          raise Exception("lb_lbfluid_print_vtk_boundary error")   
#  
#    property print_velocity:
#      def __set__(self, char* _filename):
#        if lb_lbfluid_print_velocity(_filename):
#          raise Exception("lb_lbfluid_print_vtk_velocity error")
#              
#    property print_boundary:
#      def __set__(self, char* _filename):
#        if lb_lbfluid_print_boundary(_filename):
#          raise Exception("lb_lbfluid_print_vtk_boundary error")
#          
#    property checkpoint:
#      def __set__(self, char* checkpoint_filename):
#        self.checkpoint_filename=checkpoint_filename
#        if lb_lbfluid_save_checkpoint(checkpoint_filename, self.checkpoint_binary):
#          raise Exception("lb_lbfluid_save_checkpoint error")
#      def __get__(self):
#        if lb_lbfluid_load_checkpoint(self.checkpoint_filename, self.checkpoint_binary):
#          raise Exception("lb_lbfluid_load_checkpoint error")  
#  
#    property checkpoint_style:
#      def __set__(self, int _binary):
#        self.checkpoint_binary=_binary
#      def __get__(self):
#        return self.checkpoint_binary
    def _activateMethod(self):

      self._setParamsInEsCore()
              
  class DeviceList:
    def __getitem__(self, _dev):
      return _dev
      
    
#    property gamma_odd:
#      def __set__(self, double _gamma_odd):
#        if lb_lbfluid_set_gamma_odd(_gamma_odd):
#          raise Exception("lb_lbfluid_set_gamma_odd error")
#      def __get__(self):
#        cdef double _p_gamma_odd
#        if lb_lbfluid_get_gamma_odd(&_p_gamma_odd):
#          raise Exception("lb_lbfluid_get_gamma_odd error")
#        return _p_gamma_odd
#        
#    property gamma_even:
#      def __set__(self, double _gamma_even):
#        if lb_lbfluid_set_gamma_even(_gamma_even):
#          raise Exception("lb_lbfluid_set_gamma_even error")
#      def __get__(self):
#        cdef double _p_gamma_even
#        if lb_lbfluid_get_gamma_even(&_p_gamma_even):
#          raise Exception("lb_lbfluid_get_gamma_even error")
#        return _p_gamma_even
