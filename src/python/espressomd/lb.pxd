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


IF LB_GPU == 1:

    cdef extern from "../src/config.hpp":
        pass

    # cdef extern from "../src/lattice.hpp":
    #  int lattice_switch

IF LB_GPU or LB:

###############################################
#
# Wrapper-functions for access to C-pointer
#
###############################################
  cdef inline python_lbfluid_set_density(p_dens):

    cdef double c_dens[2]
    #get pointers
    if isinstance(p_dens,float):
      c_dens[0]=p_dens
      c_dens[1]=p_dens
    else:
      return 1 
    #call c-function
    if(lb_lbfluid_set_density(c_dens)):
      raise Exception("lb_fluid_set_density error at C-level interface")

    return 0

##############################################
#
#extern functions and structs
#
##############################################

  cdef extern from "lb.hpp":

##############################################
#
#Python LB struct clone of C-struct
#
##############################################
    ctypedef struct LB_parameters:
      double rho[2]
      double viscosity[2]
      double bulk_viscosity[2]
      double agrid
      double tau
      double friction[2]
      double ext_force[3]
      double rho_lb_units[2]
      double gamma_odd[2]
      double gamma_even[2]
      int resent_halo
###############################################
#
# init struct
#
###############################################
    ctypedef struct LB_parameters:
      LB_parameters lb_params

##############################################
#
# exported C-functions from lb.hpp
#
##############################################
    int lb_lbfluid_set_tau(double c_tau)
    int lb_lbfluid_set_density(double* c_dens)
    int lb_lbfluid_get_density(double* c_dens)
    int lb_lbfluid_set_visc(double* c_visc)
    int lb_lbfluid_get_visc(double* c_visc)
    int lb_lbfluid_set_agrid(double c_agrid)
    int lb_lbfluid_get_agrid(double* c_agrid)
    int lb_lbfluid_set_friction(double* c_friction)
    int lb_lbfluid_get_friction(double* c_friction)
    int lb_lbfluid_set_gamma_odd(double* c_gamma_odd)
    int lb_lbfluid_get_gamma_odd(double* c_gamma_odd)
    int lb_lbfluid_set_gamma_even(double* c_gamma_even)
    int lb_lbfluid_get_gamma_even(double* c_gamma_even)
    int lb_lbfluid_set_ext_force(int component, double c_fx, double c_fy, double c_fz)
    int lb_lbfluid_get_ext_force(int component, double* c_fx, double* c_fy, double* c_fz)
    int lb_lbfluid_set_bulk_visc(double* c_bulk_visc)
    int lb_lbfluid_get_bulk_visc(double* c_bulk_visc)
    int lb_lbfluid_print_vtk_velocity(char* filename)
    int lb_lbfluid_print_vtk_boundary(char* filename)
    int lb_lbfluid_print_velocity(char* filename)
    int lb_lbfluid_print_boundary(char* filename)
    int lb_lbfluid_save_checkpoint(char* filename, int binary)
    int lb_lbfluid_load_checkpoint(char* filename, int binary)
#    void cython_lb_init(int switch)
