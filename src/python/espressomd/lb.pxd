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


# IF LB_GPU == 1:

# cdef extern from "../src/config.hpp":
#    pass

# cdef extern from "../src/lattice.hpp":
#  int lattice_switch

IF LB_GPU or LB:

    ###############################################
    #
    # Wrapper-functions for access to C-pointer: Set params
    #
    ###############################################
    cdef inline python_lbfluid_set_density(p_dens):

        IF SHANCHEN:
            cdef double c_dens[2]
        ELSE:
            cdef double c_dens[1]

        # get pointers
        if isinstance(p_dens, float) or isinstance(p_dens, int):
            c_dens[0] = <float > p_dens
        else:
            c_dens = p_dens
        # call c-function
        if(lb_lbfluid_set_density(c_dens)):
            raise Exception("lb_fluid_set_density error at C-level interface")

        return 0

###############################################

    cdef inline python_lbfluid_set_tau(p_tau):

        cdef double c_tau
        # get pointers
        c_tau = p_tau
        # call c-function
        if(lb_lbfluid_set_tau(c_tau)):
            raise Exception("lb_fluid_set_tau error at C-level interface")

        return 0

###############################################

    cdef inline python_lbfluid_set_visc(p_visc):

        IF SHANCHEN:
            cdef double c_visc[2]
        ELSE:
            cdef double c_visc[1]
        # get pointers
        if isinstance(p_visc, float) or isinstance(p_visc, int):
            c_visc[0] = <float > p_visc
        else:
            c_visc = p_visc
        # call c-function
        if(lb_lbfluid_set_visc(c_visc)):
            raise Exception("lb_fluid_set_visc error at C-level interface")

        return 0

###############################################

    cdef inline python_lbfluid_set_agrid(p_agrid):

        cdef double c_agrid
        # get pointers
        c_agrid = p_agrid
        # call c-function
        if(lb_lbfluid_set_tau(c_agrid)):
            raise Exception("lb_fluid_set_agrid error at C-level interface")

        return 0

###############################################

    cdef inline python_lbfluid_set_bulk_visc(p_bvisc):

        IF SHANCHEN:
            cdef double c_bvisc[2]
        ELSE:
            cdef double c_bvisc[1]
        # get pointers
        if isinstance(p_bvisc, float) or isinstance(p_bvisc, int):
            c_bvisc[0] = <float > p_bvisc
        else:
            c_bvisc = p_bvisc
        # call c-function
        if(lb_lbfluid_set_bulk_visc(c_bvisc)):
            raise Exception(
                "lb_fluid_set_bulk_visc error at C-level interface")

        return 0

###############################################

    cdef inline python_lbfluid_set_friction(p_friction):

        IF SHANCHEN:
            cdef double c_friction[2]
        ELSE:
            cdef double c_friction[1]
        # get pointers
        if isinstance(p_friction, float) or isinstance(p_friction, int):
            c_friction[0] = <float > p_friction
        else:
            c_friction = p_friction
        # call c-function
        if(lb_lbfluid_set_friction(c_friction)):
            raise Exception("lb_fluid_set_friction error at C-level interface")

        return 0

###############################################

    cdef inline python_lbfluid_set_ext_force(p_ext_force):

        cdef double c_ext_force[3]
        # get pointers
        c_ext_force = p_ext_force
        # call c-function
        if(lb_lbfluid_set_ext_force(1, c_ext_force[0], c_ext_force[1], c_ext_force[2])):
            raise Exception(
                "lb_fluid_set_ext_force error at C-level interface")

        return 0

###############################################


###############################################
#
# Wrapper-functions for access to C-pointer: Get params
#
###############################################
    cdef inline python_lbfluid_get_density(p_dens):

        IF SHANCHEN:
            cdef double c_dens[2]
        ELSE:
            cdef double c_dens[1]
        # call c-function
        if(lb_lbfluid_get_density(c_dens)):
            raise Exception("lb_fluid_get_density error at C-level interface")
        if isinstance(p_dens, float) or isinstance(p_dens, int):
            p_dens = <double > c_dens[0]
        else:
            p_dens = c_dens

        return 0

###############################################
    cdef inline python_lbfluid_get_tau(p_tau):

        cdef double c_tau[1]
        # call c-function
        if(lb_lbfluid_get_tau(c_tau)):
            raise Exception("lb_fluid_get_tau error at C-level interface")
        p_tau = <double > c_tau[0]

        return 0

###############################################
    cdef inline python_lbfluid_get_visc(p_visc):

        IF SHANCHEN:
            cdef double c_visc[2]
        ELSE:
            cdef double c_visc[1]
        # call c-function
        if(lb_lbfluid_get_visc(c_visc)):
            raise Exception(
                "lb_fluid_get_viscosity error at C-level interface")
        if isinstance(p_visc, float) or isinstance(p_visc, int):
            p_visc = <double > c_visc[0]
        else:
            p_visc = c_visc

        return 0

###############################################
    cdef inline python_lbfluid_get_agrid(p_agrid):

        cdef double c_agrid[1]
        # call c-function
        if(lb_lbfluid_get_agrid(c_agrid)):
            raise Exception("lb_fluid_get_agrid error at C-level interface")
        p_agrid = <double > c_agrid[0]

        return 0

###############################################
    cdef inline python_lbfluid_get_bulk_visc(p_bvisc):

        IF SHANCHEN:
            cdef double c_bvisc[2]
        ELSE:
            cdef double c_bvisc[1]
        # call c-function
        if(lb_lbfluid_get_bulk_visc(c_bvisc)):
            raise Exception(
                "lb_fluid_get_bulk_viscosity error at C-level interface")
        if isinstance(p_bvisc, float) or isinstance(p_bvisc, int):
            p_bvisc = <double > c_bvisc[0]
        else:
            p_bvisc = c_bvisc

        return 0

###############################################
    cdef inline python_lbfluid_get_friction(p_friction):

        IF SHANCHEN:
            cdef double c_friction[2]
        ELSE:
            cdef double c_friction[1]
        # call c-function
        if(lb_lbfluid_get_friction(c_friction)):
            raise Exception("lb_fluid_get_friction error at C-level interface")
        if isinstance(p_friction, float) or isinstance(p_friction, int):
            p_fricition = <double > c_friction[0]
        else:
            p_friction = c_friction

        return 0

###############################################

    cdef inline python_lbfluid_get_ext_force(p_ext_force):

        cdef double c_ext_force[3]
        # call c-function
        if(lb_lbfluid_get_ext_force(c_ext_force)):
            raise Exception(
                "lb_fluid_get_ext_force error at C-level interface")
        p_ext_force = c_ext_force

        return 0

###############################################


##############################################
#
# extern functions and structs
#
##############################################

    cdef extern from "lb.hpp":

        ##############################################
        #
        # Python LB struct clone of C-struct
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
        int lb_lbfluid_get_tau(double * c_tau)
        int lb_lbfluid_set_density(double * c_dens)
        int lb_lbfluid_get_density(double * c_dens)
        int lb_lbfluid_set_visc(double * c_visc)
        int lb_lbfluid_get_visc(double * c_visc)
        int lb_lbfluid_set_agrid(double c_agrid)
        int lb_lbfluid_get_agrid(double * c_agrid)
        int lb_lbfluid_set_friction(double * c_friction)
        int lb_lbfluid_get_friction(double * c_friction)
        int lb_lbfluid_set_gamma_odd(double * c_gamma_odd)
        int lb_lbfluid_get_gamma_odd(double * c_gamma_odd)
        int lb_lbfluid_set_gamma_even(double * c_gamma_even)
        int lb_lbfluid_get_gamma_even(double * c_gamma_even)
        int lb_lbfluid_set_ext_force(int component, double c_fx, double c_fy, double c_fz)
        int lb_lbfluid_get_ext_force(double * c_f)
        int lb_lbfluid_set_bulk_visc(double * c_bulk_visc)
        int lb_lbfluid_get_bulk_visc(double * c_bulk_visc)
        int lb_lbfluid_print_vtk_velocity(char * filename)
        int lb_lbfluid_print_vtk_boundary(char * filename)
        int lb_lbfluid_print_velocity(char * filename)
        int lb_lbfluid_print_boundary(char * filename)
        int lb_lbfluid_save_checkpoint(char * filename, int binary)
        int lb_lbfluid_load_checkpoint(char * filename, int binary)
        int lb_set_lattice_switch(int py_switch)
        int lb_get_lattice_switch(int * py_switch)
