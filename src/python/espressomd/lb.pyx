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
from actors cimport Actor
cimport cuda_init
import cuda_init
# cimport lb
from globals cimport *
from copy import deepcopy

####################################################
#
# Actor class
#
####################################################
cdef class HydrodynamicInteraction(Actor):
    def _lb_init(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction must define the _lb_init() method.")

####################################################
#
# LB_FLUID main class
#
####################################################
IF LB_GPU or LB:
    cdef class LB_FLUID(HydrodynamicInteraction):

        ####################################################
        #
        # validate the given parameters on actor initalization
        #
        ####################################################
        def validateParams(self):
            default_params = self.defaultParams()

            IF SHANCHEN:
                if not (self._params["dens"][0] > 0.0 and self._params["dens"][1] > 0.0):
                    raise ValueError(
                        "Density must be two positive double (ShanChen)")
            ELSE:
                if self._params["dens"] == default_params["dens"]:
                    raise Exception("LB_FLUID density not set")
                else:
                    if not (self._params["dens"] > 0.0 and (isinstance(self._params["dens"], float) or isinstance(self._params["dens"], int))):
                        raise ValueError("Density must be one positive double")

        ####################################################
        #
        # list of valid keys for parameters
        #
        ####################################################
        def validKeys(self):
            return "gpu", "agrid", "dens", "fric", "ext_force", "visc", "tau"

        ####################################################
        #
        # list of esential keys required for the fluid
        #
        ####################################################
        def requiredKeys(self):
            init = 0
            if init:
                init = 1
                return ["dens", "agrid", "visc", "tau"]
            else:
                return []

        ####################################################
        #
        # list of default parameters
        #
        ####################################################
        def defaultParams(self):
            IF SHANCHEN:
                return {"gpu": "yes",
                        "agrid": -1.0,
                        "dens": [-1.0, -1.0],
                        "fric": [-1.0, -1.0],
                        "ext_force": [0.0, 0.0, 0.0],
                        "visc": [-1.0, -1.0],
                        "bulk_visc": [-1.0, -1.0],
                        "tau": -1.0}
            ELSE:
                return {"gpu": "no",
                        "agrid": -1.0,
                        "dens": -1.0,
                        "fric": -1.0,
                        "ext_force": [0.0, 0.0, 0.0],
                        "visc": -1.0,
                        "bulk_visc": -1.0,
                        "tau": -1.0}

        ####################################################
        #
        # function that calls wrapper functions which set the parameters at C-Level
        #
        ####################################################
        def _lb_init(self):

            # Check if GPU code should be used
            if (self._params["gpu"] == "yes"):
                py_lattice_switch = 2
                print "Using LB GPU code"
            else:
                py_lattice_switch = 1
                print "Using LB CPU code"

            if lb_set_lattice_switch(py_lattice_switch):
                raise Exception("lb_set_lattice_switch error")

        def _setParamsInEsCore(self):
            default_params = self.defaultParams()

            if python_lbfluid_set_density(self._params["dens"]):
                raise Exception("lb_lbfluid_set_density error")

            if python_lbfluid_set_tau(self._params["tau"]):
                raise Exception("lb_lbfluid_set_tau error")

            if python_lbfluid_set_visc(self._params["visc"]):
                raise Exception("lb_lbfluid_set_visc error")

            if not self._params["bulk_visc"] == default_params["bulk_visc"]:
                if python_lbfluid_set_bulk_visc(self._params["bulk_visc"]):
                    raise Exception("lb_lbfluid_set_bulk_visc error")

            if python_lbfluid_set_agrid(self._params["agrid"]):
                raise Exception("lb_lbfluid_set_agrid error")

            if not self._params["fric"] == default_params["fric"]:
                if python_lbfluid_set_friction(self._params["fric"]):
                    raise Exception("lb_lbfluid_set_friction error")

            if not self._params["ext_force"] == default_params["ext_force"]:
                if python_lbfluid_set_ext_force(self._params["ext_force"]):
                    raise Exception("lb_lbfluid_set_ext_force error")

        ####################################################
        #
        # function that calls wrapper functions which get the parameters from C-Level
        #
        ####################################################
        def _getParamsFromEsCore(self):
            default_params = self.defaultParams()
            params = self._params

            if python_lbfluid_get_density(self._params["dens"]):
                raise Exception("lb_lbfluid_get_density error")

            if python_lbfluid_get_tau(self._params["tau"]):
                raise Exception("lb_lbfluid_set_tau error")

            if python_lbfluid_get_visc(self._params["visc"]):
                raise Exception("lb_lbfluid_set_visc error")

            if not self._params["bulk_visc"] == default_params["bulk_visc"]:
                if python_lbfluid_get_bulk_visc(self._params["bulk_visc"]):
                    raise Exception("lb_lbfluid_set_bulk_visc error")

            if python_lbfluid_get_agrid(self._params["agrid"]):
                raise Exception("lb_lbfluid_set_agrid error")

            if not self._params["fric"] == default_params["fric"]:
                if python_lbfluid_get_friction(self._params["fric"]):
                    raise Exception("lb_lbfluid_set_friction error")

            if not self._params["ext_force"] == default_params["ext_force"]:
                if python_lbfluid_get_ext_force(self._params["ext_force"]):
                    raise Exception("lb_lbfluid_set_ext_force error")

            return params

        ####################################################
        #
        # redef setParams for lb use
        #
        ####################################################

        def setParams(self, **_params):
            """Update parameters. Only given """
            # Check, if any key was passed, which is not known
            for k in _params.keys():
                if k not in self.validKeys():
                    raise ValueError(
                        "Only the following keys are supported: " + self.validKeys().__str__())

            print "_params", _params

            self.validateParams()

            if "dens" in _params:
                if python_lbfluid_set_density(_params["dens"]):
                    raise Exception("lb_lbfluid_set_density error")

            if "tau" in _params:
                if python_lbfluid_set_tau(_params["tau"]):
                    raise Exception("lb_lbfluid_set_tau error")

            if "visc" in _params:
                if python_lbfluid_set_visc(_params["visc"]):
                    raise Exception("lb_lbfluid_set_visc error")

            if "bulk_visc" in _params:
                if python_lbfluid_set_bulk_visc(_params["bulk_visc"]):
                    raise Exception("lb_lbfluid_set_bulk_visc error")

            if "agrid" in _params:
                if python_lbfluid_set_agrid(_params["agrid"]):
                    raise Exception("lb_lbfluid_set_agrid error")

            if "fric" in _params:
                if python_lbfluid_set_friction(_params["fric"]):
                    raise Exception("lb_lbfluid_set_friction error")

            if "ext_force" in _params:
                if python_lbfluid_set_ext_force(_params["ext_force"]):
                    raise Exception("lb_lbfluid_set_ext_force error")

            if "gpu" in _params:
                # Check if GPU code should be used
                if (_params["gpu"] == "yes"):
                    py_lattice_switch = 2
                    print "Using LB GPU code"
                else:
                    py_lattice_switch = 1
                    print "Using LB CPU code"

                if lb_set_lattice_switch(py_lattice_switch):
                    raise Exception("lb_set_lattice_switch error")

            # update paras
            self._params.update(_params)

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
    #  def _activateMethod(self):

    #    self._setParamsInEsCore()

    # class DeviceList:
    #  def __getitem__(self, _dev):
    #    return _dev


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

        ####################################################
        #
        # Activate Actor
        #
        ####################################################
        def _activateMethod(self):

            self._setParamsInEsCore()

        ####################################################
        #
        # Deactivate Actor
        #
        ####################################################
        # def _deactivateMethod(self):

        #   self._setParamsInEsCore()
