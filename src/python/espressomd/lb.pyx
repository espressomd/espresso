#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
from globals cimport *
from copy import deepcopy

# Actor class
####################################################
cdef class HydrodynamicInteraction(Actor):
    def _lb_init(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction must define the _lb_init() method.")

# LBFluid main class
####################################################
IF LB_GPU or LB:
    cdef class LBFluid(HydrodynamicInteraction):

        # validate the given parameters on actor initalization
        ####################################################
        def validate_params(self):
            default_params = self.default_params()

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

        # list of valid keys for parameters
        ####################################################
        def valid_keys(self):
            return "gpu", "agrid", "dens", "fric", "ext_force", "visc", "tau"

        # list of esential keys required for the fluid
        ####################################################
        def required_keys(self):
            return ["dens", "agrid", "visc", "tau"]

        # list of default parameters
        ####################################################
        def default_params(self):
            IF SHANCHEN:
                return {"agrid": -1.0,
                        "dens": [-1.0, -1.0],
                        "fric": [-1.0, -1.0],
                        "ext_force": [0.0, 0.0, 0.0],
                        "visc": [-1.0, -1.0],
                        "bulk_visc": [-1.0, -1.0],
                        "tau": -1.0}
            ELSE:
                return {"agrid": -1.0,
                        "dens": -1.0,
                        "fric": -1.0,
                        "ext_force": [0.0, 0.0, 0.0],
                        "visc": -1.0,
                        "bulk_visc": -1.0,
                        "tau": -1.0}

        # function that calls wrapper functions which set the parameters at C-Level
        ####################################################
        def _lb_init(self):
            py_lattice_switch = 1
            if lb_set_lattice_switch(py_lattice_switch):
                raise Exception("lb_set_lattice_switch error")

        def _set_params_in_es_core(self):
            default_params = self.default_params()

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

        # function that calls wrapper functions which get the parameters from C-Level
        ####################################################
        def _get_params_from_es_core(self):
            default_params = self.default_params()

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

            return self.params

        # input/output function wrappers for whole LB fields
        ####################################################
        def print_vtk_velocity(self, path):
            lb_lbfluid_print_vtk_velocity(path)
        def print_vtk_boundary(self, path):
            lb_lbfluid_print_vtk_boundary(path)
        def print_velocity(self, path):
            lb_lbfluid_print_velocity(path)
        def print_boundary(self, path):
            lb_lbfluid_print_boundary(path)
        def save_checkpoint(self, path, binary):
            lb_lbfluid_save_checkpoint(path, binary)
        def load_checkpoint(self, path, binary):
            lb_lbfluid_load_checkpoint(path, binary)

        # Activate Actor
        ####################################################
        def _activate_method(self):
            self.validate_params()
            self._set_params_in_es_core()
            self._lb_init()




# IF LB_GPU:
#     cdef class LBFluid_GPU(LBFluid):
#         def _lb_init(self):
#             py_lattice_switch = 2
#             if lb_set_lattice_switch(py_lattice_switch):
#                 raise Exception("lb_set_lattice_switch error")
