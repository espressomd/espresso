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

from __future__ import print_function, absolute_import
include "myconfig.pxi"
from . cimport utils
from espressomd.utils cimport *
from .utils import *
from . cimport polymer

cdef class Polymer(object):
    cdef object _params
    def __init__(self):
        self._params = self.default_params()
        pass

    def __call__(self, *args, **kwargs):
        if len(args) == 0 and len(kwargs) == 0:
            raise ValueError("Required Arguments missing " + self.required_keys().__str__() )

        for k in self.required_keys():
            if k not in kwargs:
                print(k)
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__())

        for k in kwargs:
            self._params[k] = kwargs[k]

        print(self._params)
        self.validate_params()

        cdef double start_pos[3];
        cdef double start_pos2[3];
        for i in range(3):
            start_pos[i] = self._params["start_pos"][i]
            start_pos2[i] = self._params["pos2"][i]

        return polymerC(self._params["N_P"], self._params["MPC"], self._params["bond_length"], self._params["start_id"], \
                 start_pos, self._params["mode"], self._params["shield"], self._params["max_tries"], \
                 self._params["val_poly"], self._params["charge_distance"], self._params["type_poly_neutral"], \
                 self._params["type_poly_charged"], self._params["bond_id"], self._params["angle"], \
                 self._params["angle2"], start_pos2, self._params["constraints"] )


    def validate_params(self):
        default = self.default_params()
        if self._params["N_P"] < 0 and  default["N_P"] != self._params["N_P"] :
            raise ValueError(
                    "N_P has to be a positive Integer" )
        if self._params["MPC"] < 0 and default["N_P"] != self._params["N_P"]:
            raise ValueError(
                    "MPC has to be a positive Integer" )
        if self._params["bond_length"] < 0 :
            raise ValueError(
                    "bond_length has to be a positive float" )
        if self._params["bond_id"] < 0 and default["N_P"] != self._params["N_P"]:
            raise ValueError(
                    "bond_id has to be an existing bonded interaction id" )
        if self._params["start_id"] < 0 and default["N_P"] != self._params["N_P"]:
            raise ValueError(
                    "start_id has to be a positive Integer")
        if not isinstance(self._params["start_pos"], np.ndarray) or len(self._params["start_pos"]) != 3:
            raise ValueError(
                    "start_pos has to be an numpy array with 3 Elements" )
        if not isinstance(self._params["mode"], int) and default["N_P"] != self._params["N_P"]:
            raise ValueError(
                    "mode has to be a positive Integer" )
        if self._params["shield"] < 0 and default["shield"] != self._params["shield"]:
            raise ValueError(
                    "shield has to be a positive float")
        if self._params["max_tries"] < 0 and default["max_tries"] != self._params["max_tries"]:
            raise ValueError(
                    "max_tries has to be a positive Integer")
        if not isinstance(self._params["val_poly"], float) and default["val_poly"] != self._params["val_poly"]:
            raise ValueError(
                    "val_poly has to be a float")
        if self._params["charge_distance"] < 0:
            raise ValueError(
                    "charge_distance has to be a positive Integer")
        if self._params["type_poly_neutral"] < 0:
            raise ValueError(
                    "type_poly_neutral has to be a nonnegative Integer")
        if self._params["type_poly_charged"] < 0:
            raise ValueError(
                    "type_poly_charged has to be a nonnegative Integer")
        if self._params["angle"] < 0 and default["angle"] != self._params["angle"]:
            raise ValueError(
                    "angle has to be a positive float")
        if self._params["angle2"] < 0 and default["angle2"] != self._params["angle2"]:
            raise ValueError(
                    "angle2 has to be a positive float")
        if self._params["constraints"] < 0 :
            raise ValueError(
                    "constraint has to be either 0 or 1" )
        

    def required_keys(self):
        return "N_P", "MPC", "bond_length", "bond_id"

    def valid_keys(self):
        return "N_P", "MPC", "bond_length", "bond_id", "start_id", "start_pos", "mode", "shield", "max_tries", "val_poly", "charge_distance", "type_poly_neutral", "type_poly_charged", "angle", "angle2", "constraints"

    def default_params(self):
        para=dict()
        para["N_P"] = 0 
        para["MPC"] = 0
        para["bond_length"] = 0 
        para["start_id"] = 0
        para["start_pos"] = np.array([0, 0, 0])
        para["mode"] = 1 
        para["shield"] = 0
        para["max_tries"] = 1000
        para["val_poly" ] = 0.0
        para["charge_distance"] = 1
        para["type_poly_neutral"] = 0 
        para["type_poly_charged"] = 1
        para["bond_id"] = 0
        para["angle"] = -1.0
        para["angle2"] = -1.0
        para["constraints"]=0 
        para["pos2"] = np.array([0, 0, 0])
        return para



