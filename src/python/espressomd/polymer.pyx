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
from . cimport polymer
import numpy as np

def validate_params(_params, default):
    if _params["N_P"] <= 0:
        raise ValueError(
                "N_P has to be a positive Integer" )
    if _params["MPC"] <= 1:
        raise ValueError(
                "MPC has to be a positive Integer larger than 1" )
    if _params["bond_length"] < 0 :
        raise ValueError(
                "bond_length has to be a positive float" )
    if _params["start_id"] < 0:
        raise ValueError(
                "start_id has to be a positive Integer")
    if not isinstance(_params["start_pos"], np.ndarray) or len(_params["start_pos"]) != 3:
        raise ValueError(
                "start_pos has to be an numpy array with 3 Elements" )
    if not isinstance(_params["mode"], int):
        raise ValueError(
                "mode has to be a positive Integer" )
    if _params["shield"] < 0 and default["shield"] != _params["shield"]:
        raise ValueError(
                "shield has to be a positive float")
    if _params["max_tries"] < 0 and default["max_tries"] != _params["max_tries"]:
        raise ValueError(
                "max_tries has to be a positive Integer")
    if not isinstance(_params["val_poly"], float) and default["val_poly"] != _params["val_poly"]:
        raise ValueError(
                "val_poly has to be a float")
    if _params["charge_distance"] < 0:
        raise ValueError(
                "charge_distance has to be a positive Integer")
    if _params["type_poly_neutral"] < 0:
        raise ValueError(
                "type_poly_neutral has to be a nonnegative Integer")
    if _params["type_poly_charged"] < 0:
        raise ValueError(
                "type_poly_charged has to be a nonnegative Integer")
    if _params["angle"] < 0 and default["angle"] != _params["angle"]:
        raise ValueError(
                "angle has to be a positive float")
    if _params["angle2"] < 0 and default["angle2"] != _params["angle2"]:
        raise ValueError(
                "angle2 has to be a positive float")
    if _params["constraints"] < 0 :
        raise ValueError(
                "constraint has to be either 0 or 1" )

# wrapper function to expose to the user interface
def create_polymer(**kwargs):
    params=dict()
    default_params=dict()
    default_params["N_P"] = 0 
    default_params["MPC"] = 0
    default_params["bond_length"] = 0 
    default_params["start_id"] = 0
    default_params["start_pos"] = np.array([0, 0, 0])
    default_params["mode"] = 1 
    default_params["shield"] = 0
    default_params["max_tries"] = 1000
    default_params["val_poly" ] = 0.0
    default_params["charge_distance"] = 1
    default_params["type_poly_neutral"] = 0 
    default_params["type_poly_charged"] = 1
    default_params["angle"] = -1.0
    default_params["angle2"] = -1.0
    default_params["constraints"]=0 
    default_params["pos2"] = np.array([0, 0, 0])

    params = default_params 

    valid_keys=["N_P", "MPC", "bond_length", "bond", "start_id", "start_pos", "mode", "shield", "max_tries", "val_poly", "charge_distance", "type_poly_neutral", "type_poly_charged", "angle", "angle2", "constraints"]
        
    required_keys=["N_P", "MPC", "bond_length", "bond"]

    for k in kwargs:
        if not k in valid_keys:
            raise ValueError("Unknown parameter '%s'" % k)
        params[k] = kwargs[k]

    for k in required_keys:
        if k not in kwargs:
            print(k)
            raise ValueError(
                "At least the following keys have to be given as keyword arguments: " + required_keys.__str__())

    validate_params(params, default_params)
    
    bond_id = params["bond"]._bond_id

    cdef double start_pos[3];
    cdef double start_pos2[3];
    for i in range(3):
        start_pos[i] = params["start_pos"][i]
        start_pos2[i] =params["pos2"][i]

    polymerC(params["N_P"], params["MPC"], params["bond_length"], params["start_id"], \
             start_pos, params["mode"], params["shield"], params["max_tries"], \
             params["val_poly"], params["charge_distance"], params["type_poly_neutral"], \
             params["type_poly_charged"], bond_id, params["angle"], \
             params["angle2"], start_pos2, params["constraints"])
    
    #poly=Polymer(**kwargs)




