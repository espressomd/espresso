# Copyright (C) 2010-2019 The ESPResSo project
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
include "myconfig.pxi"
from .interactions cimport bonded_ia_params
from espressomd.utils cimport handle_errors
from espressomd.utils import is_valid_type

cdef class Diamond:
    """
    Class to create a diamond-like polymer network.

    Parameters
    ----------
    a : :obj:`float`
        Size of the unit cell.
    bond_length : :obj:`float`
        Distance between adjacent monomers in the chains.
    MPC : :obj:`int`
        Monomers per chain.
    cM_dist : :obj:`int`, optional
        Distance between charged monomers.
    N_CI : :obj:`int`, optional
        Number of counter ions.
    val_nodes : :obj:`float`, optional
        Charge valency of the 8 node particles (crosslinker).
    val_cM : :obj:`float`, optional
        Valency of the charge bearing monomers.
    val_CI : :obj:`float`, optional
        Valency of the counterions.
    nonet : :obj:`bool`, optional
        False creates network, True does not crosslink the individual polymers.

    """

    def __init__(self, *args, **kwargs):
        self._params = self.default_params()
        for k in self.required_keys():
            if k not in kwargs:
                raise ValueError("At least the following keys have to be given as keyword arguments: " +
                                 self.required_keys().__str__() + " got " + kwargs.__str__())
        for k in kwargs:
            if k in self.valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a valid key" % k)

        self.validate_params()
        self._set_params_in_es_core()

    def default_params(self):
        return {"a": 0.0, "bond_length": 0.0, "MPC": 0, "N_CI": 0,
                "val_nodes": 0.0, "val_cM": 0.0, "val_CI": 0.0, "cM_dist": 1, "nonet": False}

    def required_keys(self):
        return "a", "bond_length", "MPC"

    def valid_keys(self):
        return "a", "bond_length", "MPC", "N_CI", "val_nodes", "val_cM", "val_CI", "cM_dist", "nonet"

    def validate_params(self):
        valid_keys = self.valid_keys()
        if(self._params[valid_keys[1]] < 0):
            raise ValueError("Bond length must be positive got: ",
                             self._params[valid_keys[1]])
        if(self._params[valid_keys[2]] < 0 or not is_valid_type(self._params[valid_keys[2]], int)):
            raise ValueError(
                "Monomers per chain must be positive integer got: ", self._params[valid_keys[2]])
        if(self._params[valid_keys[3]] < 0 or not is_valid_type(self._params[valid_keys[3]], int)):
            raise ValueError(
                "The number of counterions must be integer got:", self._params[valid_keys[3]])
        if(not is_valid_type(self._params[valid_keys[4]], float)):
            raise ValueError(
                "The charge of the nodes must be double got", self._params[valid_keys[4]])
        if(not is_valid_type(self._params[valid_keys[5]], float)):
            raise ValueError(
                "The charge of the monomers must be double got", self._params[valid_keys[5]])
        if(not is_valid_type(self._params[valid_keys[6]], float)):
            raise ValueError(
                "The charge of the counterions must be double got", self._params[valid_keys[6]])
        if(not is_valid_type(self._params[valid_keys[7]], int)):
            raise ValueError(
                "The distance between two charged monomers' indices must be integer ", self._params[valid_keys[7]])
        if(self._params[valid_keys[7]] == "nonet"):
            self._params[valid_keys[7]] = True
        if(bonded_ia_params.size() == 0):
            raise ValueError(
                "Please define a bonded interaction [0] before setting up polymers!")

    def __set_params_in_es_core(self):
        return create_diamond(
            partCfg(), self._params["a"], self._params["bond_length"],
            self._params["MPC"], self._params["N_CI"],
            self._params["val_nodes"], self._params["val_cM"],
            self._params["val_CI"], self._params["cM_dist"],
            int(self._params["nonet"]))

    def _set_params_in_es_core(self):
        error_code = self.__set_params_in_es_core()
        handle_errors("Failed changing bonds in create_diamond")
        if error_code == 0:
            pass
        elif error_code == -3:
            raise Exception(
                "Failed upon creating one of the monomers in Espresso!\nAborting...\n")
        else:
            raise Exception("Unknown error code: {}".format(error_code))
